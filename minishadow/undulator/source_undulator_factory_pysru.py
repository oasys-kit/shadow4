__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "12/01/2017"

#
# SHADOW Undulator preprocessors implemented in python
#
# this code replaces SHADOW's undul_phot
#
# It calculates the undulator radiation as a function of energy, theta and phi. Phi is the polar angle.
#
# It uses pySRU
#
# Available public function:
#
#     undul_phot_pysru() : like undul_phot of SHADOW but using pySRU
#
#



import numpy as np

# needed by pySRU
# try:
from pySRU.ElectronBeam import ElectronBeam as PysruElectronBeam
from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as PysruUndulator
from pySRU.Simulation import create_simulation
from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC
from pySRU.RadiationFactory import RADIATION_METHOD_APPROX_FARFIELD
# except:
#     print("Failed to import pySRU")

class SourceUndulatorFactoryPysru(object):

    @staticmethod
    def undul_phot(E_ENERGY,INTENSITY,LAMBDAU,NPERIODS,K,EMIN,EMAX,NG_E,MAXANGLE,NG_T,NG_P):

        myelectronbeam = PysruElectronBeam(Electron_energy=E_ENERGY, I_current=INTENSITY)
        myundulator = PysruUndulator(K=K, period_length=LAMBDAU, length=LAMBDAU*NPERIODS)

        #
        # polar grid matrix
        #
        photon_energy = np.linspace(EMIN,EMAX,NG_E,dtype=float)

        intens = np.zeros((NG_E,NG_T,NG_P))
        pol_deg = np.zeros_like(intens)
        theta = np.linspace(0,MAXANGLE,NG_T,dtype=float)
        phi = np.linspace(0,np.pi/2,NG_P,dtype=float)

        D = 100.0 # placed far away (100 m)

        THETA = np.outer(theta,np.ones_like(phi))
        PHI = np.outer(np.ones_like(theta),phi)

        X = (D / np.cos(THETA)) * np.sin(THETA) * np.cos(PHI)
        Y = (D / np.cos(THETA)) * np.sin(THETA) * np.sin(PHI)

        for ie,e in enumerate(photon_energy):
            print("Calculating energy %g eV (%d of %d)"%(e,ie+1,photon_energy.size))
            simulation_test = create_simulation(magnetic_structure=myundulator,electron_beam=myelectronbeam,
                                                magnetic_field=None, photon_energy=e,
                                                traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
                                                rad_method=RADIATION_METHOD_APPROX_FARFIELD,initial_condition=None,
                                                distance=D,
                                                X=X.flatten(),Y=Y.flatten(),XY_are_list=True)

            # TODO: this is not nice: I redo the calculations because I need the electric vectors to get polarization
            #       this should be avoided after refactoring pySRU to include electric field in simulations!!
            electric_field = simulation_test.radiation_fact.calculate_electrical_field(
                simulation_test.trajectory, simulation_test.source, X.flatten(), Y.flatten(), D)


            E = electric_field._electrical_field
            # pol_deg1 = (np.abs(E[:,0])**2 / (np.abs(E[:,0])**2 + np.abs(E[:,1])**2)).flatten()
            pol_deg1 = (np.abs(E[:,0]) / (np.abs(E[:,0]) + np.abs(E[:,1]))).flatten() # SHADOW definition!!

            intens1 = simulation_test.radiation.intensity.copy()
            intens1.shape = (theta.size,phi.size)
            pol_deg1.shape = (theta.size,phi.size)

            #  Conversion from pySRU units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
            intens1 *= (D*1e3)**2 # photons/mm^2 -> photons/rad^2
            intens1 /= 1e-3 * e # photons/o.1%bw -> photons/eV

            intens[ie] = intens1
            pol_deg[ie] = pol_deg1

            T0 = simulation_test.trajectory
            T = np.vstack((T0.t,T0.x,T0.y,T0.z,T0.v_x,T0.v_y,T0.v_z,T0.a_x,T0.a_y,T0.a_z))

        return {'radiation':intens,'polarization':pol_deg,'photon_energy':photon_energy,'theta':theta,'phi':phi,'trajectory':T}


