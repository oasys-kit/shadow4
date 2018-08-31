__authors__ = ["M Sanchez del Rio - ESRF ISDD Advanced Analysis and Modelling"]
__license__ = "MIT"
__date__ = "12/01/2017"

#
# SHADOW Undulator preprocessors implemented in python
#
# this code replaces SHADOW's undul_phot and undul_cdf
#
# It calculates the undulator radiation as a function of energy, theta and phi. Phi is the polar angle.
#
# It uses internal code (no dependency) hacked from pySRU
#  (see SourceUndulatorFactorySrw.py and SourceUndulatorFactoryPysru.py for SRW and native pySRU backends, respectively)
#
#
# Available public methods:
#
#     undul_phot()       : like undul_phot of SHADOW but written in python with internal code hacked from pySRU
#     undul_cdf          : like undul_cdf in SHADOW written internally in python
#
#



import numpy
import numpy as np

# scipy
import scipy.constants as codata
import scipy.integrate

class SourceUndulatorFactory(object):

    #
    # private code for internal undul_phot hacked from pySRU
    #

    # calculate a theorical trajectory in an undulator
    # adapted from pySRU: analytical_trajectory_plane_undulator():
    @staticmethod
    def _pysru_analytical_trajectory_plane_undulator(K=1.87 , gamma=2544.03131115, lambda_u=0.020, Nb_period=10, Nb_point=10, Beta_et=0.99993):


        N = Nb_period * Nb_point + 1
        ku = 2.0 * np.pi / lambda_u
        omega_u = Beta_et * codata.c * ku

        # t
        t = np.linspace(-(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), (lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), N)

        ## x and z
        z = Beta_et * t + ((K / gamma) ** 2) * (1.0 / (8.0 * omega_u)) * np.sin( 2.0 * omega_u*t)
        x = (-(K / (gamma * omega_u)) * np.cos(omega_u*t))
        # # Vx and Vz
        v_z = Beta_et + ((K / gamma) ** 2) * (1.0 / 4.0) * np.cos(2.0 *omega_u*t)
        v_x= (K / (gamma )) * np.sin(omega_u*t)
        # # Ax and Az
        a_z=-omega_u *(K / gamma) ** 2 * 0.5 * np.sin( 2.0 * omega_u*t)
        a_x= (K / (gamma )) * (omega_u ) * np.cos(omega_u*t)
        # y
        y=0.0*t
        v_y=y
        a_y=y
        return np.vstack((t,x,y,z,v_x,v_y,v_z,a_x,a_y,a_z))


    # adapted from pySRU:  energy_radiated_approximation_and_farfield()
    @staticmethod
    def _pysru_energy_radiated_approximation_and_farfield(omega=2.53465927101*10**17,electron_current=1.0,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0, D=None):

        c6 = codata.e * electron_current * 1e-9 / (8.0 * np.pi ** 2 * codata.epsilon_0 * codata.c * codata.h)

        if D is not None:
            c6 /= D**2

        N = trajectory.shape[1]
        # N = trajectory.nb_points()
        if D == None:
            # in radian :
            n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
            X = np.sqrt(x ** 2 + y ** 2 )#TODO a changer
        #in meters :
        else :
            X = np.sqrt(x ** 2 + y ** 2 + D ** 2)
            n_chap = np.array([x, y, D]) / X


        trajectory_t   = trajectory[0]
        trajectory_x   = trajectory[1]
        trajectory_y   = trajectory[2]
        trajectory_z   = trajectory[3]
        trajectory_v_x = trajectory[4]
        trajectory_v_y = trajectory[5]
        trajectory_v_z = trajectory[6]
        # trajectory_a_x = trajectory[7]
        # trajectory_a_y = trajectory[8]
        # trajectory_a_z = trajectory[9]

        E = np.zeros((3,), dtype=np.complex)
        integrand = np.zeros((3, N), dtype=np.complex)
        A1 = (n_chap[1] * trajectory_v_z - n_chap[2] * trajectory_v_y)
        A2 = (-n_chap[0] * trajectory_v_z + n_chap[2] * trajectory_v_x)
        A3 = (n_chap[0] * trajectory_v_y - n_chap[1] * trajectory_v_x)
        Alpha2 = np.exp(
            0. + 1j * omega * (trajectory_t + X / codata.c - n_chap[0] * trajectory_x
                                               - n_chap[1] * trajectory_y - n_chap[2] * trajectory_z))


        integrand[0] -= ( n_chap[1]*A3 - n_chap[2]*A2) * Alpha2
        integrand[1] -= (- n_chap[0]*A3 + n_chap[2]*A1) * Alpha2
        integrand[2] -= ( n_chap[0]*A2 - n_chap[1]*A1) * Alpha2

        for k in range(3):
            # E[k] = np.trapz(integrand[k], self.trajectory.t)
            E[k] = np.trapz(integrand[k], trajectory_t)
        E *= omega * 1j

        terme_bord = np.full((3), 0. + 1j * 0., dtype=np.complex)
        Alpha_1 = (1.0 / (1.0 - n_chap[0] * trajectory_v_x[-1]
                          - n_chap[1] * trajectory_v_y[-1] - n_chap[2] * trajectory_v_z[-1]))
        Alpha_0 = (1.0 / (1.0 - n_chap[0] * trajectory_v_x[0]
                          - n_chap[1] * trajectory_v_y[0] - n_chap[2] * trajectory_v_z[0]))

        terme_bord += ((n_chap[1] * A3[-1] - n_chap[2] * A2[-1]) * Alpha_1 *
                       Alpha2[-1])
        terme_bord -= ((n_chap[1] * A3[0] - n_chap[2] * A2[0]) * Alpha_0 *
                       Alpha2[0])
        E += terme_bord
        E *= c6**0.5
        return E


    #
    # now, the different versions of undul_phot
    #
    @staticmethod
    def undul_phot(E_ENERGY,INTENSITY,LAMBDAU,NPERIODS,K,EMIN,EMAX,NG_E,MAXANGLE,NG_T,NG_P,
                   number_of_trajectory_points=20):


        #
        # calculate trajectory
        #
        angstroms_to_eV = codata.h*codata.c/codata.e*1e10
        gamma = E_ENERGY * 1e9 / 0.511e6
        Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
        Beta_et = Beta * (1.0 - (K / (2.0 * gamma)) ** 2)


        E = np.linspace(EMIN,EMAX,NG_E,dtype=float)
        wavelength_array_in_A = angstroms_to_eV / E
        omega_array = 2*np.pi * codata.c / (wavelength_array_in_A * 1e-10)

        T = SourceUndulatorFactory._pysru_analytical_trajectory_plane_undulator(K=K, gamma=gamma, lambda_u=LAMBDAU, Nb_period=NPERIODS,
                                            Nb_point=number_of_trajectory_points,Beta_et=Beta_et)

        #
        # polar grid
        #
        D = 100.0 # placed far away (100 m)
        theta = np.linspace(0,MAXANGLE,NG_T,dtype=float)
        phi = np.linspace(0,np.pi/2,NG_P,dtype=float)

        Z2 = np.zeros((omega_array.size,theta.size,phi.size))
        POL_DEG = np.zeros_like(Z2)
        for o in range(omega_array.size):
            print("Calculating energy %8.3f eV (%d of %d)"%(E[o],o+1,omega_array.size))
            for t in range(theta.size):
                for p in range(phi.size):
                    R = D / np.cos(theta[t])
                    r = R * np.sin(theta[t])
                    X = r * np.cos(phi[p])
                    Y = r * np.sin(phi[p])
                    ElecField = SourceUndulatorFactory._pysru_energy_radiated_approximation_and_farfield(omega=omega_array[o],electron_current=INTENSITY,trajectory=T , x=X , y=Y, D=D )

                    # pol_deg = np.abs(ElecField[0])**2 / (np.abs(ElecField[0])**2 + np.abs(ElecField[1])**2)
                    pol_deg = np.abs(ElecField[0]) / (np.abs(ElecField[0]) + np.abs(ElecField[1])) # SHADOW definition
                    intensity =  (np.abs(ElecField[0]) ** 2 + np.abs(ElecField[1])** 2 + np.abs(ElecField[2])** 2)


                    #  Conversion from pySRU units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
                    intensity *= (D*1e3)**2 # photons/mm^2 -> photons/rad^2
                    intensity /= 1e-3 * E[o] # photons/o.1%bw -> photons/eV

                    Z2[o,t,p] = intensity
                    POL_DEG[o,t,p] = pol_deg

        return {'radiation':Z2,'polarization':POL_DEG,'photon_energy':E,'theta':theta,'phi':phi,'trajectory':T}

    #
    # undul_cdf
    #
    @staticmethod
    def undul_cdf(undul_phot_dict,method='trapz'):
        #
        # takes the output of undul_phot and calculate cumulative distribution functions
        #

        RN0     = undul_phot_dict['radiation']
        POL_DEG = undul_phot_dict['polarization']
        E       = undul_phot_dict['photon_energy']
        T       = undul_phot_dict['theta']
        P       = undul_phot_dict['phi']

        NG_E,NG_T,NG_P = RN0.shape
        print("undul_cdf: _NG_E,_NG_T,_NG_P, %d  %d %d \n"%(NG_E,NG_T,NG_P))

        # coordinates are polar: multiply by sin(theta) to allow dS= r^2 sin(Theta) dTheta dPhi
        YRN0 = numpy.zeros_like(RN0)
        for e in numpy.arange(NG_E):
            for t in numpy.arange(NG_T):
                for p in numpy.arange(NG_P):
                    YRN0[e,t,p] = RN0[e,t,p] * numpy.sin(T[t])


        if method == "sum":
            RN1 = YRN0.sum(axis=2) * (P[1] - P[0])             # RN1(e,t)
            RN2 = RN1.sum(axis=1)  * (T[1] - T[0])             # RN2(e)
            ZERO  = numpy.cumsum(RN0,axis=2)   * (P[1] - P[0]) # CDF(e,t,p)
            ONE   = numpy.cumsum(RN1,axis=1)   * (T[1] - T[0]) # CDF(e,t)
            if NG_E > 1:
                TWO   = numpy.cumsum(RN2)          * (E[1] - E[0]) # CDF(e)
            else:
                TWO = numpy.array([0.0])

        else:
            RN1 = numpy.trapz(YRN0,axis=2) * (P[1]-P[0])                            # RN1(e,t)
            RN2 = numpy.trapz(RN1,axis=1)  * (T[1]-T[0])                            # RN2(e)
            ZERO  = scipy.integrate.cumtrapz(RN0,initial=0,axis=2)  * (P[1] - P[0]) # CDF(e,t,p)
            ONE   = scipy.integrate.cumtrapz(RN1,initial=0,axis=1)  * (T[1] - T[0]) # CDF(e,t)
            if NG_E > 1:
                TWO   = scipy.integrate.cumtrapz(RN2,initial=0)         * (E[1] - E[0]) # CDF(e)
            else:
                TWO = numpy.array([0.0])

        print("undul_cdf: Shadow ZERO,ONE,TWO: ",ZERO.shape,ONE.shape,TWO.shape)

        if NG_E > 1:
            print("undul_cdf: Total Power emitted in the specified angles is: %g Watts."%( (RN2*E).sum()*(E[1]-E[0])*codata.e) )
        else:
            print("undul_cdf: Total Power emitted in the specified angles is: %g Watts."%( (RN2*E)*codata.e) )

        return {'cdf_EnergyThetaPhi':TWO,
                'cdf_EnergyTheta':ONE,
                'cdf_Energy':ZERO,
                'energy':E,
                'theta':T,
                'phi':P,
                'polarization':POL_DEG}
