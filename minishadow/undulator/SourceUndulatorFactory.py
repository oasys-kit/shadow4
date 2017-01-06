#
# Working version for Undulator preprocessors in python
#
# this code replaces SHADOW's undul_phot.
#
# It calculates the undulator radiation as a function of energy, theta and phi. Phi is the polar angle.
#
# Inputs: from xshundul.json (written by ShadowVui)
#
#
# Available public functions:
#
#     undul_phot()       : like undul_phot of SHADOW but written in python with internal code hacked from pySRU
#     undul_phot_pysru() : like undul_phot of SHADOW but using pySRU
#     undul_phot_srw()   : like undul_phot of SHADOW but using SRW
#     undul_cdf          : like undul_cdf in SHADOW written internally in python
#
# TODO:
#    make integral of flux
#    Calculate degree of polarization (set to one now)
#
#



import numpy
import numpy as np
import sys

# scipy
import scipy.constants as codata
from scipy import interpolate
import scipy.integrate

# needed by srw
import srwlib as sl
import array


# input/output
from SourceUndulatorInputOutput import load_uphot_dot_dat,write_uphot_dot_dat
from SourceUndulatorInputOutput import load_xshundul_dot_sha,write_xshundul_dot_sha
import json

angstroms_to_eV = codata.h*codata.c/codata.e*1e10



#
# private code for internal undul_phot hacked from pySRU
#

def _pysru_trajectory_undulator_reference(K=1.87 , gamma=2544.03131115, lambda_u=0.020, Nb_period=10, Nb_point=10, Beta_et=0.99993):

    # if Nb_point == None:
    #     Nb_point = choose_nb_pts_trajectory(K, lambda_u*Nb_period, gamma, omega, alpha=0.01,)
    #     print("Number of trajectory points (calculated automatically)",Nb_point)

    N = Nb_period * Nb_point + 1
    ku= 2.0*np.pi/lambda_u
    omega_u = Beta_et*codata.c *ku
    # trajectory =
    #   [t........]
    # 	[ X/c......]
    # 	[ Y/c ......]
    #	[ Z/c ......]
    # 	[ Vx/c .....]
    # 	[ Vy/c .....]
    # 	[ Vz/c .....]
    # 	[ Ax/c .....]
    # 	[ Ay/c .....]
    # 	[ Az/c .....]
    trajectory = np.zeros((10, N))
    # t
    trajectory[0] = np.linspace(-(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), (lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), N)
    #trajectory[0] = np.linspace(0.0,(lambda_u / (c * Beta_et)) * (Nb_period), N)
    # X et Z en fct de t
    trajectory[3] = Beta_et*trajectory[0] - ((K/gamma)**2) * (1.0/(8.0*ku*codata.c))* np.sin(2.0*omega_u * trajectory[0])
    trajectory[1] = (K/(gamma*ku*codata.c))* np.sin(omega_u * trajectory[0])
    # Vx et Vz en fct de t
    trajectory[6] = Beta_et - ((K/gamma)**2) * (2.0*omega_u/(8.0*ku*codata.c ))*np.cos(2.0*omega_u * trajectory[0])
    trajectory[4] = (K/(gamma*ku*codata.c))*omega_u* np.cos(omega_u * trajectory[0])
    # Ax et Az en fct de t
    trajectory[9] = ((2.0*omega_u*K/gamma)**2) * (1.0/(8.0*ku*codata.c))*np.sin(2.0*omega_u * trajectory[0])
    trajectory[7] = -(K/(gamma*ku*codata.c))*(omega_u**2)* np.sin(omega_u * trajectory[0])
    # trajectory *= codata.c
    # trajectory[0] *= (1.0/codata.c)
    return trajectory



def _pysru_energy_radiated(omega=2.53465927101*10**17,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0, D=None):
    N = trajectory.shape[1]

    if D == None:
        # in radian :
        n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
        X = np.sqrt(x ** 2 + y ** 2 )#TODO a changer
    else:
        # in meters :
        X = np.sqrt(x ** 2 + y ** 2 + D ** 2)
        n_chap = np.array([x, y, D]) / X


    E = np.zeros((3,), dtype=np.complex)
    integrand = np.zeros((3, N), dtype=np.complex)

    A1 = ( n_chap[1] * trajectory[6] - n_chap[2] * trajectory[5])
    A2 = (-n_chap[0] * trajectory[6] + n_chap[2] * trajectory[4])
    A3 = ( n_chap[0] * trajectory[5] - n_chap[1] * trajectory[4])
    Alpha2 = np.exp(
        0. + 1j * omega * (trajectory[0] + X / codata.c - n_chap[0] * trajectory[1]
                                           - n_chap[1] * trajectory[2] - n_chap[2] * trajectory[3]))

    integrand[0] -= ( n_chap[1]*A3 - n_chap[2]*A2) * Alpha2
    integrand[1] -= (- n_chap[0]*A3 + n_chap[2]*A1) * Alpha2
    integrand[2] -= ( n_chap[0]*A2 - n_chap[1]*A1) * Alpha2



    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
    E *= omega * 1j

    # terme_bord = np.full((3), 0. + 1j * 0., dtype=np.complex)
    # Alpha_1 = (1.0 / (1.0 - n_chap[0] * trajectory[4][-1]
    #                   - n_chap[1] * trajectory[5][-1] - n_chap[2] * trajectory[6][-1]))
    # Alpha_0 = (1.0 / (1.0 - n_chap[0] * trajectory[4][0]
    #                   - n_chap[1] * trajectory[5][0] - n_chap[2] * trajectory[6][0]))
    #
    # terme_bord += ((n_chap[1] * A3[-1] - n_chap[2] * A2[-1]) * Alpha_1 *
    #                Alpha2[-1])
    # terme_bord -= ((n_chap[1] * A3[0] - n_chap[2] * A2[0]) * Alpha_0 *
    #                Alpha2[0])
    # E += terme_bord

    return E
    # pol_deg = (np.abs(E[0])**2 / (np.abs(E[0])**2 + np.abs(E[1])**2 ) )
    # return (np.abs(E[0]) ** 2 + np.abs(E[1])** 2 + np.abs(E[2])** 2), pol_deg


#
# private code used in undul_phot_srw
#
def _srw_electron_beam(x=0., y=0., z=0., xp=0., yp=0., e=6.04, Iavg=0.2, sigX=345e-6*1.e-20, sigY=23e-6*1.e-20, mixX=0.0, mixY=0.0, sigXp=4.e-9*1.e-20/345e-6, sigYp=4.e-11*1.e-20/23e-6, sigE = 1.e-4):
  el_rest = 0.51099890221e-03
  eBeam = sl.SRWLPartBeam()
  eBeam.Iavg = Iavg
  eBeam.partStatMom1.x     =  x
  eBeam.partStatMom1.y     =  y
  eBeam.partStatMom1.z     =  z
  eBeam.partStatMom1.xp    =  xp
  eBeam.partStatMom1.yp    =  yp
  eBeam.partStatMom1.gamma =  e/el_rest
  eBeam.partStatMom1.relE0 =  1.0
  eBeam.partStatMom1.nq    = -1
  eBeam.arStatMom2[ 0] = sigX**2  #from here it is not necessary for Single Electron calculation, obviously....
  eBeam.arStatMom2[ 1] = mixX
  eBeam.arStatMom2[ 2] = sigXp**2
  eBeam.arStatMom2[ 3] = sigY**2
  eBeam.arStatMom2[ 4] = mixY
  eBeam.arStatMom2[ 5] = sigYp**2
  eBeam.arStatMom2[10] = sigE**2
  return eBeam


def _srw_drift_electron_beam(eBeam, und ):
  if isinstance(und, float):
    length = und
  elif isinstance(und, sl.SRWLMagFldU):    # Always defined in (0., 0., 0.) move the electron beam before the magnetic field.
    length = 0.0-0.55*und.nPer*und.per-eBeam.partStatMom1.z
  elif isinstance(und, sl.SRWLMagFldC):
    if isinstance(und.arMagFld[0], sl.SRWLMagFldU):
      length = und.arZc[0]-0.55*und.arMagFld[0].nPer*und.arMagFld[0].per-eBeam.partStatMom1.z
    else: raise NameError
  else: raise NameError
  eBeam.partStatMom1.z += length
  eBeam.arStatMom2[0]  += 2*length*eBeam.arStatMom2[1]+length**2*eBeam.arStatMom2[2]
  eBeam.arStatMom2[1]  += length*eBeam.arStatMom2[2]
  eBeam.arStatMom2[3]  += 2*length*eBeam.arStatMom2[4]+length**2*eBeam.arStatMom2[5]
  eBeam.arStatMom2[4]  += length*eBeam.arStatMom2[5]
  eBeam.moved = length
  return eBeam


def _srw_simple_undulator(nPer=72, per=0.0228, B=0.120215, n=1, h_or_v='v'):
  harmB = sl.SRWLMagFldH(n, h_or_v, B)
  und = sl.SRWLMagFldU([harmB], per, nPer)
  return und


def _srw_undulators(und, Xc, Yc, Zc):#for the moment only one works
  cnt = sl.SRWLMagFldC([und], array.array('d', [Xc]), array.array('d', [Yc]), array.array('d', [Zc]))
  return cnt

def _srw_single_electron_source(eBeam, cnt, mesh=sl.SRWLRadMesh(14718.4-1, 14718.4+1., 101, -15.e-6*50*3, 15e-6*50*3, 61, -15e-6*50*3, 15e-6*50*3, 61, 50.),  params=[1, 0.01, 0., 0., 20000, 1, 0]):
  wfr = sl.SRWLWfr()
  wfr.mesh = mesh
  wfr.partBeam = eBeam
  wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
  eBeam = _srw_drift_electron_beam(eBeam, cnt)
  sl.srwl.CalcElecFieldSR(wfr, 0, cnt, params)
  stk = sl.SRWLStokes()
  stk.mesh = mesh
  stk.allocate(mesh.ne, mesh.nx, mesh.ny)
  eBeam = _srw_drift_electron_beam(eBeam, -eBeam.moved)
  wfr.calc_stokes(stk)
  return stk, eBeam


def _srw_multi_electron_source(eBeam, und, mesh=sl.SRWLRadMesh(14718.4, 14718.4, 1, -15.e-6*50, 15e-6*50, 81, -15e-6*50, 15e-6*50, 81, 50.),  params=[1, 9, 1.5, 1.5, 2]):
  stk = sl.SRWLStokes()
  stk.mesh = mesh
  stk.allocate(mesh.ne, mesh.nx, mesh.ny)
  sl.srwl.CalcStokesUR(stk, eBeam, und, params)
  return stk, eBeam

def _srw_stokes0_to_arrays(stk):
  Shape = (4,stk.mesh.ny,stk.mesh.nx,stk.mesh.ne)
  data = numpy.ndarray(buffer=stk.arS, shape=Shape,dtype=stk.arS.typecode)
  data0 = data[0]
  data1 = data[1]
  x = numpy.linspace(stk.mesh.xStart,stk.mesh.xFin,stk.mesh.nx)
  y = numpy.linspace(stk.mesh.yStart,stk.mesh.yFin,stk.mesh.ny)
  e = numpy.linspace(stk.mesh.eStart,stk.mesh.eFin,stk.mesh.ne)
  Z2 = numpy.zeros((e.size,x.size,y.size))
  POL_DEG = numpy.zeros((e.size,x.size,y.size))
  for ie in range(e.size):
      for ix in range(x.size):
          for iy in range(y.size):
            Z2[ie,ix,iy] = data0[iy,ix,ie]
            Ex = numpy.sqrt(numpy.abs(0.5*(data0[iy,ix,ie]+data1[iy,ix,ie])))
            Ey = numpy.sqrt(numpy.abs(0.5*(data0[iy,ix,ie]-data1[iy,ix,ie])))
            # POL_DEG[ie,ix,iy] =  0.5*(data0[iy,ix,ie]+data1[iy,ix,ie]) / data0[iy,ix,ie]
            POL_DEG[ie,ix,iy] =  Ex / (Ex + Ey)
  return Z2,POL_DEG,e,x,y

def _srw_stokes0_to_spec(stk, fname="srw_xshundul.spec"):
  #
  # writes emission in a SPEC file (cartesian grid)
  #
  Shape = (4,stk.mesh.ny,stk.mesh.nx,stk.mesh.ne)
  data = numpy.ndarray(buffer=stk.arS, shape=Shape,dtype=stk.arS.typecode)
  data0 = data[0]
  x = numpy.linspace(stk.mesh.xStart,stk.mesh.xFin,stk.mesh.nx)
  y = numpy.linspace(stk.mesh.yStart,stk.mesh.yFin,stk.mesh.ny)
  e = numpy.linspace(stk.mesh.eStart,stk.mesh.eFin,stk.mesh.ne)
  f = open(fname,"w")
  for k in range(len(e)):
    f.write("#S %d intensity E= %f\n"%(k+1,e[k]))
    f.write("#N 3\n")
    f.write("#L X[m]  Y[m]  Intensity\n")
    for i in range(len(x)):
      for j in range(len(y)):
        f.write( "%e   %e   %e\n"%(x[i], y[j], data0[j,i,k]))
  f.close()
  sys.stdout.write('  file written: srw_xshundul.spec\n')

def _srw_interpol_object(x,y,z):
    #2d interpolation
    if numpy.iscomplex(z[0,0]):
        tck_real = interpolate.RectBivariateSpline(x,y,numpy.real(z))
        tck_imag = interpolate.RectBivariateSpline(x,y,numpy.imag(z))
        return tck_real,tck_imag
    else:
        tck = interpolate.RectBivariateSpline(x,y,z)
        return tck

def _srw_interpol(x,y,z,x1,y1):
    #2d interpolation
    if numpy.iscomplex(z[0,0]):
        tck_real,tck_imag = _srw_interpol_object(x,y,z)
        z1_real = tck_real(numpy.real(x1),numpy.real(y1))
        z1_imag = tck_imag(numpy.imag(x1),numpy.imag(y1))
        return (z1_real+1j*z1_imag)
    else:
        tck = _srw_interpol_object(x,y,z)
        z1 = tck(x1,y1)
        return z1

#
# the different versions of undul_sru
#

def undul_phot(E_ENERGY,INTENSITY,LAMBDAU,NPERIODS,K,EMIN,EMAX,NG_E,MAXANGLE,NG_T,NG_P):

    #
    # calculate trajectory
    #
    gamma = E_ENERGY * 1e9 / 0.511e6
    print("GAMMA:",gamma)
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = Beta * (1.0 - (K / (2.0 * gamma)) ** 2)



    E = np.linspace(EMIN,EMAX,NG_E,dtype=float)
    wavelength_array_in_A = angstroms_to_eV / E
    omega_array = 2*np.pi * codata.c / (wavelength_array_in_A * 1e-10)

    T = _pysru_trajectory_undulator_reference(K=K, gamma=gamma, lambda_u=LAMBDAU, Nb_period=NPERIODS,
                                        Nb_point=20,Beta_et=Beta_et)

    #
    # polar grid
    #
    D = 100.0 # placed far away (100 m)
    theta = np.linspace(0,MAXANGLE*1e-3,NG_T,dtype=float)
    phi = np.linspace(0,np.pi/2,NG_P,dtype=float)

    c6 = codata.e * INTENSITY * 1e-9 / (8.0 * np.pi ** 2 * codata.epsilon_0 * codata.c * codata.h)

    if D is not None:
        c6 /= D**2

    Z2 = np.zeros((omega_array.size,theta.size,phi.size))
    POL_DEG = np.zeros_like(Z2)
    for o in range(omega_array.size):
        print("Calculating energy %g eV (%d of %d)"%(E[o],o+1,omega_array.size))
        for t in range(theta.size):
            for p in range(phi.size):
                R = D / np.cos(theta[t])
                r = R * np.sin(theta[t])
                X = r * np.cos(phi[p])
                Y = r * np.sin(phi[p])
                ElecField = _pysru_energy_radiated(omega=omega_array[o],trajectory=T , x=X , y=Y, D=D )

                pol_deg = np.abs(ElecField[0])**2 / (np.abs(ElecField[0])**2 + np.abs(ElecField[1])**2)
                intensity =  (np.abs(ElecField[0]) ** 2 + np.abs(ElecField[1])** 2 + np.abs(ElecField[2])** 2)


                # TO DO check this conversion
                intensity *= (D*1e3)**2 # photons/mm^2 -> photons/rad^2
                intensity /= 1e-3 * E[o] # photons/o.1%bw -> photons/eV

                Z2[o,t,p] = c6*intensity
                # from photons/mm^2 to photons/rad^2

                POL_DEG[o,t,p] = pol_deg



    return {'radiation':Z2,'polarization':POL_DEG,'photon_energy':E,'theta':theta,'phi':phi}

def undul_phot_pysru(E_ENERGY,INTENSITY,LAMBDAU,NPERIODS,K,EMIN,EMAX,NG_E,MAXANGLE,NG_T,NG_P):

    from pySRU.ElectronBeam import ElectronBeam
    from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
    from pySRU.Simulation import create_simulation
    from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC
    from pySRU.RadiationFactory import RADIATION_METHOD_APPROX_FARFIELD

    beam_ESRF = ElectronBeam(Electron_energy=E_ENERGY, I_current=INTENSITY)
    ESRF18 = Undulator(K=K, period_length=LAMBDAU, length=LAMBDAU*NPERIODS)

    #
    # polar grid matrix
    #
    E = np.linspace(EMIN,EMAX,NG_E,dtype=float)

    intens = np.zeros((NG_E,NG_T,NG_P))
    pol_deg = np.zeros_like(intens)
    theta = np.linspace(0,MAXANGLE*1e-3,NG_T,dtype=float)
    phi = np.linspace(0,np.pi/2,NG_P,dtype=float)

    D = 100.0 # placed far away (100 m)

    THETA = np.outer(theta,np.ones_like(phi))
    PHI = np.outer(np.ones_like(theta),phi)

    X = (D / np.cos(THETA)) * np.sin(THETA) * np.cos(PHI)
    Y = (D / np.cos(THETA)) * np.sin(THETA) * np.sin(PHI)

    for ie,e in enumerate(E):
        print("Calculating energy %g eV (%d of %d)"%(e,ie+1,E.size))
        simulation_test = create_simulation(magnetic_structure=ESRF18,electron_beam=beam_ESRF,
                                            magnetic_field=None, photon_energy=e,
                                            traj_method=TRAJECTORY_METHOD_ANALYTIC,Nb_pts_trajectory=None,
                                            rad_method=RADIATION_METHOD_APPROX_FARFIELD,initial_condition=None,
                                            distance=D,
                                            X=X.flatten(),Y=Y.flatten(),XY_are_list=True)

        # TODO: this is not nice: I redo the calculations because I need the electric vectors to get polarization
        # TODO: this should be evoided after refactoring pySRU to include electric field in simulations!!


        electric_field = simulation_test.radiation_fact.calculate_electrical_field(
            simulation_test.trajectory, simulation_test.source, X.flatten(), Y.flatten(), D)

        # print("<><><><><><><E,X,Y",electric_field._electrical_field.shape,
        #       electric_field._X.shape,electric_field._Y.shape,
        #       X.flatten().shape, Y.flatten().shape)


        E = electric_field._electrical_field
        pol_deg1 = (np.abs(E[:,0])**2 / (np.abs(E[:,0])**2 + np.abs(E[:,1])**2)).flatten()
        # simulation_test.print_parameters()
        intens1 = simulation_test.radiation.intensity.copy()
        intens1.shape = (theta.size,phi.size)
        pol_deg1.shape = (theta.size,phi.size)
        # TO DO check this conversion
        intens1 *= (D*1e3)**2 # photons/mm^2 -> photons/rad^2
        intens1 /= 1e-3 * e # photons/o.1%bw -> photons/eV
        intens[ie] = intens1
        pol_deg[ie] = pol_deg1

    return {'radiation':intens,'polarization':pol_deg,'photon_energy':E,'theta':theta,'phi':phi}

def undul_phot_srw(E_ENERGY,INTENSITY,LAMBDAU,NPERIODS,K,EMIN,EMAX,NG_E,MAXANGLE,NG_T,NG_P):

    lambdau = LAMBDAU
    k = K
    e_energy = E_ENERGY
    nperiods = NPERIODS
    emin = EMIN
    emax = EMAX
    intensity = INTENSITY
    maxangle = MAXANGLE
    sx = 0.0 # h['SX']   #  do not use emittance at this stage
    sz = 0.0 # h['SZ']   #  do not use emittance at this stage
    ex = 0.0 # h['EX']   #  do not use emittance at this stage
    ez = 0.0 # h['EZ']   #  do not use emittance at this stage
    # nrays = h['NRAYS']
    nx = 2*NG_T - 1
    nz = nx
    ne = NG_E # int(ne)

    print ("lambdau = ",lambdau)
    print ("k = ",k)
    print ("e_energy = ",e_energy)
    print ("nperiods = ",nperiods)
    print ("intensity = ",intensity)
    print ("maxangle=%d mrad, (%d x %d points) "%(maxangle,nx,nz))
    print ("sx = ",sx)
    print ("sz = ",sz)
    print ("ex = ",ex)
    print ("ez = ",ez)
    print ("emin =%g, emax=%g, ne=%d "%(emin,emax,ne))


    #
    # define additional parameters needed by SRW
    #
    B = k/93.4/lambdau
    slit_distance = 100.0
    method = "SE" # single-electron  "ME" multi-electron
    sE = 1e-9 # 0.89e-3

    #
    # prepare inputs
    #

    # convert cm to m

    sx *= 1.0e-2
    sz *= 1.0e-2
    ex *= 1.0e-2
    ez *= 1.0e-2

    if sx != 0.0:
      sxp = ex/sx
    else:
      sxp = 0.0

    if sz != 0.0:
      szp = ez/sz
    else:
      szp = 0.0

    xxp = 0.0
    zzp = 0.0

    paramSE = [1, 0.01, 0, 0, 50000, 1, 0]
    paramME = [1, 9, 1.5, 1.5, 2]

    #
    #
    if nx==1 and nz==1: paramME[4] = 1
    params = paramSE if method=="SE" else paramME

    slit_xmin = -maxangle*1.0e-3*slit_distance
    slit_xmax =  maxangle*1.0e-3*slit_distance
    slit_zmin = -maxangle*1.0e-3*slit_distance
    slit_zmax =  maxangle*1.0e-3*slit_distance

    #
    # calculations
    #
    print("nperiods: %d, lambdau: %f, B: %f)"%(nperiods,lambdau,B))

    und = _srw_simple_undulator(nperiods,lambdau,B)
    print("e=%f,Iavg=%f,sigX=%f,sigY=%f,mixX=%f,mixY=%f,sigXp=%f,sigYp=%f,sigE=%f"%(e_energy,intensity,sx,sz,xxp,zzp,sxp,szp,sE) )
    eBeam = _srw_electron_beam(e=e_energy,Iavg=intensity,sigX=sx,sigY=sz,mixX=xxp,mixY=zzp,sigXp=sxp,sigYp=szp,sigE=sE)


    cnt = _srw_undulators(und, 0., 0., 0.)
    sys.stdout.flush()

    mesh = sl.SRWLRadMesh(emin,emax,ne,slit_xmin,slit_xmax,nx,slit_zmin,slit_zmax,nz,slit_distance)
    if (method == 'SE'):
        print ("Calculating SE...")
        stkSE, eBeam = _srw_single_electron_source(eBeam, cnt, mesh, params)
        sys.stdout.write('  done\n')
        sys.stdout.write('  saving SE Stokes...'); sys.stdout.flush()
        stk = stkSE
    else:
        print ("Calculating ME...")
        stkME, eBeam = _srw_multi_electron_source(eBeam, und) # cnt, mesh, params)
        sys.stdout.write('  done\n')
        sys.stdout.write('  saving SE Stokes...'); sys.stdout.flush()
        stk = stkME

    #
    # dump file with radiation on cartesian grid
    #
    _srw_stokes0_to_spec(stk,fname="srw_xshundul.spec")

    #
    # interpolate for polar grid
    #

    # polar grid
    theta = numpy.linspace(0,MAXANGLE*1e-3,NG_T)
    phi = numpy.linspace(0,numpy.pi/2,NG_P)
    Z2 = numpy.zeros((NG_E,NG_T,NG_P))
    POL_DEG = numpy.zeros((NG_E,NG_T,NG_P))

    # interpolate on polar grid
    radiation,pol_deg,e,x,y = _srw_stokes0_to_arrays(stk)
    for ie in range(e.size):
      tck = _srw_interpol_object(x,y,radiation[ie])
      tck_pol_deg = _srw_interpol_object(x,y,pol_deg[ie])
      for itheta in range(theta.size):
        for iphi in range(phi.size):
          R = slit_distance / numpy.cos(theta[itheta])
          r = R * numpy.sin(theta[itheta])
          X = r * numpy.cos(phi[iphi])
          Y = r * numpy.sin(phi[iphi])
          tmp = tck(X,Y)

          # TO DO check this conversion
          tmp *= (slit_distance*1e3)**2 # photons/mm^2 -> photons/rad^2
          tmp /= 1e-3 * e[ie] # photons/o.1%bw -> photons/eV

          Z2[ie,itheta,iphi] = tmp



          POL_DEG[ie,itheta,iphi] = tck_pol_deg(X,Y)

    # !C SHADOW defines the degree of polarization by |E| instead of |E|^2
    # !C i.e.  P = |Ex|/(|Ex|+|Ey|)   instead of   |Ex|^2/(|Ex|^2+|Ey|^2)
    # POL_DEG = numpy.sqrt(POL_DEG2)/(numpy.sqrt(POL_DEG2)+numpy.sqrt(1.0-POL_DEG2))

    return {'radiation':Z2,'polarization':POL_DEG,'photon_energy':e,'theta':theta,'phi':phi}


#
# undul_cdf
#
def undul_cdf(uphot_dot_dat_dict,method='trapz',do_plot=False):
    #
    # read uphot.dat file (like in SHADOW undul_phot_dump)
    #
    # from SourceUndulator import load_uphot_dot_dat
    # uphot_dot_dat_dict = load_uphot_dot_dat(file_in=file_in,do_plot=do_plot)


    RN0     = uphot_dot_dat_dict['radiation']
    POL_DEG = uphot_dot_dat_dict['polarization']
    E       = uphot_dot_dat_dict['photon_energy']
    T       = uphot_dot_dat_dict['theta']
    P       = uphot_dot_dat_dict['phi']

    NG_E,NG_T,NG_P = RN0.shape
    print("undul_cdf: NG_E,NG_T,NG_P, %d  %d %d \n"%(NG_E,NG_T,NG_P))


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
        TWO   = numpy.cumsum(RN2)          * (E[1] - E[0]) # CDF(e)
    else:
        RN1 = numpy.trapz(YRN0,axis=2) * (P[1]-P[0])                            # RN1(e,t)
        RN2 = numpy.trapz(RN1,axis=1)  * (T[1]-T[0])                            # RN2(e)
        ZERO  = scipy.integrate.cumtrapz(RN0,initial=0,axis=2)  * (P[1] - P[0]) # CDF(e,t,p)
        ONE   = scipy.integrate.cumtrapz(RN1,initial=0,axis=1)  * (T[1] - T[0]) # CDF(e,t)
        TWO   = scipy.integrate.cumtrapz(RN2,initial=0)         * (E[1] - E[0]) # CDF(e)



    print("undul_cdf: Shadow ZERO,ONE,TWO: ",ZERO.shape,ONE.shape,TWO.shape)
    print("undul_cdf: Total Power emitted in the specified angles is: %g Watts."%( (RN2*E).sum()*(E[1]-E[0])*codata.e) )


    if do_plot:
        plot(E,TWO,title="NEW TWO",xtitle="E",ytitle="TWO",show=0)
        plot_image(ONE,numpy.arange(NG_E),numpy.arange(NG_T),title="NEW ONE",xtitle="index Energy",ytitle="index Theta",show=0)
        plot_image(ZERO[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="NEW ZERO[0]",xtitle="index Theta",ytitle="index Phi",show=0)
        #plot_image(POL_DEGREE[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="POL_DEGREE[0]",xtitle="index Theta",ytitle="index Phi",show=0)


    return {'cdf_EnergyThetaPhi':TWO,'cdf_EnergyTheta':ONE,'cdf_Energy':ZERO,'energy':E,'theta':T,'phi':P,'polarization':POL_DEG}


#
# Tests
#

# for test purposes only
def run_shadow3_using_preprocessors(jsn):
    from SourceUndulator import SourceUndulator
    u = SourceUndulator()
    u.load_json_shadowvui_dictionary(jsn)
    u.run_using_preprocessors()

def test_undul_phot(h,do_plot=True):
    undul_phot_dict = undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                    LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                    EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                    MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                    NG_P = h["NG_P"])

    write_uphot_dot_dat(undul_phot_dict,file_out="uphot_minishadow.dat")

    undul_phot_srw_dict = undul_phot_srw(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                    LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                    EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                    MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                    NG_P = h["NG_P"])

    write_uphot_dot_dat(undul_phot_srw_dict,file_out="uphot_srw.dat")

    #
    # compare radiation
    #
    d0 = load_uphot_dot_dat("uphot.dat",do_plot=do_plot,show=False)
    d1 = load_uphot_dot_dat("uphot_minishadow.dat",do_plot=do_plot,show=True)
    d2 = load_uphot_dot_dat("uphot_srw.dat",do_plot=do_plot,show=False)

    rad0 = d0['radiation']
    rad1 = d1['radiation']
    rad2 = d2['radiation']

    rad0 /= rad0.max()
    rad1 /= rad1.max()
    rad2 /= rad2.max()


    tmp = numpy.where(rad0 > 0.1)
    print(">>> test_undul_phot: minishadow/shadow3: %4.2f %% "%(numpy.average( 100*numpy.abs((rad1[tmp]-rad0[tmp])/rad0[tmp]) )))
    print(">>> test_undul_phot:        srw/shadow3: %4.2f %% "%(numpy.average( 100*numpy.abs((rad2[tmp]-rad0[tmp])/rad0[tmp]) )))
    print(">>> test_undul_phot: minishadow/srw    : %4.2f %% "%(numpy.average( 100*numpy.abs((rad1[tmp]-rad2[tmp])/rad2[tmp]) )))



def test_undul_cdf(do_plot=True):
    #
    # uphot.dat must exist
    #

    radiation = load_uphot_dot_dat(file_in="uphot.dat")

    cdf2 = undul_cdf(radiation,method='sum',do_plot=False)
    write_xshundul_dot_sha(cdf2,file_out="xshundul2.sha")


    cdf3 = undul_cdf(radiation,method='trapz',do_plot=False)
    write_xshundul_dot_sha(cdf3,file_out="xshundul3.sha")

    cdf1 = load_xshundul_dot_sha(file_in="xshundul.sha", do_plot=do_plot,show=False)
    cdf2 = load_xshundul_dot_sha(file_in="xshundul2.sha",do_plot=do_plot,show=False)
    cdf3 = load_xshundul_dot_sha(file_in="xshundul3.sha",do_plot=do_plot,show=True)

    ZERO1 = cdf1['cdf_Energy']
    ONE1 = cdf1['cdf_EnergyTheta']
    TWO1 = cdf1['cdf_EnergyThetaPhi']

    ZERO2 = cdf2['cdf_Energy']
    ONE2 = cdf2['cdf_EnergyTheta']
    TWO2 = cdf2['cdf_EnergyThetaPhi']

    ZERO3 = cdf3['cdf_Energy']
    ONE3 = cdf3['cdf_EnergyTheta']
    TWO3 = cdf3['cdf_EnergyThetaPhi']

    tmp = numpy.where(ZERO1 > 0.1*ZERO1.max())
    print("test_undul_cdf: ZERO:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ZERO2[tmp]-ZERO1[tmp])/ZERO1[tmp]) )))
    print("test_undul_cdf: ZERO: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ZERO3[tmp]-ZERO1[tmp])/ZERO1[tmp]) )))

    tmp = numpy.where(ONE1 > 0.1*ONE1.max())
    print(r"test_undul_cdf: ONE:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ONE2[tmp]-ONE1[tmp])/ONE1[tmp]) )))
    print(r"test_undul_cdf: ONE: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ONE3[tmp]-ONE1[tmp])/ONE1[tmp]) )))

    tmp = numpy.where(TWO1 > 0.1*TWO1.max())
    print("test_undul_cdf: TWO:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((TWO2[tmp]-TWO1[tmp])/TWO1[tmp]) )))
    print("test_undul_cdf: TWO: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((TWO3[tmp]-TWO1[tmp])/TWO1[tmp]) )))



if __name__ == "__main__":

    from srxraylib.plot.gol import plot_image,plot, plot_show

        # "EMIN":       10500.0000,
        # "EMAX":       10550.0000,
        # "INTENSITY":      0.200000003,
        # "EMIN":       10200.0000,
        # "EMAX":       10650.0000,

    tmp = \
        """
        {
        "LAMBDAU":     0.0320000015,
        "K":      0.250000000,
        "E_ENERGY":       6.03999996,
        "E_ENERGY_SPREAD":    0.00100000005,
        "NPERIODS": 50,
        "EMIN":       10200.0000,
        "EMAX":       10650.0000,
        "INTENSITY":      1.0,
        "MAXANGLE":     0.0149999997,
        "NG_E": 11,
        "NG_T": 51,
        "NG_P": 11,
        "NG_PLOT(1)":"1",
        "NG_PLOT(2)":"No",
        "NG_PLOT(3)":"Yes",
        "UNDUL_PHOT_FLAG(1)":"4",
        "UNDUL_PHOT_FLAG(2)":"Shadow code",
        "UNDUL_PHOT_FLAG(3)":"Urgent code",
        "UNDUL_PHOT_FLAG(4)":"SRW code",
        "UNDUL_PHOT_FLAG(5)":"Gaussian Approx",
        "UNDUL_PHOT_FLAG(6)":"python code by Sophie",
        "SEED": 36255,
        "SX":     0.0399999991,
        "SZ":    0.00100000005,
        "EX":   4.00000005E-07,
        "EZ":   3.99999989E-09,
        "FLAG_EMITTANCE(1)":"1",
        "FLAG_EMITTANCE(2)":"No",
        "FLAG_EMITTANCE(3)":"Yes",
        "NRAYS": 15000,
        "F_BOUND_SOUR": 0,
        "FILE_BOUND":"NONESPECIFIED",
        "SLIT_DISTANCE":       1000.00000,
        "SLIT_XMIN":      -1.00000000,
        "SLIT_XMAX":       1.00000000,
        "SLIT_ZMIN":      -1.00000000,
        "SLIT_ZMAX":       1.00000000,
        "NTOTALPOINT": 10000000,
        "JUNK4JSON":0
        }
        """
    h = json.loads(tmp)


    # test_undul_phot(h,do_plot=True)
    # test_undul_cdf(do_plot=True)

    do_plot_intensity = 1
    do_plot_polarization = 1
    #
    run_shadow3_using_preprocessors(h)
    undul_phot_preprocessor_dict = load_uphot_dot_dat("uphot.dat")

    if do_plot_intensity: plot_image(undul_phot_preprocessor_dict['radiation'][0,:,:],undul_phot_preprocessor_dict['theta']*1e6,undul_phot_preprocessor_dict['phi']*180/numpy.pi,
               title="INTENS UNDUL_PHOT_PREPROCESSOR: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

    if do_plot_polarization: plot_image(undul_phot_preprocessor_dict['polarization'][0,:,:],undul_phot_preprocessor_dict['theta']*1e6,undul_phot_preprocessor_dict['phi']*180/numpy.pi,
               title="POL_DEG UNDUL_PHOT_PREPROCESSOR: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)


    # internal code
    undul_phot_dict = undul_phot(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                    LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                    EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                    MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                    NG_P = h["NG_P"])

    if do_plot_intensity: plot_image(undul_phot_dict['radiation'][0,:,:],undul_phot_dict['theta']*1e6,undul_phot_dict['phi']*180/numpy.pi,
               title="INTENS UNDUL_PHOT: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
    if do_plot_polarization: plot_image(undul_phot_dict['polarization'][0,:,:],undul_phot_dict['theta']*1e6,undul_phot_dict['phi']*180/numpy.pi,
               title="POL_DEG UNDUL_PHOT: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

    # pySRU
    undul_phot_pysru_dict = undul_phot_pysru(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                    LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                    EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                    MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                    NG_P = h["NG_P"])
    if do_plot_intensity: plot_image(undul_phot_pysru_dict['radiation'][0,:,:],undul_phot_pysru_dict['theta']*1e6,undul_phot_pysru_dict['phi']*180/numpy.pi,
               title="INTENS UNDUL_PHOT_PYSRU: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
    if do_plot_polarization: plot_image(undul_phot_pysru_dict['polarization'][0,:,:],undul_phot_pysru_dict['theta']*1e6,undul_phot_pysru_dict['phi']*180/numpy.pi,
               title="POL_DEG UNDUL_PHOT_PYSRU: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

    # srw
    undul_phot_srw_dict = undul_phot_srw(E_ENERGY = h["E_ENERGY"],INTENSITY = h["INTENSITY"],
                                    LAMBDAU = h["LAMBDAU"],NPERIODS = h["NPERIODS"],K = h["K"],
                                    EMIN = h["EMIN"],EMAX = h["EMAX"],NG_E = h["NG_E"],
                                    MAXANGLE = h["MAXANGLE"],NG_T = h["NG_T"],
                                    NG_P = h["NG_P"])
    if do_plot_intensity: plot_image(undul_phot_srw_dict['radiation'][0,:,:],undul_phot_srw_dict['theta']*1e6,undul_phot_srw_dict['phi']*180/numpy.pi,
               title="INTENS UNDUL_PHOT_SRW: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
    if do_plot_polarization: plot_image(undul_phot_srw_dict['polarization'][0,:,:],undul_phot_srw_dict['theta']*1e6,undul_phot_srw_dict['phi']*180/numpy.pi,
               title="POL_DEG UNDUL_PHOT_SRW: RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)


    undul_cdf(undul_phot_srw_dict)

    x = undul_phot_srw_dict["photon_energy"]
    y0 = (undul_phot_preprocessor_dict["radiation"]).sum(axis=2).sum(axis=1)
    y1 = (undul_phot_dict["radiation"]).sum(axis=2).sum(axis=1)
    y2 = (undul_phot_pysru_dict["radiation"]).sum(axis=2).sum(axis=1)
    y3 = (undul_phot_srw_dict["radiation"]).sum(axis=2).sum(axis=1)
    # yrange=[y0.min()*0.9,y3.max()*1.1]
    if do_plot_intensity: plot(x,y0,x,y1,x,y2,x,y3,xtitle="Photon energy [eV]",ytitle="Flux[photons/s/eV/rad^2]",legend=["preprocessor","internal","pySRU","SRW"])
    if do_plot_polarization or do_plot_intensity: plot_show()