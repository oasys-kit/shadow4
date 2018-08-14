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
# It can use three backends:
#      - the internal code (no dependency) hacked from pySRU
#      - pySRU
#      - SRW
#
#
# Available public functions:
#
#     undul_phot()       : like undul_phot of SHADOW but written in python with internal code hacked from pySRU
#     undul_phot_pysru() : like undul_phot of SHADOW but using pySRU
#     undul_phot_srw()   : like undul_phot of SHADOW but using SRW
#     undul_cdf          : like undul_cdf in SHADOW written internally in python
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
try:
    import srwlib as sl
    import array
except:
    print("Failed to import SRW")

# needed by pySRU
try:
    from pySRU.ElectronBeam import ElectronBeam
    from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as Undulator
    from pySRU.Simulation import create_simulation
    from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC
    from pySRU.RadiationFactory import RADIATION_METHOD_APPROX_FARFIELD
except:
    print("Failed to import pySRU")

angstroms_to_eV = codata.h*codata.c/codata.e*1e10


#
# private code for internal undul_phot hacked from pySRU
#

# calculate a theorical trajectory in an undulator
# adapted from pySRU: analytical_trajectory_plane_undulator():
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
# now, the different versions of undul_phot
#

def undul_phot(E_ENERGY,INTENSITY,LAMBDAU,NPERIODS,K,EMIN,EMAX,NG_E,MAXANGLE,NG_T,NG_P,
               number_of_trajectory_points=20):


    #
    # calculate trajectory
    #
    gamma = E_ENERGY * 1e9 / 0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = Beta * (1.0 - (K / (2.0 * gamma)) ** 2)


    E = np.linspace(EMIN,EMAX,NG_E,dtype=float)
    wavelength_array_in_A = angstroms_to_eV / E
    omega_array = 2*np.pi * codata.c / (wavelength_array_in_A * 1e-10)

    T = _pysru_analytical_trajectory_plane_undulator(K=K, gamma=gamma, lambda_u=LAMBDAU, Nb_period=NPERIODS,
                                        Nb_point=number_of_trajectory_points,Beta_et=Beta_et)

    #
    # polar grid
    #
    D = 100.0 # placed far away (100 m)
    theta = np.linspace(0,MAXANGLE*1e-3,NG_T,dtype=float)
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
                ElecField = _pysru_energy_radiated_approximation_and_farfield(omega=omega_array[o],electron_current=INTENSITY,trajectory=T , x=X , y=Y, D=D )

                pol_deg = np.abs(ElecField[0])**2 / (np.abs(ElecField[0])**2 + np.abs(ElecField[1])**2)
                intensity =  (np.abs(ElecField[0]) ** 2 + np.abs(ElecField[1])** 2 + np.abs(ElecField[2])** 2)


                #  Conversion from pySRU units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
                intensity *= (D*1e3)**2 # photons/mm^2 -> photons/rad^2
                intensity /= 1e-3 * E[o] # photons/o.1%bw -> photons/eV

                Z2[o,t,p] = intensity
                POL_DEG[o,t,p] = pol_deg

    return {'radiation':Z2,'polarization':POL_DEG,'photon_energy':E,'theta':theta,'phi':phi,'trajectory':T}

def undul_phot_pysru(E_ENERGY,INTENSITY,LAMBDAU,NPERIODS,K,EMIN,EMAX,NG_E,MAXANGLE,NG_T,NG_P):


    myelectronbeam = ElectronBeam(Electron_energy=E_ENERGY, I_current=INTENSITY)
    myundulator = Undulator(K=K, period_length=LAMBDAU, length=LAMBDAU*NPERIODS)

    #
    # polar grid matrix
    #
    photon_energy = np.linspace(EMIN,EMAX,NG_E,dtype=float)

    intens = np.zeros((NG_E,NG_T,NG_P))
    pol_deg = np.zeros_like(intens)
    theta = np.linspace(0,MAXANGLE*1e-3,NG_T,dtype=float)
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
        pol_deg1 = (np.abs(E[:,0])**2 / (np.abs(E[:,0])**2 + np.abs(E[:,1])**2)).flatten()

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
    # _srw_stokes0_to_spec(stk,fname="srw_xshundul.spec")

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

          #  Conversion from SRW units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
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
