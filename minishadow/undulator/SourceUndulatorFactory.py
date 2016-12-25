
import numpy
import scipy.constants as codata
# import json
# import os

# from numpy.testing import assert_equal, assert_almost_equal

import scipy.integrate

# import Shadow

# import platform
# if platform.system() == "Linux":
#     SHADOW3_BINARY = "/users/srio/OASYS_VE/shadow3/shadow3"
# else:
#     SHADOW3_BINARY = "/Users/srio/Oasys/OASYS_VE/shadow3/shadow3"

angstroms_to_eV = codata.h*codata.c/codata.e*1e10

# for undul_phot_srw
# import json
# import numpy
import srwlib as sl
import array
# import sys
from scipy import interpolate




#
# Working version for Undulator preprocessors in python
#
# this code replaces SHADOW's undul_phot.
#
# It calculates the undulator radiation as a function of energy, theta and phi. Phi is the polar angle.
#
# Inputs: from xshundul.json (written by ShadowVui)
# Output: file uphot.dat, like in shadow
#
#
# TODO:
#    Follow changes in Sophie's code
#    Calculate degree of polarization (set to one now)
#
#
import numpy as np
import scipy.constants as codata
import json
import sys

angstroms_to_eV = codata.h*codata.c/codata.e*1e10


#
# undul_phot
# COPY OF SOPHIE'S CODE
#

def trajectory_undulator_reference(K=1.87 , gamma=2544.03131115, lambda_u=0.020, Nb_period=10, Nb_point=10, Beta_et=0.99993):
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

def energy_radiated(omega=2.53465927101*10**17,trajectory=np.zeros((11,10)) , x=0.00 , y=0.0, D=None):
    N = trajectory.shape[1]

    if D == None:
        # in radian :
        n_chap = np.array([x, y, 1.0 - 0.5 * (x ** 2 + y ** 2)])
    else:
        # in meters :
        R = np.sqrt(x ** 2 + y ** 2 + D ** 2)
        n_chap = np.array([x, y, D]) / R

    E = np.full((3,), 0. + 1j * 0., dtype=np.complex)
    integrand = np.full((3,N), 0. + 1j * 0., dtype=np.complex)
    Alpha=trajectory[7]*(n_chap[2]-trajectory[6]) - (n_chap[0]-trajectory[4])*trajectory[9]
    Alpha2=np.exp(0. + 1j * omega * (trajectory[0] - n_chap[0]*trajectory[1]-n_chap[2]*trajectory[3]))
    integrand[0] += (-(n_chap[1]**2)*trajectory[7]-n_chap[2]*Alpha)*Alpha2
    integrand[1] += n_chap[1]*(n_chap[0]*trajectory[7]+n_chap[2]*trajectory[9])*Alpha2
    integrand[2] += (-(n_chap[1]**2)*trajectory[9]+n_chap[0]*Alpha)*Alpha2
    integrand *= (1.0 / (1.0 - n_chap[0]*trajectory[4]-n_chap[2]*trajectory[6])) ** 2
    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory[0])
        #E[k] = integrate.simps(integrand[k], trajectory[0])
    # np.linalg.norm
    pol_deg = (np.abs(E[0])**2 / (np.abs(E[0])**2 + np.abs(E[1])**2 ) )
    return (np.abs(E[0]) ** 2 + np.abs(E[1])** 2 + np.abs(E[2])** 2), pol_deg
#
# END COPY OF SOPHIE'S CODE
#

def undul_phot(myinput):
    #
    # read inputs from a file created by ShadowVUI ----------------------------
    #
    if isinstance(myinput,str):
        inFileTxt = myinput # "xshundul.json"
        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)

    elif isinstance(myinput,dict):
        h = myinput
    else:
        raise Exception("Unknown input")


    # list all non-empty keywords
    print ("-----------------------------------------------------")
    for i,j in h.items():
        if (j != None):
            print ("%s = %s" % (i,j))
    print ("-----------------------------------------------------")
    print ("k: ",h['K'])

    #
    # calculate trajectory
    #
    gamma = h["E_ENERGY"] * 1e9 / 0.511e6
    print("GAMMA:",gamma)
    Beta = np.sqrt(1.0 - (1.0 / gamma ** 2))
    Beta_et = Beta * (1.0 - (h['K'] / (2.0 * gamma)) ** 2)
    T = trajectory_undulator_reference(K=h['K'], gamma=gamma, lambda_u=h["LAMBDAU"], Nb_period=h["NPERIODS"],
                                        Nb_point=20,Beta_et=Beta_et)



    E = np.linspace(h["EMIN"],h["EMAX"],h["NG_E"],dtype=float)
    wavelength_array_in_A = angstroms_to_eV / E
    omega_array = 2*np.pi * codata.c / (wavelength_array_in_A * 1e-10)

    #
    # cartesian grid
    #
    gridding = 1 # 0=cartesian, 1=polar

    if gridding == 0:
        D = None
        X = np.linspace(0.0,h["MAXANGLE"]*1e-3,h["NG_T"],dtype=float)
        Y = np.linspace(0.0,h["MAXANGLE"]*1e-3,h["NG_T"],dtype=float)

    else:
        D = 100.0 # placed far away (100 m)
        theta = np.linspace(0,h["MAXANGLE"]*1e-3,h["NG_T"],dtype=float)
        phi = np.linspace(0,np.pi/2,h["NG_P"],dtype=float)

    c6= codata.e*1e-10/(8.0*np.pi**2*codata.epsilon_0*codata.c*codata.h)

    if gridding == 0:
        Z2 = np.zeros((E.size,X.size,Y.size))
        POL_DEG = np.zeros_like(Z2)
        for o in range(omega_array.size):
            print("Calculating energy %g eV (%d of %d)"%(E[o],o+1,omega_array.size))
            for i in range(X.size):
                for j in range(Y.size):
                    intensity, pol_deg=energy_radiated(omega=omega_array[o],trajectory=T , x=X[i] , y=Y[j], D=D )
                    Z2[o,i,j] = c6*intensity
                    POL_DEG[o,i,j] = pol_deg
    elif gridding == 1:
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
                    intensity, pol_deg = energy_radiated(omega=omega_array[o],trajectory=T , x=X , y=Y, D=D )
                    Z2[o,t,p] = c6*intensity
                    POL_DEG[o,t,p] = pol_deg


    #
    # create uphot.dat file (like in SHADOW undul_phot)
    #
    # file_out = "uphot.dat"
    # f = open(file_out,'w')
    # f.write("%d  %d  %d \n"%(h["NG_E"],h["NG_T"],h["NG_P"]))
    # for e in range(h["NG_E"]):
    #     f.write("%20.10f \n"%(E[e]))
    #
    # for e in range(h["NG_E"]):
    #     for t in range(h["NG_T"]):
    #         f.write("%20.10f \n"%(theta[t]))
    #
    # for e in range(h["NG_E"]):
    #     for t in range(h["NG_T"]):
    #         for p in range(h["NG_P"]):
    #             f.write("%20.10f \n"%(phi[p]))
    #
    #
    # for e in range(h["NG_E"]):
    #     for t in range(h["NG_T"]):
    #         for p in range(h["NG_P"]):
    #             f.write("%20.10f \n"%Z2[e,t,p])
    #
    # for e in range(h["NG_E"]):
    #     for t in range(h["NG_T"]):
    #         for p in range(h["NG_P"]):
    #             f.write("%20.10f \n"%POL_DEG[e,t,p])
    #
    # f.close()
    # print("File written to disk: %s"%file_out)

    return {'radiation':Z2,'polarization':POL_DEG,'photon_energy':E,'theta':theta,'phi':phi}


#
# undul_phot_srw
#



def ElectronBeam(x=0., y=0., z=0., xp=0., yp=0., e=6.04, Iavg=0.2, sigX=345e-6*1.e-20, sigY=23e-6*1.e-20, mixX=0.0, mixY=0.0, sigXp=4.e-9*1.e-20/345e-6, sigYp=4.e-11*1.e-20/23e-6, sigE = 1.e-4):
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


def DriftElectronBeam(eBeam, und ):
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


def SimpleUndulator(nPer=72, per=0.0228, B=0.120215, n=1, h_or_v='v'):
  harmB = sl.SRWLMagFldH(n, h_or_v, B)
  und = sl.SRWLMagFldU([harmB], per, nPer)
  return und


# def Undulator(nPer=72, per=0.0228, B=[0.120215], n=[1], h_or_v=['v']):
#   assert (len(B)==len(n)), "Wrong length of input arrays"
#   assert (len(B)==len(h_or_v)), "Wrong length of input arrays"
#   harms = [ sl.SRWLMagFldH(n[i], h_or_v[i], B[i]) for i in range(len(B)) ]
#   und = sl.SRWLMagFldU(harms, per, nPer)
#   return und
#
#
def Undulators(und, Xc, Yc, Zc):#for the moment only one works
  cnt = sl.SRWLMagFldC([und], array.array('d', [Xc]), array.array('d', [Yc]), array.array('d', [Zc]))
  return cnt


#def SrwSESource(eBeam, cnt, mesh=sl.SRWLRadMesh(12000., 16000., 101, -15.e-6*50*3, 15e-6*50*3, 61, -15e-6*50*3, 15e-6*50*3, 61, 50.),  params=[1, 0.01, 0., 0., 20000, 1, 0]):
def SrwSESource(eBeam, cnt, mesh=sl.SRWLRadMesh(14718.4-1, 14718.4+1., 101, -15.e-6*50*3, 15e-6*50*3, 61, -15e-6*50*3, 15e-6*50*3, 61, 50.),  params=[1, 0.01, 0., 0., 20000, 1, 0]):
  wfr = sl.SRWLWfr()
  wfr.mesh = mesh
  wfr.partBeam = eBeam
  wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
  eBeam = DriftElectronBeam(eBeam, cnt)
  sl.srwl.CalcElecFieldSR(wfr, 0, cnt, params)
  stk = sl.SRWLStokes()
  stk.mesh = mesh
  stk.allocate(mesh.ne, mesh.nx, mesh.ny)
  eBeam = DriftElectronBeam(eBeam, -eBeam.moved)
  wfr.calc_stokes(stk)
  return stk, eBeam


def SrwMESource(eBeam, und, mesh=sl.SRWLRadMesh(14718.4, 14718.4, 1, -15.e-6*50, 15e-6*50, 81, -15e-6*50, 15e-6*50, 81, 50.),  params=[1, 9, 1.5, 1.5, 2]):
#def SrwMESource(eBeam, und, mesh=sl.SRWLRadMesh(1000., 21000., 10001, -15.e-6*50, 15e-6*50, 1, -15e-6*50, 15e-6*50, 1, 50.),  params=[1, 21, 1.5, 1.5, 1]):
  stk = sl.SRWLStokes()
  stk.mesh = mesh
  stk.allocate(mesh.ne, mesh.nx, mesh.ny)
  sl.srwl.CalcStokesUR(stk, eBeam, und, params)
  return stk, eBeam

#
# def save(stk, eBeam, fname="SrwStokes"):
#   pickle.dump( stk, open(fname+"_stk.dat", "wb") )
#   pickle.dump( eBeam, open(fname+"_ebeam.dat", "wb") )

def Stokes0ToArrays(stk):
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

def Stokes0ToSpec(stk, fname="srw_xshundul.spec"):
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

#
#
#
def myinterpol_object(x,y,z):
    #2d interpolation
    if numpy.iscomplex(z[0,0]):
        tck_real = interpolate.RectBivariateSpline(x,y,numpy.real(z))
        tck_imag = interpolate.RectBivariateSpline(x,y,numpy.imag(z))
        return tck_real,tck_imag
    else:
        tck = interpolate.RectBivariateSpline(x,y,z)
        return tck

def myinterpol(x,y,z,x1,y1):
    #2d interpolation
    if numpy.iscomplex(z[0,0]):
        tck_real,tck_imag = myinterpol_object(x,y,z)
        z1_real = tck_real(numpy.real(x1),numpy.real(y1))
        z1_imag = tck_imag(numpy.imag(x1),numpy.imag(y1))
        return (z1_real+1j*z1_imag)
    else:
        tck = myinterpol_object(x,y,z)
        z1 = tck(x1,y1)
        return z1


def undul_phot_srw(myinput):
    #
    # read inputs from a file created by ShadowVUI ----------------------------
    #
    if isinstance(myinput,str):
        inFileTxt = myinput # "xshundul.json"
        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)

    elif isinstance(myinput,dict):
        h = myinput
    else:
        raise Exception("Unknown input")

    print ("-----------------------------------------------------")
    for i,j in h.items():
        if (j != None):
            print ("%s = %s" % (i,j))
    print ("-----------------------------------------------------")
    print ("k: ",h['K'])


    lambdau = h['LAMBDAU']
    k = h['K']
    e_energy = h['E_ENERGY']
    nperiods = h['NPERIODS']
    emin = h['EMIN']
    emax = h['EMAX']
    intensity = h['INTENSITY']
    maxangle = h['MAXANGLE']
    sx = 0.0 # h['SX']   #  do not use emittance at this stage
    sz = 0.0 # h['SZ']   #  do not use emittance at this stage
    ex = 0.0 # h['EX']   #  do not use emittance at this stage
    ez = 0.0 # h['EZ']   #  do not use emittance at this stage
    nrays = h['NRAYS']
    nx = 2*h['NG_T'] - 1
    nz = nx
    ne = h["NG_E"] # int(ne)

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
    print ("nrays = ",nrays)
    print ("emin =%g, emax=%g, ne=%d "%(emin,emax,ne))


    #
    # define additional parameters needed by SRW
    #
    B = k/93.4/lambdau
    slit_distance = 50.0
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

    und = SimpleUndulator(nperiods,lambdau,B)
    print("e=%f,Iavg=%f,sigX=%f,sigY=%f,mixX=%f,mixY=%f,sigXp=%f,sigYp=%f,sigE=%f"%(e_energy,intensity,sx,sz,xxp,zzp,sxp,szp,sE) )
    eBeam = ElectronBeam(e=e_energy,Iavg=intensity,sigX=sx,sigY=sz,mixX=xxp,mixY=zzp,sigXp=sxp,sigYp=szp,sigE=sE)


    cnt = Undulators(und, 0., 0., 0.)
    sys.stdout.flush()

    mesh = sl.SRWLRadMesh(emin,emax,ne,slit_xmin,slit_xmax,nx,slit_zmin,slit_zmax,nz,slit_distance)
    if (method == 'SE'):
        print ("Calculating SE...")
        stkSE, eBeam = SrwSESource(eBeam, cnt, mesh, params)
        sys.stdout.write('  done\n')
        sys.stdout.write('  saving SE Stokes...'); sys.stdout.flush()
        stk = stkSE
    else:
        print ("Calculating ME...")
        stkME, eBeam = SrwMESource(eBeam, und) # cnt, mesh, params)
        sys.stdout.write('  done\n')
        sys.stdout.write('  saving SE Stokes...'); sys.stdout.flush()
        stk = stkME

    #
    # dump file with radiation on cartesian grid
    #
    Stokes0ToSpec(stk,fname="srw_xshundul.spec")

    #
    # interpolate for polar grid
    #

    # polar grid
    theta = numpy.linspace(0,h["MAXANGLE"]*1e-3,h["NG_T"])
    phi = numpy.linspace(0,numpy.pi/2,h["NG_P"])
    Z2 = numpy.zeros((h["NG_E"],h["NG_T"],h["NG_P"]))
    POL_DEG = numpy.zeros((h["NG_E"],h["NG_T"],h["NG_P"]))

    # interpolate on polar grid
    radiation,pol_deg,e,x,y = Stokes0ToArrays(stk)
    for ie in range(e.size):
      tck = myinterpol_object(x,y,radiation[ie])
      tck_pol_deg = myinterpol_object(x,y,pol_deg[ie])
      for itheta in range(theta.size):
        for iphi in range(phi.size):
          R = slit_distance / numpy.cos(theta[itheta])
          r = R * numpy.sin(theta[itheta])
          X = r * numpy.cos(phi[iphi])
          Y = r * numpy.sin(phi[iphi])
          tmp = tck(X,Y)
          Z2[ie,itheta,iphi] = tmp
          POL_DEG[ie,itheta,iphi] = tck_pol_deg(X,Y)

    # !C SHADOW defines the degree of polarization by |E| instead of |E|^2
    # !C i.e.  P = |Ex|/(|Ex|+|Ey|)   instead of   |Ex|^2/(|Ex|^2+|Ey|^2)
    # POL_DEG = numpy.sqrt(POL_DEG2)/(numpy.sqrt(POL_DEG2)+numpy.sqrt(1.0-POL_DEG2))

    #
    # plot at first energy point
    #
    # title = r"K=%4.2f, N=%d, E=%3.2f GeV, $\lambda_u$=%2.0f mm, $E_{ph}$=%7.2f eV, Max photons=%g"%(
    #         h['K'],h["NPERIODS"],h["E_ENERGY"],1e3*h["LAMBDAU"],e[0],Z2[0].max())
    # xtitle = "theta [rad]"
    # ytitle = "phi [rad]"
    # from lipt import plot_surface
    # plot_surface(Z2[0]/Z2[0].max(),theta,phi,title=title,xtitle=xtitle,ytitle=ytitle)

    #
    # create uphot.dat file (like in SHADOW undul_phot)
    #
    # file_out = "uphot.dat"
    # f = open(file_out,'w')
    # f.write("%d  %d  %d \n"%(h["NG_E"],h["NG_T"],h["NG_P"]))
    # for ie in range(h["NG_E"]):
    #     f.write("%20.10f \n"%(e[ie]))
    #
    # for ie in range(h["NG_E"]):
    #     for t in range(h["NG_T"]):
    #         f.write("%20.10f \n"%(theta[t]))
    #
    # for ie in range(h["NG_E"]):
    #     for t in range(h["NG_T"]):
    #         for p in range(h["NG_P"]):
    #             f.write("%20.10f \n"%(phi[p]))
    #
    #
    # for ie in range(h["NG_E"]):
    #     for t in range(h["NG_T"]):
    #         for p in range(h["NG_P"]):
    #             f.write("%20.10f \n"%Z2[ie,t,p])
    #
    # for ie in range(h["NG_E"]):
    #     for t in range(h["NG_T"]):
    #         for p in range(h["NG_P"]):
    #             f.write("%20.10f \n"%(POL_DEG[ie,t,p]))
    #
    # f.close()
    # print("File written to disk: %s"%file_out)

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


    # corrdinates are polar: multiply by sin(theta) to allow dS= r^2 sin(Theta) dTheta dPhi
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
def test_undul_phot():
    tmp = \
        """
        {
        "LAMBDAU":     0.0320000015,
        "K":      0.250000000,
        "E_ENERGY":       6.03999996,
        "E_ENERGY_SPREAD":    0.00100000005,
        "NPERIODS": 50,
        "EMIN":       10500.0000,
        "EMAX":       10550.0000,
        "INTENSITY":      0.200000003,
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
    undul_phot(h)
    tmp = np.loadtxt("uphot.dat",skiprows=1)
    print("Obtained result[700]: %g (comparing to 6.09766e+16)"%tmp[7000])
    assert( np.abs(tmp[7000] - 6.09766e+16) < 1e13)


def test_undul_phot_srw():
    tmp = \
        """
        {
        "LAMBDAU":     0.0320000015,
        "K":      0.250000000,
        "E_ENERGY":       6.03999996,
        "E_ENERGY_SPREAD":    0.00100000005,
        "NPERIODS": 50,
        "EMIN":       10500.0000,
        "EMAX":       10550.0000,
        "INTENSITY":      0.200000003,
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
    undul_phot_srw(h)
    tmp = numpy.loadtxt("uphot.dat",skiprows=1)
    print("Obtained result[700]: %g (comparing to 6.09766e+16)"%tmp[7000])
    # assert( np.abs(tmp[7000] - 6.09766e+16) < 1e13)


def test_undul_cdf(do_plot=True):
    #
    # undul_phot.dat must exist
    #
    from SourceUndulator import shadow3_commands, load_uphot_dot_dat, load_xshundun_dot_sha, write_xshundun_dot_sha

    shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
                    input_file="shadow3_undul_cdf.inp")


    radiation = load_uphot_dot_dat(file_in="uphot.dat")

    cdf2 = undul_cdf(radiation,method='sum',do_plot=False)
    write_xshundun_dot_sha(cdf2,file_out="xshundul2.sha")


    cdf3 = undul_cdf(radiation,method='trapz',do_plot=False)
    write_xshundun_dot_sha(cdf3,file_out="xshundul3.sha")

    TWO1,ONE1,ZERO1,E1,T1,P1,POL_DEGREE1 = load_xshundun_dot_sha(file_in="xshundul.sha", do_plot=do_plot)
    TWO2,ONE2,ZERO2,E2,T2,P2,POL_DEGREE2 = load_xshundun_dot_sha(file_in="xshundul2.sha",do_plot=do_plot)
    TWO3,ONE3,ZERO3,E3,T3,P3,POL_DEGREE3 = load_xshundun_dot_sha(file_in="xshundul3.sha",do_plot=do_plot)


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

    # undul_phot("xshundul.json")
    # test_undul_phot()

    undul_phot_srw("xshundul.json")
    # test_undul_phot_srw()

    # test_undul_cdf(do_plot=False)
