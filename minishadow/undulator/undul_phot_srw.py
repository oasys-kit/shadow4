#
# script for running SRW to create a SHADOW source
#
import json
import numpy
import srwlib as sl
import array
import sys
from scipy import interpolate

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
  x = numpy.linspace(stk.mesh.xStart,stk.mesh.xFin,stk.mesh.nx)
  y = numpy.linspace(stk.mesh.yStart,stk.mesh.yFin,stk.mesh.ny)
  e = numpy.linspace(stk.mesh.eStart,stk.mesh.eFin,stk.mesh.ne)
  Z2 = numpy.zeros((e.size,x.size,y.size))
  for ie in range(e.size):
      for ix in range(x.size):
          for iy in range(y.size):
            Z2[ie,ix,iy] = data0[iy,ix,ie]
  return Z2,e,x,y

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


    # interpolate on polar grid
    radiation,e,x,y = Stokes0ToArrays(stk)
    for ie in range(e.size):
      tck = myinterpol_object(x,y,radiation[ie])
      for itheta in range(theta.size):
        for iphi in range(phi.size):
          R = slit_distance / numpy.cos(theta[itheta])
          r = R * numpy.sin(theta[itheta])
          X = r * numpy.cos(phi[iphi])
          Y = r * numpy.sin(phi[iphi])
          tmp = tck(X,Y)
          Z2[ie,itheta,iphi] = tmp


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
    file_out = "uphot.dat"
    f = open(file_out,'w')
    f.write("%d  %d  %d \n"%(h["NG_E"],h["NG_T"],h["NG_P"]))
    for ie in range(h["NG_E"]):
        f.write("%20.10f \n"%(e[ie]))

    for ie in range(h["NG_E"]):
        for t in range(h["NG_T"]):
            f.write("%20.10f \n"%(theta[t]))

    for ie in range(h["NG_E"]):
        for t in range(h["NG_T"]):
            for p in range(h["NG_P"]):
                f.write("%20.10f \n"%(phi[p]))


    for ie in range(h["NG_E"]):
        for t in range(h["NG_T"]):
            for p in range(h["NG_P"]):
                f.write("%20.10f \n"%Z2[ie,t,p])

    for ie in range(h["NG_E"]):
        for t in range(h["NG_T"]):
            for p in range(h["NG_P"]):
                f.write("%20.10f \n"%(1.0)) # POL_DEG[ie,t,p])

    f.close()
    print("File written to disk: %s"%file_out)


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






if __name__=="__main__":
  undul_phot_srw("xshundul.json")
  # test_undul_phot_srw()


