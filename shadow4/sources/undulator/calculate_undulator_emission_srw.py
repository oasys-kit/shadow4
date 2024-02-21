"""
Calculation of undulator far field and backpropagation to the center of the undulator
SRW implementation



Available public functions:

    calculate_undulator_emission_srw()


"""

try:
    import oasys_srw.srwlib as sl
except:
    try:
        import srwlib as sl
    except:
        raise ImportError("SRW not imported")

from copy import deepcopy
import numpy
import sys
from scipy import interpolate

def _srw_stokes0_to_arrays(stk):
  Shape = (4,stk.mesh.ny,stk.mesh.nx,stk.mesh.ne)
  data = numpy.ndarray(buffer=stk.arS, shape=Shape, dtype=stk.arS.typecode)
  data0 = data[0]
  data1 = data[1]
  x = numpy.linspace(stk.mesh.xStart, stk.mesh.xFin, stk.mesh.nx)
  y = numpy.linspace(stk.mesh.yStart, stk.mesh.yFin, stk.mesh.ny)
  e = numpy.linspace(stk.mesh.eStart, stk.mesh.eFin, stk.mesh.ne)
  Z2 = numpy.zeros((e.size, x.size, y.size))
  POL_DEG = numpy.zeros((e.size, x.size, y.size))
  for ie in range(e.size):
      for ix in range(x.size):
          for iy in range(y.size):
            Z2[ie, ix, iy] = data0[iy, ix, ie]
            # this is shadow definition, that uses POL_DEG = |Ex|/(|Ex|+|Ey|)
            Ex = numpy.sqrt(numpy.abs(0.5 * (data0[iy, ix, ie] + data1[iy, ix, ie])))
            Ey = numpy.sqrt(numpy.abs(0.5 * (data0[iy, ix, ie] - data1[iy, ix, ie])))
            POL_DEG[ie,ix,iy] = Ex / (Ex + Ey)
  return Z2, POL_DEG, e, x, y

def _srw_stokes0_to_spec(stk, fname="srw_xshundul.spec"):
  #
  # writes emission in a SPEC file (cartesian grid)
  #
  Shape = (4, stk.mesh.ny, stk.mesh.nx, stk.mesh.ne)
  data = numpy.ndarray(buffer=stk.arS, shape=Shape, dtype=stk.arS.typecode)
  data0 = data[0]
  x = numpy.linspace(stk.mesh.xStart, stk.mesh.xFin, stk.mesh.nx)
  y = numpy.linspace(stk.mesh.yStart, stk.mesh.yFin, stk.mesh.ny)
  e = numpy.linspace(stk.mesh.eStart, stk.mesh.eFin, stk.mesh.ne)
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

def _SRWEFieldAsNumpy(srwwf):
    """
    Extracts electrical field from a SRWWavefront
    :param srw_wavefront: SRWWavefront to extract electrical field from.
    :return: 4D numpy array: [energy, horizontal, vertical, polarisation={0:horizontal, 1: vertical}]
    """
    dim_x = srwwf.mesh.nx
    dim_y = srwwf.mesh.ny
    number_energies = srwwf.mesh.ne

    x_polarization = _SRWArrayToNumpyComplexArray(srwwf.arEx, dim_x, dim_y, number_energies)
    y_polarization = _SRWArrayToNumpyComplexArray(srwwf.arEy, dim_x, dim_y, number_energies)

    return [x_polarization, y_polarization]

def _SRWArrayToNumpyComplexArray(srw_array, dim_x, dim_y, number_energies, polarized=True):
    """
    Converts a SRW array to a numpy.array.
    :param srw_array: SRW array
    :param dim_x: size of horizontal dimension
    :param dim_y: size of vertical dimension
    :param number_energies: Size of energy dimension
    :return: 4D numpy array: [energy, horizontal, vertical, polarisation={0:horizontal, 1: vertical}]
    """
    re = numpy.array(srw_array[::2], dtype=float)
    im = numpy.array(srw_array[1::2], dtype=float)
    return _reshape(re + 1j * im, dim_x, dim_y, number_energies, polarized)

def _reshape(numpy_array, dim_x, dim_y, number_energies, polarized=True):
    if polarized: numpy_array = numpy_array.reshape((dim_y, dim_x, number_energies, 1))
    else:         numpy_array.reshape((dim_y, dim_x, number_energies))
    numpy_array = numpy_array.swapaxes(0, 2)
    return numpy_array.copy()

def _srw_interpol_object(x, y, z):
    #2d interpolation
    if numpy.iscomplex(z[0, 0]):
        tck_real = interpolate.RectBivariateSpline(x, y, numpy.real(z))
        tck_imag = interpolate.RectBivariateSpline(x, y, numpy.imag(z))
        return tck_real,tck_imag
    else:
        tck = interpolate.RectBivariateSpline(x, y, z)
        return tck

def _srw_interpol(x, y, z, x1, y1):
    #2d interpolation
    if numpy.iscomplex(z[0, 0]):
        tck_real, tck_imag = _srw_interpol_object(x, y, z)
        z1_real = tck_real(numpy.real(x1), numpy.real(y1))
        z1_imag = tck_imag(numpy.imag(x1), numpy.imag(y1))
        return (z1_real + 1j * z1_imag)
    else:
        tck = _srw_interpol_object(x, y, z)
        z1 = tck(x1, y1)
        return z1

def _undul_phot_SRW(E_ENERGY, INTENSITY, LAMBDAU, NPERIODS, K, EMIN, EMAX, NG_E, MAXANGLE, NG_T, NG_P,
                    number_of_trajectory_points=100,
                    flag_size=0,
                    distance=100.0,
                    srw_range=0.05,       # **** [5]: Horizontal Range modification factor at Resizing (1. means no modification)
                    srw_resolution=50.0,  # **** [6]: Horizontal Resolution modification factor at Resizing
                    srw_semianalytical=0, # **** [3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
                    ):

    lambdau = LAMBDAU
    k = K
    e_energy = E_ENERGY
    nperiods = NPERIODS
    emin = EMIN
    emax = EMAX
    intensity = INTENSITY
    maxangle = MAXANGLE

    nx = NG_T
    nz = nx

    ne = NG_E
    u_length = nperiods * lambdau

    print("SRW calculation: lambdau = ",lambdau)
    print("SRW calculation: k = ",k)
    print("SRW calculation: e_energy = ",e_energy)
    print("SRW calculation: nperiods = ",nperiods)
    print("SRW calculation: intensity = ",intensity)
    print("SRW calculation: maxangle=%g rad, (%d x %d points) "%(maxangle, nx, nz))
    print("SRW calculation: emin =%g, emax=%g, ne=%d "%(emin,emax,ne))
    print("SRW calculation: UNDULATOR LENGTH: = ", u_length, -0.5 * lambdau * (nperiods + 4))

    #
    # define additional parameters needed by SRW
    #
    B = k / 93.4 / lambdau

    #
    # prepare inputs
    #
    slit_xmin = -maxangle * distance
    slit_xmax =  maxangle * distance
    slit_zmin = -maxangle * distance
    slit_zmax =  maxangle * distance

    print("SRW calculation: H slit: %f, %f, %f" % (slit_xmin, slit_xmax, nx))
    print("SRW calculation: V slit: %f, %f, %f" % (slit_zmin, slit_zmax, nz))
    print("SRW calculation: nperiods: %d, lambdau: %f, B: %f)" % (nperiods, lambdau, B))
    print("SRW calculation: e=%f, Iavg=%f " % (e_energy, intensity) )

    #
    # calculations
    #
    ####################################################
    # LIGHT SOURCE

    z_start = -0.5 * lambdau * (nperiods + 4)

    part_beam = sl.SRWLPartBeam()
    part_beam.Iavg = intensity
    part_beam.partStatMom1.x = 0.0
    part_beam.partStatMom1.y = 0.0
    part_beam.partStatMom1.z = z_start #  -2.45 ##########################
    part_beam.partStatMom1.xp = 0.0
    part_beam.partStatMom1.yp = 0.0
    part_beam.partStatMom1.gamma = e_energy / 0.51099890221e-03 # 11741.70710144324
    part_beam.partStatMom1.relE0 = 1.0 #
    part_beam.partStatMom1.nq = -1 #
    part_beam.arStatMom2[0] = 0.0
    part_beam.arStatMom2[1] = 0.0
    part_beam.arStatMom2[2] = 0.0
    part_beam.arStatMom2[3] = 0.0
    part_beam.arStatMom2[4] = 0.0
    part_beam.arStatMom2[5] = 0.0
    part_beam.arStatMom2[10] = 0.0


    magnetic_fields = []
    magnetic_fields.append(sl.SRWLMagFldH(1, 'v',
                                       _B=B,
                                       _ph=0.0,
                                       _s=-1,
                                       _a=1.0))
    magnetic_structure = sl.SRWLMagFldU(_arHarm=magnetic_fields, _per=lambdau, _nPer=nperiods)
    magnetic_field_container = sl.SRWLMagFldC(_arMagFld=[magnetic_structure],
                                           _arXc=sl.array('d', [0.0]),
                                           _arYc=sl.array('d', [0.0]),
                                           _arZc=sl.array('d', [0.0]))


    sys.stdout.flush()

    mesh = sl.SRWLRadMesh(_eStart=emin,
                       _eFin=emax,
                       _ne=ne,
                       _xStart=slit_xmin,
                       _xFin=slit_xmax,
                       _nx=nx,
                       _yStart=slit_zmin,
                       _yFin=slit_zmax,
                       _ny=nz,
                       _zStart=distance)

    print ("Calculating source (SRW)...", mesh.ne, mesh.eStart, mesh.eFin, mesh.xStart, mesh.xFin, mesh.yStart, mesh.yFin)

    wfr = sl.SRWLWfr()
    wfr.allocate(mesh.ne, mesh.nx, mesh.ny)
    wfr.mesh = mesh
    wfr.partBeam = part_beam # part_beam
    wfr.unitElFld = 1 #Electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)

    sl.srwl.CalcElecFieldSR(wfr, 0, magnetic_field_container, [1, 0.01, 0, 0, number_of_trajectory_points * nperiods, 1, 0])
    ###
    # part_beam = _srw_drift_electron_beam(part_beam, -part_beam.moved)
    ###

    sys.stdout.write('  done\n')
    sys.stdout.flush()

    #
    # use wavefront
    #

    # stk = sl.SRWLStokes()
    # stk.mesh = mesh
    # stk.allocate(mesh.ne, mesh.nx, mesh.ny)
    # wfr.calc_stokes(stk)
    #
    # radiation, pol_deg, e, x, y = _srw_stokes0_to_arrays(stk)
    wSigma, wPi = _SRWEFieldAsNumpy(wfr)
    wModSigma = numpy.abs(wSigma[:, :, :, 0])
    wModPi = numpy.abs(wPi[:, :, :, 0])
    wAngleSigma = numpy.angle(wSigma[:, :, :, 0])
    wAnglePi = numpy.angle(wPi[:, :, :, 0])

    radiation = wModSigma**2 + wModPi**2
    # !C SHADOW defines the degree of polarization by |E| instead of |E|^2
    # !C i.e.  P = |Ex|/(|Ex|+|Ey|)   instead of   |Ex|^2/(|Ex|^2+|Ey|^2)
    pol_deg = wModSigma / (wModSigma + wModPi)
    x = numpy.linspace(mesh.xStart, mesh.xFin, mesh.nx)
    y = numpy.linspace(mesh.yStart, mesh.yFin, mesh.ny)
    e = numpy.linspace(mesh.eStart, mesh.eFin, mesh.ne)

    #
    # interpolate for polar grid
    #

    # polar grid
    theta          = numpy.linspace(0, MAXANGLE, NG_T)
    phi            = numpy.linspace(0, numpy.pi / 2, NG_P)
    Z2             = numpy.zeros((NG_E, NG_T, NG_P))
    POL_DEG        = numpy.zeros((NG_E, NG_T, NG_P))
    efield_s_mod   = numpy.zeros_like(Z2)
    efield_s_angle = numpy.zeros_like(Z2)
    efield_p_mod   = numpy.zeros_like(Z2)
    efield_p_angle = numpy.zeros_like(Z2)

    for ie in range(e.size):
      tck = _srw_interpol_object(x, y, radiation[ie])
      tck_pol_deg = _srw_interpol_object(x, y, pol_deg[ie])
      tck_efield_s_mod = _srw_interpol_object(x, y, wModSigma[ie])
      tck_efield_s_angle = _srw_interpol_object(x, y, wAngleSigma[ie])
      tck_efield_p_mod = _srw_interpol_object(x, y, wModPi[ie])
      tck_efield_p_angle = _srw_interpol_object(x, y, wAnglePi[ie])
      for itheta in range(theta.size):
        for iphi in range(phi.size):
          R = distance / numpy.cos(theta[itheta])
          r = R * numpy.sin(theta[itheta])
          X = r * numpy.cos(phi[iphi])
          Y = r * numpy.sin(phi[iphi])
          tmp = tck(X, Y)

          #  Conversion from SRW units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
          tmp *= (distance * 1e3)**2 # photons/mm^2 -> photons/rad^2
          tmp /= 1e-3 * e[ie] # photons/o.1%bw -> photons/eV

          Z2[ie, itheta, iphi] = tmp
          POL_DEG[ie, itheta, iphi] = tck_pol_deg(X, Y)
          efield_s_mod[ie, itheta, iphi] = tck_efield_s_mod(X, Y)
          efield_s_angle[ie, itheta, iphi] = tck_efield_s_angle(X, Y)
          efield_p_mod[ie, itheta, iphi] = tck_efield_p_mod(X, Y)
          efield_p_angle[ie, itheta, iphi] = tck_efield_p_angle(X, Y)

    out = { 'radiation'             : Z2,
            'polarization'          : POL_DEG,
            'photon_energy'         : e,
            'theta'                 : theta,
            'phi'                   : phi,
            'trajectory'            : None,
            'e_amplitude_sigma'     : efield_s_mod * numpy.exp(1j * efield_s_angle),
            'e_amplitude_pi'        : efield_p_mod * numpy.exp(1j * efield_p_angle),
            'CART_radiation'        : radiation,
            'CART_polarizartion'    : pol_deg,
            'CART_x'                : x / distance, # angle in rad
            'CART_y'                : y / distance, # angle in rad
            'CART_e_amplitude_sigma': wSigma[:, :, :, 0],
            'CART_e_amplitude_pi'   : wPi[:, :, :, 0],
            }

    #
    # backpropagation
    #
    if flag_size == 2:
        ####################################################
        # BEAMLINE

        srw_oe_array = []
        srw_pp_array = []

        drift_before_oe_0 = sl.SRWLOptD(-distance)
        # Wavefront Propagation Parameters:
        # [0]: Auto-Resize (1) or not (0) Before propagation
        # [1]: Auto-Resize (1) or not (0) After propagation
        # [2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
        # **** [3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
        # [4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
        # **** [5]: Horizontal Range modification factor at Resizing (1. means no modification)
        # **** [6]: Horizontal Resolution modification factor at Resizing
        # **** [7]: Vertical Range modification factor at Resizing
        # **** [8]: Vertical Resolution modification factor at Resizing
        # [9]: Type of wavefront Shift before Resizing (not yet implemented)
        # [10]: New Horizontal wavefront Center position after Shift (not yet implemented)
        # [11]: New Vertical wavefront Center position after Shift (not yet implemented)
        # pp_drift_before_oe_0 = [0, 0, 1.0, 0, 0, 0.05, 50.0, 0.05, 50.0, 0, 0.0, 0.0]
        # srw_h_range = 0.05, srw_h_resolution = 50.0, srw_v_range = 0.05, srw_v_resolution = 50.0,
        pp_drift_before_oe_0 = [0, 0, 1.0, srw_semianalytical, 0, srw_range, srw_resolution, srw_range, srw_resolution, 0, 0.0, 0.0]
        print('pp_drift_before_oe_0: ', pp_drift_before_oe_0)

        srw_oe_array.append(drift_before_oe_0)
        srw_pp_array.append(pp_drift_before_oe_0)

        ####################################################
        # PROPAGATION

        optBL = sl.SRWLOptC(srw_oe_array, srw_pp_array)

        print("Backpropagating source (SRW)...")
        sl.srwl.PropagElecField(wfr, optBL)

        BPwSigma, BPwPi = _SRWEFieldAsNumpy(wfr)

        mesh1 = deepcopy(wfr.mesh)
        print("Backpropagated mesh...", mesh1.ne, mesh1.eStart, mesh1.eFin, mesh1.ne, mesh1.xStart, mesh1.xFin, mesh1.nx, mesh1.yStart, mesh1.yFin, mesh1.ny)
        BPx = numpy.linspace(mesh1.xStart, mesh1.xFin, mesh1.nx)
        BPy = numpy.linspace(mesh1.yStart, mesh1.yFin, mesh1.ny)

        out['CART_BACKPROPAGATED_e_amplitude_sigma']= BPwSigma[:, :, :, 0]
        out['CART_BACKPROPAGATED_e_amplitude_pi']= BPwPi[:, :, :, 0]
        out['CART_BACKPROPAGATED_radiation']= numpy.abs(BPwSigma[:, :, :, 0])**2 + numpy.abs(BPwPi[:, :, :, 0])**2
        out['CART_BACKPROPAGATED_x']= BPx # length in m
        out['CART_BACKPROPAGATED_y']= BPy # length in m

    return out

def calculate_undulator_emission_srw(
                     electron_energy              = 6.0,
                     electron_current             = 0.2,
                     undulator_period             = 0.018,
                     undulator_nperiods           = 100,
                     K                            = 1.0,
                     photon_energy                = 2000.0,
                     EMAX                         = 20000.0,
                     NG_E                         = 10,
                     MAXANGLE                     = 0.1,
                     number_of_points             = 100,
                     NG_P                         = 100,
                     number_of_trajectory_points  = 100,
                     flag_size                    = 2,
                     distance                     = 100.0,
                     srw_range                    = 0.05,
                     srw_resolution               = 50.0,
                     srw_semianalytical           = 0,
                     ):
    """
    Calculate undulator emission (far field) and backpropagation to the center of the undulator using srw.

    Parameters
    ----------
    electron_energy: float, optional
        The electron energy in GeV.
    electron_current: float, optional
        The electron beam current in A.
    undulator_period: float, optional
        The undulator period in m.
    undulator_nperiods: float, optional
        The undulator number of periods.
    K: float, optional
        The undulator K factor (vertical).
    photon_energy: float, optional
        The photon energy (or starting energy for arrays) in eV.
    EMAX: float, optional
        The end photon energy (for arrays).
    NG_E: int, optional
        The number of points in photon energy (>1 for array).
    MAXANGLE: float, optional
        The maximum half-angle for delimiting the calculation in rad.
    number_of_points: int, optional
        The number of points in theta (elevation angle) or x and y.
    NG_P: int, optional
        The number of points in psi (azimuthal angle).
    number_of_trajectory_points: int, optional
        The number of trajectory points (per period).
    flag_size: int, optional
        A flag to control the model used for sampling the radiation size at the center of the undulator:
        0=point, 1=Gaussian, 2=backpropagated the far field.
    distance: float, optional
        The distance to place the image plane where far field is calculated in m.
    srw_range: float, optional
        The SRW range factor (affecting both the H and V directions).
    srw_resolution: float, optional
        The SRW resolution factor (affecting both the H and V directions).
    srw_semianalytical: int, optional
        A flag to indicate if SRW uses the semianalytical treatment of the phases (0=No, 1=Yes).

    Returns
    -------
    dict
        A dictionary with calculated values and arrays.

    """
    return _undul_phot_SRW(
                        electron_energy,
                        electron_current,
                        undulator_period,
                        undulator_nperiods,
                        K,
                        photon_energy,
                        EMAX,
                        NG_E,
                        MAXANGLE,
                        number_of_points,
                        NG_P,
                        number_of_trajectory_points=number_of_trajectory_points,
                        flag_size=flag_size,
                        distance=distance,
                        srw_range=srw_range,
                        srw_resolution=srw_resolution,
                        srw_semianalytical=srw_semianalytical,
                        )
#
#
#
if __name__ == "__main__":
    dict1 = calculate_undulator_emission_srw(
                     electron_energy              = 6.0,
                     electron_current             = 0.2,
                     undulator_period             = 0.025,
                     undulator_nperiods           = 188.0,
                     K                            = 1.681183,
                     photon_energy                = 5591.0,
                     EMAX                         = 5700.0,
                     NG_E                         = 11,
                     MAXANGLE                     = 3 * 2e-5,
                     number_of_points             = 251,
                     NG_P                         = 11,
                     number_of_trajectory_points  = 51,
                     flag_size                    = 2,
                     distance                     = 100.0,
                     srw_range                    = 0.05,
                     srw_resolution               = 50.0,
                     srw_semianalytical           = 1,
    )

    from srxraylib.plot.gol import plot_image, plot
    if False:
        plot_image(dict1['radiation'][0], dict1['theta'], dict1['phi'], aspect='auto',  title="first", show=0)
        plot_image(dict1['radiation'][-1], dict1['theta'], dict1['phi'], aspect='auto', title="last", show=1)

        plot_image(dict1['CART_radiation'][0], 1e6 * dict1['CART_x'], 1e6 * dict1['CART_y'], aspect='auto', title="first", show=0)
        plot_image(dict1['CART_radiation'][-1], 1e6 * dict1['CART_x'], 1e6 * dict1['CART_y'], aspect='auto', title="last", show=1)

    if True:
        plot_image(numpy.abs(dict1['CART_e_amplitude_sigma'][0])**2  + numpy.abs(dict1['CART_e_amplitude_pi'][0])**2 ,
                   100 * 1e6 * dict1['CART_x'], 100 * 1e6 * dict1['CART_y'], aspect='auto', title="first", show=0)
        plot_image(numpy.abs(dict1['CART_e_amplitude_sigma'][-1])**2 + numpy.abs(dict1['CART_e_amplitude_pi'][-1])**2,
                   100 * 1e6 * dict1['CART_x'], 100 * 1e6 * dict1['CART_y'], aspect='auto', title="last", show=0)

        i_prop0 =  dict1['CART_BACKPROPAGATED_radiation'][0]  # numpy.abs(dict1['CART_BACKPROPAGATED_e_amplitude_sigma'][0])**2  + numpy.abs(dict1['CART_BACKPROPAGATED_e_amplitude_pi'][0])**2
        i_prop1 =  dict1['CART_BACKPROPAGATED_radiation'][-1] # numpy.abs(dict1['CART_BACKPROPAGATED_e_amplitude_sigma'][-1])**2 + numpy.abs(dict1['CART_BACKPROPAGATED_e_amplitude_pi'][-1])**2
        x_um = 1e6 * dict1['CART_BACKPROPAGATED_x']
        y_um = 1e6 * dict1['CART_BACKPROPAGATED_y']
        plot_image( i_prop0, x_um, y_um, aspect='auto', title="PROPAGATED first", show=0)
        plot_image( i_prop1, x_um, y_um, aspect='auto', title="PROPAGATED last", show=1)
        plot(x_um, i_prop0[x_um.size//2,:],
             y_um, i_prop0[:, y_um.size//2])
