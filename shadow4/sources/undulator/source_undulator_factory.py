#
# SHADOW3 Undulator preprocessors implemented in python
#
# this code replaces SHADOW3's undul_phot and undul_cdf
#
# It calculates the undulator radiation as a function of energy, theta and phi. Phi is the polar angle.
#
# It uses internal code (no dependency) hacked from pySRU
#  (see SourceUndulatorFactorySrw.py and SourceUndulatorFactoryPysru.py for SRW and native pySRU backends, respectively)
#
#
# Available functions:
#
#     calculate_undulator_emission: interface for undul_phot (like undul_phot of SHADOW3 but written in python with
#                                   internal code hacked from pySRU).
#     undul_cdf                   : like undul_cdf in SHADOW3.
#
#
import numpy
import numpy as np
from scipy import interpolate

import scipy.constants as codata
import scipy.integrate

#
# calculate a theorical trajectory in an undulator
# adapted from pySRU: analytical_trajectory_plane_undulator()
#
def _pysru_analytical_trajectory_plane_undulator(K=1.87, gamma=2544.03131115, lambda_u=0.020,
                                                 Nb_period=10, Nb_point=10, Beta_et=0.99993):
    N = int(Nb_period * Nb_point + 1)
    ku = 2.0 * np.pi / lambda_u
    omega_u = Beta_et * codata.c * ku

    # t
    t = np.linspace(-(lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), (lambda_u / (codata.c * Beta_et)) * (Nb_period / 2), N)
    ## x and z
    z = Beta_et * t + ((K / gamma)**2) * (1.0 / (8.0 * omega_u)) * np.sin( 2.0 * omega_u * t)
    x = (-(K / (gamma * omega_u)) * np.cos(omega_u * t))
    # # Vx and Vz
    v_z = Beta_et + ((K / gamma)**2) * (1.0 / 4.0) * np.cos(2.0 * omega_u * t)
    v_x= (K / (gamma )) * np.sin(omega_u * t)
    # # Ax and Az
    a_z=-omega_u * (K / gamma)**2 * 0.5 * np.sin(2.0 * omega_u * t)
    a_x= (K / (gamma )) * (omega_u ) * np.cos(omega_u * t)
    # y
    y=0.0*t
    v_y=y
    a_y=y
    # note that x,y,z are in units of c: X[m] = c * x
    # note that v_x,v_y,v_z are beta (in units of c): V_v[m  s**-1] = c * v_v
    return np.vstack((t, x, y, z, v_x, v_y, v_z, a_x, a_y, a_z))

# adapted from pySRU:  energy_radiated_approximation_and_farfield()
def _pysru_energy_radiated_approximation_and_farfield(omega=2.53465927101e17, electron_current=1.0,
                                                      trajectory=None, x=0.00 , y=0.0, D=None):

    if trajectory is None: trajectory=np.zeros((11,10))

    c6 = codata.e * electron_current * 1e-9 / (8.0 * np.pi ** 2 * codata.epsilon_0 * codata.c * codata.h)

    if D is not None:
        c6 /= D**2

    N = trajectory.shape[1]
    if D == None: # in radian
        n_chap = np.array([x, y, 1.0 - 0.5 * (x**2 + y**2)])
        X = np.sqrt(x ** 2 + y ** 2 )#TODO a changer
    else :  #in meters
        X = np.sqrt(x**2 + y**2 + D**2)
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

    E = np.zeros((3,), dtype=complex)
    integrand = np.zeros((3, N), dtype=complex)
    A1 = ( n_chap[1] * trajectory_v_z - n_chap[2] * trajectory_v_y)
    A2 = (-n_chap[0] * trajectory_v_z + n_chap[2] * trajectory_v_x)
    A3 = ( n_chap[0] * trajectory_v_y - n_chap[1] * trajectory_v_x)
    Alpha2 = np.exp(
        0. + 1j * omega * (trajectory_t + X / codata.c - n_chap[0] * trajectory_x
                                           - n_chap[1] * trajectory_y - n_chap[2] * trajectory_z))


    integrand[0] -= ( n_chap[1]*A3 - n_chap[2]*A2) * Alpha2
    integrand[1] -= (-n_chap[0]*A3 + n_chap[2]*A1) * Alpha2
    integrand[2] -= ( n_chap[0]*A2 - n_chap[1]*A1) * Alpha2

    for k in range(3):
        E[k] = np.trapz(integrand[k], trajectory_t)
    E *= omega * 1j

    terme_bord = np.full((3), 0. + 1j * 0., dtype=complex)
    Alpha_1 = (1.0 / (1.0 - n_chap[0] * trajectory_v_x[-1]
                      - n_chap[1] * trajectory_v_y[-1] - n_chap[2] * trajectory_v_z[-1]))
    Alpha_0 = (1.0 / (1.0 - n_chap[0] * trajectory_v_x[0]
                      - n_chap[1] * trajectory_v_y[0] - n_chap[2] * trajectory_v_z[0]))

    terme_bord += ((n_chap[1] * A3[-1] - n_chap[2] * A2[-1]) * Alpha_1 * Alpha2[-1])
    terme_bord -= ((n_chap[1] * A3[0]  - n_chap[2] * A2[0])  * Alpha_0 * Alpha2[0])
    E += terme_bord
    E *= c6**0.5
    return E


def _get_radiation_interpolated_cartesian(radiation, photon_energy, thetabm, phi,
                                          npointsx=100, npointsz=100, thetamax=None):
    """
    Interpolates the radiation array (in polar coordinates) to cartesian coordinates.

    Parameters
    ----------
    npointsx : int, optional
        The number of points in X.
    npointsz : int, optional
        The number of points in Z.
    thetamax : None or float, optional
        Maximum value of theta. By default (None) it uses the maximum theta.

    Returns
    -------
    tuple
        (radiation, array_x, array_y) in units of W/rad^2.
    """

    # radiation, photon_energy, thetabm, phi = self.get_result_radiation_polar()

    if thetamax is None:
        thetamax = thetabm.max()

    vx = numpy.linspace(-1.1 * thetamax, 1.1 * thetamax, npointsx)
    vz = numpy.linspace(-1.1 * thetamax, 1.1 * thetamax, npointsz)
    VX = numpy.outer(vx, numpy.ones_like(vz))
    VZ = numpy.outer(numpy.ones_like(vx), vz)
    VY = numpy.sqrt(1 - VX**2 - VZ**2)

    IN_thetabm = numpy.outer(thetabm, numpy.ones_like(phi))
    IN_phi = numpy.outer(numpy.ones_like(thetabm), phi)
    IN = numpy.zeros( (IN_thetabm.size, 2))
    IN[:,0] = IN_thetabm.flatten()
    IN[:,1] = IN_phi.flatten()

    # THETA = numpy.abs(numpy.arctan(numpy.sqrt(VX**2 + VZ**2) / VY))
    THETA = numpy.arctan(numpy.sqrt(VX**2 + VZ**2) / VY)
    # PHI = numpy.arctan2(numpy.abs(VZ), numpy.abs(VX))
    PHI = numpy.arctan2(VZ, VX)

    print(">> interpolating: THETA min, max: ", THETA.min(), THETA.max())
    print(">> interpolating: PHI min, max: ", PHI.min(), PHI.max())

    radiation_interpolated = numpy.zeros((radiation.shape[0], npointsx, npointsz))

    for i in range(radiation.shape[0]):
        interpolator_value = interpolate.RectBivariateSpline(thetabm, phi, radiation[i])
        radiation_interpolated[i] = interpolator_value.ev(THETA, PHI)

        # print(">>>>>>>>>>>>>>", IN_thetabm.shape, IN_phi.shape, radiation[i].shape)
        # interpolator_value = interpolate.NearestNDInterpolator(IN, radiation[i].flatten(), rescale=1)
        # print(">>>>>>>>>>>>>>>>", interpolator_value(THETA, PHI).shape)
        # radiation_interpolated[i] = interpolator_value(THETA, PHI)

    return radiation_interpolated, photon_energy, vx, vz

def _propagate_complex_amplitude_from_arrays(
                                            amplitude_flatten,
                                            X_flatten,
                                            Y_flatten,
                                            det_X_flatten=None,
                                            det_Y_flatten=None,
                                            wavelength=1e-10,
                                            propagation_distance=100):
    #
    # Fresnel-Kirchhoff integral (neglecting inclination factor)
    #
    if det_X_flatten is None: det_X_flatten = X_flatten
    if det_Y_flatten is None: det_Y_flatten = Y_flatten

    ngood = det_X_flatten.size

    fla_complex_amplitude_propagated = numpy.zeros(ngood, dtype=complex)

    Propagation_distance = numpy.ones_like(X_flatten) * propagation_distance

    wavenumber = 2 * numpy.pi / wavelength

    for i in range(ngood):
        r = numpy.sqrt( (X_flatten - det_X_flatten[i])**2 +
                        (Y_flatten - det_Y_flatten[i])**2 +
                        Propagation_distance**2 )

        fla_complex_amplitude_propagated[i] = (amplitude_flatten / r * numpy.exp(1.j * wavenumber *  r)).sum()


    # added srio@esrf.eu 2018-03-23 to conserve energy - TODO: review method!
    i0 = numpy.abs(amplitude_flatten)**2
    i1 = numpy.abs(fla_complex_amplitude_propagated)**2

    fla_complex_amplitude_propagated *= i0.sum() / i1.sum()

    return fla_complex_amplitude_propagated

# now, undul_phot
def _undul_phot(E_ENERGY,INTENSITY, LAMBDAU, NPERIODS, K, EMIN, EMAX, NG_E, MAXANGLE, NG_T, NG_P,
               number_of_trajectory_points=20, flag_size=0, distance=100.0, magnification=0.01, ):
    #
    # calculate trajectory
    #
    angstroms_to_eV = codata.h * codata.c / codata.e * 1e10
    gamma = E_ENERGY * 1e9 / 0.511e6
    Beta = np.sqrt(1.0 - (1.0 / gamma**2))
    Beta_et = Beta * (1.0 - (K / (2.0 * gamma))**2)

    E = np.linspace(EMIN, EMAX, NG_E, dtype=float)
    wavelength_array_in_A = angstroms_to_eV / E
    omega_array = 2 * np.pi * codata.c / (wavelength_array_in_A * 1e-10)


    #
    # calculate source in polar coordinates
    #
    Nb_pts_trajectory = int(number_of_trajectory_points * NPERIODS)
    print("Calculating trajectory...")
    T = _pysru_analytical_trajectory_plane_undulator(K=K, gamma=gamma, lambda_u=LAMBDAU, Nb_period=NPERIODS,
                                        Nb_point=Nb_pts_trajectory,
                                        Beta_et=Beta_et)
    print("Done calculating trajectory...")

    #
    # polar grid
    #
    D = distance # placed far away (100 m)

    #
    # for divergences sample angle in [0,MAXANGLE]
    #
    # theta = np.linspace(-MAXANGLE, MAXANGLE, NG_T, dtype=float)
    # phi = (numpy.arange(NG_P) / NG_P - 0.5) * numpy.pi # np.linspace(-numpy.pi/2, numpy.pi/2, NG_P, dtype=float)
    theta = np.linspace(0, MAXANGLE, NG_T, dtype=float)
    phi = numpy.linspace(0, numpy.pi / 2, NG_P)
    # phi = (numpy.arange(NG_P) / NG_P - 0.5) * 2 * numpy.pi # np.linspace(-numpy.pi/2, numpy.pi/2, NG_P, dtype=float)

    THETA = numpy.outer(theta, numpy.ones_like(phi))
    PHI = numpy.outer(numpy.ones_like(theta), phi)
    r = distance / numpy.cos(THETA) * numpy.sin(THETA)
    xx = r * numpy.cos(PHI)
    yy = r * numpy.sin(PHI)

    POL_DEG = np.zeros((omega_array.size, theta.size, phi.size))
    EFIELD_X = np.zeros_like(POL_DEG, dtype=complex)
    EFIELD_Y = np.zeros_like(POL_DEG, dtype=complex)

    print("Calculating radiation...", POL_DEG.shape, theta.shape, phi.shape)
    # for o in range(omega_array.size):
    #     print("Calculating energy %8.3f eV (%d of %d)"%(E[o],o+1,omega_array.size))
    #     for t in range(theta.size):
    #         for p in range(phi.size):
    #             R = D / np.cos(theta[t])
    #             r = R * np.sin(theta[t])
    #             X = r * np.cos(phi[p])
    #             Y = r * np.sin(phi[p])
    #             ElecField = _pysru_energy_radiated_approximation_and_farfield(omega=omega_array[o],
    #                                                                           electron_current=INTENSITY,
    #                                                                           trajectory=T ,
    #                                                                           x=X ,
    #                                                                           y=Y,
    #                                                                           D=D)
    #
    #             # if numpy.abs(phi[p]) > numpy.pi / 3: ElecField *= 0
    #
    #             pol_deg = np.abs(ElecField[0]) / (np.abs(ElecField[0]) + np.abs(ElecField[1])) # SHADOW definition
    #             intensity =  (np.abs(ElecField[0]) ** 2 + np.abs(ElecField[1])** 2 + np.abs(ElecField[2])** 2)
    #
    #             #  Conversion from pySRU units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
    #             coeff = (D * 1e3) ** 2  # photons/mm^2 -> photons/rad^2
    #             coeff /= 1e-3 * E[o]  # photons/o.1%bw -> photons/eV
    #             intensity *= coeff
    #
    #             Z2[o, t, p] = intensity
    #             POL_DEG[o, t, p] = pol_deg
    #             EFIELD_X[o, t, p] = ElecField[0] * numpy.sqrt(coeff)
    #             EFIELD_Y[o, t, p] = ElecField[1] * numpy.sqrt(coeff)
    #             EFIELD_Z[o, t, p] = ElecField[2] * numpy.sqrt(coeff)


    from pySRU.Simulation import create_simulation
    from pySRU.ElectronBeam import ElectronBeam as PysruElectronBeam
    from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as PysruUndulator
    from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC
    from pySRU.RadiationFactory import RADIATION_METHOD_APPROX_FARFIELD

    for ie, e in enumerate(E):
        print("Calculating undulator source at energy %g eV (%d of %d) (%d x %d points) %d traj."%(e,
                                            ie+1, E.size, NG_T, NG_P, Nb_pts_trajectory))
        simulation_test = create_simulation(magnetic_structure=PysruUndulator(
                                                                        K=K,
                                                                        period_length=LAMBDAU,
                                                                        length=LAMBDAU*NPERIODS),
                                            electron_beam=PysruElectronBeam(
                                                                        Electron_energy=E_ENERGY,
                                                                        I_current=INTENSITY),
                                            magnetic_field=None,
                                            photon_energy=e,
                                            traj_method=TRAJECTORY_METHOD_ANALYTIC,
                                            Nb_pts_trajectory=Nb_pts_trajectory, #None,
                                            rad_method=RADIATION_METHOD_APPROX_FARFIELD,
                                            initial_condition=None,
                                            distance=D,
                                            X=xx.flatten(),
                                            Y=yy.flatten(),
                                            XY_are_list=True)

        Efield = simulation_test.electric_field._electrical_field
        tmp_x = Efield[:, 0].copy()
        tmp_y = Efield[:, 1].copy()
        tmp_x.shape = (theta.size, phi.size)
        tmp_y.shape = (theta.size, phi.size)

        EFIELD_X[ie, :, :] = tmp_x
        EFIELD_Y[ie, :, :] = tmp_y

        POL_DEG[ie, :, :] = numpy.abs(tmp_x) / (numpy.abs(tmp_x) + numpy.abs(tmp_y)) # SHADOW definition
    print("Done calculating radiation...")

    radiation = numpy.abs(EFIELD_X)**2 + numpy.abs(EFIELD_Y)**2
    out = {'radiation'         : radiation,
            'polarization'     : POL_DEG,
            'photon_energy'    : E,
            'theta'            : theta,
            'phi'              : phi,
            'trajectory'       : T,
            'e_amplitude_sigma': EFIELD_X, # todo: verify!
            'e_amplitude_pi'   : EFIELD_Y, # todo: verify!
            }

    # from srxraylib.plot.gol import plot
    # plot(out['theta'], out['radiation'].sum(axis=2).sum(axis=0), title='accumulated intensity far field')
    #
    # interpolate to cartesian grid
    #
    npointsx = theta.size
    npointsz = theta.size
    thetamax = None
    print("Interpolating to cartesian grid (%d x %d points)" % (npointsx, npointsz))
    # interpolate intensity
    radiation_interpolated, photon_energy, vx, vz = _get_radiation_interpolated_cartesian(
        radiation, E, theta, phi, npointsx=npointsx, npointsz=npointsz, thetamax=thetamax)

    out['CART_radiation'] = radiation_interpolated
    out['CART_x'] = vx  # angle in rad
    out['CART_y'] = vz  # angle in rad

    ####################################################

    #
    # back propagation
    #
    if flag_size == 2:
        #
        # recalculate....
        #
        if True:
            #
            # for sizes sample angle in [-MAXANGLE,MAXANGLE]
            #
            theta = np.linspace(-MAXANGLE, MAXANGLE, NG_T, dtype=float)
            # NG_P *= 2
            phi = (numpy.arange(NG_P) / NG_P - 0.5) * numpy.pi # np.linspace(-numpy.pi/2, numpy.pi/2, NG_P, dtype=float)
            THETA = numpy.outer(theta, numpy.ones_like(phi))
            PHI = numpy.outer(numpy.ones_like(theta), phi)
            r = distance / numpy.cos(THETA) * numpy.sin(THETA)
            xx = r * numpy.cos(PHI)
            yy = r * numpy.sin(PHI)

            EFIELD_X = np.zeros((omega_array.size, theta.size, phi.size), dtype=complex)

            print("RE-Calculating radiation...", POL_DEG.shape, theta.shape, phi.shape)

            from pySRU.Simulation import create_simulation
            from pySRU.ElectronBeam import ElectronBeam as PysruElectronBeam
            from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as PysruUndulator
            from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC
            from pySRU.RadiationFactory import RADIATION_METHOD_APPROX_FARFIELD

            for ie, e in enumerate(E):
                print("Calculating undulator source at energy %g eV (%d of %d) (%d x %d points) %d traj." % (e,
                                                                                                             ie + 1, E.size,
                                                                                                             NG_T, NG_P,
                                                                                                             Nb_pts_trajectory))
                simulation_test = create_simulation(magnetic_structure=PysruUndulator(
                    K=K,
                    period_length=LAMBDAU,
                    length=LAMBDAU * NPERIODS),
                    electron_beam=PysruElectronBeam(
                        Electron_energy=E_ENERGY,
                        I_current=INTENSITY),
                    magnetic_field=None,
                    photon_energy=e,
                    traj_method=TRAJECTORY_METHOD_ANALYTIC,
                    Nb_pts_trajectory=Nb_pts_trajectory,  # None,
                    rad_method=RADIATION_METHOD_APPROX_FARFIELD,
                    initial_condition=None,
                    distance=D,
                    X=xx.flatten(),
                    Y=yy.flatten(),
                    XY_are_list=True)

                Efield = simulation_test.electric_field._electrical_field
                tmp_x = Efield[:, 0].copy()
                tmp_x.shape = (theta.size, phi.size)
                EFIELD_X[ie, :, :] = tmp_x

            print("Done calculating radiation...")
            ca = EFIELD_X
        else:
            ca = out['e_amplitude_sigma']
        # x_in = input_wavefront.get_coordinate_x()
        # y_in = input_wavefront.get_coordinate_y()
        # print(x_in.min(), x_in.max(), y_in.min(), y_in.max())

        det_half = MAXANGLE * distance * magnification
        x_axis = numpy.linspace(-det_half, det_half, 1000)
        y_axis = numpy.linspace(-det_half, det_half, 1000)

        xx_fla = xx.flatten()
        yy_fla = yy.flatten()

        # backpropagated_efield_sigma = numpy.zeros((omega_array.size, theta.size, phi.size), dtype=complex)
        CART_BACKPROPAGATED_radiation= numpy.zeros((omega_array.size, x_axis.size), dtype=float)

        print(">>>>>ca: ", ca.shape, E.shape)
        for ie, e in enumerate(E):
            fla_complex_amplitude_propagated_x_axis = _propagate_complex_amplitude_from_arrays(
                ca[ie].flatten(),
                xx_fla,
                yy_fla,
                det_X_flatten=x_axis,
                det_Y_flatten=y_axis * 0,
                wavelength=wavelength_array_in_A[ie] * 1e-10,
                propagation_distance=-distance)

            fla_complex_amplitude_propagated_y_axis = _propagate_complex_amplitude_from_arrays(
                ca[ie].flatten(),
                xx_fla,
                yy_fla,
                det_X_flatten=x_axis * 0,
                det_Y_flatten=y_axis,
                wavelength=wavelength_array_in_A[ie] * 1e-10,
                propagation_distance=-distance)

            # area = (x_axis[1] - x_axis[0]) * numpy.abs(x_axis) * 2 * numpy.pi
            # x_norm = 1.0 / (0.025**2 * area)
            # print(x_norm)
            # plot(x_axis * 0.025, x_norm)
            ii_x = numpy.abs(fla_complex_amplitude_propagated_x_axis) ** 2
            ii_y = numpy.abs(fla_complex_amplitude_propagated_y_axis) ** 2
            ii_x /= ii_x.max()
            ii_y /= ii_y.max()
            # from srxraylib.plot.gol import plot
            # if ie < 4:
            #     plot(x_axis, ii_x,
            #          y_axis, ii_y,
            #          legend=['x axis', 'y axis'])

            CART_BACKPROPAGATED_radiation[ie, :] = ii_x

        out['BACKPROPAGATED_radiation'] = CART_BACKPROPAGATED_radiation
        out['BACKPROPAGATED_r'] = x_axis  # size in m


        #
        # plot(x_axis * 0.025, numpy.abs(fla_complex_amplitude_propagated) ** 2,
        #      y_axis * 0.025, numpy.abs(fla_complex_amplitude_propagated) ** 2,
        #      legend=['x', 'y'])
        #

        #
        # tmp = numpy.zeros((x.size, y.size), dtype=complex)
        # tmp[:, y.size // 2] = fla_complex_amplitude_propagated[0:x.size]
        # tmp[x.size // 2, :] = fla_complex_amplitude_propagated[x.size:]
        # output_wavefront = GenericWavefront2D.initialize_wavefront_from_arrays(x * 0.025, y * 0.025,
        #                                                                        tmp,
        #                                                                        wavelength=wavelength)








        # out['CART_BACKPROPAGATED_radiation'] = radiation_interpolated
        # out['CART_BACKPROPAGATED_x'] = vx * D  # size in m
        # out['CART_BACKPROPAGATED_y'] = vz * D  # size in m


    #####################################################
    return out


def calculate_undulator_emission(electron_energy              = 6.0,
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
                                 magnification                = 0.01,
                                 ):
    return _undul_phot(electron_energy,
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
                       magnification=magnification,
                       )

#
# undul_cdf
#
def undul_cdf(undul_phot_dict, method='trapz'):
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
        ZERO  = numpy.cumsum(RN0,axis=2) * (P[1] - P[0])   # CDF(e,t,p)
        ONE   = numpy.cumsum(RN1,axis=1) * (T[1] - T[0])   # CDF(e,t)
        if NG_E > 1:
            TWO = numpy.cumsum(RN2)      * (E[1] - E[0]) # CDF(e)
        else:
            TWO = numpy.array([0.0])
    else:
        RN1 = numpy.trapz(YRN0,axis=2) * (P[1]-P[0])                            # RN1(e,t)
        RN2 = numpy.trapz(RN1,axis=1)  * (T[1]-T[0])                            # RN2(e)
        ZERO  = scipy.integrate.cumtrapz(RN0,initial=0,axis=2)  * (P[1] - P[0]) # CDF(e,t,p)
        ONE   = scipy.integrate.cumtrapz(RN1,initial=0,axis=1)  * (T[1] - T[0]) # CDF(e,t)
        if NG_E > 1:
            TWO = scipy.integrate.cumtrapz(RN2,initial=0)       * (E[1] - E[0]) # CDF(e)
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

if __name__ == "__main__":

    dict1 = calculate_undulator_emission(
                     electron_energy              = 6.0,
                     electron_current             = 0.2,
                     undulator_period             = 0.025,
                     undulator_nperiods           = 188.0,
                     K                            = 1.681183,
                     photon_energy                = 5591.0,
                     EMAX                         = 5700.0,
                     NG_E                         = 1,
                     MAXANGLE                     = 30e-6,
                     number_of_points             = 151,
                     NG_P                         = 11,
                     number_of_trajectory_points  = 20,
                     flag_size                    = 2,
                     distance                     = 100.0,
                     magnification                = 0.0125,
    )

    from srxraylib.plot.gol import plot_image, plot
    # plot_image(dict1['radiation'][0], dict1['theta'], dict1['phi'], aspect='auto',  title="first", show=0)
    plot_image(dict1['radiation'][-1], dict1['theta'], numpy.degrees(dict1['phi']), aspect='auto',
               title="last", xtitle="theta", ytitle="phi [deg]", show=0)

    profiles = numpy.abs(dict1['e_amplitude_sigma'][-1])**2
    plot(100 * dict1['theta'], profiles[:, 0],
         100 * dict1['theta'], profiles[:, 2],
         100 * dict1['theta'], profiles[:, 4],
         100 * dict1['theta'], profiles[:, 6],
         100 * dict1['theta'], profiles[:, 8],
         100 * dict1['theta'], profiles[:, 10],
         xtitle="r [m] (%d points)" % (dict1['theta'].size),
         ytitle="intensity (%d points)" % (profiles.shape[0]),
         )

    ii = numpy.abs(dict1['e_amplitude_sigma'])**2 + numpy.abs(dict1['e_amplitude_pi'])**2
    plot_image(ii[-1], dict1['theta'], numpy.degrees(dict1['phi']), aspect='auto',
               title="last (from amplitudes)", xtitle="theta", ytitle="phi [deg]", show=1)

    if True:
        plot_image(dict1['CART_radiation'][0],
                   100 * 1e6 * dict1['CART_x'], 100 * 1e6 * dict1['CART_y'], aspect='auto', title="first", show=0)
        plot_image(dict1['CART_radiation'][-1],
                   100 * 1e6 * dict1['CART_x'], 100 * 1e6 * dict1['CART_y'], aspect='auto', title="last", show=1)


    # BACK...

    ii = dict1['BACKPROPAGATED_radiation']
    plot(1e6 * dict1['BACKPROPAGATED_r'], ii[-1],
               title="backpropagated last", xtitle="x [um]", ytitle="intensity", show=1)

    # plot_image(dict1['CART_BACKPROPAGATED_radiation'][-1],
    #            1e6 * dict1['CART_BACKPROPAGATED_x'], 1e6 * dict1['CART_BACKPROPAGATED_y'], aspect='auto',
    #            title="last CART BACKPROPAGATED", xtitle="x [um]", ytitle="y [um]", show=1)

    if False:
        plot_image(dict1['CART_radiation'][0],
                   100 * 1e6 * dict1['CART_x'], 100 * 1e6 * dict1['CART_y'], aspect='auto', title="first", show=0)
        plot_image(dict1['CART_radiation'][-1],
                   100 * 1e6 * dict1['CART_x'], 100 * 1e6 * dict1['CART_y'], aspect='auto', title="last", show=1)