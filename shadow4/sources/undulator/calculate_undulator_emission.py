"""
Calculation of undulator far field and backpropagation to the center of the undulator
Implementation using internal code (no dependency) hacked from pySRU for source simulation
and Wofry2D for backpropagation.


Available public functions:

    calculate_undulator_emission()

    undul_cdf() [like undul_cdf in SHADOW3].

"""
import numpy
import numpy as np
from scipy import interpolate
import time

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
                                          npointsx=100, npointsz=100, thetamax=None,
                                          distance=100):
    if thetamax is None:
        thetamax = thetabm.max()

    distancemax = 1.0 * thetamax * distance

    #v are coordinates at the image plane
    vx = numpy.linspace(-distancemax, distancemax, npointsx)
    vz = numpy.linspace(-distancemax, distancemax, npointsz)
    VX = numpy.outer(vx, numpy.ones_like(vz))
    VZ = numpy.outer(numpy.ones_like(vx), vz)
    VY = numpy.sqrt(distance ** 2 +  VX ** 2 + VZ ** 2)

    THETA = numpy.arctan(numpy.sqrt(VX ** 2 + VZ ** 2) / VY)
    # abs to be sure that PHI is in [0,90] deg (like the polar data)
    PHI = numpy.arctan2(numpy.abs(VZ), numpy.abs(VX))


    # from srxraylib.plot.gol import plot_scatter
    # plot_scatter(THETA, PHI, xtitle="THETA", ytitle="PHI", show=0)
    # plot_scatter(numpy.outer(thetabm, numpy.ones_like(phi)),
    #              numpy.outer(numpy.ones_like(thetabm), phi ),
    #              xtitle="thetabm", ytitle="phi", show=1)

    radiation_interpolated = numpy.zeros((radiation.shape[0], npointsx, npointsz))

    for i in range(radiation.shape[0]):
        interpolator_value = interpolate.RectBivariateSpline(thetabm, phi, radiation[i])
        radiation_interpolated[i] = interpolator_value.ev(THETA, PHI)

    return radiation_interpolated, photon_energy, vx, vz

# hacked from Wofry (2D propagator integral)
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
               number_of_trajectory_points=20, flag_size=0, distance=100.0, magnification=0.01,
                flag_backprop_recalculate_source=0,
                flag_backprop_weight=0,
                weight_ratio=0.5,
                ):
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
    print("Calculating trajectory...")
    T = _pysru_analytical_trajectory_plane_undulator(K=K, gamma=gamma, lambda_u=LAMBDAU, Nb_period=NPERIODS,
                                        Nb_point=number_of_trajectory_points,
                                        Beta_et=Beta_et)
    print("Done calculating trajectory...")

    #
    # polar grid
    #

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

    t0 = time.time()
    print("Calculating radiation...", POL_DEG.shape, theta.shape, phi.shape)

    Nb_pts_trajectory = T.shape[1]

    use_pysru = 0 # todo: delete pysru code...

    if use_pysru:
        from pySRU.Simulation import create_simulation
        from pySRU.ElectronBeam import ElectronBeam as PysruElectronBeam
        from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as PysruUndulator
        from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC
        from pySRU.RadiationFactory import RADIATION_METHOD_APPROX_FARFIELD

        for ie, e in enumerate(E):
            print("   Calculating undulator source at energy %g eV (%d of %d) (%d x %d points) %d traj."%(e,
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
                                                distance=distance,
                                                X=xx.flatten(),
                                                Y=yy.flatten(),
                                                XY_are_list=True)

            Efield = simulation_test.electric_field._electrical_field

            # Conversion from pySRU units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
            coeff = (distance * 1e3) ** 2  # photons/mm^2 -> photons/rad^2
            coeff /= 1e-3 * e  # photons/o.1%bw -> photons/eV

            tmp_x = Efield[:, 0].copy()
            tmp_y = Efield[:, 1].copy()
            tmp_x.shape = (theta.size, phi.size)
            tmp_y.shape = (theta.size, phi.size)

            EFIELD_X[ie, :, :] = tmp_x * numpy.sqrt(coeff)
            EFIELD_Y[ie, :, :] = tmp_y * numpy.sqrt(coeff)

            POL_DEG[ie, :, :] = numpy.abs(tmp_x) / (numpy.abs(tmp_x) + numpy.abs(tmp_y))  # SHADOW definition
    else:
        for o in range(omega_array.size):
            print("   Calculating energy %8.3f eV (%d of %d)"%(E[o],o+1,omega_array.size))
            for t in range(theta.size):
                for p in range(phi.size):
                    R = distance / np.cos(theta[t])
                    r = R * np.sin(theta[t])
                    X = r * np.cos(phi[p])
                    Y = r * np.sin(phi[p])
                    Efield = _pysru_energy_radiated_approximation_and_farfield(omega=omega_array[o],
                                                                                  electron_current=INTENSITY,
                                                                                  trajectory=T,
                                                                                  x=X ,
                                                                                  y=Y,
                                                                                  D=distance)

                    # Conversion from pySRU units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
                    coeff = (distance * 1e3) ** 2  # photons/mm^2 -> photons/rad^2
                    coeff /= 1e-3 * E[o]  # photons/o.1%bw -> photons/eV

                    EFIELD_X[o, t, p] = Efield[0] * numpy.sqrt(coeff)
                    EFIELD_Y[o, t, p] = Efield[1] * numpy.sqrt(coeff)
                    POL_DEG [o, t, p] = numpy.abs(Efield[0]) / (numpy.abs(Efield[0]) + numpy.abs(Efield[1]))  # SHADOW definition


    print("Done calculating radiation... (%f s)" % (time.time() - t0))

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
        radiation, E, theta, phi, npointsx=npointsx, npointsz=npointsz, thetamax=thetamax,
        distance=distance)

    out['CART_radiation'] = radiation_interpolated
    out['CART_x'] = vx / distance # angle in rad
    out['CART_y'] = vz / distance # angle in rad

    ####################################################

    #
    # back propagation
    #

    # todo: p=polarization
    if flag_size == 2:
        #
        # recalculate or reuse the radiation....
        #
        if flag_backprop_recalculate_source:
            #
            # Redefine grids
            # for sizes sample angle in [-MAXANGLE,MAXANGLE]
            #
            theta = np.linspace(-MAXANGLE, MAXANGLE, NG_T, dtype=float)
            phi = (numpy.arange(NG_P) / NG_P - 0.5) * numpy.pi # np.linspace(-numpy.pi/2, numpy.pi/2, NG_P, dtype=float)
            THETA = numpy.outer(theta, numpy.ones_like(phi))
            PHI = numpy.outer(numpy.ones_like(theta), phi)
            r = distance / numpy.cos(THETA) * numpy.sin(THETA)
            xx = r * numpy.cos(PHI)
            yy = r * numpy.sin(PHI)

            EFIELD_X = np.zeros((omega_array.size, theta.size, phi.size), dtype=complex)
            t0 = time.time()
            print("RE-Calculating radiation... ", POL_DEG.shape, theta.shape, phi.shape)

            if use_pysru:
                from pySRU.Simulation import create_simulation
                from pySRU.ElectronBeam import ElectronBeam as PysruElectronBeam
                from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as PysruUndulator
                from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC
                from pySRU.RadiationFactory import RADIATION_METHOD_APPROX_FARFIELD

                for ie, e in enumerate(E):
                    print("  Calculating undulator source at energy %g eV (%d of %d) (%d x %d points) %d traj." % (e,
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
                        distance=distance,
                        X=xx.flatten(),
                        Y=yy.flatten(),
                        XY_are_list=True)

                    Efield = simulation_test.electric_field._electrical_field
                    tmp_x = Efield[:, 0].copy()
                    tmp_x.shape = (theta.size, phi.size)
                    EFIELD_X[ie, :, :] = tmp_x
            else:
                for o in range(omega_array.size):
                    print("   Calculating energy %8.3f eV (%d of %d)" % (E[o], o + 1, omega_array.size))
                    for t in range(theta.size):
                        for p in range(phi.size):
                            R = distance / np.cos(theta[t])
                            r = R * np.sin(theta[t])
                            X = r * np.cos(phi[p])
                            Y = r * np.sin(phi[p])
                            Efield = _pysru_energy_radiated_approximation_and_farfield(omega=omega_array[o],
                                                                                       electron_current=INTENSITY,
                                                                                       trajectory=T,
                                                                                       x=X,
                                                                                       y=Y,
                                                                                       D=distance)

                            # Conversion from pySRU units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
                            coeff = (distance * 1e3) ** 2  # photons/mm^2 -> photons/rad^2
                            coeff /= 1e-3 * E[o]  # photons/o.1%bw -> photons/eV

                            EFIELD_X[o, t, p] = Efield[0] * numpy.sqrt(coeff)

            print("Done re-calculating radiation... (%f s)" % (time.time() - t0))
            ca = EFIELD_X
            do_concatenate = 0
        else: # reusing xx, yy
            print("Re using radiation...", POL_DEG.shape, theta.shape, phi.shape)
            ca = out['e_amplitude_sigma']
            do_concatenate = 1

        # detector grid (calculate only the X axis)
        det_half = MAXANGLE * distance * magnification
        x_axis = numpy.linspace(-det_half, det_half, 1000)
        y_axis = numpy.linspace(-det_half, det_half, 1000)

        xx_fla = xx.flatten()
        yy_fla = yy.flatten()

        CART_BACKPROPAGATED_radiation= numpy.zeros((omega_array.size, x_axis.size), dtype=float)

        for ie, e in enumerate(E):
            ca_fla = ca[ie].flatten()

            if do_concatenate: # copy data from first quadrant to the others
                ca_fla = numpy.concatenate((ca_fla,  ca_fla,  ca_fla,  ca_fla))
                xx_fla = numpy.concatenate((xx_fla, -xx_fla, -xx_fla,  xx_fla))
                yy_fla = numpy.concatenate((yy_fla,  yy_fla, -yy_fla, -yy_fla))

            fla_complex_amplitude_propagated_x_axis = _propagate_complex_amplitude_from_arrays(
                ca_fla,
                xx_fla,
                yy_fla,
                det_X_flatten=x_axis,
                det_Y_flatten=numpy.zeros_like(y_axis),
                wavelength=wavelength_array_in_A[ie] * 1e-10,
                propagation_distance=-distance)
            ii_x = numpy.abs(fla_complex_amplitude_propagated_x_axis) ** 2
            ii_x /= ii_x.max()

            if False:
                fla_complex_amplitude_propagated_y_axis = _propagate_complex_amplitude_from_arrays(
                    ca_fla,
                    xx_fla,
                    yy_fla,
                    det_X_flatten=numpy.zeros_like(x_axis),
                    det_Y_flatten=y_axis,
                    wavelength=wavelength_array_in_A[ie] * 1e-10,
                    propagation_distance=-distance)

                ii_y = numpy.abs(fla_complex_amplitude_propagated_y_axis) ** 2
                ii_y /= ii_y.max()
                from srxraylib.plot.gol import plot
                if ie < 4:
                    plot(x_axis, ii_x,
                         y_axis, ii_y,
                         legend=['x axis', 'y axis'])

            #
            # add a Gaussian weight (the Gaussian affects the amplitude, thus squared for intensities)
            #
            if flag_backprop_weight:
                xmax = x_axis.max()
                sigma = weight_ratio * xmax
                weight = numpy.exp(-(x_axis - x_axis.mean()) ** 2 / (2 * sigma ** 2))
            else:
                weight = numpy.ones_like(ii_x)

            CART_BACKPROPAGATED_radiation[ie, :] = ii_x * weight**2

        out['BACKPROPAGATED_radiation'] = CART_BACKPROPAGATED_radiation
        out['BACKPROPAGATED_r'] = x_axis  # size in m

    return out


def calculate_undulator_emission(electron_energy                  = 6.0,
                                 electron_current                 = 0.2,
                                 undulator_period                 = 0.018,
                                 undulator_nperiods               = 100,
                                 K                                = 1.0,
                                 photon_energy                    = 2000.0,
                                 EMAX                             = 20000.0,
                                 NG_E                             = 10,
                                 MAXANGLE                         = 0.1,
                                 number_of_points                 = 100,
                                 NG_P                             = 100,
                                 number_of_trajectory_points      = 100,
                                 flag_size                        = 2,
                                 distance                         = 100.0,
                                 magnification                    = 0.01,
                                 flag_backprop_recalculate_source = 0,
                                 flag_backprop_weight             = 0,
                                 weight_ratio                     = 0.5,
                                 ):
    """
    Calculate undulator emission (far field) and backpropagation to the center of the undulator using internal code.

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
    magnification: float, optional
        The magnification for backpropagation.
    flag_backprop_recalculate_source: int, optional
        for code_undul_phot in ["internal", "pysru"] and flag_size=2: apply Gaussian weight to backprop amplitudes (0=No, 1=Yes).
    flag_backprop_weight: int, optional
        for code_undul_phot in ["internal", "pysru"] and flag_size=2: apply Gaussian weight to backprop amplitudes.
    weight_ratio: float, optional
        for flag_backprop_weight=1: the Gaussian sigma in units of r.max().

    Returns
    -------
    dict
        A dictionary with calculated values and arrays.

    """
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
                       flag_backprop_recalculate_source=flag_backprop_recalculate_source,
                       flag_backprop_weight=flag_backprop_weight,
                       weight_ratio=weight_ratio,
                       )

if __name__ == "__main__":

    dict1 = calculate_undulator_emission(
                     electron_energy                  = 6.0,
                     electron_current                 = 0.2,
                     undulator_period                 = 0.025,
                     undulator_nperiods               = 188.0,
                     K                                = 1.681183,
                     photon_energy                    = 5591.0,
                     EMAX                             = 5700.0,
                     NG_E                             = 1,
                     MAXANGLE                         = 30e-6,
                     number_of_points                 = 151,
                     NG_P                             = 11,
                     number_of_trajectory_points      = 20,
                     flag_size                        = 2,
                     distance                         = 100.0,
                     magnification                    = 0.0125,
                     flag_backprop_recalculate_source = 1,
                     flag_backprop_weight             = 1,
                     weight_ratio                     = 0.5,
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

    if False:
        plot_image(dict1['CART_radiation'][0],
                   1e6 * dict1['CART_x'], 1e6 * dict1['CART_y'], aspect='auto', title="first",
                   xtitle="theta_x [um]", ytitle="theta_y [um]", show=0)
        plot_image(dict1['CART_radiation'][-1],
                   1e6 * dict1['CART_x'], 1e6 * dict1['CART_y'], aspect='auto', title="last",
                   xtitle="theta_x [um]", ytitle="theta_y [um]", show=1)

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