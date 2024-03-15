"""
Calculation of undulator far field and backpropagation to the center of the undulator
Implementation using ppySRU for source simulation and Wofry2D for backpropagation.


Available public functions:

    calculate_undulator_emission_pysru()


"""

try:
    from pySRU.ElectronBeam import ElectronBeam as PysruElectronBeam
    from pySRU.MagneticStructureUndulatorPlane import MagneticStructureUndulatorPlane as PysruUndulator
    from pySRU.Simulation import create_simulation
    from pySRU.TrajectoryFactory import TRAJECTORY_METHOD_ANALYTIC
    from pySRU.RadiationFactory import RADIATION_METHOD_APPROX_FARFIELD
except:
    raise ImportError("pySRU not imported")

# for propagation, etc.
from syned.beamline.beamline_element import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates

import scipy.constants as codata
from srxraylib.plot.gol import plot, plot_image

try:

    from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters

    from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D

    from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
    from wofryimpl.propagator.propagators2D.fresnel import Fresnel2D
    from wofryimpl.propagator.propagators2D.fresnel_convolution import FresnelConvolution2D
    from wofryimpl.propagator.propagators2D.fraunhofer import Fraunhofer2D
    from wofryimpl.propagator.propagators2D.integral import Integral2D
    from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D

    # from orangecontrib.esrf.wofry.util.light_source import WOPySRULightSource  # TODO: from wofryimpl...
    from wofryimpl.propagator.light_source_pysru import WOPySRULightSource
    from wofryimpl.beamline.beamline import WOBeamline
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen

    from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters

    from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
    from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen

    from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D

except:
    raise ImportError("wofry or wofryimpl not imported")

import numpy
from scipy import interpolate
from shadow4.sources.undulator.calculate_undulator_emission import _get_radiation_interpolated_cartesian

def _pysru_wofry_2D_run(photon_energy,
                        energy_in_GeV=6,
                        current=0.2,
                        K_vertical=1.68118,
                        period_length=0.025,
                        number_of_periods=188,
                        distance=100,
                        gapH=0.006,
                        gapV=0.006,
                        h_slit_points=51,
                        v_slit_points=51,
                        number_of_trajectory_points=5000,
                        traj_method=1,
                        rad_method=2,
                        magnification=0.025,
                        ):
    #
    # Import section
    #
    # import numpy

    # from syned.beamline.beamline_element import BeamlineElement
    # from syned.beamline.element_coordinates import ElementCoordinates
    # from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
    #
    # from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
    #
    # from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
    # from wofryimpl.propagator.propagators2D.fresnel import Fresnel2D
    # from wofryimpl.propagator.propagators2D.fresnel_convolution import FresnelConvolution2D
    # from wofryimpl.propagator.propagators2D.fraunhofer import Fraunhofer2D
    # from wofryimpl.propagator.propagators2D.integral import Integral2D
    # from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
    #
    # from srxraylib.plot.gol import plot, plot_image
    plot_from_oe = 1000  # set to a large number to avoid plots

    for i in range(photon_energy.size):
        ##########  SOURCE ##########

        #
        # create output_wavefront
        #
        #
        # from orangecontrib.esrf.wofry.util.light_source import WOPySRULightSource  # TODO: from wofryimpl...
        light_source = WOPySRULightSource.initialize_from_keywords(
            energy_in_GeV=energy_in_GeV,
            current=current,
            K_vertical=K_vertical,
            period_length=period_length,
            number_of_periods=number_of_periods,
            distance=distance,
            gapH=gapH,
            gapV=gapV,
            photon_energy=photon_energy[i],
            h_slit_points=h_slit_points,
            v_slit_points=v_slit_points,
            number_of_trajectory_points=number_of_trajectory_points,
            traj_method=traj_method,
            rad_method=rad_method,
        )

        output_wavefront = light_source.get_wavefront()

        if plot_from_oe <= 0: plot_image(output_wavefront.get_intensity(), output_wavefront.get_coordinate_x(),
                                         output_wavefront.get_coordinate_y(), aspect='auto', title='SOURCE')

        if i == 0:
            tmp_s = output_wavefront.get_complex_amplitude(polarization=0).copy()
            tmp_p = output_wavefront.get_complex_amplitude(polarization=1).copy()
            CART_e_amplitude_sigma = numpy.zeros((photon_energy.size, tmp_s.shape[0], tmp_s.shape[1]), dtype=complex)
            CART_e_amplitude_pi = numpy.zeros_like(CART_e_amplitude_sigma)
            CART_e_amplitude_sigma[0,:,:] = tmp_s
            CART_e_amplitude_pi[0,:,:] = tmp_p
        else:
            CART_e_amplitude_sigma[i,:,:] = output_wavefront.get_complex_amplitude(polarization=0).copy()
            CART_e_amplitude_pi[i,:,:] =    output_wavefront.get_complex_amplitude(polarization=1).copy()

            ##########  OPTICAL SYSTEM ##########

        ##########  OPTICAL ELEMENT NUMBER 1 ##########

        input_wavefront = output_wavefront.duplicate()
        # from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen

        optical_element = WOScreen()

        # drift_before -100 m
        #
        # propagating
        #
        #
        propagation_elements = PropagationElements()
        beamline_element = BeamlineElement(optical_element=optical_element,
                                           coordinates=ElementCoordinates(p=-distance,
                                                                          q=0.000000,
                                                                          angle_radial=numpy.radians(0.000000),
                                                                          angle_azimuthal=numpy.radians(0.000000)))
        propagation_elements.add_beamline_element(beamline_element)
        propagation_parameters = PropagationParameters(wavefront=input_wavefront, propagation_elements=propagation_elements)
        # self.set_additional_parameters(propagation_parameters)
        #
        propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
        propagation_parameters.set_additional_parameters('magnification_x', magnification)
        propagation_parameters.set_additional_parameters('magnification_y', magnification)
        #
        propagator = PropagationManager.Instance()
        try:
            propagator.add_propagator(FresnelZoomXY2D())
        except:
            pass
        output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                     handler_name='FRESNEL_ZOOM_XY_2D')

        if i == 0:
            tmp_s = output_wavefront.get_complex_amplitude(polarization=0).copy()
            tmp_p = 0 # output_wavefront.get_complex_amplitude(polarization=1).copy() # todo: add polarization in wofry propagation
            CART_BACKPROPAGATED_radiation = numpy.zeros((photon_energy.size, tmp_s.shape[0], tmp_s.shape[1]), dtype=float)
            CART_BACKPROPAGATED_radiation[0, :, :] = numpy.abs(tmp_s)**2 + numpy.abs(tmp_p)**2
            CART_BACKPROPAGATED_x = output_wavefront.get_coordinate_x()
            CART_BACKPROPAGATED_y = output_wavefront.get_coordinate_y()
        else:
            tmp_s = output_wavefront.get_complex_amplitude(polarization=0).copy()
            tmp_p = output_wavefront.get_complex_amplitude(polarization=1).copy()
            CART_BACKPROPAGATED_radiation[i, :, :] = numpy.abs(tmp_s)**2 + numpy.abs(tmp_p)**2

        #
        # ---- plots -----
        #
        if plot_from_oe <= 1: plot_image(output_wavefront.get_intensity(), output_wavefront.get_coordinate_x(),
                                         output_wavefront.get_coordinate_y(), aspect='auto', title='OPTICAL ELEMENT NR 1')

        # for completeness, define python code
        if i == 0:
            # from wofryimpl.beamline.beamline import WOBeamline
            bl = WOBeamline(light_source=light_source)
            propagator_info = {
                "propagator_class_name": "FresnelZoomXY2D",
                "propagator_handler_name": FresnelZoomXY2D.HANDLER_NAME,
                "propagator_additional_parameters_names": ['shift_half_pixel', 'magnification_x','magnification_y'],
                "propagator_additional_parameters_values": [1, magnification, magnification]}
            bl.append_beamline_element(beamline_element, propagator_info)

            print("\n\n################## PYTHON CODE #########################\n")
            print(bl.to_python_code())
            print("\n################## END PYTHON CODE #########################\n")

        return CART_e_amplitude_sigma, CART_e_amplitude_pi, CART_BACKPROPAGATED_radiation,\
               CART_BACKPROPAGATED_x, CART_BACKPROPAGATED_y

def _undul_phot_pySRU(E_ENERGY, INTENSITY, LAMBDAU, NPERIODS, K, EMIN, EMAX, NG_E, MAXANGLE, NG_T, NG_P,
                      number_of_trajectory_points=20, flag_size=2,
                      distance=100.0,
                      magnification=0.05,
                      flag_backprop_recalculate_source=0, # 0=interpolate, 1=recalculate
                      flag_backprop_weight=0,
                      weight_ratio=0.5,
                      ):
    myelectronbeam = PysruElectronBeam(Electron_energy=E_ENERGY, I_current=INTENSITY)
    myundulator = PysruUndulator(K=K, period_length=LAMBDAU, length=LAMBDAU*NPERIODS)
    D = distance
    #
    # polar grid matrix
    #
    photon_energy = numpy.linspace(EMIN, EMAX, NG_E, dtype=float)

    intens = numpy.zeros((NG_E, NG_T, NG_P))
    pol_deg = numpy.zeros_like(intens)
    theta = numpy.linspace(0, MAXANGLE, NG_T, dtype=float)
    phi = numpy.linspace(0, numpy.pi / 2, NG_P, dtype=float)
    efield_x = numpy.zeros_like(intens, dtype=complex)
    efield_y = numpy.zeros_like(intens, dtype=complex)
    efield_z = numpy.zeros_like(intens, dtype=complex)


    THETA = numpy.outer(theta, numpy.ones_like(phi))
    PHI = numpy.outer(numpy.ones_like(theta), phi)

    X = (D / numpy.cos(THETA)) * numpy.sin(THETA) * numpy.cos(PHI)
    Y = (D / numpy.cos(THETA)) * numpy.sin(THETA) * numpy.sin(PHI)
    Nb_pts_trajectory = int(number_of_trajectory_points * NPERIODS)

    for ie, e in enumerate(photon_energy):
        print("Calculating undulator source at energy %g eV (%d of %d) (%d x %d points) %d traj."%(e,
                                            ie+1, photon_energy.size, NG_T, NG_P, Nb_pts_trajectory))
        simulation_test = create_simulation(magnetic_structure=myundulator,
                                            electron_beam=myelectronbeam,
                                            magnetic_field=None,
                                            photon_energy=e,
                                            traj_method=TRAJECTORY_METHOD_ANALYTIC,
                                            Nb_pts_trajectory=Nb_pts_trajectory, #None,
                                            rad_method=RADIATION_METHOD_APPROX_FARFIELD,
                                            initial_condition=None,
                                            distance=D,
                                            X=X.flatten(),
                                            Y=Y.flatten(),
                                            XY_are_list=True)

        E = simulation_test.electric_field._electrical_field
        EEx = E[:, 0].copy()
        EEy = E[:, 1].copy()
        EEz = E[:, 2].copy()
        EEx.shape = (theta.size, phi.size)
        EEy.shape = (theta.size, phi.size)
        EEz.shape = (theta.size, phi.size)

        # pol_deg1 = (numpy.abs(E[:,0])**2 / (numpy.abs(E[:,0])**2 + numpy.abs(E[:,1])**2)).flatten()
        pol_deg1 = (numpy.abs(E[:, 0]) / (numpy.abs(E[:, 0]) + numpy.abs(E[:, 1]))).flatten() # SHADOW definition!!
        pol_deg1.shape = (theta.size, phi.size)

        intens1 = simulation_test.radiation.intensity.copy()
        intens1.shape = (theta.size, phi.size)
        #  Conversion from pySRU units (photons/mm^2/0.1%bw) to SHADOW units (photons/rad^2/eV)
        coeff = (D * 1e3)**2 # photons/mm^2 -> photons/rad^2
        coeff /= 1e-3 * e # photons/o.1%bw -> photons/eV
        intens1 *= coeff

        # store
        intens[ie] = intens1
        pol_deg[ie] = pol_deg1
        efield_x[ie] = EEx * numpy.sqrt(coeff)
        efield_y[ie] = EEy * numpy.sqrt(coeff)
        efield_z[ie] = EEz * numpy.sqrt(coeff)

    T0 = simulation_test.trajectory
    T = numpy.vstack((T0.t, T0.x, T0.y, T0.z, T0.v_x, T0.v_y, T0.v_z, T0.a_x, T0.a_y, T0.a_z))

    # interpolate to cartesian grid
    print("Interpolating to cartesian grid (%d x %d points)" % (theta.size, theta.size))
    # interpolate intensity
    radiation_interpolated, photon_energy, vx, vz = _get_radiation_interpolated_cartesian(
        intens, photon_energy, theta, phi, npointsx=theta.size, npointsz=theta.size, thetamax=None,
        distance=distance)

    out = {'radiation'         : intens,
            'polarization'     : pol_deg,
            'photon_energy'    : photon_energy,
            'theta'            : theta,
            'phi'              : phi,
            'trajectory'       : T,
            'e_amplitude_sigma': efield_x,  # todo: verify!
            'e_amplitude_pi'   : efield_y,  # todo: verify!
            'CART_radiation'   : radiation_interpolated,
            'CART_x'           : vx / distance,  # angle in rad
            'CART_y'           : vz / distance,  # angle in rad
            }

    # for backpropagation we need the electric field (interpolated or recalculated)
    # efield_method = 'recalculate' # 'interpolate' # 'recalculate'
    if flag_size == 2:
        npointsx = theta.size
        npointsz = npointsx
        thetamax = None
        if flag_backprop_recalculate_source == 0: # 'interpolate'

            efield_x_mod_interpolated, photon_energy, vx, vz = _get_radiation_interpolated_cartesian(
                numpy.abs(efield_x), photon_energy, theta, phi, npointsx=npointsx, npointsz=npointsz,
                thetamax=thetamax, distance=distance)
            efield_x_angle_interpolated, photon_energy, vx, vz = _get_radiation_interpolated_cartesian(
                numpy.angle(efield_x), photon_energy, theta, phi, npointsx=npointsx, npointsz=npointsz,
                thetamax=thetamax, distance=distance)
            efield_y_mod_interpolated, photon_energy, vx, vz = _get_radiation_interpolated_cartesian(
                numpy.abs(efield_y), photon_energy, theta, phi, npointsx=npointsx, npointsz=npointsz,
                thetamax=thetamax, distance=distance)
            efield_y_angle_interpolated, photon_energy, vx, vz = _get_radiation_interpolated_cartesian(
                numpy.angle(efield_y), photon_energy, theta, phi, npointsx=npointsx, npointsz=npointsz,
                thetamax=thetamax, distance=distance)
            efield_x_interpolated = efield_x_mod_interpolated * numpy.exp(1j * efield_x_angle_interpolated)
            efield_y_interpolated = efield_y_mod_interpolated * numpy.exp(1j * efield_y_angle_interpolated)
            out['CART_e_amplitude_sigma'] = efield_x_interpolated
            out['CART_e_amplitude_pi'] = efield_y_interpolated
            ############################################################################################################
            ##########  OPTICAL ELEMENT NUMBER 1 ##########

            # from syned.beamline.beamline_element import BeamlineElement
            # from syned.beamline.element_coordinates import ElementCoordinates
            # from wofry.propagator.propagator import PropagationManager, PropagationElements, PropagationParameters
            #
            # from wofry.propagator.wavefront2D.generic_wavefront import GenericWavefront2D
            # from wofryimpl.beamline.optical_elements.ideal_elements.screen import WOScreen
            #
            # from wofryimpl.propagator.propagators2D.fresnel_zoom_xy import FresnelZoomXY2D
            # import scipy.constants as codata

            backpropagated_intensity = numpy.zeros((photon_energy.size, npointsx, npointsz))
            for ie, e in enumerate(photon_energy):
                print("Backpropagating energy %g eV (%d of %d)" % (e, ie + 1, photon_energy.size))
                print("wavelength: ", codata.h * codata.c / codata.e / e)

                input_wavefront = GenericWavefront2D.initialize_wavefront_from_arrays(
                    vx, vz, efield_x_interpolated[ie],
                    z_array_pi=efield_y_interpolated[ie],
                    wavelength=codata.h * codata.c / codata.e / e)

                optical_element = WOScreen()

                # drift_before -100 m
                #
                # propagating
                #
                #
                propagation_elements = PropagationElements()
                beamline_element = BeamlineElement(optical_element=optical_element,
                                                   coordinates=ElementCoordinates(p=-D, q=0.000000,
                                                                                  angle_radial=numpy.radians(0.000000),
                                                                                  angle_azimuthal=numpy.radians(
                                                                                      0.000000)))
                propagation_elements.add_beamline_element(beamline_element)
                propagation_parameters = PropagationParameters(wavefront=input_wavefront,
                                                               propagation_elements=propagation_elements)
                # self.set_additional_parameters(propagation_parameters)
                #
                propagation_parameters.set_additional_parameters('shift_half_pixel', 1)
                propagation_parameters.set_additional_parameters('magnification_x', magnification)
                propagation_parameters.set_additional_parameters('magnification_y', magnification)
                #
                propagator = PropagationManager.Instance()
                try:
                    propagator.add_propagator(FresnelZoomXY2D())
                except:
                    pass
                output_wavefront = propagator.do_propagation(propagation_parameters=propagation_parameters,
                                                             handler_name='FRESNEL_ZOOM_XY_2D')

                backpropagated_intensity[ie, :, :] = output_wavefront.get_intensity(polarization=0)

                print(beamline_element.info())

            out['CART_BACKPROPAGATED_radiation'] = backpropagated_intensity
            out['CART_BACKPROPAGATED_x'] = output_wavefront.get_coordinate_x()  # length in m
            out['CART_BACKPROPAGATED_y'] = output_wavefront.get_coordinate_y()  # length in m


        elif flag_backprop_recalculate_source == 1: # 'recalculate':
            CART_e_amplitude_sigma, CART_e_amplitude_pi, CART_BACKPROPAGATED_radiation,\
                CART_BACKPROPAGATED_x, CART_BACKPROPAGATED_y = _pysru_wofry_2D_run(photon_energy,
                                                                                   energy_in_GeV=E_ENERGY,
                                                                                   current=INTENSITY,
                                                                                   K_vertical=K,
                                                                                   period_length=LAMBDAU,
                                                                                   number_of_periods=NPERIODS,
                                                                                   distance=distance,
                                                                                   gapH=2 * MAXANGLE * distance,
                                                                                   gapV=2 * MAXANGLE * distance,
                                                                                   # photon_energy=photon_energy[i],
                                                                                   h_slit_points=NG_T,
                                                                                   v_slit_points=NG_T,
                                                                                   number_of_trajectory_points=int(number_of_trajectory_points * NPERIODS),
                                                                                   traj_method=1,
                                                                                   rad_method=2,
                                                                                   magnification=magnification,
                                                                                   )
            out['CART_e_amplitude_sigma']        = CART_e_amplitude_sigma
            out['CART_e_amplitude_pi']           = CART_e_amplitude_pi

            #
            # add a Gaussian weight (the Gaussian affects the amplitude, thus squared for intensities)
            #
            if flag_backprop_weight:
                x_axis = CART_BACKPROPAGATED_x
                y_axis = CART_BACKPROPAGATED_y
                xmax = x_axis.max()
                sigma = weight_ratio * xmax
                weight = numpy.exp(-(x_axis - x_axis.mean()) ** 2 / (2 * sigma ** 2))
                weight = weight**2
                WEIGHT = numpy.outer(weight, weight)
                ne = CART_BACKPROPAGATED_radiation.shape[0]
                for i in range(ne):
                    CART_BACKPROPAGATED_radiation[i, :, :] = CART_BACKPROPAGATED_radiation[i, :, :] * WEIGHT

            out['CART_BACKPROPAGATED_radiation'] = CART_BACKPROPAGATED_radiation
            out['CART_BACKPROPAGATED_x']         = CART_BACKPROPAGATED_x  # length in m
            out['CART_BACKPROPAGATED_y']         = CART_BACKPROPAGATED_y  # length in m

        else:
            raise Exception("Invalid flag_backprop_recalculate_source: %d" % flag_backprop_recalculate_source)

    return out


def calculate_undulator_emission_pysru(
                     electron_energy                  = 6.0,
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
                     magnification                    = 0.05,
                     flag_backprop_recalculate_source = 0,
                     flag_backprop_weight             = 0,
                     weight_ratio                     = 0.5,
                     ):
    """
    Calculate undulator emission (far field) and backpropagation to the center of the undulator using pysru+wofry.

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
    return _undul_phot_pySRU(
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
                        magnification=magnification,
                        flag_backprop_recalculate_source=flag_backprop_recalculate_source,
                        flag_backprop_weight=flag_backprop_weight,
                        weight_ratio=weight_ratio,
                        )

if __name__ == "__main__":
    import numpy
    dict1 = calculate_undulator_emission_pysru(
                     electron_energy                  = 6.0,
                     electron_current                 = 0.2,
                     undulator_period                 = 0.025,
                     undulator_nperiods               = 188.0,
                     K                                = 1.681183,
                     photon_energy                    = 5591.0,
                     EMAX                             = 5700.0,
                     NG_E                             = 1,
                     MAXANGLE                         = 30e-6,
                     number_of_points                 = 31,
                     NG_P                             = 11,
                     number_of_trajectory_points      = 20,
                     flag_size                        = 2,
                     distance                         = 100.0,
                     magnification                    = 0.05,
                     flag_backprop_recalculate_source = 1,
                     flag_backprop_weight             = 1,
                     weight_ratio                     = 0.5,
    )
    from srxraylib.plot.gol import plot_image, plot
    if True:

        plot_image(dict1['radiation'][0], dict1['theta'], dict1['phi'], aspect='auto',  title="first", show=0)
        plot_image(dict1['radiation'][-1], dict1['theta'], dict1['phi'], aspect='auto', title="last", show=1)

        plot_image(dict1['CART_radiation'][0], 1e6 * dict1['CART_x'], 1e6 * dict1['CART_y'], aspect='auto', title="first", show=0)
        plot_image(dict1['CART_radiation'][-1], 1e6 * dict1['CART_x'], 1e6 * dict1['CART_y'], aspect='auto', title="last", show=1)




    if True:
        plot_image(numpy.abs(dict1['CART_e_amplitude_sigma'][0])**2  + numpy.abs(dict1['CART_e_amplitude_pi'][0])**2 ,
                   100 * 1e6 * dict1['CART_x'], 100 * 1e6 * dict1['CART_y'], aspect='auto', title="first", show=0)
        plot_image(numpy.abs(dict1['CART_e_amplitude_sigma'][-1])**2 + numpy.abs(dict1['CART_e_amplitude_pi'][-1])**2,
                   100 * 1e6 * dict1['CART_x'], 100 * 1e6 * dict1['CART_y'], aspect='auto', title="last", show=1)

        i_prop0 =  dict1['CART_BACKPROPAGATED_radiation'][0]  # numpy.abs(dict1['CART_BACKPROPAGATED_e_amplitude_sigma'][0])**2  + numpy.abs(dict1['CART_BACKPROPAGATED_e_amplitude_pi'][0])**2
        i_prop1 =  dict1['CART_BACKPROPAGATED_radiation'][-1] # numpy.abs(dict1['CART_BACKPROPAGATED_e_amplitude_sigma'][-1])**2 + numpy.abs(dict1['CART_BACKPROPAGATED_e_amplitude_pi'][-1])**2
        x_um = 1e6 * dict1['CART_BACKPROPAGATED_x']
        y_um = 1e6 * dict1['CART_BACKPROPAGATED_y']
        plot_image( i_prop0, x_um, y_um, aspect='auto', title="PROPAGATED first", show=0)
        plot_image( i_prop1, x_um, y_um, aspect='auto', title="PROPAGATED last", show=1)
        plot(x_um, i_prop0[x_um.size//2,:],
             y_um, i_prop0[:, y_um.size//2])