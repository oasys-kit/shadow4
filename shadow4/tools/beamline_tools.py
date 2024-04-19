import numpy
import scipy.constants as codata


def beamline_get_source_beam(beamline):

    n = beamline.get_beamline_elements_number()
    print(">>>>>> n: ", n)
    if n == 0:
        light_source = beamline.get_light_source()
        beam0 = light_source.get_beam()
    else:
        beam0 = beamline.get_beamline_element_at(0).get_input_beam()
        if beam0 is None: beam0 = light_source.get_beam()

    return beam0

def beamline_get_last_beam(beamline):

    n = beamline.get_beamline_elements_number()
    print(">>>>>> n: ", n)
    if n == 0:
        light_source = beamline.get_light_source()
        beam1 = light_source.get_beam()
    else:
        last_bel = beamline.get_beamline_element_at(n-1)
        beam0 = last_bel.get_input_beam()
        if beam0 is None:
            beam1, _ = beamline.run_beamline()
        else:
            beam1, _ = last_bel.trace_beam()

    return beam1

def flux_summary(beamline, spectrum_energy=None, spectrum_flux=None, e_min=None, e_max=None, nbins=102):

    # beamline = get_beamline()
    beam0 = beamline_get_source_beam(beamline)
    if spectrum_energy is None or spectrum_flux is None:
        spectrum_flux, spectrum_energy = beamline.get_light_source().get_flux()

    beam1 = beamline_get_last_beam(beamline)
    photon_energy = beam1.get_column(26, nolost=1)

    #
    # SHADOW data
    #
    txt = "\n"

    txt += "# SHADOW SOURCE --------- \n"
    txt += "\n"
    txt += " Source Central Energy: %.3f eV \n" % (0.5 * (photon_energy.max() - photon_energy.min()) + photon_energy.min())
    txt += " Source Energy Range  : %.3f to %.3f eV \n" % (photon_energy.min(), photon_energy.max())
    txt += " Source ΔE: %f eV \n" % (round(photon_energy.max() - photon_energy.min(), 2))

    #
    #
    #
    if e_min is None: e_min = photon_energy.min()
    if e_max is None: e_max = photon_energy.max()

    if e_min < photon_energy.min() : e_min = photon_energy.min()
    if e_max > photon_energy.max(): e_max = photon_energy.max()


    is_monochromatic = beamline.get_light_source().get_magnetic_structure().is_monochromatic()

    if is_monochromatic: raise Exception("Cannot calculate flux and power for monochromatic beams")

    if beamline.get_light_source().get_magnetic_structure().is_monochromatic(): # monochromatic
        raise NotImplementedError()
    else: # polychromatic
        ticket0 = beam0.histo1(26, nbins=nbins, xrange=[e_min, e_max], nolost=1, ref=23)
        ticket1 = beam1.histo1(26, nbins=nbins, xrange=[e_min, e_max], nolost=1, ref=23)

    e_delta = e_max - e_min

    if e_delta < (4 * ticket1['fwhm']):
        raise ValueError("Source \u0394E (" + str(round(e_delta, 2)) + " eV) should be at least 4 times bigger than the bandwidth (" + str(round(ticket1['fwhm'], 3)) + " eV)")


    txt += "\n"
    txt += "# SHADOW BEAMLINE --------- \n"
    txt += " \n"
    txt += " Shadow Intensity (Initial): %f (from histogram: %f)\n" % (beam0.intensity(nolost=1), numpy.array(ticket0['histogram']).sum())
    txt += " Shadow Intensity (Final)  : %f (from histogram: %f)\n" % (beam1.intensity(nolost=1), numpy.array(ticket1['histogram']).sum())
    txt += " Efficiency: %f %% \n" % (100 * beam1.intensity(nolost=1) / beam0.intensity(nolost=1))

    txt += "\n"
    txt += " Bandwidth (at the Image Plane): %.3f eV \n" % (ticket1['fwhm'])
    T_eV = beam1.intensity(nolost=1) / (ticket1['fwhm']) / \
           (beam0.intensity(nolost=1) / (photon_energy.max() - photon_energy.min()))
    txt += " Flux Transmitivity per eV: %.3g %% \n" % (100 * T_eV)

    #
    # SOURCE soectrum
    #
    txt += "\n"
    txt += "# SOURCE SPECTRUM --------- \n"
    txt += "\n"

    peak_i = numpy.argmax(spectrum_flux)
    peak_f = spectrum_flux[peak_i]
    peak_e = spectrum_energy[peak_i]
    txt += " Peak Flux from Source (at %.3f eV): %.3g ph/s/0.1%%bw = %.3g ph/s/eV \n" % (peak_e, peak_f, peak_f / (1e3 * peak_e))
    txt += " Averaged Flux from Source: %.3g ph/s/0.1%%bw = %.3g ph/s/eV \n" % (spectrum_flux.mean(), spectrum_flux.mean() / (1e3 * spectrum_energy.mean()))
    initial_flux = numpy.trapz(spectrum_flux / (1e-3 * spectrum_energy), spectrum_energy)
    averaged_flux = initial_flux / (spectrum_energy.max() - spectrum_energy.min())
    txt += " Integrated Flux from Source: Total: %.3g ph/s = %.3g ph/s/eV\n" % (initial_flux, averaged_flux)

    #
    # Calculation using an approximated method (an interpolated value of the source flux, supposed almost constants at the source)
    #
    txt += "\n"
    txt += "# FLUX CALCULATION (APPROXIMATED)--------- \n"
    txt += "\n"

    interpolated_flux = numpy.interp(ticket0["bin_center"],
                                     spectrum_energy,
                                     spectrum_flux,
                                     left=spectrum_flux[0],
                                     right=spectrum_flux[-1])

    interpolated_flux_per_ev = interpolated_flux / (1e-3 * ticket0["bin_center"])

    txt += " Initial Flux from Source (interpolated at E=%f eV]): %g  ph/s/0.1%%bw = %g ph/s/eV" % (ticket0["bin_center"][0],
                                                                            interpolated_flux[0],
                                                                            interpolated_flux_per_ev[0])

    txt += "\n"
    flux_at_sample = interpolated_flux_per_ev * T_eV * ticket1['fwhm']
    txt += " ---> Integrated Flux at image: %.3g ph/s \n" % ( flux_at_sample[0] )

    ticket = beam1.histo2(1, 3, nbins=100, nolost=1, ref=23)

    dx = ticket['fwhm_v'] * 1e3 # mm
    dy = ticket['fwhm_h'] * 1e3 # mm

    txt += " ---> Flux Density  : %.3g ph/s/mm^2 \n" % (flux_at_sample[0] / (dx * dy))

    power_at_sample = flux_at_sample[0] * ticket0["bin_center"][0] * codata.e
    txt += " ---> Integrated Power at image: %.3g W\n" % (power_at_sample)
    txt += " ---> Power Density  : %.3g W/mm^2 (over %f x %f um2) \n" % (power_at_sample / (dx * dy), 1e3 * dx, 1e3 * dy)


    #
    # Calculation using an exact method (via calibrated histograms)
    #
    txt += "\n"
    txt += "# FLUX CALCULATION (EXACT)--------- \n"
    txt += "\n"

    txt += " Initial Flux from Source (integrated over histogram): %g ph/s" % (
                                        numpy.trapz(interpolated_flux_per_ev, ticket0['bin_center']))

    txt += "\n"
    flux_at_sample = interpolated_flux_per_ev * ticket1['histogram'] / ticket0['histogram']
    flux_at_sample_integrated = numpy.trapz(flux_at_sample, ticket1['bin_center'])

    txt += " ---> Integrated Flux at image: %.3g ph/s \n" % (flux_at_sample_integrated)
    txt += " ---> Flux Density  : %.3g ph/s/mm^2  (over %f x %f um2) \n" % (flux_at_sample_integrated / (dx * dy), 1e3 * dx, 1e3 * dy)

    power_at_sample = flux_at_sample * ticket0["bin_center"] * codata.e
    power_at_sample_integrated = numpy.trapz(power_at_sample, ticket1['bin_center'])
    step = ticket1['bin_center'][1] - ticket1['bin_center'][0]
    txt += " ---> Integrated Power at image: %.3g W = %g\n" % (power_at_sample_integrated, power_at_sample.sum() * step)
    txt += " ---> Power Density  : %.3g W/mm^2 (over %f x %f um2) \n" % (power_at_sample_integrated / (dx * dy), 1e3 * dx, 1e3 * dy)


    return(txt)


if __name__ == "__main__":

    def get_beamline():
        from shadow4.beamline.s4_beamline import S4Beamline

        beamline = S4Beamline()

        # electron beam
        from shadow4.sources.s4_electron_beam import S4ElectronBeam
        electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
        electron_beam.set_sigmas_all(sigma_x=3.01836e-05, sigma_y=4.36821e-06, sigma_xp=3.63641e-06,
                                     sigma_yp=1.37498e-06)

        # magnetic structure
        from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian
        source = S4UndulatorGaussian(
            period_length=0.042,  # syned Undulator parameter (length in m)
            number_of_periods=38.571,  # syned Undulator parameter
            photon_energy=5000.0,  # Photon energy (in eV)
            delta_e=10.0,  # Photon energy width (in eV)
            ng_e=100,  # Photon energy scan number of points
            flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
            flag_energy_spread=0,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
            harmonic_number=1,  # harmonic number
            flag_autoset_flux_central_cone=1,  # value to set the flux peak
            flux_central_cone=681709040139326.4,  # value to set the flux peak
        )

        # light source
        from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource
        light_source = S4UndulatorGaussianLightSource(name='GaussianUndulator', electron_beam=electron_beam,
                                                      magnetic_structure=source, nrays=15000, seed=5676561)
        beam = light_source.get_beam()

        beamline.set_light_source(light_source)

        # optical element number XX
        from syned.beamline.shape import Rectangle
        boundary_shape = Rectangle(x_left=-0.001, x_right=0.001, y_bottom=-0.001, y_top=0.001)

        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
        optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator', boundary_shape=boundary_shape,
                                   i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                                   i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=27.2, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
        beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX
        boundary_shape = None

        from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror
        optical_element = S4PlaneMirror(name='Plane Mirror', boundary_shape=boundary_shape,
                                        f_reflec=1, f_refl=5, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Ni', coating_density=8.902, coating_roughness=0)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=2.7, q=0, angle_radial=1.563796327, angle_azimuthal=1.570796327,
                                         angle_radial_out=1.563796327)
        movements = None
        from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement
        beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

        beam, mirr = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX
        boundary_shape = None

        from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror
        optical_element = S4PlaneMirror(name='Plane Mirror', boundary_shape=boundary_shape,
                                        f_reflec=1, f_refl=5, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Ni', coating_density=8.902, coating_roughness=0)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=0.825, q=0, angle_radial=1.563796327, angle_azimuthal=3.141592654,
                                         angle_radial_out=1.563796327)
        movements = None
        from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirrorElement
        beamline_element = S4PlaneMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                movements=movements, input_beam=beam)

        beam, mirr = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX

        from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
        optical_element = S4Empty(name='Empty Element')

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=4.71238898,
                                         angle_radial_out=3.141592654)
        from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
        beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

        beam, mirr = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX
        from syned.beamline.shape import Rectangle
        boundary_shape = Rectangle(x_left=-0.0015, x_right=0.0015, y_bottom=-0.0015, y_top=0.0015)

        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
        optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator', boundary_shape=boundary_shape,
                                   i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                                   i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=5.475, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
        beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX
        from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
        optical_element = S4PlaneCrystal(name='Plane Crystal',
                                         boundary_shape=None, material='Si',
                                         miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                         f_bragg_a=False, asymmetry_angle=0.0,
                                         is_thick=1, thickness=0.001,
                                         f_central=1, f_phot_cent=0, phot_cent=5000.0,
                                         file_refl='bragg.dat',
                                         f_ext=0,
                                         material_constants_library_flag=0,
                                         # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                         )
        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=1.9, q=0, angle_radial=1.164204344, angle_azimuthal=0,
                                         angle_radial_out=1.164204344)
        movements = None
        from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
        beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                                 movements=movements, input_beam=beam)

        beam, mirr = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX
        from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystal
        optical_element = S4PlaneCrystal(name='Plane Crystal',
                                         boundary_shape=None, material='Si',
                                         miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                         f_bragg_a=False, asymmetry_angle=0.0,
                                         is_thick=1, thickness=0.001,
                                         f_central=1, f_phot_cent=0, phot_cent=5000.0,
                                         file_refl='bragg.dat',
                                         f_ext=0,
                                         material_constants_library_flag=0,
                                         # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                         )
        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=0.012, q=0, angle_radial=1.164204344, angle_azimuthal=3.141592654,
                                         angle_radial_out=1.164204344)
        movements = None
        from shadow4.beamline.optical_elements.crystals.s4_plane_crystal import S4PlaneCrystalElement
        beamline_element = S4PlaneCrystalElement(optical_element=optical_element, coordinates=coordinates,
                                                 movements=movements, input_beam=beam)

        beam, mirr = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX
        from syned.beamline.shape import Rectangle
        boundary_shape = Rectangle(x_left=-0.0005, x_right=0.0005, y_bottom=-0.0005, y_top=0.0005)

        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
        optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator', boundary_shape=boundary_shape,
                                   i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                                   i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=12.888, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
        beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX
        boundary_shape = None

        from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
        optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                            surface_calculation=0,
                                            min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                            p_focus=51.500000, q_focus=0.150000, grazing_angle=0.006000,
                                            is_cylinder=1, cylinder_direction=0, convexity=1,
                                            f_reflec=1, f_refl=5, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                            coating_material='Ni', coating_density=8.902, coating_roughness=0)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=0.5, q=0.0575, angle_radial=1.564796327, angle_azimuthal=0,
                                         angle_radial_out=1.564796327)
        movements = None
        from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
        beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                    movements=movements, input_beam=beam)

        beam, mirr = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX
        boundary_shape = None

        from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror
        optical_element = S4EllipsoidMirror(name='Ellipsoid Mirror', boundary_shape=boundary_shape,
                                            surface_calculation=0,
                                            min_axis=2.000000, maj_axis=2.000000, pole_to_focus=1.000000,
                                            p_focus=51.590000, q_focus=0.060000, grazing_angle=0.006000,
                                            is_cylinder=1, cylinder_direction=0, convexity=1,
                                            f_reflec=1, f_refl=5, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                            coating_material='Ni', coating_density=8.902, coating_roughness=0)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=0.0325, q=0.06, angle_radial=1.564796327, angle_azimuthal=1.570796327,
                                         angle_radial_out=1.564796327)
        movements = None
        from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirrorElement
        beamline_element = S4EllipsoidMirrorElement(optical_element=optical_element, coordinates=coordinates,
                                                    movements=movements, input_beam=beam)

        beam, mirr = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX

        from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4Empty
        optical_element = S4Empty(name='Empty Element')

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=4.71238898,
                                         angle_radial_out=3.141592654)
        from shadow4.beamline.optical_elements.ideal_elements.s4_empty import S4EmptyElement
        beamline_element = S4EmptyElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

        beam, mirr = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # test plot
        if 0:
            from srxraylib.plot.gol import plot_scatter
            plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
                         title='(Intensity,Photon Energy)', plot_histograms=0)
            plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1),
                         title='(X,Z) in microns')


        return beamline


    beamline = get_beamline()
    a = numpy.loadtxt("/home/srio/Oasys/spectrum.dat.csv")
    print(a.shape)
    print(flux_summary(beamline, spectrum_energy=a[:,0], spectrum_flux=a[:,1]))

    if 0:
        from srxraylib.plot.gol import plot_scatter, plot

        # plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
        #              title='(Intensity,Photon Energy)', plot_histograms=0)
        plot_scatter(1e6 * beam0.get_column(1, nolost=1), 1e6 * beam0.get_column(3, nolost=1), title='SOURCE RESTORED (X,Z) in microns')

        plot(spectrum_energy, spectrum_flux)

    if 0:
        from srxraylib.plot.gol import plot_scatter, plot

        plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
                     title='(Intensity,Photon Energy)', plot_histograms=0)
        plot_scatter(1e6 * beam1.get_column(1, nolost=1), 1e6 * beam1.get_column(3, nolost=1), title='LAST RESTORED (X,Z) in microns')

        plot(spectrum_energy, spectrum_flux)