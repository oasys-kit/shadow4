import numpy
import scipy.constants as codata
from shadow4.sources.s4_light_source_base import S4LightSourceBase

def beamline_get_source_beam(beamline):

    n = beamline.get_beamline_elements_number()
    if n == 0:
        light_source = beamline.get_light_source()
        beam0 = light_source.get_beam()
    else:
        beam0 = beamline.get_beamline_element_at(0).get_input_beam()
        if beam0 is None: beam0 = light_source.get_beam()

    return beam0

def beamline_get_last_beam(beamline):

    n = beamline.get_beamline_elements_number()
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


    beam1 = beamline_get_last_beam(beamline)
    photon_energy = beam1.get_column(26, nolost=1)

    if isinstance(beamline.get_light_source(), S4LightSourceBase):
        return "*** ERROR *** To be implemented for geometrical sources"

    is_monochromatic = beamline.get_light_source().get_magnetic_structure().is_monochromatic()
    if spectrum_energy is None or spectrum_flux is None:
        spectrum_energy, spectrum_flux, _ = beamline.get_light_source().calculate_spectrum()

    #
    # SHADOW data
    #
    txt = "\n"

    txt += "# SHADOW SOURCE --------- \n"
    txt += "\n"
    try:
        if is_monochromatic:
            txt += " Source Energy (monochromatic): %.3f eV \n" % (photon_energy[0])

            txt += "\n"
            txt += "# SHADOW BEAMLINE --------- \n"
            txt += " \n"
            txt += " Shadow Intensity (Initial): %f \n" % (beam0.intensity(nolost=1))
            txt += " Shadow Intensity (Final)  : %f \n" % (beam1.intensity(nolost=1))
            txt += " Efficiency: %f %% \n" % (100 * beam1.intensity(nolost=1) / beam0.intensity(nolost=1))

            txt += "\n"
            txt += " Bandwidth considered : 1 eV \n"
            T_eV = beam1.intensity(nolost=1) / beam0.intensity(nolost=1)
            txt += " Flux Transmitivity per eV: %.3g %% \n" % (100 * T_eV)
        else:
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




            # if is_monochromatic: raise Exception("Cannot calculate flux and power for monochromatic beams")

            if beamline.get_light_source().get_magnetic_structure().is_monochromatic(): # monochromatic
                raise NotImplementedError()
            else: # polychromatic
                ticket0 = beam0.histo1(26, nbins=nbins, xrange=[e_min, e_max], nolost=1, ref=23)
                ticket1 = beam1.histo1(26, nbins=nbins, xrange=[e_min, e_max], nolost=1, ref=23)

            e_delta = e_max - e_min

            if e_delta < (4 * ticket1['fwhm']):
                txt += ("\n\n\n*** WARNING *** Source \u0394E (" + str(round(e_delta, 2)) + " eV) should be at least 4 times bigger than the bandwidth (" + str(round(ticket1['fwhm'], 3)) + " eV)\n\n\n")

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

    except:
        txt += "*** Error calculating SHADOW SOURCE ***"
        return (txt)

    #
    # SOURCE spectrum
    #
    txt += "\n"
    txt += "# SOURCE SPECTRUM --------- \n"
    txt += "\n"

    try:

        if is_monochromatic:
            peak_i = spectrum_flux.size // 2
            peak_f = spectrum_flux[peak_i]
            peak_e = spectrum_energy[peak_i]
            txt += "Flux from Source (at %.3f eV): %.3g ph/s/0.1%%bw = %.3g ph/s/eV \n" % (
            (peak_e, peak_f, peak_f / (1e-3 * peak_e)))
        else:
            peak_i = numpy.argmax(spectrum_flux)
            peak_f = spectrum_flux[peak_i]
            peak_e = spectrum_energy[peak_i]
            txt += " Peak Flux from Source (at %.3f eV): %.3g ph/s/0.1%%bw = %.3g ph/s/eV \n" % (peak_e, peak_f, peak_f / (1e-3 * peak_e))
            txt += " Averaged Flux from Source: %.3g ph/s/0.1%%bw = %.3g ph/s/eV \n" % (spectrum_flux.mean(), spectrum_flux.mean() / (1e-3 * spectrum_energy.mean()))
            initial_flux = numpy.trapz(spectrum_flux / (1e-3 * spectrum_energy), spectrum_energy)
            averaged_flux = initial_flux / (spectrum_energy.max() - spectrum_energy.min())
            txt += " Integrated Flux from Source: Total: %.3g ph/s = %.3g ph/s/eV\n" % (initial_flux, averaged_flux)
    except:
        txt += "*** Error calculating SHADOW SPECTRUM ***"
        return (txt)

    #
    # Calculation using an approximated method (an interpolated value of the source flux, supposed almost constants at the source)
    #
    txt += "\n"
    txt += "# FLUX CALCULATION (APPROXIMATED)--------- \n"
    txt += "\n"

    try:
        if is_monochromatic:
            txt += "\n"
            flux_at_sample = peak_f * T_eV
            txt += " ---> Flux at image: %.3g ph/s/eV \n" % (flux_at_sample)

            ticket = beam1.histo2(1, 3, nbins=100, nolost=1, ref=23)

            dx = ticket['fwhm_v'] * 1e3  # mm
            dy = ticket['fwhm_h'] * 1e3  # mm

            txt += " ---> Flux Density  : %.3g ph/s/eV/mm^2 \n" % (flux_at_sample / (dx * dy))

            power_at_sample = flux_at_sample * peak_e * codata.e
            txt += " ---> Power at image: %.3g W/eV\n" % (power_at_sample)
            txt += " ---> Power Density  : %.3g W/eV/mm^2 (over %f x %f um2) \n" % (
            power_at_sample / (dx * dy), 1e3 * dx, 1e3 * dy)
        else:

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
    except:
        txt += "*** Error calculating FLUX CALCULATION (APPROXIMATED) ***"
        return (txt)

    #
    # Calculation using an exact method (via calibrated histograms)
    #
    try:
        if not is_monochromatic:
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
    except:
        txt += "*** Error calculating FLUX CALCULATION (EXACT) ***"
        return (txt)

    return(txt)


if __name__ == "__main__":

    def get_beamline():

        #
        #
        #
        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
        light_source = SourceGeometrical(name='SourceGeometrical', nrays=5000, seed=5676561)
        light_source.set_spatial_type_rectangle(width=0.100000, height=0.200000)
        light_source.set_depth_distribution_off()
        light_source.set_angular_distribution_flat(hdiv1=-0.000000, hdiv2=0.000000, vdiv1=-0.000005, vdiv2=0.000005)
        light_source.set_energy_distribution_singleline(1000.000000, unit='eV')
        light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
        beam = light_source.get_beam()

        # test plot
        from srxraylib.plot.gol import plot_scatter
        rays = beam.get_rays()
        plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')

        from shadow4.beamline.s4_beamline import S4Beamline
        beamline = S4Beamline()
        beamline.set_light_source(light_source)

        return beamline


    beamline = get_beamline()
    a = numpy.loadtxt("/home/srio/Oasys/spectrum.dat.csv")
    print(a.shape)
    # print(flux_summary(beamline, spectrum_energy=a[:,0], spectrum_flux=a[:,1]))
    print(flux_summary(beamline))

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