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

#
#focnew
#
# def focnew_summary(beamline):
#     return focnew(beamline)['text']

def focnew(beamline=None, beam=None, nolost=1, mode=0, center=[0.0,0.0]):
    """
    FocNew tool that finds the best focus looking at the evolution of the beam calculated after the standard deviation
    of some beam variables.

    It locates the waist of the beam by minimizing the variance in real space with respect to the optical axis in the
    X, Z planes and for the circle of least confusion.
    For many purposes, the focus can be defined as the position where the beam spread is minimum. Thus we have to find
    the minimum of the variance.
    If x_i is the coordinate-vector of the i-th ray, and v_i is its direction vector then the variance will be:
    var_x(t) = (1/N) Sum[(x_i + v_i t)^2]. Expanding the sum, and taking the derivative equal to zero we obtain "t", or
    location of the beast focus: t = Sum[x_i v_i] / Sum[v_i^2].

    The user may also choose a center of the distribution other than the optical axis. This utility is very important
    in the study of monochromators and spectrographs, since it allows the user to verify the real location of the
    focal position.


    Parameters
    ----------
    beamline : instance of S4beamline or None, optional
        The beamline instance. Note that either beamline or beam should be defined.
    beam : instance of S4Beam or None
        The beam instance. Note that either beamline or beam should be defined.
    nolost : int, optional
        * 0=uses all rays,
        * 1=uses only good rays (non-lost rays),
        * 2=uses only lost rays.
    mode : int, optional
        A flag to define the center:
        * 0 = center at origin,
        * 1 = Center at barycenter (coordinate mean),
        * 2 = External center.
    center : list or tuple, optional
        The (x,z) coordinates of the center. Used if mode=2.

    Returns
    -------
    dict
        A dictionary (ticket) with the calculation results:
        * ticket['nolost']         # input flag
        * ticket['mode']           # input flag
        * ticket['center_at']      #  text of mode: 'Origin','Baricenter' or 'External'
        * ticket['AX']             # \
        * ticket['AZ']             #  focnew coefficients (to be used by focnew_scan)
        * ticket['AT']             # /
        * ticket['x_waist']        # position of waist X
        * ticket['z_waist']        # position of waist Z
        * ticket['t_waist']        # position of waist T (averaged)
        * ticket['text'] = txt     # a text with focnew info

    Notes
    -----
    One must be careful in using FOCNEW and verify with PLOTXY that the calculations are correct. This is because with
    strongly asymmetric or multimodel distribution the standard deviation is not a good measure of the image “size”.
    Also, if the beam are not Gaussian, the relation between the FWHM and sigma is not defined (It is FWHM=2.355 sigma
    for beams that follow Gaussian distributions).

    """
    if beam is not None:
        beam1 = beam
    else:
        if beamline is not None:
            beam1 = beamline_get_last_beam(beamline)
        else:
            raise Exception("Empty beam. Please define either beamline or beam.")

    AX, AZ, AT = beam1.focnew_coeffs(nolost=nolost, mode=mode, center=center)


    # store versors
    ZBAR = AZ[3]
    VZBAR = AZ[5]
    #
    XBAR = AX[3]
    VXBAR = AX[5]
    #
    TBAR = ZBAR + XBAR
    VTBAR = VZBAR + VXBAR

    # #reset coeffs
    # if mode != 1:
    #     AZ[3] = 0.0
    #     AZ[4] = 0.0
    #     AZ[5] = 0.0
    #
    #     AX[3] = 0.0
    #     AX[4] = 0.0
    #     AX[5] = 0.0
    #
    #     AT[3] = 0.0
    #     AT[4] = 0.0
    #     AT[5] = 0.0

    #get Y coordinate of the three waists

    if numpy.abs(AZ[0]-AZ[5]) > 1e-30:
        TPARZ = (AZ[4] - AZ[1]) / (AZ[0] - AZ[5])
    else:
        TPARZ = 0.0

    if numpy.abs(AX[0]-AX[5]) > 1e-30:
        TPARX = (AX[4] - AX[1]) / (AX[0] - AX[5])
    else:
        TPARX = 0.0

    if numpy.abs(AT[0]-AX[5]) > 1e-30:
        TPART = (AT[4] - AT[1]) / (AT[0] - AT[5])
    else:
        TPART = 0.0

    #prepare text output

    NMODE = ['Origin','Baricenter','External']
    txt = ""
    txt += '-----------------------------------------------------------------------------\n'
    txt += 'Center at : %s\n'%(NMODE[mode])
    if mode == 2: txt += 'X = %f    Z = %f\n'%(center[0],center[1])
    txt += '-----------------------------------------------------------------------------\n'

    SIGX = numpy.sqrt(numpy.abs( AX[0] * TPARX**2 + 2.0 * AX[1] * TPARX + AX[2] - ( AX[3] + 2.0 * AX[4] * TPARX + AX[5] * TPARX**2)))
    SIGZ = numpy.sqrt(numpy.abs( AZ[0] * TPARZ**2 + 2.0 * AZ[1] * TPARZ + AZ[2] - ( AZ[3] + 2.0 * AZ[4] * TPARZ + AZ[5] * TPARZ**2)))
    SIGT = numpy.sqrt(numpy.abs( AT[0] * TPART**2 + 2.0 * AT[1] * TPART + AT[2] - ( AT[3] + 2.0 * AT[4] * TPART + AT[5] * TPART**2)))

    SIGX0 = numpy.sqrt(numpy.abs(AX[2] - AX[3]))
    SIGZ0 = numpy.sqrt(numpy.abs(AZ[2] - AZ[3]))
    SIGT0 = numpy.sqrt(numpy.abs(AT[2] - AT[3]))

    txt += '.............   X AXIS (column 1)    ............\n'
    txt += 'X coefficients :   %g %g %g\n'%(AX[0],AX[1],AX[2])
    txt += 'Center : %g   Average versor : %g\n'%(numpy.sqrt(numpy.abs(XBAR)),numpy.sqrt(numpy.abs(VXBAR)))
    txt += 'Focus along X at       :  %g\n'%(TPARX)
    txt += 'Waist size at best focus (rms)	:  %g\n'%(SIGX)
    txt += 'Waist size at origin                :  %g\n'%(SIGX0)

    txt += '.............   Z AXIS (column 3)    ............\n'
    txt += 'Z coefficients :   %g %g %g\n'%(AZ[0],AZ[1],AZ[2])
    txt += 'Center : %g   Average versor : %g\n'%(numpy.sqrt(numpy.abs(ZBAR)),numpy.sqrt(numpy.abs(VZBAR)))
    txt += 'Focus along Z at       :  %g\n'%(TPARZ)
    txt += 'Waist size at best focus (rms)	:  %g\n'%(SIGZ)
    txt += 'Waist size at origin                :  %g\n'%(SIGZ0)

    txt += '.............  L E A S T  C O N F U S I O N  ...............\n'
    txt += 'XZ coefficients :   %g %g %g\n'%(AT[0],AT[1],AT[2])
    txt += 'Center : %g   Average versor : %g\n'%(numpy.sqrt(numpy.abs(TBAR)),numpy.sqrt(numpy.abs(VTBAR)))
    txt += 'Circle of least confusion :  %g\n'%(TPART)
    txt += 'Waist size at best focus (rms)	:  %g\n'%(SIGT)
    txt += 'Waist size at origin                :  %g\n'%(SIGT0)

    #store all outputs
    ticket = {}
    # copy the inputs
    ticket['nolost'] = nolost
    ticket['mode'] = mode
    ticket['center_at'] = NMODE[mode]
    # coefficients
    ticket['AX'] = AX
    ticket['AZ'] = AZ
    ticket['AT'] = AT
    # position of waists
    ticket['x_waist'] = TPARX
    ticket['z_waist'] = TPARZ
    ticket['t_waist'] = TPART
    # text output
    ticket['text'] = txt

    return ticket



def focnew_scan(A, x):
    """
    Scans the RMS of the beam size along the optical axis using the focnew coefficients.

    Parameters
    ----------
    A : numpy array or list or tuple
        array of 6 coefficients as a result of S4beam.focnew_coeffs().
    x : numpy array
        the abscissas array.

    Returns
    -------
    numpy array
        the array with RMS values.
    """
    x1 = numpy.array(x)
    y = numpy.sqrt(numpy.abs( A[0] * x1**2 + 2.0 * A[1] * x1 + A[2] - (A[3] + 2.0 * A[4] * x1 + A[5] * x1**2)))
    return y


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

        return beamline, beam


    beamline, beam = get_beamline()
    a = numpy.loadtxt("/home/srio/Oasys/spectrum.dat.csv")
    print(a.shape)
    # print(flux_summary(beamline, spectrum_energy=a[:,0], spectrum_flux=a[:,1]))
    print(flux_summary(beamline))

    # print(focnew_summary(beamline))


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

    print(beamline.get_beamline_elements_number())
    # print(focnew(beam=beam))
    # print(focnew(beam=beamline_get_last_beam(beamline)))
    print(focnew(beamline=beamline))