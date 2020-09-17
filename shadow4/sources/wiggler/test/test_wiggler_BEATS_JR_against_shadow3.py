#
# This is a comparison shadow3-shadow4 for a short-wiggler using the ALBA magnetic field map in
# the SESAME ring.
#
# The goal is to check whether the vertical divergences match.
#
# in shadow4 the sampling of the vertical divergences are calculated using the Ephoton/Ecritical_max,
# but it turns out that the use of Ecritical_max is overestimated in some cases, in particular this
# ALBA magnetic map.
#


#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy
from srxraylib.sources import srfunc


# shadow4
from syned.storage_ring.electron_beam import ElectronBeam
from shadow4.sources.wiggler.s4_wiggler import S4Wiggler
from shadow4.sources.wiggler.s4_wiggler_light_source import S4WigglerLightSource
from shadow4.beam.beam import Beam
from shadow4.compatibility.beam3 import Beam3


def none_conv(result):
    """Just in case of a None result, to prevent script crashing"""
    if result != None:
        return result
    else:
        new_result = 0.0
        return new_result


def calculate_spectrum(emin, emax, npoints=500, do_plot=False):
    #
    # script to run the wiggler preprocessor (created by ShadowOui:Wiggler)
    #

    (traj, pars) = srfunc.wiggler_trajectory(
        b_from=1,
        inData="Bz_Alba_rev3.dat",
        nPer=1,
        nTrajPoints=501,
        ener_gev=2.5,
        per=0.755,
        kValue=211.0,
        trajFile="tmp.traj",
        shift_x_flag=4,
        shift_x_value=0.0,
        shift_betax_flag=4,
        shift_betax_value=0.0)

    #
    # calculate cdf and write file for Shadow/Source
    #

    srfunc.wiggler_cdf(traj,
                       enerMin=emin,
                       enerMax=emax,
                       enerPoints=1001,
                       outFile=b'xshwig.sha',
                       elliptical=False)

    calculate_spectrum = True

    if calculate_spectrum:
        e, f, w = srfunc.wiggler_spectrum(traj,
                                          enerMin=emin,
                                          enerMax=emax,
                                          nPoints=npoints,
                                          electronCurrent=400 * 1e-3,
                                          outFile="spectrum.dat",
                                          elliptical=False)

        if do_plot:
            from srxraylib.plot.gol import plot
            plot(e, f, xlog=False, ylog=False, show=False,
                 xtitle="Photon energy [eV]", ytitle="Flux [Photons/s/0.1%bw]", title="Flux")
            plot(e, w, xlog=False, ylog=False, show=False,
                 xtitle="Photon energy [eV]", ytitle="Spectral Power [W/eV]", title="Spectral Power")
    #
    # end script
    #
    return e, w


def run_shadow4(photon_energy, n_rays=5e5, emittance=True):

    use_emittances = emittance


    if use_emittances:
        # in mm
        EPSI_X = 2.574e-05
        EPSI_Z = 2.574e-07
        SIGMAX = 0.8208
        SIGMAZ = 0.0142

        sigmax = SIGMAX * 1e-3
        sigmaz = SIGMAZ * 1e-3
        sigmaxprime = EPSI_X * 1e-3 / sigmax
        sigmazprime = EPSI_Z * 1e-3 / sigmaz
    else:
        sigmax = 0.0
        sigmaz = 0.0
        sigmaxprime = 0.0
        sigmazprime = 0.0

    e_min = photon_energy
    e_max = photon_energy
    NRAYS = n_rays


    nTrajPoints = 501
    shift_x_flag = 4
    shift_x_value = 0.0
    shift_betax_flag = 4
    shift_betax_value = 0.0


    #
    # syned
    #

    syned_electron_beam = ElectronBeam(energy_in_GeV=2.5,current=0.4,
                                       moment_xx=(sigmax)**2,
                                       moment_xpxp=(sigmaxprime)**2,
                                       moment_yy=(sigmaz)**2,
                                       moment_ypyp=(sigmazprime)**2,
                                       )


    # wiggler: B from file
    # syned_wiggler = MagneticStructure1DField.initialize_from_file("Bz_Alba_rev3.dat")


    if e_min == e_max:
        ng_e = 1
    else:
        ng_e = 10

    sourcewiggler = S4Wiggler(
                    magnetic_field_periodic=0,
                    file_with_magnetic_field="Bz_Alba_rev3.dat",
                    flag_emittance=use_emittances,
                    emin=e_min,
                    emax=e_max,
                    ng_e=ng_e,
                    ng_j=nTrajPoints)


    # sourcewiggler.set_electron_initial_conditions_by_label(position_label="value_at_zero",
    #                                                        velocity_label="value_at_zero")

    sourcewiggler.set_electron_initial_conditions(shift_x_flag=shift_x_flag,
                                                  shift_x_value=shift_x_value,
                                                  shift_betax_flag=shift_betax_flag,
                                                  shift_betax_value=shift_betax_value)

    ls = S4WigglerLightSource(name="", electron_beam=syned_electron_beam, wiggler_magnetic_structure=sourcewiggler)
    print(ls.info())

    beam = ls.get_beam(NRAYS=NRAYS)

    return beam


def run_shadow3(photon_energy, n_rays=5e5, emittance=True):
    #
    # script to run the wiggler preprocessor (created by ShadowOui:Wiggler)
    #

    (traj, pars) = srfunc.wiggler_trajectory(
        b_from=1,
        inData="Bz_Alba_rev3.dat",
        nPer=1,
        nTrajPoints=501,
        ener_gev=2.5,
        per=0.755,
        kValue=211.0,
        trajFile="tmp.traj",
        shift_x_flag=4,
        shift_x_value=0.0,
        shift_betax_flag=4,
        shift_betax_value=0.0)

    #
    # calculate cdf and write file for Shadow/Source
    #

    srfunc.wiggler_cdf(traj,
                       enerMin=photon_energy,  # 5000.0,
                       enerMax=photon_energy + 1.0,
                       enerPoints=1001,
                       outFile=b'xshwig.sha',
                       elliptical=False)

    #
    # end script
    #

    #
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.BENER = 2.5
    oe0.CONV_FACT = 1000.0
    # oe0.EPSI_X = 2.574e-05
    # oe0.EPSI_Z = 2.574e-07
    if emittance:
        oe0.EPSI_X = 2.574e-05
        oe0.EPSI_Z = 2.574e-07
    oe0.FDISTR = 0
    oe0.FILE_TRAJ = b'xshwig.sha'
    oe0.FSOUR = 0
    oe0.FSOURCE_DEPTH = 0
    oe0.F_COLOR = 0
    oe0.F_PHOT = 0
    oe0.F_WIGGLER = 1
    oe0.HDIV1 = 1.0
    oe0.HDIV2 = 1.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 420024
    oe0.NCOL = 0
    oe0.NPOINT = n_rays
    oe0.N_COLOR = 0
    oe0.POL_DEG = 0.0
    # oe0.PH1 = 40000.0
    # oe0.PH2 = 40001.0
    oe0.SIGMAX = 0.8208
    oe0.SIGMAZ = 0.014
    if emittance:
        oe0.SIGMAX = 0.8208
        oe0.SIGMAZ = 0.0142
    else:
        oe0.SIGMAX = 0.0
        oe0.SIGMAZ = 0.0

    oe0.SIGMAY = 0.0
    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0

    oe1.DUMMY = 0.1
    oe1.FWRITE = 3
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.N_SCREEN = 1
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 30000.0

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")


    return beam, oe0


if __name__ == "__main__":
    from Shadow.ShadowTools import plotxy


    from srxraylib.plot.gol import set_qt
    set_qt()

    n_rays = 15000

    # photon_energies = numpy.linspace(100, 100000, 1665)
    photon_energies = numpy.linspace(100, 20000, 200)

    fwhm_v_3 = numpy.zeros_like(photon_energies)
    fwhm_v_4 = numpy.zeros_like(photon_energies)

    for i,photon_energy in enumerate(photon_energies):
        #
        # shadow3
        #
        beam3, oe0 = run_shadow3(photon_energy, n_rays=n_rays, emittance=False)


        #
        # shadow4
        #
        beam4 = run_shadow4(photon_energy, n_rays=n_rays, emittance=False)



        tkt3 = beam3.histo2(4, 6, nolost=1, ref=23, nbins=201)
        beam4_3 = Beam3.initialize_from_shadow4_beam(beam4)
        tkt4 = beam4_3.histo2(4, 6, nolost=1, ref=23, nbins=201)

        if i==0:
            plotxy(beam3, 4, 6, nolost=1, ref=23, nbins=201, title="shadow3")
            plotxy(beam4_3, 4, 6, nolost=1, ref=23, nbins=201, title="shadow4")


        print(">>>>>", photon_energy, tkt3["fwhm_v"], tkt4["fwhm_v"], tkt4["fwhm_v"].shape)
        fwhm_v_3[i] = tkt3["fwhm_v"]
        fwhm_v_4[i] = tkt4["fwhm_v"]


    f = open("tmp.dat",'w')
    for i in range(photon_energies.size):
        f.write("%f  %f  %f %g \n" % (photon_energies[i], fwhm_v_3[i], fwhm_v_4[i], (fwhm_v_3[i] - fwhm_v_4[i])))
    f.close()
    print("File tmp.dat written to disk.") # moved to ./test_wiggler_BEATS_JR_against_shadow3.dat