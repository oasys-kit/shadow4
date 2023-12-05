#
# examples for SourceUndulator to be used in ShadowOui
#
import numpy
from syned.storage_ring.electron_beam import ElectronBeam

from shadow4.sources.undulator.s4_undulator import S4Undulator
from shadow4.sources.undulator.s4_undulator_light_source import S4UndulatorLightSource

if __name__ == "__main__":


    from srxraylib.plot.gol import plot, set_qt
    set_qt()

    do_plots = True


    ebeam = ElectronBeam(energy_in_GeV=6.04,
                 energy_spread = 0.0,
                 current = 0.2,
                 moment_xx=(400e-6)**2,
                 moment_xxp=0.0,
                 moment_xpxp=(10e-6)**2,
                 moment_yy=(10e-6)**2,
                 moment_yyp=0.0,
                 moment_ypyp=(4e-6)**2 )

    und = S4Undulator(
        K_vertical=0.25,  # syned Undulator parameter
        period_length=0.032,  # syned Undulator parameter
        number_of_periods=50,  # syned Undulator parameter
        emin=10490.0,  # Photon energy scan from energy (in eV)
        emax=10510.0,  # Photon energy scan to energy (in eV)
        ng_e=3,  # Photon energy scan number of points
        maxangle=0.015,  # Maximum radiation semiaperture in RADIANS
        ng_t=100,  # Number of points in angle theta
        ng_p=11,  # Number of points in angle phi
        ng_j=20,  # Number of points in electron trajectory (per period) for internal calculation only
        code_undul_phot="internal",  # internal, pysru, srw
        flag_emittance=0,  # when sampling rays: Use emittance (0=No, 1=Yes)
        flag_size=2,  # when sampling rays: 0=point,1=Gaussian,2=FT(Divergences)
        )


    print("gamma: ", ebeam.gamma())
    print("resonance: ", und.resonance_energy(ebeam.gamma()))
    und.set_energy_monochromatic(und.resonance_energy(ebeam.gamma()))
    # sourceundulator._MAXANGLE *= 1.2

    # print(und.info())

    ls = S4UndulatorLightSource(name="", electron_beam=ebeam, magnetic_structure=und,
                                nrays=15000,seed=5655452)
    beam = ls.get_beam()

    print(ls.info())
    #
    # plot
    #
    if do_plots:
        from srxraylib.plot.gol import plot_image, plot_scatter

        radiation,photon_energy, theta,phi = ls.get_radiation_polar()
        plot_image(radiation[0],1e6*theta,phi,aspect='auto',title="intensity",xtitle="theta [urad]",ytitle="phi [rad]")

        radiation_interpolated,photon_energy, vx,vz = ls.get_radiation_interpolated_cartesian()
        plot_image(radiation_interpolated[0],vx,vz,aspect='auto',title="intensity interpolated in cartesian grid",xtitle="vx",ytitle="vy")

        polarization = ls.get_result_polarization()
        plot_image(polarization[0],1e6*theta,phi,aspect='auto',title="polarization",xtitle="theta [urad]",ytitle="phi [rad]")



    print("Beam intensity: ",beam.get_column(23).sum())
    print("Beam intensity s-pol: ",beam.get_column(24).sum())
    print("Beam intensity: p-pol",beam.get_column(25).sum())

    #
    # plot
    #
    if do_plots:
        plot_scatter(1e6*beam.rays[:,0],1e6*beam.rays[:,2],title="real space",xtitle="X [um]",ytitle="Z [um]",show=False)
        plot_scatter(1e6*beam.rays[:,3],1e6*beam.rays[:,5],title="divergence space",xtitle="X [urad]",ytitle="Z [urad]",show=True)

        plot(ls.get_photon_size_distribution()[0]*1e6,
             ls.get_photon_size_distribution()[1],
             title="Photon size distribution",xtitle="R [um]",ytitle="Intensity [a.u.]")

    # check the correct size sampling (values must agree for FLAG_SIZE=1!!!)
    x_photon = beam.rays[:,0]
    z_photon = beam.rays[:,2]
    R = numpy.sqrt(x_photon**2 + z_photon**2)
    print(">> s_phot, Std R", ls.get_result_photon_size_sigma(), numpy.sqrt((R ** 2).sum() / (R.size - 1)))
    print(">> s_phot, Std X", ls.get_result_photon_size_sigma(), numpy.sqrt((x_photon ** 2).sum() / (x_photon.size - 1)))
    print(">> s_phot, Std Z", ls.get_result_photon_size_sigma(), numpy.sqrt((z_photon ** 2).sum() / (z_photon.size - 1)))
