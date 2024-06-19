
# electron beam
from shadow4.sources.s4_electron_beam import S4ElectronBeam
electron_beam = S4ElectronBeam(energy_in_GeV=2,energy_spread=0,current=0.4)
electron_beam.set_sigmas_all(sigma_x=0,sigma_y=0,sigma_xp=0,sigma_yp=0)

#magnetic structure
from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet
import numpy

Dtheta = 10e-3
magnetic_field = -0.87
radius = numpy.abs(S4BendingMagnet.calculate_magnetic_radius(-0.87, electron_beam.energy())) # radius = 7.668140119497748
length = Dtheta * radius # 0.07668140119497747
sample_emission_cone_in_horizontal = 1

source = S4BendingMagnet(
                 radius=radius, # from syned BM, can be obtained as numpy.abs(S4BendingMagnet.calculate_magnetic_radius(-0.87, electron_beam.energy()))
                 magnetic_field=magnetic_field, # from syned BM
                 length=length, # from syned BM = abs(BM divergence * magnetic_radius)
                 emin=282.0,     # Photon energy scan from energy (in eV)
                 emax=282.0,     # Photon energy scan to energy (in eV)
                 ng_e=100,     # Photon energy scan number of points
                 flag_emittance=0, # when sampling rays: Use emittance (0=No, 1=Yes)
                 epsi_dx=0.0,  # position of X waist [m]
                 epsi_dz=0.0 , # position of Z waist [m]
                 )



#light source
from shadow4.sources.bending_magnet.s4_bending_magnet_light_source import S4BendingMagnetLightSource
light_source = S4BendingMagnetLightSource(name='BendingMagnet', electron_beam=electron_beam, magnetic_structure=source, nrays=50000, seed=5676561)
beam = light_source.get_beam(sample_emission_cone_in_horizontal=sample_emission_cone_in_horizontal)

# # test plot
# from srxraylib.plot.gol import plot_scatter
# rays = beam.get_rays()
# plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 3], title="(X,X')")

#
# x(y) scatter phase space
#
if True:
    import matplotlib.pylab as plt
    from srxraylib.plot.gol import plot_scatter, plot, plot_show

    theta = numpy.linspace(-0.5 * Dtheta, 0.5 * Dtheta, 100)
    x = radius * theta**2 / 2
    y = radius * theta
    xp = theta
    # plot(x, xp, color='r', show=0)

    fig, axScatter, axHistx, axHisty = plot_scatter(beam.get_column(2), beam.get_column(1), title="trajectory  (X,Y)", show=0)
    axScatter.plot(y, x, color='r')

    fig, axScatter, axHistx, axHisty = plot_scatter(1e6 * beam.get_column(1), 1e6 * beam.get_column(4), title="phase space (X,X')",
            xtitle=r"X [$\mu$m]", ytitle=r"X' [$\mu$rad]", show=0)

    if sample_emission_cone_in_horizontal:
        axScatter.plot(1e6 * x, 1e6 * xp, color='r')
        plt.savefig("alsu_phasespaceH.png")
    else:
        plt.savefig("alsu_phasespaceHnoemissioncone.png")

    fig, axScatter, axHistx, axHisty = plot_scatter(1e6 * beam.get_column(4), 1e6 * beam.get_column(6), title="phase space (Z,Z')",
            xtitle=r"Z [$\mu$m]", ytitle=r"Z' [$\mu$rad]", show=0)

    if sample_emission_cone_in_horizontal:
        plt.savefig("alsu_phasespaceV.png")

    plot_show()