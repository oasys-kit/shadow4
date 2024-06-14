import numpy
import matplotlib.pylab as plt
from srxraylib.plot.gol import plot_scatter, plot, plot_show


# electron beam
from shadow4.sources.s4_electron_beam import S4ElectronBeam
electron_beam = S4ElectronBeam(energy_in_GeV=6.04,energy_spread=0.001,current=0.2)
electron_beam.set_sigmas_all(sigma_x=7.8e-05,sigma_y=4.87179e-05,sigma_xp=3.6e-05,sigma_yp=1.05556e-06)

#magnetic structure
from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet

photon_energy = 8000.0

R = numpy.abs(S4BendingMagnet.calculate_magnetic_radius(-0.8, electron_beam.energy()))
Dtheta = 1e-3

print(">>>> radius: ", R)
source = S4BendingMagnet(
                 radius=R, # from syned BM, can be obtained as numpy.abs(S4BendingMagnet.calculate_magnetic_radius(0.8, electron_beam.energy()))
                 magnetic_field=-0.8, # from syned BM
                 length=Dtheta * R, # from syned BM = abs(BM divergence * magnetic_field)
                 emin=photon_energy,     # Photon energy scan from energy (in eV)
                 emax=photon_energy,     # Photon energy scan to energy (in eV)
                 ng_e=100,               # Photon energy scan number of points
                 flag_emittance=0, # when sampling rays: Use emittance (0=No, 1=Yes)
                 )
#light source
from shadow4.sources.bending_magnet.s4_bending_magnet_light_source import S4BendingMagnetLightSource
light_source = S4BendingMagnetLightSource(name='BendingMagnet',
                                          electron_beam=electron_beam,
                                          magnetic_structure=source,
                                          nrays=25000,
                                          seed=5676561)


beam = light_source.get_beam(sample_emission_cone_in_horizontal=1)
rays = beam.get_rays()
# test plot
# plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')


#
# Stokes P3 histogram
#
# tktP3 = beam.histo1(6, ref=33)
# plot(
#     tktP3["bin_path"], tktP3["histogram_path"], )


#
# x(y) scatter phase space
#
if False:
    theta = numpy.linspace(-0.5 * Dtheta, 0.5 * Dtheta, 100)
    x = R * theta**2 / 2
    y = R * theta
    xp = theta

    fig, axScatter, axHistx, axHisty = plot_scatter(beam.get_column(2), beam.get_column(1), title="trajectory  (X,Y)", show=0)
    axScatter.plot(y, x, color='r')

    fig, axScatter, axHistx, axHisty = plot_scatter(beam.get_column(1), beam.get_column(4), title="phase space (X,X')", show=0)
    axScatter.plot(x, xp, color='r')

    plot_show()


#
# psi histogram
#
tktI = beam.histo1(6, ref=23)
tktS = beam.histo1(6, ref=24)
tktP = beam.histo1(6, ref=25)

plot(
    tktI["bin_path"], tktI["histogram_path"],
    tktS["bin_path"], tktS["histogram_path"],
    tktP["bin_path"], tktP["histogram_path"],
    xtitle=r"$\psi$ [rad]", ytitle="Intensity [arbitrary units] at E=%d eV" % (photon_energy),
)


# normalization to the area
I0 = numpy.trapz(tktI["histogram"], tktI["bin_center"] * 1e3)
T0 = numpy.trapz(light_source.tot, light_source.angle_array_mrad)

plot(
    light_source.angle_array_mrad,   light_source.tot,
    light_source.angle_array_mrad,     light_source.s,
    light_source.angle_array_mrad,     light_source.p,
    tktI["bin_path"] * 1e3, tktI["histogram_path"] * T0 / I0,
    tktS["bin_path"] * 1e3, tktS["histogram_path"] * T0 / I0,
    tktP["bin_path"] * 1e3, tktP["histogram_path"] * T0 / I0,
    xtitle=r"$\psi$ [mrad]", ytitle="Flux [photons/s/eV] at E=%d eV" % (photon_energy),
    linestyle=['--','--','--', None, None, None],
    color=['b','g','r', 'b','g','r'],
    # legend=['Flux', r'$\sigma$-polarized', r'$\pi$-polarized',
    #         'Histogram', r'$\sigma$-polarized', r'$\pi$-polarized'],

    show=0)

plt.savefig('esrf1_histo_psi.pdf')
plt.show()
