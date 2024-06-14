import numpy
from srxraylib.plot.gol import plot_scatter, plot
import matplotlib.pylab as plt

nbins = 101
# electron beam
from shadow4.sources.s4_electron_beam import S4ElectronBeam
electron_beam = S4ElectronBeam(energy_in_GeV=6.04,energy_spread=0.001,current=0.2)
electron_beam.set_sigmas_all(sigma_x=7.8e-05,sigma_y=4.87179e-05,sigma_xp=3.6e-05,sigma_yp=1.05556e-06)

#magnetic structure
from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet
source = S4BendingMagnet(
                 radius=25.18408918746048, # from syned BM, can be obtained as S4BendingMagnet.calculate_magnetic_radius(0.8, electron_beam.energy())
                 magnetic_field=-0.8, # from syned BM
                 length=0.02518408918746048, # from syned BM = abs(BM divergence * magnetic_field)
                 emin=100.0,        # Photon energy scan from energy (in eV)
                 emax=100000.0,     # Photon energy scan to energy (in eV)
                 ng_e=nbins,        # Photon energy scan number of points
                 flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
                 )


#light source
from shadow4.sources.bending_magnet.s4_bending_magnet_light_source import S4BendingMagnetLightSource
light_source = S4BendingMagnetLightSource(name='BendingMagnet',
                                          electron_beam=electron_beam,
                                          magnetic_structure=source,
                                          nrays=25000,
                                          seed=5676561)

beam = light_source.get_beam()
rays = beam.get_rays()
# test plot
# plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')


tktE = beam.histo1(-11, ref=23, nbins=nbins)
plot(
    tktE["bin_path"], tktE["histogram_path"], )

e, f, w = light_source.calculate_spectrum(shift_half_interval=1)
ff = f / (e * 1e-3)

n_total = numpy.trapz(ff, e)
h_total = numpy.trapz(tktE["histogram"], tktE["bin_center"])
plot(
     e, ff / n_total,
     tktE["bin_path"], tktE["histogram_path"] / h_total,
     xlog=0, ylog=0, marker=['+','+'], show=0)

#
# calibrated histogram
#

e, f, w = light_source.calculate_spectrum(tktE["bin_center"])
ff = f / (e * 1e-3) # in ph/eV

n_total = numpy.trapz(ff, e)
h_total = numpy.trapz(tktE["histogram"], tktE["bin_center"])

# lin
fig, ax = plot(
     e, ff,
     tktE["bin_path"], tktE["histogram_path"] * n_total / h_total,
     xlog=0, ylog=0, marker=['',''], xtitle='photon energy [eV]', ytitle='Flux [photons/s/eV]', show=0)

# ax.errorbar(tktE["bin_center"],
#             tktE["histogram"] * n_total / h_total, linestyle='', marker='',
#             yerr = tktE['histogram_sigma'][0] * n_total / h_total * 0.5)

plt.savefig('esrf1_histo_Elin.pdf')


# log
e, f, w = light_source.calculate_spectrum(numpy.linspace(100,100000,1000))
ff = f / (e * 1e-3) # in ph/eV

fig, ax = plot(
     e, ff,
     tktE["bin_path"], tktE["histogram_path"] * n_total / h_total,
     xlog=1, ylog=1, marker=['',''], xtitle='photon energy [eV]', ytitle='Flux [photons/s/eV]', show=0)

# ax.errorbar(tktE["bin_center"],
#             tktE["histogram"] * n_total / h_total, linestyle='', marker='',
#             yerr = tktE['histogram_sigma'][0] * n_total / h_total * 0.5)

plt.savefig('esrf1_histo_Elog.pdf')
plt.show()