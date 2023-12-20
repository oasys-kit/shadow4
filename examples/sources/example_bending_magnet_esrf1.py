

# electron beam
from shadow4.sources.s4_electron_beam import S4ElectronBeam
electron_beam = S4ElectronBeam(energy_in_GeV=6.04,energy_spread=0.001,current=0.2)
electron_beam.set_sigmas_all(sigma_x=7.8e-05,sigma_y=4.87179e-05,sigma_xp=3.6e-05,sigma_yp=1.05556e-06)

#magnetic structure
from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet
source = S4BendingMagnet(
                 radius=25.18408918746048, # from syned BM, can be obtained as S4BendingMagnet.calculate_magnetic_radius(0.8, electron_beam.energy())
                 magnetic_field=0.8, # from syned BM
                 length=0.02518408918746048, # from syned BM = abs(BM divergence * magnetic_field)
                 emin=8000.0,     # Photon energy scan from energy (in eV)
                 emax=8000.0001,     # Photon energy scan to energy (in eV)
                 ng_e=100,     # Photon energy scan number of points
                 flag_emittance=1, # when sampling rays: Use emittance (0=No, 1=Yes)
                 )


#light source
from shadow4.sources.bending_magnet.s4_bending_magnet_light_source import S4BendingMagnetLightSource
light_source = S4BendingMagnetLightSource(name='BendingMagnet',
                                          electron_beam=electron_beam,
                                          magnetic_structure=source,
                                          nrays=25000,
                                          seed=5676561)

import time
t0 = time.time()
beam = light_source.get_beam()
print(">> time: ", time.time()-t0)
# test plot
from srxraylib.plot.gol import plot_scatter, plot
rays = beam.get_rays()
# plot_scatter(1e6 * rays[:, 0], 1e6 * rays[:, 2], title='(X,Z) in microns')


tktI = beam.histo1(6, ref=23)
tktS = beam.histo1(6, ref=24)
tktP = beam.histo1(6, ref=25)
plot(
    tktI["bin_path"], tktI["histogram_path"],
    tktS["bin_path"], tktS["histogram_path"],
    tktP["bin_path"], tktP["histogram_path"],
)

tktP3 = beam.histo1(6, ref=33)
plot(
    tktP3["bin_path"], tktP3["histogram_path"], )

import numpy
plot(
    light_source.angle_array_mrad,   light_source.fm[:, 0]/ numpy.max(light_source.fm[:, 0]),
    light_source.angle_array_mrad, light_source.fm_s[:, 0]/ numpy.max(light_source.fm[:, 0]),
    light_source.angle_array_mrad, light_source.fm_p[:, 0]/ numpy.max(light_source.fm[:, 0]),
    tktI["bin_path"] * 1e3, tktI["histogram_path"] / numpy.max(tktI["histogram_path"]),
    tktS["bin_path"] * 1e3, tktS["histogram_path"] / numpy.max(tktI["histogram_path"]),
    tktP["bin_path"] * 1e3, tktP["histogram_path"] / numpy.max(tktI["histogram_path"]),
    xtitle="Angle / mrad", ytitle="Flux at Emin=%f eV" % (8000))