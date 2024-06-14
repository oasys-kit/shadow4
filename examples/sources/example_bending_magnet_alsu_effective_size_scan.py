import numpy
from srxraylib.plot.gol import plot

def run_source(Dtheta = 10e-3):
    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam
    electron_beam = S4ElectronBeam(energy_in_GeV=2,energy_spread=0,current=0.4)
    electron_beam.set_sigmas_all(sigma_x=0,sigma_y=0,sigma_xp=0,sigma_yp=0)

    #magnetic structure
    from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet


    magnetic_field = -0.87
    radius = numpy.abs(S4BendingMagnet.calculate_magnetic_radius(-0.87, electron_beam.energy())) # radius = 7.668140119497748
    length = Dtheta * radius # 0.07668140119497747

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
    beam = light_source.get_beam(sample_emission_cone_in_horizontal=0)

    return beam


#
# main
#

npoint = 10
FWHM = numpy.zeros(npoint)
SIGMA = numpy.zeros(npoint)
DTHETA = numpy.linspace(1, 20, npoint) * 1e-3

for i in range(npoint):
    print('Calculating index %d if %d' % (i, npoint))
    beam = run_source(Dtheta=DTHETA[i])
    beam.retrace(0.0)
    tkt = beam.histo1(1, nbins=100, calculate_widths=1)

    if 0: plot(1e6 * tkt["bin_path"], tkt["histogram_path"],
        xtitle=r"x [$\mu$m]", ytitle="Intensity [arbitrary units]",
        title="FWHM: %f, StDev: %f" % (1e6 * FWHM[i], 1e6 * SIGMA[i] )
        )

    FWHM[i] = tkt['fwhm']
    SIGMA[i] = beam.get_column(1).std()

    if numpy.isnan(FWHM[i]):   is

plot(DTHETA, FWHM,
     DTHETA, SIGMA,
     legend=['FWHM', 'StDev'])

print(FWHM)