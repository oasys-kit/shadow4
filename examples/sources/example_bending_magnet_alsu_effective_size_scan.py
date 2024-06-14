import numpy
from srxraylib.plot.gol import plot
import matplotlib.pylab as plt

def run_source(Dtheta = 10e-3, all_terms=0):
    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam
    electron_beam = S4ElectronBeam(energy_in_GeV=2,energy_spread=0,current=0.4)
    # if all_terms == 0:
    #     electron_beam.set_sigmas_all(sigma_x=0,sigma_y=0,sigma_xp=0,sigma_yp=0)
    # else:
    electron_beam.set_sigmas_all(sigma_x=5e-6, sigma_y=5e-6, sigma_xp=70e-12/5e-6, sigma_yp=70e-12/5e-6)

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
                     flag_emittance=all_terms, # when sampling rays: Use emittance (0=No, 1=Yes)
                     epsi_dx=0.0,  # position of X waist [m]
                     epsi_dz=0.0 , # position of Z waist [m]
                     )



    #light source
    from shadow4.sources.bending_magnet.s4_bending_magnet_light_source import S4BendingMagnetLightSource
    light_source = S4BendingMagnetLightSource(name='BendingMagnet', electron_beam=electron_beam, magnetic_structure=source, nrays=50000, seed=5676561)
    beam = light_source.get_beam(sample_emission_cone_in_horizontal=all_terms)

    return beam, radius


#
# main
#

npoint = 10
FWHM = numpy.zeros(npoint)
SIGMA = numpy.zeros(npoint)
SIGMA_ALLTERMS = numpy.zeros(npoint)
DTHETA = numpy.linspace(1, 10, npoint) * 1e-3


for i in range(npoint):
    print('Calculating index %d if %d' % (i, npoint))
    beam, radius = run_source(Dtheta=DTHETA[i])
    beam.retrace(0.0)

    if 0:
        tkt = beam.histo1(1, nbins=100, calculate_widths=1)
        plot(1e6 * tkt["bin_path"], tkt["histogram_path"],
        xtitle=r"x [$\mu$m]", ytitle="Intensity [arbitrary units]",
        title="FWHM: %f, StDev: %f" % (1e6 * FWHM[i], 1e6 * SIGMA[i] )
        )
        FWHM[i] = tkt['fwhm']

    SIGMA[i] = beam.get_column(1).std()

for i in range(npoint):
    print('Calculating index %d if %d' % (i, npoint))
    beam, radius = run_source(Dtheta=DTHETA[i], all_terms=1)
    beam.retrace(0.0)
    tkt = beam.histo1(1, nbins=100, calculate_widths=1)
    FWHM[i] = tkt['fwhm']
    print(">>>>>>>>>>>>>", tkt['fwhm'])
    # plot(tkt['bin_path'],tkt['histogram_path'])
    SIGMA_ALLTERMS[i] = beam.get_column(1).std()

    # if numpy.isnan(FWHM[i]):   is


DTHETA_NUCARA = numpy.linspace(1, 10, npoint * 10) * 1e-3
SIGMA_NUCARA = radius**2 * (3/2 + numpy.sin(DTHETA_NUCARA) / 2 / DTHETA_NUCARA - 4 * numpy.sin(DTHETA_NUCARA / 2) / DTHETA_NUCARA - \
                            (1 - 2 / DTHETA_NUCARA * numpy.sin(DTHETA_NUCARA / 2))**2 )
SIGMA_NUCARA = numpy.sqrt(SIGMA_NUCARA)

# plot(DTHETA_NUCARA, SIGMA_NUCARA, show=0 )

plot(
     # 1e3 * DTHETA, 1e6 * FWHM,
     1e3 * DTHETA, 1e6 * SIGMA,
     1e3 * DTHETA, 1e6 * SIGMA_ALLTERMS,
     1e3 * DTHETA_NUCARA, 1e6 * SIGMA_NUCARA,
     legend=['StDev from rays (only trajectory)', 'StDev from rays (emittance+emission)', 'StDev analytical (geometrical)'],
     xtitle="Horizontal acceptance [mrad]", ytitle=r"Effective H size [$\mu$m]",
     linestyle=['', '', None], marker=['+', 'x', None], show=0)

plt.savefig("alsu_std.pdf")

SIGMA_APPROX = numpy.sqrt(SIGMA_NUCARA**2 + (5e-6)**2 + (163e-6 * DTHETA_NUCARA*radius)**2 )
plot(
     1e3 * DTHETA, 1e6 * SIGMA_ALLTERMS,
     1e3 * DTHETA_NUCARA, 1e6 * SIGMA_APPROX,
     1e3 * DTHETA, 1e6 * FWHM,
     1e3 * DTHETA_NUCARA, 1e6 * 2.355 * SIGMA_APPROX,
     legend=['StDev from rays', 'StDev analytical',
             'FWHM from rays', '2.355*StDev analytical'],
     xtitle="Horizontal acceptance [mrad]", ytitle=r"Effective H size [$\mu$m]",
     linestyle=['', None, '', None], marker=['x', None, 'x', None], show=0)

plt.savefig("alsu_std_fwhm.pdf")

plt.show()