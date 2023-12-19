import numpy
from syned.storage_ring.electron_beam import ElectronBeam

from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet
from shadow4.sources.bending_magnet.s4_bending_magnet_light_source import S4BendingMagnetLightSource


if __name__ == "__main__":
    from srxraylib.plot.gol import plot_scatter, set_qt

    set_qt()

    flag_emittance = True        # when sampling rays: Use emittance (0=No, 1=Yes)

    electron_beam = ElectronBeam(energy_in_GeV=1.9, current=0.4,
                                 moment_xx   = (39e-6)**2,
                                 moment_xpxp = (2000e-12 / 51e-6)**2,
                                 moment_yy   = (31e-6)**2,
                                 moment_ypyp = (30e-12 / 31e-6)**2,
                                 )

    emin = 1000.0  # Photon energy scan from energy (in eV)
    emax = 1001.0  # Photon energy scan to energy (in eV)
    ng_e = 200     # Photon energy scan number of points

    bm = S4BendingMagnet.initialize_from_magnetic_field_divergence_and_electron_energy(magnetic_field=-1.26754,
                                                                                       divergence=69e-3,
                                                                                       electron_energy_in_GeV=1.9,
                                                                                       emin=emin,  # Photon energy scan from energy (in eV)
                                                                                       emax=emax,  # Photon energy scan to energy (in eV)
                                                                                       ng_e=ng_e,  # Photon energy scan number of points
                                                                                       flag_emittance=flag_emittance,  # when sampling rays: Use emittance (0=No, 1=Yes)
                                                                                       )

    print(bm.info())

    light_source = S4BendingMagnetLightSource(electron_beam=electron_beam, magnetic_structure=bm,
                                              nrays=25000, seed=123456)

    beam = light_source.get_beam(F_COHER=0, EPSI_DX=0.0, EPSI_DZ=0.0, verbose=False)

    rays = beam.rays

    plot_scatter(rays[:,0]*1e6, rays[:,2]*1e6,xtitle="X um",ytitle="Z um")
    plot_scatter(rays[:,1], rays[:,0]*1e6,xtitle="Y m",ytitle="X um")
    plot_scatter(rays[:,1], rays[:,2]*1e6,xtitle="Y m",ytitle="Z um")
    plot_scatter(rays[:,3]*1e6, rays[:,5]*1e6,xtitle="X' urad",ytitle="Z' urad")

    from srxraylib.plot.gol import plot
    tktI = beam.histo1(6, ref=23)
    tktS = beam.histo1(6, ref=24)
    tktP = beam.histo1(6, ref=25)
    plot(
        light_source.angle_array_mrad * 1e-3, light_source.fm[:, 0]/ numpy.max(light_source.fm[:, 0]),
        light_source.angle_array_mrad * 1e-3, light_source.fm_s[:, 0]/ numpy.max(light_source.fm[:, 0]),
        light_source.angle_array_mrad * 1e-3, light_source.fm_p[:, 0]/ numpy.max(light_source.fm[:, 0]),
        tktI["bin_path"], tktI["histogram_path"] / numpy.max(tktI["histogram_path"]),
        tktS["bin_path"], tktS["histogram_path"] / numpy.max(tktI["histogram_path"]),
        tktP["bin_path"], tktP["histogram_path"] / numpy.max(tktI["histogram_path"]),
         )

    tktP3 = beam.histo1(6, ref=33)
    plot(
        tktP3["bin_path"], tktP3["histogram_path"],)
