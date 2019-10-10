
from syned.storage_ring.electron_beam import ElectronBeam
from syned.storage_ring.magnetic_structures.bending_magnet import BendingMagnet

from shadow4.sources.bending_magnet.source_bending_magnet import SourceBendingMagnet



if __name__ == "__main__":
    from srxraylib.plot.gol import plot_scatter, set_qt
    set_qt()

    flag_emittance = True        # when sampling rays: Use emittance (0=No, 1=Yes)

    syned_electron_beam = ElectronBeam(energy_in_GeV=1.9,current=0.4,
                                       moment_xx   = (39e-6)**2,
                                       moment_xpxp = (2000e-12 / 51e-6)**2,
                                       moment_yy   = (31e-6)**2,
                                       moment_ypyp = (30e-12 / 31e-6)**2,
                                       )

    # syned_bending_magnet = BendingMagnet(radius=-5.0,magnetic_field=-1.26754,length=5.0*0.069)

    syned_bending_magnet = BendingMagnet.initialize_from_magnetic_field_divergence_and_electron_energy(
        magnetic_field=-1.26754,divergence=69e-3,electron_energy_in_GeV=1.9)
    #
    # syned_bending_magnet = BendingMagnet.initialize_from_magnetic_radius_divergence_and_electron_energy(
    #     magnetic_radius=-5.0,divergence=69e-3,electron_energy_in_GeV=1.9)


    emin = 1000.0                # Photon energy scan from energy (in eV)
    emax = 1001.0                # Photon energy scan to energy (in eV)
    ng_e = 200                # Photon energy scan number of points
    ng_j = 100                # Number of points in electron trajectory (per period) for internal calculation only




    bm = SourceBendingMagnet(syned_electron_beam=syned_electron_beam,
                 syned_bending_magnet=syned_bending_magnet,
                 emin=emin,               # Photon energy scan from energy (in eV)
                 emax=emax,               # Photon energy scan to energy (in eV)
                 ng_e=ng_e,                    # Photon energy scan number of points
                 ng_j=ng_j,                    # Number of points in electron trajectory (per period) for internal calculation only
                 flag_emittance=flag_emittance,           # when sampling rays: Use emittance (0=No, 1=Yes)
                )

    print(bm.info())

    rays = bm.calculate_rays(F_COHER=0,NRAYS=5000,SEED=123456,EPSI_DX=0.0,EPSI_DZ=0.0,verbose=False)


    plot_scatter(rays[:,0]*1e6,rays[:,2]*1e6,xtitle="X um",ytitle="Z um")
    plot_scatter(rays[:,1],rays[:,0]*1e6,xtitle="Y m",ytitle="X um")
    plot_scatter(rays[:,1],rays[:,2]*1e6,xtitle="Y m",ytitle="Z um")
    plot_scatter(rays[:,3]*1e6,rays[:,5]*1e6,xtitle="X' urad",ytitle="Z' urad")
