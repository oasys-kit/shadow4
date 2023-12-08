import numpy
from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian

if __name__ == "__main__":
    from syned.util.json_tools import load_from_json_file
    from srxraylib.plot.gol import plot,plot_scatter, set_qt
    set_qt()

    do_plot = 1

    # Gaussian
    if False:
        gs = SourceGaussian(
                     nrays=10000,
                     sigmaX=1.0e-6,
                     sigmaY=0.0,
                     sigmaZ=1.0e-6,
                     sigmaXprime=0.0002,
                     sigmaZprime=0.0002,
                     # real_space_center     = numpy.array([0.0,0.0,0.0]),
                     # direction_space_center= numpy.array([0.0,0.0]),
                                     )

        beam = gs.get_beam()

        if do_plot:
            plot_scatter(1e6*beam.rays[:,0], 1e6*beam.rays[:,2], xrange=[-6,6], yrange=[-3,3], title="Gaussian")

        print(gs.info())

        file_name = 'example_sources_file_io_gaussian.json'
        gs.to_json(file_name=file_name)

        gs_new = load_from_json_file(file_name,
                                     exec_commands=
                                     "from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian")

        print(gs_new)
        beam2 = gs_new.get_beam()
        if do_plot:
            plot_scatter(1e6*beam2.rays[:,0], 1e6*beam2.rays[:,2], xrange=[-6,6], yrange=[-3,3], title="Gaussian **LOADED FROM FILE**")


    # rectangle
    if False:
        gs = SourceGeometrical(nrays=1000)
        gs.set_spatial_type_rectangle(4,2)
        beam1 = gs.get_beam()

        if do_plot:
            plot_scatter(beam1.rays[:,0], beam1.rays[:,2], xrange=[-3,3], yrange=[-3,3], title="Rectangle")

        print(gs.info())

        file_name = 'example_sources_file_io_rectangle.json'
        gs.to_json(file_name=file_name)

        gs_new = load_from_json_file(file_name,
                                     exec_commands=
                                     "from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical")
        print(gs_new)
        beam2 = gs_new.get_beam()
        if do_plot:
            plot_scatter(beam2.rays[:,0], beam2.rays[:,2], xrange=[-3,3], yrange=[-3,3], title="Rectangle **LOADED FROM FILE**")

        print(gs_new.get_info())
        print(gs_new.to_python_code())

    # electron beam
    from shadow4.sources.s4_electron_beam import S4ElectronBeam

    electron_beam = S4ElectronBeam(energy_in_GeV=6.04, energy_spread=0.001, current=0.2)
    electron_beam.set_sigmas_all(sigma_x=7.8e-05, sigma_y=4.87179e-05, sigma_xp=3.6e-05, sigma_yp=1.05556e-06)

    # Bending magnet
    from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet

    bm = S4BendingMagnet(
                 radius=25.18408918746048, # from syned BM, can be obtained as S4BendingMagnet.calculate_magnetic_radius(0.8, electron_beam.energy())
                 magnetic_field=0.8, # from syned BM
                 length=0.02518408918746048, # from syned BM = abs(BM divergence * magnetic_field)
                 emin=8000.0,     # Photon energy scan from energy (in eV)
                 emax=8100.0,     # Photon energy scan to energy (in eV)
                 ng_e=100,     # Photon energy scan number of points
                 ng_j=100,     # Number of points in electron trajectory (per period) for internal calculation only
                 flag_emittance=1, # when sampling rays: Use emittance (0=No, 1=Yes)
                 )
    print(bm.info())

    # light source
    from shadow4.sources.bending_magnet.s4_bending_magnet_light_source import S4BendingMagnetLightSource

    light_source = S4BendingMagnetLightSource(name='BendingMagnet',
                                              electron_beam=electron_beam,
                                              magnetic_structure=bm,
                                              nrays=10000,
                                              seed=5676561)
    print(light_source.info())
    file_name = 'example_sources_file_io_bending_magnet.json'
    light_source.to_json(file_name=file_name)

    light_source_new = load_from_json_file(file_name,
                                 exec_commands=
                                 ["from shadow4.sources.bending_magnet.s4_bending_magnet import S4BendingMagnet",
                                  "from shadow4.sources.s4_electron_beam import S4ElectronBeam",
                                  "from shadow4.sources.bending_magnet.s4_bending_magnet_light_source import S4BendingMagnetLightSource"])
    print(light_source_new.info())


    import time

    if do_plot:
        t0 = time.time()
        beam1 = light_source.get_beam(verbose=0)
        t1 = time.time()
        plot_scatter(beam1.rays[:, 0], beam1.rays[:, 2], title="BM ", show=0)
        plot_scatter(beam1.rays[:, 3], beam1.rays[:, 5], title="BM DIVERGENCES", show=1)

        t2 = time.time()
        beam2 = light_source_new.get_beam()
        t3 = time.time()
        plot_scatter(beam2.rays[:, 0], beam2.rays[:, 2], title="BM **LOADED FROM FILE**", show=0)
        plot_scatter(beam2.rays[:, 3], beam2.rays[:, 5], title="BM DIVERGENCES **LOADED FROM FILE**")


        print("Running time 1: ", t1 - t0)
        print("Running time 2: ", t3 - t2)