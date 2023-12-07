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
    if True:
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
