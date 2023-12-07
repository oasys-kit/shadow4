from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian
from shadow4.tools.graphics import plotxy
from srxraylib.plot.gol import plot_scatter

def example_point_source():

    a = SourceGaussian.initialize_point_source(
                 sigmaXprime=3e-6,
                 sigmaZprime=1e-6,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0],
                 nrays=10000,)
    print(a.info())

    x,y,z = a._get_arrays_real_space()
    print("x.shape:",x.shape)
    print("y.shape:",y.shape)
    print("z.shape:",z.shape)

    xp,zp = a._get_arrays_direction_space()
    print("xp.shape:",xp.shape)
    print("zp.shape:",zp.shape)


    VP = a._get_volume_divergences()
    print("VP",VP.shape,VP.size)


    Vx = a._get_volume_real_space()
    print("Vx: ",Vx.shape)

    V = a._get_volume()
    print("V: ",V.shape)


    beam = a.get_beam()

    plotxy(beam,4,6,title="point source")

def example_collimated_source():
    #
    a = SourceGaussian.initialize_collimated_source(
                 sigmaX=1.0,
                 sigmaY=0.0,
                 sigmaZ=1.0,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0],
                 nrays=10000,)

    print(a.info())
    beam = a.get_beam()

    plotxy(beam,1,3,title="collimated source")

    plot_scatter(beam.get_column(1),beam.get_column(3),xtitle="col 1",ytitle="col 3",title="collimated source")
    print(beam.info())

def example_double_gaussian():
    #
    #
    a = SourceGaussian(
                                                sigmaX=2.0e-6,
                                                sigmaY=1.0e-3,
                                                sigmaZ=1.0e-6,
                                                sigmaXprime=3e-6,
                                                sigmaZprime=3e-6,
                                                real_space_center=[0.0, 0.0, 0.0],
                                                direction_space_center=[0.0, 0.0],
                                                nrays=10000,)

    print(a.info())

    beam = a.get_beam()
    plot_scatter(beam.get_column(2), beam.get_column(1), title="double Gaussian", xtitle="col 2", ytitle="col 1")
    plot_scatter(beam.get_column(1), beam.get_column(3), title="double Gaussian", xtitle="col 1", ytitle="col 3")
    plot_scatter(beam.get_column(4), beam.get_column(6), title="double Gaussian", xtitle="col 4", ytitle="col 6")
    print(beam.info())

    print(isinstance(a, SourceGaussian))

if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()
    example_point_source()
    example_collimated_source()
    example_double_gaussian()