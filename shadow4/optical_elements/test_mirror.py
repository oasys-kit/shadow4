import numpy

from syned.beamline.beamline import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.mirrors.mirror import Mirror as SyMirror

from shadow4.beam.beam import Beam
from shadow4.sources.source_geometrical.gaussian import SourceGaussian
from syned.beamline.shape import Plane
from syned.beamline.shape import Rectangle, Ellipse
from shadow4.optical_elements.mirror import Mirror

if __name__ == "__main__":

    from shadow4.compatibility.beam3 import Beam3
    from Shadow.ShadowTools import plotxy


    do_plot = False
    #
    # source
    #
    # beam0 = Beam.initialize_as_pencil(N=500)
    source = SourceGaussian.initialize_from_keywords(number_of_rays=100000,
                 sigmaX=0.0,
                 sigmaY=0.0,
                 sigmaZ=0.0,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,)
    beam0 = Beam()
    beam0.genSource(source)
    print(beam0.info())

    if do_plot:
        beam0s3 = Beam3.initialize_from_shadow4_beam(beam0)
        plotxy(beam0s3, 4, 6, title="Image 0", nbins=201)


    #
    # syned definitopns
    #
    surface_shape = Plane() # SurfaceShape()
    rlen1 = 5e-05
    rlen2 = 5e-05
    rwidx1 = 2e-05
    rwidx2 = 2e-05
    boundary_shape = Rectangle(x_left=-rwidx2,x_right=rwidx1,y_bottom=-rlen2,y_top=rlen1) # None

    symirror1 = SyMirror(
                name="M1",
                surface_shape=surface_shape,
                boundary_shape=boundary_shape,
                coating=None,
                coating_thickness=None)

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 6.0,
                                           angle_radial = 88.8 * numpy.pi / 180,)

    beamline_element_syned = BeamlineElement(optical_element=symirror1, coordinates=coordinates_syned)

    #
    # shadow definitions
    #
    mirror1 = Mirror(beamline_element_syned=beamline_element_syned)
    # print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam(beam0)
    print(mirr1.info())

    #
    # check
    #

    if True: #do_plot:
        # beam1s3 = Beam3.initialize_from_shadow4_beam(beam1)
        # plotxy(beam1s3, 1, 3, title="Image 1", nbins=101)
        mirr1s3 = Beam3.initialize_from_shadow4_beam(mirr1)
        plotxy(mirr1s3, 2, 1, title="Footprint 1", nbins=101, nolost=1)

    # #
    # # M2
    # #
    # mirror2 = Mirror()
    # mirror2.set_positions(10,100,3e-3)
    # mirror2.set_surface_conic([0,0,0,0,0,0,0,0,-1,0])
    #
    # print(mirror2.info())
    # beam2, mirr2 = mirror2.trace_beam(beam1)
    #
    # beam2s3 = Beam3.initialize_from_shadow4_beam(beam2)
    # plotxy(beam2s3, 1, 3, title="Image 2", nbins=101)
    # mirr2s3 = Beam3.initialize_from_shadow4_beam(mirr2)
    # plotxy(mirr2s3, 2, 1, title="Footprint 2", nbins=101)
