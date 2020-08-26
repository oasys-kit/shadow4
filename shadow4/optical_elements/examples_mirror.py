import numpy

from syned.beamline.beamline import BeamlineElement
from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.mirrors.mirror import Mirror as SyMirror

from shadow4.beam.beam import Beam
from shadow4.sources.source_geometrical.gaussian import SourceGaussian

from shadow4.syned.shape import Rectangle, Ellipse, TwoEllipses # TODO from syned.beamline.shape
from shadow4.syned.shape import Toroidal, Conic, NumericalMesh # TODO from syned.beamline.shape
from shadow4.syned.shape import Plane, Sphere, Ellipsoid, Paraboloid, Hyperboloid # TODO from syned.beamline.shape
from shadow4.syned.shape import SphericalCylinder # TODO from syned.beamline.shape


from shadow4.optical_elements.mirror import Mirror

from shadow4.compatibility.beam3 import Beam3
from Shadow.ShadowTools import plotxy


def test_branch_1(do_plot=True):
    #
    # source
    #
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

    # surface shape

    surface_shape = Conic(conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0])
    # surface_shape = Toroidal(min_radius=5.0, maj_radius=2.0) # Plane() # SurfaceShape()


    # boundaries
    # boundary_shape = None
    rlen1 = 5e-05
    rlen2 = 5e-05
    rwidx1 = 2e-05
    rwidx2 = 2e-05
    boundary_shape = Rectangle(x_left=-rwidx2,x_right=rwidx1,y_bottom=-rlen2,y_top=rlen1)
    # boundary_shape = Rectangle(x_left=-1e-05, x_right=2e-05, y_bottom=-5e-04, y_top=7e-04)
    # boundary_shape = Ellipse(a_axis_min=-rwidx2/2, a_axis_max=rwidx2/2, b_axis_min=-rlen2/2, b_axis_max=rlen2/2)
    # boundary_shape = Ellipse(a_axis_min=-0e-05, a_axis_max=1e-05, b_axis_min=-1.5e-05, b_axis_max=2.5e-05)
    # boundary_shape = Ellipse(a_axis_min=0, a_axis_max=1e-05, b_axis_min=-0.0005, b_axis_max=0)

    # rlen1 = -2.5e-05
    # rlen2 = 5e-05
    # rwidx1 = 1e-05
    # rwidx2 = 2e-05
    #
    # boundary_shape = TwoEllipses(
    #     a1_axis_min=-rwidx1 / 2, a1_axis_max=rwidx1 / 2, b1_axis_min=-rlen1 / 2, b1_axis_max=rlen1 / 2,
    #     a2_axis_min=-rwidx2 / 2, a2_axis_max=rwidx2 / 2, b2_axis_min=-rlen2 / 2, b2_axis_max=rlen2 / 2)

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
    print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam(beam0)
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        beam1s3 = Beam3.initialize_from_shadow4_beam(beam1)
        plotxy(beam1s3, 1, 3, title="Image 1", nbins=101, nolost=1)
        mirr1s3 = Beam3.initialize_from_shadow4_beam(mirr1)
        plotxy(mirr1s3, 2, 1, title="Footprint 1", nbins=101, nolost=1)

    #
    # M2
    #
    mirror2 = Mirror()

    mirror2.set_positions(10,100,3e-3)
    mirror2.set_surface_conic([0,0,0,0,0,0,0,0,-1,0])
    mirror2.set_boundaries_rectangle(-100e-6,100e-6,-150e-6,150e-6)
    print(mirror2.info())
    #
    #
    if do_plot:
        beam2, mirr2 = mirror2.trace_beam(beam1)
        beam2s3 = Beam3.initialize_from_shadow4_beam(beam2)
        plotxy(beam2s3, 1, 3, title="Image 2", nbins=101, nolost=1)
        mirr2s3 = Beam3.initialize_from_shadow4_beam(mirr2)
        plotxy(mirr2s3, 2, 1, title="Footprint 2", nbins=101, nolost=1)

def test_branch_2(do_plot=True):
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

    # surface shape
    surface_shape = Toroidal(min_radius=0.157068, maj_radius=358.124803-0.157068) # Plane() # SurfaceShape()


    # boundaries
    boundary_shape = None #Rectangle(x_left=-rwidx2,x_right=rwidx1,y_bottom=-rlen2,y_top=rlen1)

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
    print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam(beam0)
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        beam1s3 = Beam3.initialize_from_shadow4_beam(beam1)
        plotxy(beam1s3, 1, 3, title="Image 1", nbins=101, nolost=1)
        mirr1s3 = Beam3.initialize_from_shadow4_beam(mirr1)
        plotxy(mirr1s3, 2, 1, title="Footprint 1", nbins=101, nolost=1)

def test_branch_3(do_plot=True):
    #
    # source
    #
    # beam0 = Beam.initialize_as_pencil(N=500)
    source = SourceGaussian.initialize_from_keywords(number_of_rays=100000,
                 sigmaX=0.0,
                 sigmaY=0.0,
                 sigmaZ=0.0,
                 sigmaXprime=1e-4,
                 sigmaZprime=1e-4,)
    beam0 = Beam()
    beam0.genSource(source)
    print(beam0.info())

    if do_plot:
        beam0s3 = Beam3.initialize_from_shadow4_beam(beam0)
        plotxy(beam0s3, 4, 6, title="Image 0", nbins=201)

    #
    # syned definitopns
    #

    # surface shape

    surface_shape = NumericalMesh("%s/test_shadow4.hdf5" % OASYS_HOME)


    # boundaries
    rlen1 = 0.6
    rlen2 = 0.6
    rwidx1 = 0.05
    rwidx2 = 0.05
    boundary_shape = Rectangle(x_left=-rwidx2,x_right=rwidx1,y_bottom=-rlen2,y_top=rlen1)


    symirror1 = SyMirror(
                name="M1",
                surface_shape=surface_shape,
                boundary_shape=boundary_shape,
                coating=None,
                coating_thickness=None)

    coordinates_syned = ElementCoordinates(p = 100.0,
                                           q = 1000.0,
                                           angle_radial = 88.8 * numpy.pi / 180,)

    beamline_element_syned = BeamlineElement(optical_element=symirror1, coordinates=coordinates_syned)

    #
    # shadow definitions
    #
    mirror1 = Mirror(beamline_element_syned=beamline_element_syned)
    print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam(beam0)
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        beam1s3 = Beam3.initialize_from_shadow4_beam(beam1)
        plotxy(beam1s3, 1, 3, title="Image 1", nbins=101, nolost=1)
        mirr1s3 = Beam3.initialize_from_shadow4_beam(mirr1)
        plotxy(mirr1s3, 2, 1, title="Footprint 1", nbins=101, nolost=1)


def test_branch_4(do_plot=True):


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
    beam0.set_photon_wavelength(5e-10)
    print(beam0.info())

    if do_plot:
        beam0s3 = Beam3.initialize_from_shadow4_beam(beam0)
        plotxy(beam0s3, 4, 6, title="Image 0", nbins=201)

    #
    # syned definitopns
    #

    # surface shape
    surface_shape = Conic(conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0])

    # boundaries
    rlen1 = 5e-05
    rlen2 = 5e-05
    rwidx1 = 2e-05
    rwidx2 = 2e-05
    boundary_shape = Rectangle(x_left=-rwidx2,x_right=rwidx1,y_bottom=-rlen2,y_top=rlen1)

    symirror1 = SyMirror(
                name="M1",
                surface_shape=surface_shape,
                boundary_shape=boundary_shape,
                coating="%s/SiC.dat" % OASYS_HOME,
                coating_thickness=None)

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 6.0,
                                           angle_radial = 88.8 * numpy.pi / 180,)

    beamline_element_syned = BeamlineElement(optical_element=symirror1, coordinates=coordinates_syned)

    #
    # shadow definitions
    #
    mirror1 = Mirror(beamline_element_syned=beamline_element_syned)
    print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam(beam0)
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        beam1s3 = Beam3.initialize_from_shadow4_beam(beam1)
        plotxy(beam1s3, 1, 3, title="Image 1", nbins=101, nolost=1)
        mirr1s3 = Beam3.initialize_from_shadow4_beam(mirr1)
        plotxy(mirr1s3, 2, 1, title="Footprint 1", nbins=101, nolost=1)

def test_branch_5(surface_type, do_plot=True):
    #
    # source
    #
    source = SourceGaussian.initialize_from_keywords(number_of_rays=100000,
                 sigmaX=0.0,
                 sigmaY=0.0,
                 sigmaZ=0.0,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,)
    beam0 = Beam()
    beam0.genSource(source)
    # print(beam0.info())



    #
    # syned definitopns
    #

    # surface shape

    if surface_type == "plane":
        surface_shape = Plane()
    elif surface_type == "sphere":
        surface_shape = Sphere()
        surface_shape.initialize_from_p_q(10.0, 10.0, grazing_angle=(90.0 - 88.8) * numpy.pi / 180)
        print(">>><<<<>>><<<", surface_shape.info())
    elif surface_type == "spherical_cylinder_tangential":
        surface_shape = SphericalCylinder()
        surface_shape.set_direction_tangential()
        surface_shape.initialize_from_p_q(10.0, 10.0, grazing_angle=(90.0 - 88.8) * numpy.pi / 180)
    elif surface_type == "spherical_cylinder_sagittal":
        surface_shape = SphericalCylinder()
        surface_shape.set_direction_sagittal()
        surface_shape.initialize_from_p_q(10.0, 10.0, grazing_angle=(90.0 - 88.8) * numpy.pi / 180)
    elif surface_type == "ellipsoid":
        surface_shape = Ellipsoid()
        #TODO: this is to be discussed in the future: the ellipsoid is only defined by the major and minor axes, which
        #is not enough to characterise the mirror surface. We need a third parameter that could be the incident
        #angle, or p. We cannot retrieve in shadow the ellipsoidal surface with only the axes, we need another
        #parameter. Shadow uses the subtended angle from the ellipsoid center to the mirror pole ELL_THE
        surface_shape.initialize_from_p_q(10.0, 10.0, grazing_angle=(90.0 - 88.8) * numpy.pi / 180)
    elif surface_type == "conic":
        surface_shape = Conic(conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0])
    else:
        raise Exception("undefined surface shape")


    # boundaries
    boundary_shape = None

    symirror1 = SyMirror(
                name="M1",
                surface_shape=surface_shape,
                boundary_shape=boundary_shape,
                coating=None,
                coating_thickness=None)

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 10.0,
                                           angle_radial = 88.8 * numpy.pi / 180,)

    beamline_element_syned = BeamlineElement(optical_element=symirror1, coordinates=coordinates_syned)

    #
    # shadow definitions
    #
    mirror1 = Mirror(beamline_element_syned=beamline_element_syned)

    #
    # run
    #
    print(mirror1.info())
    beam1, mirr1 = mirror1.trace_beam(beam0)


    #
    # check
    #

    if do_plot:
        beam1s3 = Beam3.initialize_from_shadow4_beam(beam1)
        plotxy(beam1s3, 1, 3, nbins=101, nolost=1, title=surface_type)
        # mirr1s3 = Beam3.initialize_from_shadow4_beam(mirr1)
        # plotxy(mirr1s3, 2, 1, title="Footprint 1", nbins=101, nolost=1)


if __name__ == "__main__":

    # see corresponding Oasys workspace in shadow4/oasys_workspaces

    from srxraylib.plot.gol import set_qt
    set_qt()

    OASYS_HOME = "/Users/srio/Oasys/"

    # test_branch_1(do_plot=False) # two plane mirrors
    # test_branch_2(do_plot=False) # toroid
    # test_branch_3(do_plot=False) # mesh
    # test_branch_4(do_plot=False) # prerefl

    for myconicshape in ["plane", "sphere", "spherical_cylinder_tangential", "spherical_cylinder_sagittal","ellipsoid"]:
        test_branch_5(myconicshape,do_plot=True) # conic mirrors
