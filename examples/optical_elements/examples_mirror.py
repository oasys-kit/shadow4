import numpy

from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.sources.source_geometrical.source_gaussian import SourceGaussian

from syned.beamline.shape import Rectangle, Direction, Side

from shadow4.beamline.s4_optical_element_decorators import SurfaceCalculation

from shadow4.beamline.optical_elements.mirrors.s4_conic_mirror import S4ConicMirror, S4ConicMirrorElement
from shadow4.beamline.optical_elements.mirrors.s4_toroid_mirror import S4ToroidMirror, S4ToroidMirrorElement
from shadow4.beamline.optical_elements.mirrors.s4_numerical_mesh_mirror import S4NumericalMeshMirror, S4NumericalMeshMirrorElement
from shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror, S4PlaneMirrorElement
from shadow4.beamline.optical_elements.mirrors.s4_ellipsoid_mirror import S4EllipsoidMirror, S4EllipsoidMirrorElement
from shadow4.beamline.optical_elements.mirrors.s4_hyperboloid_mirror import S4HyperboloidMirror, S4HyperboloidMirrorElement
from shadow4.beamline.optical_elements.mirrors.s4_paraboloid_mirror import S4ParaboloidMirror, S4ParaboloidMirrorElement
from shadow4.beamline.optical_elements.mirrors.s4_sphere_mirror import S4SphereMirror, S4SphereMirrorElement

from shadow4.tools.graphics import plotxy

from shadow4.physical_models.prerefl.prerefl import PreRefl

def example_branch_1(do_plot=True):
    #
    # source
    #
    source = SourceGaussian(nrays=100000,
                 sigmaX=0.0,
                 sigmaY=0.0,
                 sigmaZ=0.0,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,)
    beam0 = S4Beam()
    beam0.generate_source(source)
    print(beam0.info())

    if do_plot:
        plotxy(beam0, 4, 6, title="Image 0", nbins=201)

    #
    # syned definitopns
    #

    # surface shape


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


    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 6.0,
                                           angle_radial = numpy.radians(88.8))

    #
    # shadow definitions
    #
    mirror1 = S4ConicMirrorElement(optical_element=S4ConicMirror(name="M1",
                                                       conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
                                                       boundary_shape=boundary_shape),
                              coordinates=coordinates_syned,
                              input_beam=beam0)

    print(mirror1.info())

    
    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam()
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        plotxy(beam1, 1, 3, title="Image 1", nbins=101, nolost=1)
        plotxy(mirr1, 2, 1, title="Footprint 1", nbins=101, nolost=1)

    #
    # M2
    #
    
    mirror2 = S4ConicMirrorElement(optical_element=S4ConicMirror(conic_coefficients=[0,0,0,0,0,0,0,0,-1,0],
                                                                 boundary_shape=Rectangle(-100e-6,100e-6,-150e-6,150e-6)),
                                   coordinates=ElementCoordinates(p = 10.0, q = 100.0, angle_radial = (0.5*numpy.pi - 0.003)),
                                   input_beam=beam1)

    print(mirror2.info())
    #
    #
    if do_plot:
        beam2, mirr2 = mirror2.trace_beam()
        plotxy(beam2, 1, 3, title="Image 2", nbins=101, nolost=1)
        plotxy(mirr2, 2, 1, title="Footprint 2", nbins=101, nolost=1)

def example_branch_2(do_plot=True):
    source = SourceGaussian(nrays=100000,
                 sigmaX=0.0,
                 sigmaY=0.0,
                 sigmaZ=0.0,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,)
    beam0 = S4Beam()
    beam0.generate_source(source)
    print(beam0.info())

    if do_plot:
        plotxy(beam0, 4, 6, title="Image 0", nbins=201)

    #
    # syned definitopns
    #

    # boundaries
    boundary_shape = None #Rectangle(x_left=-rwidx2,x_right=rwidx1,y_bottom=-rlen2,y_top=rlen1)

    #
    # shadow definitions
    #
    mirror1 = S4ToroidMirrorElement(optical_element=S4ToroidMirror(name="M1",
                                                                     surface_calculation=SurfaceCalculation.EXTERNAL,
                                                                     min_radius=0.157068,
                                                                     maj_radius=358.124803 - 0.157068,
                                                                     f_torus=0,
                                                                     boundary_shape=boundary_shape),
                                      coordinates=ElementCoordinates(p = 10.0,
                                                                     q = 6.0,
                                                                     angle_radial = numpy.radians(88.8)),
                                      input_beam=beam0)


    print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam()
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        plotxy(beam1, 1, 3, title="Image 1", nbins=101, nolost=1)
        plotxy(mirr1, 2, 1, title="Footprint 1", nbins=101, nolost=1)

def example_branch_3(surface_shape_file, do_plot=True):
    #
    # source
    #
    # beam0 = Beam.initialize_as_pencil(N=500)
    source = SourceGaussian(nrays=100000,
                 sigmaX=0.0,
                 sigmaY=0.0,
                 sigmaZ=0.0,
                 sigmaXprime=1e-4,
                 sigmaZprime=1e-4,)
    beam0 = S4Beam()
    beam0.generate_source(source)
    print(beam0.info())

    if do_plot:
        plotxy(beam0, 4, 6, title="Image 0", nbins=201)

    #
    # syned definitopns
    #

    # boundaries
    rlen1 = 0.6
    rlen2 = 0.6
    rwidx1 = 0.05
    rwidx2 = 0.05

    #
    # shadow definitions
    #
    mirror1 = S4NumericalMeshMirrorElement(optical_element=S4NumericalMeshMirror(name="M1",
                                                                                 surface_data_file=surface_shape_file,
                                                                                 boundary_shape=Rectangle(x_left=-rwidx2,x_right=rwidx1,y_bottom=-rlen2,y_top=rlen1)),
                                           coordinates=ElementCoordinates(p = 100.0,
                                                                        q = 1000.0,
                                                                        angle_radial = numpy.radians(88.8)),
                                           input_beam=beam0)
    print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam()
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        plotxy(beam1, 1, 3, title="Image 1", nbins=101, nolost=1)
        plotxy(mirr1, 2, 1, title="Footprint 1", nbins=101, nolost=1)


def example_branch_4(do_plot=True, f_refl=0):


    #
    # source
    #
    # beam0 = Beam.initialize_as_pencil(N=500)
    source = SourceGaussian(nrays=100000,
                 sigmaX=0.0,
                 sigmaY=0.0,
                 sigmaZ=0.0,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,)
    beam0 = S4Beam()
    numpy.random.seed(123456)
    beam0.generate_source(source)
    beam0.set_photon_wavelength(5e-10)
    print(beam0.info())

    if do_plot:
        plotxy(beam0, 4, 6, title="Image 0", nbins=201)

    #
    # syned definitopns
    #

    # surface shape
    # boundaries
    rlen1 = 5e-05
    rlen2 = 5e-05
    rwidx1 = 2e-05
    rwidx2 = 2e-05

    boundary_shape = Rectangle(x_left=-rwidx2,x_right=rwidx1,y_bottom=-rlen2,y_top=rlen1)

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 6.0,
                                           angle_radial = numpy.radians(88.8))


    #
    # shadow definitions
    #
    if f_refl == 0: # prerefl

        PreRefl.prerefl(interactive=False, SYMBOL="SiC", DENSITY=3.217, FILE="SiC.dat",
                        E_MIN=100.0, E_MAX=20000.0, E_STEP=100.0)
        mirror1 = S4ConicMirrorElement(optical_element=S4ConicMirror(name="M1",
                                                                     conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
                                                                     boundary_shape=boundary_shape,
                                                                     f_reflec=1, f_refl=f_refl, file_refl="SiC.dat"),
                                       coordinates=coordinates_syned,
                                       input_beam=beam0)
    elif f_refl == 1: # refraction index
        import xraylib
        refraction_index = xraylib.Refractive_Index("SiC", 2.4797, 3.217)
        mirror1 = S4ConicMirrorElement(optical_element=S4ConicMirror(name="M1",
                                                                     conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
                                                                     boundary_shape=boundary_shape,
                                                                     f_reflec=1, f_refl=f_refl, file_refl="", refraction_index=refraction_index),
                                       coordinates=coordinates_syned,
                                       input_beam=beam0)
    elif f_refl == 2:  # user file: 1D  vs angle
        mirror1 = S4ConicMirrorElement(optical_element=S4ConicMirror(name="M1",
                                                                     conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
                                                                     boundary_shape=boundary_shape,
                                                                     f_reflec=1, f_refl=f_refl, file_refl="xoppy_f1f2_139980555361648.dat"),
                                       coordinates=coordinates_syned,
                                       input_beam=beam0)
    elif f_refl == 3:  # user file 1D vs energy
        mirror1 = S4ConicMirrorElement(optical_element=S4ConicMirror(name="M1",
                                                                     conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
                                                                     boundary_shape=boundary_shape,
                                                                     f_reflec=1, f_refl=f_refl, file_refl="xoppy_f1f2_139981943656272.dat"),
                                       coordinates=coordinates_syned,
                                       input_beam=beam0)
    elif f_refl == 4:  # user file
        mirror1 = S4ConicMirrorElement(optical_element=S4ConicMirror(name="M1",
                                                           conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
                                                           boundary_shape=boundary_shape,
                                                           f_reflec=1, f_refl=f_refl, file_refl="xoppy_f1f2_139980938100080.dat"),
                                       coordinates=coordinates_syned,
                                       input_beam=beam0)

    print(mirror1.info())

    #
    # run
    #
    beam1, mirr1 = mirror1.trace_beam()
    print(mirr1.info())

    #
    # check
    #

    if do_plot:
        plotxy(beam1, 1, 3, title="Image 1", nbins=101, nolost=1)
        plotxy(mirr1, 2, 1, title="Footprint 1", nbins=101, nolost=1)

def example_branch_5(surface_type, do_plot=True):
    #
    # source
    #
    source = SourceGaussian(nrays=100000,
                 sigmaX=0.0,
                 sigmaY=0.0,
                 sigmaZ=0.0,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,)
    beam0 = S4Beam()
    beam0.generate_source(source)
    # print(beam0.info())



    #
    # syned definitopns
    #
    # boundaries
    boundary_shape = None

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 10.0,
                                           angle_radial = numpy.radians(88.8))

    # surface shape

    if surface_type == "plane":
        mirror1 = S4PlaneMirrorElement(optical_element=S4PlaneMirror(name="M1",
                                                                     boundary_shape=boundary_shape),
                                       coordinates=coordinates_syned,
                                       input_beam=beam0)
    elif surface_type == "sphere":
        mirror1 = S4SphereMirrorElement(optical_element=S4SphereMirror(name="M1",
                                                                       boundary_shape=boundary_shape,
                                                                       is_cylinder=False,
                                                                       surface_calculation=SurfaceCalculation.INTERNAL,
                                                                       p_focus=10.0,
                                                                       q_focus=10.0,
                                                                       grazing_angle=numpy.radians(90.0 - 88.8)),
                                       coordinates=coordinates_syned,
                                       input_beam=beam0)
    elif surface_type == "spherical_cylinder_tangential":
        mirror1 = S4SphereMirrorElement(optical_element=S4SphereMirror(name="M1",
                                                                       boundary_shape=boundary_shape,
                                                                       is_cylinder=True,
                                                                       cylinder_direction=Direction.TANGENTIAL,
                                                                       surface_calculation=SurfaceCalculation.INTERNAL,
                                                                       p_focus=10.0,
                                                                       q_focus=10.0,
                                                                       grazing_angle=numpy.radians(90.0 - 88.8)),
                                       coordinates=coordinates_syned,
                                       input_beam=beam0)
    elif surface_type == "spherical_cylinder_sagittal":
        mirror1 = S4SphereMirrorElement(optical_element=S4SphereMirror(name="M1",
                                                                       boundary_shape=boundary_shape,
                                                                       is_cylinder=True,
                                                                       cylinder_direction=Direction.SAGITTAL,
                                                                       surface_calculation=SurfaceCalculation.INTERNAL,
                                                                       p_focus=10.0,
                                                                       q_focus=10.0,
                                                                       grazing_angle=numpy.radians(90.0 - 88.8)),
                                       coordinates=coordinates_syned,
                                       input_beam=beam0)
    elif surface_type == "ellipsoid":
        mirror1 = S4EllipsoidMirrorElement(optical_element=S4EllipsoidMirror(name="M1",
                                                                             boundary_shape=boundary_shape,
                                                                             is_cylinder=False,
                                                                             surface_calculation=SurfaceCalculation.INTERNAL,
                                                                             p_focus=20.0,
                                                                             q_focus=10.0,
                                                                             grazing_angle=0.003),
                                          coordinates=coordinates_syned,
                                          input_beam=beam0)
    elif surface_type == "elliptical_cylinder":
        mirror1 = S4EllipsoidMirrorElement(optical_element=S4EllipsoidMirror(name="M1",
                                                                             boundary_shape=boundary_shape,
                                                                             is_cylinder=True,
                                                                             cylinder_direction=Direction.TANGENTIAL,
                                                                             surface_calculation=SurfaceCalculation.INTERNAL,
                                                                             p_focus=20.0,
                                                                             q_focus=10.0,
                                                                             grazing_angle=0.003),
                                          coordinates=coordinates_syned,
                                          input_beam=beam0)
    elif surface_type == "hyperboloid":
        mirror1 = S4HyperboloidMirrorElement(optical_element=S4HyperboloidMirror(name="M1",
                                                                                 boundary_shape=boundary_shape,
                                                                                 is_cylinder=False,
                                                                                 surface_calculation=SurfaceCalculation.INTERNAL,
                                                                                 p_focus=20.0,
                                                                                 q_focus=10.0,
                                                                                 grazing_angle=0.003),
                                             coordinates=coordinates_syned,
                                             input_beam=beam0)
    elif surface_type == "hyperbolic_cylinder":
        mirror1 = S4HyperboloidMirrorElement(optical_element=S4HyperboloidMirror(name="M1",
                                                                                 boundary_shape=boundary_shape,
                                                                                 is_cylinder=True,
                                                                                 cylinder_direction=Direction.TANGENTIAL,
                                                                                 surface_calculation=SurfaceCalculation.INTERNAL,
                                                                                 p_focus=20.0,
                                                                                 q_focus=10.0,
                                                                                 grazing_angle=0.003),
                                             coordinates=coordinates_syned,
                                             input_beam=beam0)
    elif surface_type == "paraboloid":
        mirror1 = S4ParaboloidMirrorElement(optical_element=S4ParaboloidMirror(name="M1",
                                                                                 boundary_shape=boundary_shape,
                                                                                 is_cylinder=False,
                                                                                 surface_calculation=SurfaceCalculation.INTERNAL,
                                                                                 at_infinity=Side.SOURCE,
                                                                                 p_focus=20.0,
                                                                                 q_focus=10.0,
                                                                                 grazing_angle=0.003),
                                             coordinates=coordinates_syned,
                                             input_beam=beam0)
    elif surface_type == "parabolic_cylinder":
        mirror1 = S4ParaboloidMirrorElement(optical_element=S4ParaboloidMirror(name="M1",
                                                                               boundary_shape=boundary_shape,
                                                                               is_cylinder=True,
                                                                               cylinder_direction=Direction.TANGENTIAL,
                                                                               surface_calculation=SurfaceCalculation.INTERNAL,
                                                                               at_infinity=Side.SOURCE,
                                                                               p_focus=20.0,
                                                                               q_focus=10.0,
                                                                               grazing_angle=0.003),
                                            coordinates=coordinates_syned,
                                       input_beam=beam0)
    elif surface_type == "conic":
        mirror1 = S4ConicMirrorElement(optical_element=S4ConicMirror(name="M1",
                                                                     boundary_shape=boundary_shape,
                                                                     conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0]),
                                       coordinates=coordinates_syned,
                                       input_beam=beam0)
    else:
        raise Exception("undefined surface shape")


    #
    # run
    #
    print(mirror1.info())
    beam1, mirr1 = mirror1.trace_beam()


    #
    # check
    #

    if do_plot:
        plotxy(beam1, 1, 3, nbins=101, nolost=1, title=surface_type)
        # plotxy(mirr1, 2, 1, title="Footprint 1", nbins=101, nolost=1)


if __name__ == "__main__":

    # see corresponding Oasys workspace in shadow4/oasys_workspaces

    from srxraylib.plot.gol import set_qt
    set_qt()

    do_plot = True

    example_branch_1(do_plot=do_plot) # two plane mirrors
    example_branch_2(do_plot=do_plot) # toroid
    example_branch_3("test_shadow4.hdf5",do_plot=do_plot) # mesh #TODO: remote acces
    example_branch_4(do_plot=do_plot, f_refl=0) # prerefl
    example_branch_4(do_plot=do_plot, f_refl=1)  # refraction index
    example_branch_4(do_plot=do_plot, f_refl=2)  # user file 1D, angle[mrad], reflectivity
    example_branch_4(do_plot=do_plot, f_refl=3)  # user file 1D, angle[mrad], reflectivity
    example_branch_4(do_plot=do_plot, f_refl=4)  # user file 2D, energy [eV], angle[mrad], reflectivity

    for myconicshape in ["plane", \
                         "sphere", "spherical_cylinder_tangential", "spherical_cylinder_sagittal", \
                         "ellipsoid", "elliptical_cylinder", \
                         "hyperboloid","hyperbolic_cylinder",\
                         "paraboloid","parabolic_cylinder"]:
        example_branch_5(myconicshape,do_plot=do_plot) # conic mirrors

