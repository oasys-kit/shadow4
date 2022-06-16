import numpy
import os

from shadow4.syned.shape import Direction, Convexity
from shadow4.syned.shape import Sphere, SphericalCylinder
from shadow4.syned.shape import Ellipsoid, EllipticalCylinder
from shadow4.syned.shape import Hyperboloid, HyperbolicCylinder
from shadow4.syned.shape import Paraboloid, ParabolicCylinder
from shadow4.syned.shape import Toroid, Conic, SurfaceData, Plane, Side

from shadow4.optical_surfaces.s4_conic import S4Conic
from shadow4.optical_surfaces.s4_mesh import S4Mesh
from shadow4.optical_surfaces.s4_toroid import S4Toroid

class S4OpticalElement(object):

    def __init__(self):
        pass

    def info(self):
        raise NotImplementedError()

    def to_python_code(self, data=None):
        raise NotImplementedError()

    def get_surface_shape_instance(self):
        raise NotImplementedError()

    def get_optical_surface_instance(self):
        raise NotImplementedError()

class SurfaceCalculation:
    INTERNAL = 0
    EXTERNAL = 1


class S4PlaneOpticalElement(S4OpticalElement):
    def __init__(self):
        S4OpticalElement.__init__(self)

        self._plane_surface_shape = Plane()

    def get_surface_shape_instance(self):
        return self._plane_surface_shape

    def get_optical_surface_instance(self):
        return S4Conic.initialize_as_plane()

class S4CurvedOpticalElement(S4OpticalElement):

    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 ):
        S4OpticalElement.__init__(self)

        self._surface_calculation = surface_calculation
        self._is_cylinder = is_cylinder
        self._curved_surface_shape = None # auxiliary field

    def get_curved_surface_shape(self):
        return self._curved_surface_shape

    def get_surface_shape_instance(self):
        return self._curved_surface_shape

##################################################
# CURVED OPTICAL ELEMENTS
##################################################

class S4SphereOpticalElement(S4CurvedOpticalElement):
    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 radius=0.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):
        S4CurvedOpticalElement.__init__(self, surface_calculation, is_cylinder)

        if self._surface_calculation == SurfaceCalculation.EXTERNAL:
            if self._is_cylinder: self._curved_surface_shape = SphericalCylinder.create_spherical_cylinder_from_radius(radius, convexity, cylinder_direction)
            else:                 self._curved_surface_shape = Sphere.create_sphere_from_radius(radius, convexity)
        else:
            if self._is_cylinder: self._curved_surface_shape = SphericalCylinder.create_spherical_cylinder_from_p_q(p_focus, q_focus, grazing_angle, convexity, cylinder_direction)
            else:                 self._curved_surface_shape = Sphere.create_sphere_from_p_q(p_focus, q_focus, grazing_angle, convexity)

    def get_optical_surface_instance(self):
        surface_shape = self.get_surface_shape()

        switch_convexity = 0 if surface_shape.get_convexity() == Convexity.DOWNWARD else 1

        radius = surface_shape.get_radius()
        if isinstance(surface_shape, SphericalCylinder):
            cylindrical = 1
            cylangle = 0.0 if surface_shape.get_cylinder_direction() == Direction.TANGENTIAL else (0.5 * numpy.pi)

        elif isinstance(surface_shape, Sphere):
            cylindrical = 0
            cylangle    = 0.0

        print(">>>>> S4SphereOpticalElement.get_optical_surface_instance(): R, cyl, cyl_angle, optical element, ", radius, cylindrical, cylangle, surface_shape)

        return S4Conic.initialize_as_sphere_from_curvature_radius(radius, cylindrical=cylindrical, cylangle=cylangle,
                                                                  switch_convexity=switch_convexity)



class S4EllipsoidOpticalElement(S4CurvedOpticalElement):

    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 min_axis=0.0,
                 maj_axis=0.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):
        S4CurvedOpticalElement.__init__(self, surface_calculation, is_cylinder)

        if self._surface_calculation == SurfaceCalculation.EXTERNAL:
            if self._is_cylinder: self._curved_surface_shape = EllipticalCylinder.create_elliptical_cylinder_from_axes(min_axis, maj_axis, p_focus, convexity, cylinder_direction)
            else:                 self._curved_surface_shape = Ellipsoid.create_ellipsoid_from_axes(min_axis, maj_axis, p_focus, convexity)
        else:
            if self._is_cylinder: self._curved_surface_shape = EllipticalCylinder.create_elliptical_cylinder_from_p_q(p_focus, q_focus, grazing_angle, convexity, cylinder_direction)
            else:                 self._curved_surface_shape = Ellipsoid.create_ellipsoid_from_p_q(p_focus, q_focus, grazing_angle, convexity)

    def get_optical_surface_instance(self):
        surface_shape = self.get_surface_shape()

        switch_convexity = 0 if surface_shape.get_convexity() == Convexity.DOWNWARD else 1

        if isinstance(surface_shape, EllipticalCylinder):
            print(">>>>> EllipticalCylinder optical element", surface_shape)
            cylindrical = 1
            cylangle = 0.0 if surface_shape.get_cylinder_direction() == Direction.TANGENTIAL else (0.5 * numpy.pi)
        elif isinstance(surface_shape, Ellipsoid):
            print(">>>>> Ellipsoid optical element", surface_shape)
            cylindrical = 0
            cylangle    = 0.0

        return S4Conic.initialize_as_ellipsoid_from_focal_distances(surface_shape.get_p_focus(),
                                                                    surface_shape.get_q_focus(),
                                                                    surface_shape.get_grazing_angle(),
                                                                    cylindrical=cylindrical,
                                                                    cylangle=cylangle,
                                                                    switch_convexity=switch_convexity)
class S4HyperboloidOpticalElement(S4CurvedOpticalElement):
    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 min_axis=0.0,
                 maj_axis=0.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):
        S4CurvedOpticalElement.__init__(self, surface_calculation, is_cylinder)

        if self._surface_calculation == SurfaceCalculation.EXTERNAL:
            if self._is_cylinder: self._curved_surface_shape = HyperbolicCylinder.create_hyperbolic_cylinder_from_axes(min_axis, maj_axis, p_focus, convexity, cylinder_direction)
            else:                 self._curved_surface_shape = Hyperboloid.create_hyperboloid_from_axes(min_axis, maj_axis, p_focus, convexity)
        else:
            if self._is_cylinder: self._curved_surface_shape = HyperbolicCylinder.create_hyperbolic_cylinder_from_p_q(p_focus, q_focus, grazing_angle, convexity, cylinder_direction)
            else:                 self._curved_surface_shape = Hyperboloid.create_hyperboloid_from_p_q(p_focus, q_focus, grazing_angle, convexity)

    def get_optical_surface_instance(self):
        surface_shape = self.get_surface_shape()

        switch_convexity = 0 if surface_shape.get_convexity() == Convexity.DOWNWARD else 1

        if isinstance(surface_shape, HyperbolicCylinder):
            print(">>>>> HyperbolicCylinder optical element", surface_shape)
            cylindrical = 1
            cylangle = 0.0 if surface_shape.get_cylinder_direction() == Direction.TANGENTIAL else (0.5 * numpy.pi)
        elif isinstance(surface_shape, Hyperboloid):
            print(">>>>> Hyperboloid optical element", surface_shape)
            cylindrical = 0
            cylangle    = 0.0

        return S4Conic.initialize_as_hyperboloid_from_focal_distances(surface_shape.get_p_focus(),
                                                                      surface_shape.get_q_focus(),
                                                                      surface_shape.get_grazing_angle(),
                                                                      cylindrical=cylindrical,
                                                                      cylangle=cylangle,
                                                                      switch_convexity=switch_convexity)

class S4ToroidalOpticalElement(S4CurvedOpticalElement):
    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 min_radius=0.0,
                 maj_radius=0.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):
        S4CurvedOpticalElement.__init__(self, surface_calculation, False)

        if self._surface_calculation == SurfaceCalculation.EXTERNAL:
            self._curved_surface_shape = Toroid.create_toroid_from_radii(min_radius, maj_radius)
        else:
            self._curved_surface_shape = Toroid.create_toroid_from_p_q(p_focus, q_focus, grazing_angle)

    def get_optical_surface_instance(self):
        surface_shape = self.get_surface_shape()

        print(">>>>> Toroidal optical element", surface_shape._min_radius, surface_shape._maj_radius)

        toroid = S4Toroid()
        toroid.set_toroid_radii(surface_shape.get_maj_radius(), surface_shape.get_min_radius())
        return toroid

class S4ConicOpticalElement(S4CurvedOpticalElement):
    def __init__(self, conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]):
        S4CurvedOpticalElement.__init__(self, surface_calculation=SurfaceCalculation.INTERNAL, is_cylinder=False)

        self._curved_surface_shape = Conic(conic_coefficients=conic_coefficients)

    def get_optical_surface_instance(self):
        surface_shape = self.get_surface_shape()
        print(">>>>> Conic optical element")
        return S4Conic.initialize_from_coefficients(surface_shape.get_conic_coefficients())

class S4ParaboloidOpticalElement(S4CurvedOpticalElement):
    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 parabola_parameter=0.0,
                 at_infinity=Side.SOURCE,
                 pole_to_focus=None,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):
        S4CurvedOpticalElement.__init__(self, surface_calculation, is_cylinder)

        if self._surface_calculation == SurfaceCalculation.EXTERNAL:
            if self._is_cylinder: self._curved_surface_shape = ParabolicCylinder.create_parabolic_cylinder_from_parabola_parameter(parabola_parameter, at_infinity, pole_to_focus, convexity, cylinder_direction)
            else:                 self._curved_surface_shape = Paraboloid.create_paraboloid_from_parabola_parameter(parabola_parameter, at_infinity, pole_to_focus, convexity)
        else:
            if self._is_cylinder: self._curved_surface_shape = ParabolicCylinder.create_parabolic_cylinder_from_p_q(p_focus, q_focus, grazing_angle, at_infinity, convexity, cylinder_direction)
            else:                 self._curved_surface_shape = Paraboloid.create_paraboloid_from_p_q(p_focus, q_focus, grazing_angle, at_infinity, convexity)

    def get_optical_surface_instance(self):
        surface_shape = self.get_surface_shape()

        switch_convexity = 0 if surface_shape.get_convexity() == Convexity.DOWNWARD else 1

        if surface_shape.get_at_infinity() == Side.SOURCE:
            p = 1e20
            q = surface_shape.get_pole_to_focus()
        else:
            q = 1e20
            p = surface_shape.get_pole_to_focus()

        if isinstance(surface_shape, ParabolicCylinder):
            print(">>>>> ParabolicCylinder optical element", surface_shape)
            cylindrical = 1
            cylangle = 0.0 if surface_shape.get_cylinder_direction() == Direction.TANGENTIAL else (0.5 * numpy.pi)
        elif isinstance(surface_shape, Paraboloid):
            print(">>>>> Paraboloid optical element", surface_shape)
            cylindrical = 0
            cylangle    = 0.0

        return S4Conic.initialize_as_paraboloid_from_focal_distances(p, q, surface_shape.get_grazing_angle(), cylindrical=cylindrical, cylangle=cylangle, switch_convexity=switch_convexity)


class S4SurfaceDataOpticalElement(S4CurvedOpticalElement):
    def __init__(self, xx=None, yy=None, zz=None, surface_data_file=None):
        S4CurvedOpticalElement.__init__(self, surface_calculation=SurfaceCalculation.INTERNAL, is_cylinder=False)

        self._curved_surface_shape = SurfaceData(xx, yy, zz, surface_data_file)

    def get_optical_surface_instance(self):
        surface_shape = self.get_surface_shape()

        print(">>>>> SurfaceData optical element")
        num_mesh = S4Mesh()

        if surface_shape.has_surface_data():
            num_mesh.load_surface_data(surface_shape)
        elif surface_shape.has_surface_data_file():
            filename, file_extension = os.path.splitext(surface_shape._surface_data_file)

            if file_extension.lower() in [".h5", ".hdf", ".hdf5"]:
                num_mesh.load_h5file(surface_shape._surface_data_file)
            else:
                num_mesh.load_file(surface_shape._surface_data_file) # 3 columns ASCII

        return num_mesh

class S4AdditiveSurfaceDataOpticalElement(S4CurvedOpticalElement):
    def __init__(self, xx=None, yy=None, zz=None, surface_data_file=None, base_surface_function=None):
        S4CurvedOpticalElement.__init__(self, surface_calculation=SurfaceCalculation.INTERNAL, is_cylinder=False)

        self._curved_surface_shape = SurfaceData(xx, yy, zz, surface_data_file)

        if base_surface_function is None:
            self._base_surface_function = lambda x,y: x * 0.0
        else:
            self._base_surface_function = base_surface_function

    def get_optical_surface_instance(self):
        surface_shape = self.get_surface_shape()

        print(">>>>> AdditiveSurfaceData optical element")
        num_mesh = S4Mesh()

        if surface_shape.has_surface_data():
            num_mesh.load_surface_data(surface_shape)
        elif surface_shape.has_surface_data_file():
            filename, file_extension = os.path.splitext(surface_shape._surface_data_file)

            if file_extension.lower() in [".h5", ".hdf", ".hdf5"]:
                num_mesh.load_h5file(surface_shape._surface_data_file)
            else:
                num_mesh.load_file(surface_shape._surface_data_file) # 3 columns ASCII

        # TODO: add mesh with evaluated surface
        x, y = num_mesh.get_mesh_x_y()
        z = num_mesh.get_mesh_z()
        Y = numpy.outer(numpy.ones_like(x), y)
        X = numpy.outer(x, numpy.ones_like(y))

        zadd = self._base_surface_function(X, Y)
        print(">>>>", X.shape, Y.shape, zadd.shape)
        # from srxraylib.plot.gol import plot_surface
        # plot_surface(zadd, x, y, xtitle="x")
        num_mesh.add_to_mesh( zadd )
        return num_mesh