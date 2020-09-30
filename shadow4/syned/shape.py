
import numpy
from syned.syned_object import SynedObject
from collections import OrderedDict

class Convexity:
    NONE = -1
    UPWARD = 0
    DOWNWARD = 1

class Direction:
    TANGENTIAL = 0
    SAGITTAL = 1

class Side:
    SOURCE = 0
    IMAGE = 1
#
# main shape subclasses:
#      SurfaceShape to caracterize the shape (sphere etc.) of the optical element surface
#      BoundaryShape to characterize the optical element dimensions (rectangle, etc.)
#
class Shape(SynedObject):
    def __init__(self):
        SynedObject.__init__(self)

class SurfaceShape(Shape):
    def __init__(self, convexity = Convexity.UPWARD):
        Shape.__init__(self)

        self._convexity = convexity

    def get_convexity(self):
        return self._convexity

class BoundaryShape(Shape):
    def __init__(self):
        Shape.__init__(self)
        
    def get_boundaries(self):
        raise NotImplementedError()

#
# Subclasses for SurfaceShape
#

class Cylinder:
    def __init__(self, cylinder_direction=Direction.TANGENTIAL):
        self._cylinder_direction = cylinder_direction

    def get_cylinder_direction(self):
        return self._cylinder_direction

class Conic(SurfaceShape):
    def __init__(self, 
                 conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]):
        SurfaceShape.__init__(self, convexity=Convexity.NONE)

        self._conic_coefficients = conic_coefficients

    def get_conic_coefficients(self):
        return self._conic_coefficients

class Plane(SurfaceShape):
    def __init__(self):
        SurfaceShape.__init__(self, convexity=Convexity.NONE)

class Sphere(SurfaceShape):
    def __init__(self, radius=1.0, convexity=Convexity.UPWARD):
        SurfaceShape.__init__(self, convexity=convexity)
        self._radius = radius

    @classmethod
    def create_sphere_from_radius(cls, radius=0.0, convexity=Convexity.UPWARD):
        return Sphere(radius, convexity)

    @classmethod
    def create_sphere_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, convexity=Convexity.UPWARD):
        sphere = Sphere(convexity=convexity)
        sphere.initialize_from_p_q(p, q, grazing_angle)

        return sphere

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        self._radius = Sphere.get_radius_from_p_q(p, q, grazing_angle)

    @classmethod
    def get_radius_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003):
        # 1/p + 1/q = 2/(R cos(pi/2 - gr.a.))
        return (2*p*q/(p+q))/numpy.sin(grazing_angle)

    def get_radius(self):
        return self._radius


class SphericalCylinder(Sphere, Cylinder):
    def __init__(self, 
                 radius=1.0, 
                 convexity=Convexity.UPWARD, 
                 cylinder_direction=Direction.TANGENTIAL):
        Sphere.__init__(self, radius, convexity)
        Cylinder.__init__(self, cylinder_direction)

    @classmethod
    def create_spherical_cylinder_from_radius(cls, radius=0.0, convexity=Convexity.UPWARD, cylinder_direction=Direction.TANGENTIAL):
        return SphericalCylinder(radius, convexity, cylinder_direction)

    @classmethod
    def create_spherical_cylinder_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, convexity=Convexity.UPWARD, cylinder_direction=Direction.TANGENTIAL):
        spherical_cylinder = SphericalCylinder(convexity=convexity, cylinder_direction=cylinder_direction)
        spherical_cylinder.initialize_from_p_q(p, q, grazing_angle)

        return spherical_cylinder

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        if self._cylinder_direction == Direction.TANGENTIAL:
            self._radius = Sphere.get_radius_from_p_q(p, q, grazing_angle)
        elif self._cylinder_direction == Direction.SAGITTAL:
            self._radius = SphericalCylinder.get_radius_from_p_q_sagittal(p, q, grazing_angle)

    @classmethod
    def get_radius_from_p_q_sagittal(cls, p=2.0, q=1.0, grazing_angle=0.003):
        # 1/p + 1/q = 2 cos(pi/2 - gr.a.)/r
        return (2*p*q/(p+q))*numpy.sin(grazing_angle)

class Ellipsoid(SurfaceShape):
    """
    Ellipsoid: Revolution ellipsoid (rotation around major axis).
    It is defined with three parameters: axes of the ellipse and an additional parameter
    defining the position of the origin of the mirror. This additional parameter can be "p", "x0", "y0"
    or the angle beta from the ellipsoid center (tan(beta)=y0/x0). For simplicity, we store "p" in syned.
    """
    def __init__(self, min_axis=0.0, maj_axis=0.0, p_focus=0.0, convexity=Convexity.UPWARD):
        SurfaceShape.__init__(self, convexity)

        self._min_axis = min_axis
        self._maj_axis = maj_axis
        self._p_focus  = p_focus


    @classmethod
    def create_ellipsoid_from_axes(cls, min_axis=0.0, maj_axis=0.0, p_focus=0.0, convexity=Convexity.UPWARD):
        return Ellipsoid(min_axis, maj_axis, p_focus, convexity)

    @classmethod
    def create_ellipsoid_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, convexity=Convexity.UPWARD):
        ellipsoid = Ellipsoid(convexity=convexity)
        ellipsoid.initialize_from_p_q(p, q, grazing_angle)

        return ellipsoid

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        self._min_axis, self._maj_axis = Ellipsoid.get_axis_from_p_q(p, q, grazing_angle)
        self._p_focus = p

    def initialize_from_shadow_parameters(self, axmaj=2.0, axmin=1.0, ell_the=0.003, convexity=Convexity.UPWARD):
        tanbeta2 = numpy.tan(ell_the) ** 2
        y = axmaj * axmin / numpy.sqrt(axmin ** 2 + axmaj ** 2 * tanbeta2)
        z = y * numpy.tan(ell_the)
        c = numpy.sqrt(axmaj ** 2 - axmin ** 2)
        p = numpy.sqrt( (y + c)**2 + z**2)

        self.__init__(axmin, axmaj, p, convexity)

    def get_axes(self):
        return self._min_axis, self._maj_axis

    def get_p_q(self, grazing_angle=0.003):
        return Ellipsoid.get_p_q_from_axis(self._min_axis, self._maj_axis, grazing_angle)

    # semiaxes etc
    def get_a(self):
        return 0.5 * self._maj_axis

    def get_b(self):
        return 0.5 * self._min_axis

    def get_c(self):
        return numpy.sqrt(self.get_a()**2 - self.get_b()**2)

    def get_p_focus(self):
        return self._p_focus

    def get_q_focus(self):
        return 2 * self.get_a() - self.get_p_focus()

    def get_eccentricity(self):
        return self.get_c() / self.get_a()

    def get_grazing_angle(self):
        return numpy.arcsin(self.get_b() / numpy.sqrt(self.get_p_focus() * self.get_q_focus()))

    def get_mirror_center(self):
        coor_along_axis_maj = (self.get_p_focus()**2 - self.get_q_focus()**1) / (4 * self.get_c())
        coor_along_axis_min = self.get_b * numpy.sqrt(1 - (coor_along_axis_maj / self.get_a())**2)
        return coor_along_axis_maj, coor_along_axis_min

    def get_angle_pole_from_origin(self):
        x1, x2 = self.get_mirror_center()
        return numpy.arctan(x2 / x1)

    @classmethod
    def get_axis_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003):
        # see calculation of ellipse axis in shadow_kernel.f90 row 3605
        min_axis = 2*numpy.sqrt(p*q)*numpy.sin(grazing_angle)
        maj_axis = (p + q)

        return min_axis, maj_axis

    @classmethod
    def get_p_q_from_axis(cls, min_axis=2.0, maj_axis=1.0, grazing_angle=0.003):
        a = maj_axis/2
        b = min_axis/2
        p = a + numpy.sqrt(a**2 - (b/numpy.sin(grazing_angle))**2)
        q = maj_axis - p

        return p, q



class EllipticalCylinder(Ellipsoid, Cylinder):
    def __init__(self, 
                 min_axis=0.0, 
                 maj_axis=0.0, 
                 p_focus=0.0,
                 convexity=Convexity.UPWARD,
                 cylinder_direction=Direction.TANGENTIAL):
        Ellipsoid.__init__(self, min_axis, maj_axis, p_focus, convexity)
        Cylinder.__init__(self, cylinder_direction)

    @classmethod
    def create_elliptical_cylinder_from_axes(cls, min_axis=0.0, maj_axis=0.0, p_focus=0.0, convexity=Convexity.UPWARD, cylinder_direction=Direction.TANGENTIAL):
        return EllipticalCylinder(min_axis, maj_axis, p_focus, convexity, cylinder_direction)

    @classmethod
    def create_elliptical_cylinder_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, convexity=Convexity.UPWARD, cylinder_direction=Direction.TANGENTIAL):
        elliptical_cylinder = EllipticalCylinder(convexity=convexity, cylinder_direction=cylinder_direction)
        elliptical_cylinder.initialize_from_p_q(p, q, grazing_angle)

        return elliptical_cylinder

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        if self._cylinder_direction == Direction.SAGITTAL: raise NotImplementedError("Operation not possible for SAGITTAL direction")

        super().initialize_from_p_q(p, q, grazing_angle)

    def get_p_q(self, grazing_angle=0.003):
        if self._cylinder_direction == Direction.SAGITTAL: raise NotImplementedError("Operation not possible for SAGITTAL direction")

        return super().get_p_q(grazing_angle)

class Hyperboloid(SurfaceShape):
    """
    Hyperboloid: Revolution hyperboloid (two sheets: rotation around major axis).
    It is defined with three parameters: axes of the hyperbola and an additional parameter
    defining the position of the origin of the mirror. This additional parameter can be "p", "x0", "y0"
    or the angle beta from the ellipsoid center (tan(beta)=y0/x0). For simplicity, we store "p" in syned.
    """
    def __init__(self, min_axis=0.0, maj_axis=0.0, p_focus=0.0, convexity=Convexity.UPWARD):
        SurfaceShape.__init__(self, convexity)

        self._min_axis = min_axis
        self._maj_axis = maj_axis
        self._p_focus  = p_focus

    @classmethod
    def create_hyperboloid_from_axes(cls, min_axis=0.0, maj_axis=0.0, p_focus=0.0, convexity=Convexity.UPWARD):
        return Hyperboloid(min_axis, maj_axis, p_focus, convexity)

    @classmethod
    def create_hyperboloid_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, convexity=Convexity.UPWARD):
        hyperboloid = Hyperboloid(convexity=convexity)
        hyperboloid.initialize_from_p_q(p, q, grazing_angle)

        return hyperboloid

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        self._min_axis, self._maj_axis = Hyperboloid.get_axis_from_p_q(p, q, grazing_angle)
        self._p_focus = p

    # TODO:
    def initialize_from_shadow_parameters(self, axmaj=2.0, axmin=1.0, ell_the=0.003, convexity=Convexity.UPWARD):
        raise NotImplementedError("TODO")

    def get_axes(self):
        return self._min_axis, self._maj_axis

    def get_p_q(self, grazing_angle=0.003):
        return Hyperboloid.get_p_q_from_axis(self._min_axis, self._maj_axis, grazing_angle)

    # semiaxes etc
    def get_a(self):
        return 0.5 * self._maj_axis

    def get_b(self):
        return 0.5 * self._min_axis

    def get_c(self):
        return numpy.sqrt(self.get_a()**2 + self.get_b()**2)

    def get_p_focus(self):
        return self._p_focus

    def get_q_focus(self):
        return self.get_p_focus() - 2 * self.get_a()

    def get_eccentricity(self):
        return self.get_c / self.get_a()

    def get_grazing_angle(self):
        return numpy.arcsin(self.get_b() / numpy.sqrt(self.get_p_focus() * self.get_q_focus()))

    #TODO:
    def get_mirror_center(self):
        raise NotImplementedError("TODO")

    def get_angle_pole_from_origin(self):
        raise NotImplementedError("TODO")

    @classmethod
    def get_axis_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, branch_sign=+1):
        min_axis = 2*numpy.sqrt(p*q)*numpy.sin(grazing_angle)
        maj_axis = (p - q) * branch_sign

        return min_axis, maj_axis

    # TODO:
    @classmethod
    def get_p_q_from_axis(cls, min_axis=2.0, maj_axis=1.0, grazing_angle=0.003):
        raise NotImplementedError("TODO")

class HyperbolicCylinder(Hyperboloid, Cylinder):
    def __init__(self, 
                 min_axis=0.0, 
                 maj_axis=0.0, 
                 p_focus=0.0,
                 convexity=Convexity.UPWARD, 
                 cylinder_direction=Direction.TANGENTIAL):
        Hyperboloid.__init__(self, min_axis, maj_axis, p_focus, convexity)
        Cylinder.__init__(self, cylinder_direction)


    @classmethod
    def create_hyperbolic_cylinder_from_axes(cls, min_axis=0.0, maj_axis=0.0, p_focus=0.0, convexity=Convexity.UPWARD, cylinder_direction=Direction.TANGENTIAL):
        return HyperbolicCylinder(min_axis, maj_axis, p_focus, convexity, cylinder_direction)

    @classmethod
    def create_hyperbolic_cylinder_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, convexity=Convexity.UPWARD, cylinder_direction=Direction.TANGENTIAL):
        hyperbolic_cylinder = HyperbolicCylinder(convexity=convexity, cylinder_direction=cylinder_direction)
        hyperbolic_cylinder.initialize_from_p_q(p, q, grazing_angle)

        return hyperbolic_cylinder

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        if self._cylinder_direction == Direction.SAGITTAL: raise NotImplementedError("Operation not possible for SAGITTAL direction")

        super().initialize_from_p_q(p, q, grazing_angle)

    def get_p_q(self, grazing_angle=0.003):
        if self._cylinder_direction == Direction.SAGITTAL: raise NotImplementedError("Operation not possible for SAGITTAL direction")

        return super().get_p_q(grazing_angle)

class Paraboloid(SurfaceShape):
    """
    Paraboloid: Revolution paraboloid (rotation around symmetry axis).
    It is defined with three parameters: the parabola_parameter and two more parameters
    defining the position of the origin of the mirror.
    The parabola_parameter = 2 * focal_distance = - 0.5 * ccc_9 / ccc_2
    The additional parameter can be the focal distances
    ("p" or "q", one is infinity), "x0", "y0" or the grazing angle.
    Here, we selected the at_infinity and the finite focal distance p or q or distance from
    the mirror pole to focus (pole to focus).

    """
    def __init__(self,
                 parabola_parameter=0.0,
                 at_infinity=Side.SOURCE,
                 pole_to_focus=None,
                 convexity=Convexity.UPWARD):
        SurfaceShape.__init__(self, convexity)

        self._parabola_parameter = parabola_parameter
        self._at_infinity = at_infinity
        self._pole_to_focus = pole_to_focus

    @classmethod
    def create_paraboloid_from_parabola_parameter(cls, parabola_parameter=0.0, at_infinity=Side.SOURCE, pole_to_focus=None, convexity=Convexity.UPWARD):
        return Paraboloid(parabola_parameter, at_infinity=at_infinity, pole_to_focus=pole_to_focus, convexity=convexity)

    @classmethod
    def create_paraboloid_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, at_infinity=Side.SOURCE, convexity=Convexity.UPWARD):
        paraboloid = Paraboloid(convexity=convexity)
        paraboloid.initialize_from_p_q(p, q, grazing_angle=grazing_angle, at_infinity=at_infinity)

        return paraboloid

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003, at_infinity=Side.SOURCE):
        self._parabola_parameter = Paraboloid.get_parabola_parameter_from_p_q(p=p, q=q, grazing_angle=grazing_angle, at_infinity=at_infinity)
        self._at_infinity = at_infinity
        if at_infinity == Side.SOURCE:
            self._pole_to_focus = q
        elif at_infinity == Side.IMAGE:
            self._pole_to_focus = p

    @classmethod
    def get_parabola_parameter_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, at_infinity=Side.SOURCE):
        if at_infinity == Side.IMAGE:
            return 2*p*(numpy.sin(grazing_angle))**2
        elif at_infinity == Side.SOURCE:
            return 2*q*(numpy.sin(grazing_angle))**2

    def get_parabola_parameter(self):
        return self._parabola_parameter

    def get_at_infinity(self):
        return self._at_infinity

    def get_pole_to_focus(self):
        return self._pole_to_focus

    def get_grazing_angle(self):
        return numpy.arcsin( numpy.sqrt( self.get_parabola_parameter() / (2 * self.get_pole_to_focus())))


class ParabolicCylinder(Paraboloid, Cylinder):
    def __init__(self,
                 parabola_parameter=0.0,
                 at_infinity=Side.SOURCE,
                 pole_to_focus=None,
                 convexity=Convexity.UPWARD,
                 cylinder_direction=Direction.TANGENTIAL):
        Paraboloid.__init__(self, parabola_parameter=parabola_parameter, at_infinity=at_infinity,
                            pole_to_focus=pole_to_focus, convexity=convexity)
        Cylinder.__init__(self, cylinder_direction)

    @classmethod
    def create_parabolic_cylinder_from_parabola_parameter(cls, parabola_parameter=0.0, at_infinity=Side.SOURCE, pole_to_focus=None, convexity=Convexity.UPWARD, cylinder_direction=Direction.TANGENTIAL):
        return ParabolicCylinder(parabola_parameter, at_infinity, pole_to_focus, convexity, cylinder_direction)

    @classmethod
    def create_parabolic_cylinder_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, at_infinity=Side.SOURCE, convexity=Convexity.UPWARD, cylinder_direction=Direction.TANGENTIAL):
        parabolic_cylinder = ParabolicCylinder(convexity=convexity, cylinder_direction=cylinder_direction)
        parabolic_cylinder.initialize_from_p_q(p, q, grazing_angle, at_infinity)

        return parabolic_cylinder

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003, at_infinity=Side.SOURCE):
        if self._cylinder_direction == Direction.SAGITTAL:
            raise NotImplementedError("Operation not possible for SAGITTAL direction")

        return super().initialize_from_p_q(p, q, grazing_angle, at_infinity)

# TODO: consider rename to Toroid?
class Toroid(SurfaceShape):
    def __init__(self, min_radius=0.0, maj_radius=0.0):
        SurfaceShape.__init__(self, convexity=Convexity.NONE)
        
        self._min_radius = min_radius
        self._maj_radius = maj_radius

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("min_radius"         , "Minor radius r   ", "m" ),
                    ("maj_radius"         , "Major radius R (optical=R+r)", "m" ),
            ] )

    @classmethod
    def create_toroid_from_radii(cls, min_radius=0.0, maj_radius=0.0):
        return Toroid(min_radius, maj_radius)

    @classmethod
    def create_toroid_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003):
        toroid = Toroid()
        toroid.initialize_from_p_q(p, q, grazing_angle)

        return toroid

    def get_radii(self):
        return self._min_radius, self._maj_radius

    def get_min_radius(self):
        return self._min_radius

    def get_maj_radius(self):
        return self._maj_radius

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        self._maj_radius = Sphere.get_radius_from_p_q(p, q, grazing_angle)
        self._min_radius = SphericalCylinder.get_radius_from_p_q_sagittal(p, q, grazing_angle)

        # FROM SHADOW3:
        #! C
        #! C NOTE : The major radius is the in reality the radius of the torus
        #! C max. circle. The true major radius is then
        #! C
        #        R_MAJ	=   R_MAJ - R_MIN
        self._maj_radius -= self._min_radius


# This is exactly the same as OasysSurfaceData
# TODO: consider moving this definition from oasys1 to syned
# class OasysSurfaceData(object):
class SurfaceData(SurfaceShape):
    def __init__(self,
                 xx=None,
                 yy=None,
                 zz=None,
                 surface_data_file=None):
        self._xx = xx
        self._yy = yy
        self._zz = zz
        self._surface_data_file=surface_data_file

    def has_surface_data(self):
        return not (self._xx is None or self._yy is None or self._zz is None)

    def has_surface_data_file(self):
        return not self._surface_data_file is None
#
# subclasses for BoundaryShape
#


class Rectangle(BoundaryShape):
    def __init__(self, x_left=-0.010, x_right=0.010, y_bottom=-0.020, y_top=0.020):
        super().__init__()

        self._x_left   = x_left
        self._x_right  = x_right
        self._y_bottom = y_bottom
        self._y_top    = y_top

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("x_left"          , "x (width) minimum (signed)   ", "m" ),
                    ("x_right"         , "x (width) maximum (signed)   ", "m" ),
                    ("y_bottom"        , "y (length) minimum (signed)  ", "m" ),
                    ("y_top"           , "y (length) maximum (signed)  ", "m" ),
            ] )

    def get_boundaries(self):
        return self._x_left, self._x_right, self._y_bottom, self._y_top

    def set_boundaries(self,x_left=-0.010, x_right=0.010, y_bottom=-0.020, y_top=0.020):
        self._x_left = x_left
        self._x_right = x_right
        self._y_bottom = y_bottom
        self._y_top = y_top

    def set_width_and_length(self,width=10e-3,length=30e-3):
        self._x_left = -0.5 * width
        self._x_right = 0.5 * width
        self._y_bottom = -0.5 * length
        self._y_top = 0.5 * length

class Ellipse(BoundaryShape):
    def __init__(self, a_axis_min=-10e-6, a_axis_max=10e-6, b_axis_min=-5e-6, b_axis_max=5e-6):
        super().__init__()

        self._a_axis_min   = a_axis_min
        self._a_axis_max  = a_axis_max
        self._b_axis_min = b_axis_min
        self._b_axis_max    = b_axis_max
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("a_axis_min"         , "x (width) axis starts (signed)  ", "m" ),
                    ("a_axis_max"        , "x (width) axis ends (signed)    ", "m" ),
                    ("b_axis_min"       , "y (length) axis starts (signed) ", "m" ),
                    ("b_axis_max"          , "y (length) axis ends (signed)   ", "m" ),
            ] )

    def get_boundaries(self):
        return self._a_axis_min, self._a_axis_max, self._b_axis_min, self._b_axis_max

    def get_axis(self):
        return numpy.abs(self._a_axis_max - self._a_axis_min), numpy.abs(self._b_axis_max - self._b_axis_min)


class TwoEllipses(BoundaryShape):
    def __init__(self,
                 a1_axis_min=-10e-6, a1_axis_max=10e-6, b1_axis_min=-5e-6, b1_axis_max=5e-6,
                 a2_axis_min=-20e-6, a2_axis_max=20e-6, b2_axis_min=-8e-6, b2_axis_max=8e-6):
        super().__init__()

        self._a1_axis_min = a1_axis_min
        self._a1_axis_max = a1_axis_max
        self._b1_axis_min = b1_axis_min
        self._b1_axis_max = b1_axis_max
        self._a2_axis_min = a2_axis_min
        self._a2_axis_max = a2_axis_max
        self._b2_axis_min = b2_axis_min
        self._b2_axis_max = b2_axis_max
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("a1_axis_min", "x (width) axis 1 starts (signed)  ", "m" ),
                    ("a1_axis_max", "x (width) axis 1 ends (signed)    ", "m" ),
                    ("b1_axis_min", "y (length) axis 1 starts (signed) ", "m" ),
                    ("b1_axis_max", "y (length) axis 1 ends (signed)   ", "m" ),
                    ("a2_axis_min", "x (width) axis 2 starts (signed)  ", "m"),
                    ("a2_axis_max", "x (width) axis 2 ends (signed)    ", "m"),
                    ("b2_axis_min", "y (length) axis 2 starts (signed) ", "m"),
                    ("b2_axis_max", "y (length) axis 2 ends (signed)   ", "m"),
            ] )

    def get_boundaries(self):
        return \
            self._a1_axis_min, self._a1_axis_max, self._b1_axis_min, self._b1_axis_max, \
            self._a2_axis_min, self._a2_axis_max, self._b2_axis_min, self._b2_axis_max

    def get_axis(self):
        return \
            numpy.abs(self._a1_axis_max - self._a1_axis_min), numpy.abs(self._b1_axis_max - self._b2_axis_min), \
            numpy.abs(self._a2_axis_max - self._a2_axis_min), numpy.abs(self._b2_axis_max - self._b2_axis_min)


class Circle(BoundaryShape):
    def __init__(self,radius=50e-6,x_center=0.0,y_center=0.0):
        super().__init__()

        self._radius = radius
        self._x_center = x_center
        self._y_center = y_center
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("radius"              , "radius  ", "m" ),
                    ("x_center"            , "x center (signed)    ", "m" ),
                    ("y_center"            , "y center (signed)    ", "m" ),
            ] )

    def get_boundaries(self):
        return self._radius, self._x_center, self._y_center

    def set_boundaries(self, radius=1.0, x_center=0.0, y_center=0.0):
        self._radius = radius
        self._x_center = x_center
        self._y_center = y_center

    def get_radius(self):
        return self._radius

    def get_center(self):
        return [self._x_center,self._y_center]

class Polygon(BoundaryShape):
    def __init__(self,x=[],y=[]):
        super().__init__()

        self._x = numpy.array(x)
        self._y = numpy.array(y)
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("x"            , "x vertices    ", "m" ),
                    ("y"            , "y vertices    ", "m" ),
            ] )

    def get_boundaries(self):
        return self._x, self._y

    def set_boundaries(self, x, y):
        self._x = numpy.array(x)
        self._y = numpy.array(y)

    def get_number_of_vertices(self):
        n = numpy.array(self._x).size
        if (numpy.abs(self._x[0] - self._x[-1]) < 1e-10)  and (numpy.abs(self._y[0] - self._y[-1]) < 1e-10):
            # print(">>>>> same first and last point")
            n -= 1
        return n

    def get_polygon(self):
        polygon = []
        for i in range(self.get_number_of_vertices()):
            polygon.append([self._x[i], self._y[i]])

        return polygon

    def check_inside_vector(self, x0, y0):
        # see https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
        poly = self.get_polygon()
        n = len(poly)
        x = numpy.array(x0)
        y = numpy.array(y0)

        inside = numpy.zeros(x.size, numpy.bool_)
        p2x = 0.0
        p2y = 0.0
        xints = 0.0
        p1x, p1y = poly[0]

        for i in range(n + 1):
            p2x, p2y = poly[i % n]

            idx = numpy.nonzero((y > min(p1y, p2y)) & (y <= max(p1y, p2y)) & (x <= max(p1x, p2x)))[0]
            if len(idx > 0): # added intuitively by srio TODO: make some tests to compare with self.check_insize
                if p1y != p2y:
                    xints = (y[idx] - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                if p1x == p2x:
                    inside[idx] = ~inside[idx]
                else:
                    idxx = idx[x[idx] <= xints]
                    inside[idxx] = ~inside[idxx]

            p1x, p1y = p2x, p2y
        return inside

    def check_inside(self, x, y):
        return [self.check_inside_one_point(xi, yi) for xi, yi in zip(x, y)]

    def check_inside_one_point(self, x0, y0):
        # see https://stackoverflow.com/questions/36399381/whats-the-fastest-way-of-checking-if-a-point-is-inside-a-polygon-in-python
        poly = self.get_polygon()
        x = x0
        y = y0
        n = len(poly)
        inside = False
        p2x = 0.0
        p2y = 0.0
        xints = 0.0
        p1x, p1y = poly[0]
        for i in range(n + 1):
            p2x, p2y = poly[i % n]
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xints:
                            inside = not inside
            p1x, p1y = p2x, p2y

        return inside

    def check_outside(self, x0, y0):
        inside = self.check_inside(x0, y0)
        return ~inside


class MultiplePatch(BoundaryShape):
    def __init__(self, patch_list=None):
        super().__init__()

        if patch_list is None:
            self._patch_list = []
        else:
            self._patch_list = patch_list
        self._set_support_text([
                    ("multiple_patch_list",  "Multiple Patch", ""),
            ])


    # overwrites the SynedObject method for dealing with list
    def to_dictionary(self):
        dict_to_save = OrderedDict()
        dict_to_save.update({"CLASS_NAME":self.__class__.__name__})

        dict_to_save["multiple_patch_list"] = [el.to_dictionary() for el in self._patch_list]

        return dict_to_save


    def reset(self):
        self._patch_list = []

    def get_number_of_patches(self):
        return len(self._patch_list)

    def get_boundaries(self):
        boundaries_list = []
        for i in range(self.get_number_of_patches()):
            boundaries_list.extend(list(self._patch_list[i].get_boundaries()))
        return tuple(boundaries_list)

    def append_patch(self,patch=BoundaryShape()):
        self._patch_list.append(patch)

    def append_rectangle(self,x_left=-0.010,x_right=0.010,y_bottom=-0.020,y_top=0.020):
        self.append_patch(Rectangle(x_left=x_left, x_right=x_right, y_bottom=y_bottom, y_top=y_top))

    def append_circle(self,radius, x_center=0.0, y_center=0.0):
        self.append_patch(Circle(radius, x_center=x_center, y_center=y_center))

    def append_ellipse(self,a_axis_min, a_axis_max, b_axis_min, b_axis_max):
        self.append_patch(Ellipse(a_axis_min, a_axis_max, b_axis_min, b_axis_max))

    def append_polygon(self,x, y):
        self.append_patch(Polygon(x, y))

    def get_patches(self):
        return self._patch_list

    def get_patch(self,index):
        return self.get_patches()[index]

    def get_name_of_patch(self,index):
        return self._patch_list[index].__class__.__name__

class DoubleRectangle(MultiplePatch):
    def __init__(self, x_left1=-0.010, x_right1=0.0, y_bottom1=-0.020, y_top1=0.0,
                        x_left2=-0.010, x_right2=0.010, y_bottom2=-0.001, y_top2=0.020):
        super().__init__()
        self.reset()
        self.append_patch(Rectangle(x_left=x_left1, x_right=x_right1, y_bottom=y_bottom1, y_top=y_top1))
        self.append_patch(Rectangle(x_left=x_left2, x_right=x_right2, y_bottom=y_bottom2, y_top=y_top2))

    def set_boundaries(self,x_left1=-0.010, x_right1=0.0, y_bottom1=-0.020, y_top1=0.0,
                        x_left2=-0.010, x_right2=0.010, y_bottom2=-0.001, y_top2=0.020):
        self._patch_list[0].set_boundaries(x_left1, x_right1, y_bottom1, y_top1)
        self._patch_list[1].set_boundaries(x_left2, x_right2, y_bottom2, y_top2)

class DoubleEllipse(MultiplePatch):
    def __init__(self, a_axis_min1=-0.010, a_axis_max1=0.0,   b_axis_min1=-0.020, b_axis_max1=0.0,
                       a_axis_min2=-0.010, a_axis_max2=0.010, b_axis_min2=-0.001, b_axis_max2=0.020):
        super().__init__()
        self.reset()
        self.append_patch(Ellipse(a_axis_min1, a_axis_max1, b_axis_min1, b_axis_max1))
        self.append_patch(Ellipse(a_axis_min2, a_axis_max2, b_axis_min2, b_axis_max2))

    def set_boundaries(self,a_axis_min1=-0.010, a_axis_max1=0.0,   b_axis_min1=-0.020, b_axis_max1=0.0,
                            a_axis_min2=-0.010, a_axis_max2=0.010, b_axis_min2=-0.001, b_axis_max2=0.020):
        self._patch_list[0].set_boundaries(a_axis_min1,a_axis_max1,b_axis_min1,b_axis_max1)
        self._patch_list[1].set_boundaries(a_axis_min2,a_axis_max2,b_axis_min2,b_axis_max2)

class DoubleCircle(MultiplePatch):
    def __init__(self, radius1=50e-6,x_center1=0.0,y_center1=0.0,
                       radius2=50e-6,x_center2=100e-6,y_center2=100e-6):
        super().__init__()
        self.reset()
        self.append_patch(Circle(radius1,x_center1,y_center1))
        self.append_patch(Circle(radius2,x_center2,y_center2))

    def set_boundaries(self,radius1=50e-6,x_center1=0.0,y_center1=0.0,
                            radius2=50e-6,x_center2=100e-6,y_center2=100e-6):
        self._patch_list[0].set_boundaries(radius1,x_center1,y_center1)
        self._patch_list[1].set_boundaries(radius2,x_center2,y_center2)



if __name__=="__main__":

    pass


    # ell = Ellipsoid()
    # ell.initialize_from_p_q(20, 10, 0.2618)
    # print ("ellipse axes: ",ell._min_axis/2, ell._maj_axis/2)
    #
    #
    # ell = Ellipsoid(min_axis=ell._min_axis, maj_axis=ell._maj_axis)
    # print("for grazing angle 0.2618, ellipse p,q = ",ell.get_p_q(0.2618))


    # ell = Ellipsoid()
    # p = 20
    # q = 10
    # theta_graz = 0.2618
    # ell.initialize_from_p_q(p, q, theta_graz)
    # print ("ellipse p, q: ",ell.get_p(), ell.get_q())
    # print("ellipse grazing_angle: ", ell.get_grazing_angle())
    # assert (numpy.abs(p - ell.get_p()) < 1e-10 )
    # assert (numpy.abs(q - ell.get_q()) < 1e-10)
    # assert (numpy.abs(theta_graz - ell.get_grazing_angle()) < 1e-10)



    p = 20
    q = 10
    theta_graz = 0.003
    at_infinity = Side.SOURCE
    par = Paraboloid.create_paraboloid_from_p_q(p=p, q=q, grazing_angle=theta_graz, at_infinity=at_infinity, convexity=Convexity.UPWARD)
    print("inputs  p, q, theta_graz: ", p, q, theta_graz, at_infinity)
    print ("ellipse p or q: ",par.get_pole_to_focus())
    print("ellipse par: ", par.get_parabola_parameter())
    print("ellipse grazing_angle: ", par.get_grazing_angle())
    if par.get_at_infinity() == Side.SOURCE:
        assert (numpy.abs(q - par.get_pole_to_focus()) < 1e-10 )
    else:
        assert (numpy.abs(p - par.get_pole_to_focus()) < 1e-10)
    assert (numpy.abs(theta_graz - par.get_grazing_angle()) < 1e-10)


    # circle = Circle(3.0)
    #
    # print(circle.get_radius(),circle.get_center())
    # print(circle.get_boundaries())




    # patches = MultiplePatch()
    #
    # patches.append_rectangle(-0.02,-0.01,-0.001,0.001)
    # patches.append_rectangle(0.01,0.02,-0.001,0.001)
    # patches.append_polygon([-0.02,-0.02,0.02,0.02], [-0.02,0.02,0.02,-0.02])
    #
    # print(patches.get_number_of_patches(),patches.get_boundaries())
    # for patch in patches.get_patches():
    #     print(patch.info())
    # print("Patch 0 is: ",patches.get_name_of_patch(0))
    # print("Patch 1 is: ",patches.get_name_of_patch(1))
    # print(patches.get_boundaries())




    # double_rectangle = DoubleRectangle()
    # double_rectangle.set_boundaries(-0.02,-0.01,-0.001,0.001,0.01,0.02,-0.001,0.001)
    # print("Rectangle 0 is: ",double_rectangle.get_name_of_patch(0))
    # print("Rectangle 1 is: ",double_rectangle.get_name_of_patch(1))
    # print(double_rectangle.get_boundaries())


    # angle = numpy.linspace(0, 2 * numpy.pi, 5)
    # x = numpy.sin(angle) + 0.5
    # y = numpy.cos(angle) + 0.5
    # poly = Polygon(x=x, y=y)
    # print(poly.info())
    # print("vertices: ", poly.get_number_of_vertices())
    # from srxraylib.plot.gol import plot,set_qt
    # set_qt()
    # plot(x,y)
    # print(poly.get_polygon())
    # print(poly.check_inside([0.5,0],[0.5,5]))
    # print(poly.check_outside([0.5, 0], [0.5, 5]))



    # patches = MultiplePatch()
    # patches.append_polygon(numpy.array([-1,-1,1,1]),numpy.array([-1,1,1,-1]))
    # x = [-0.00166557,  0.12180897, -0.11252591, -0.12274196,  0.00586896, -0.12999401, -0.12552975, -0.0377907,  -0.01094828, -0.13689862]
    # y = [ 0.16279557, -0.00085991,  0.01349174, -0.01371226,  0.01480265, -0.04810334, 0.07198068, -0.03725407,  0.13301309, -0.00296213]
    # x = numpy.array(x)
    # y = numpy.array(y)
    # patch = patches.get_patch(0)
    # # # print(patch.check_inside(x,y))
    # for i in range(x.size):
    #     tmp = patch.check_inside_one_point(x[i], y[i])
    #     print(x[i], y[i], tmp )
    # print(patch.check_inside(x, y))
    # print(patch.check_inside_vector(x, y))
