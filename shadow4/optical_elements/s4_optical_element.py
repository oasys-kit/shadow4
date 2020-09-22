from shadow4.syned.shape import Direction, Convexity, \
    Sphere, SphericalCylinder, \
    Ellipsoid, EllipticalCylinder, \
    Hyperboloid, HyperbolicCylinder, \
    Paraboloid, ParabolicCylinder, \
    Toroid, Conic, Plane

class S4OpticalElement(object):

    def __init__(self):
        pass

    def info(self):
        raise NotImplementedError()

    def to_python_code(self, data=None):
        raise NotImplementedError()


class SurfaceCalculation:
    INTERNAL = 0
    EXTERNAL = 1

class S4CurvedOpticalElement(S4OpticalElement):

    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 ):
        S4OpticalElement.__init__(self)

        self._surface_calculation = surface_calculation
        self._is_cylinder = is_cylinder
        self._curved_surface_shape = None # auxiliary field

class S4PlaneOpticalElement(S4OpticalElement):
    def __init__(self):
        S4OpticalElement.__init__(self)
        
        self._plane_surface_shape = Plane()

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
            self._curved_surface_shape = Toroid.create_toroid_from_radii(min_radius, maj_radius, maj_radius)
        else:
            self._curved_surface_shape = Toroid.create_toroid_from_p_q(p_focus, q_focus, grazing_angle)

class S4ConicOpticalElement(S4CurvedOpticalElement):
    def __init__(self, conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]):
        S4CurvedOpticalElement.__init__(self, surface_calculation=SurfaceCalculation.INTERNAL, is_cylinder=False)

        self._conic_surface_shape = Conic(conic_coefficients=conic_coefficients)

class S4ParaboloidOpticalElement(S4CurvedOpticalElement):
    pass
