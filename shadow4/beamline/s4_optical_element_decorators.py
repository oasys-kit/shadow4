"""
Defines the functionality for shadow4 optical elements.

Basically it adds some methods:

- get_info: to be overloaded with specific info for a particular optical element.

- to_python_code: to be overloaded to create the python code for the specifoc optical element.

- get_surface_shape_instance: returns the SYNED shape instance.

- get_optical_surface_instance: returns the S4 optical surface instance.

"""


import numpy
import os

from syned.beamline.shape import Direction, Convexity
from syned.beamline.shape import Sphere, SphericalCylinder
from syned.beamline.shape import Ellipsoid, EllipticalCylinder
from syned.beamline.shape import Hyperboloid, HyperbolicCylinder
from syned.beamline.shape import Paraboloid, ParabolicCylinder
from syned.beamline.shape import Toroid, Conic, NumericalMesh, Plane, Side

from shadow4.optical_surfaces.s4_conic import S4Conic
from shadow4.optical_surfaces.s4_mesh import S4Mesh
from shadow4.optical_surfaces.s4_toroid import S4Toroid

class S4OpticalElementDecorator(object):
    """
    Base abstract class.
    """

    def __init__(self):
        pass

    def get_info(self):
        """
        Returns string containing some specific information of the shadow4 optical element.

        Returns
        -------
        str
        """
        return "\n\nSpecific information not available (get_info() method not overloaded).\n\n"

    def to_python_code(self, **kwargs):
        """
        Returns string containing python code to create shadow4 optical element.

        Raises
        ------
        NotImplementedError.
        """
        raise NotImplementedError()

    def get_surface_shape_instance(self): # return a SYNED object of type Shape
        """
        Returns a SYNED object of type Shape.

        Raises
        ------
        NotImplementedError.
        """
        raise NotImplementedError()

    def get_optical_surface_instance(self): # return a Shadow4 object of type S4OpticalSurface
        raise NotImplementedError()

class SurfaceCalculation:
    """
    class supporting the flag for external or internal calculation:

    INTERNAL = 0

    EXTERNAL = 1

    """
    INTERNAL = 0
    EXTERNAL = 1

class S4PlaneOpticalElementDecorator(S4OpticalElementDecorator):
    """
    Creates the shadow4 optical element decorator for a plane surface.
    """
    def __init__(self):

        S4OpticalElementDecorator.__init__(self)

        self._plane_surface_shape = Plane()

    def get_surface_shape_instance(self):
        """
        Returns the SYNED instance of type Shape.

        Returns
        -------
        instance of syned.beamline.shape.Plane
        """
        return self._plane_surface_shape

    def get_optical_surface_instance(self):
        """
        Returns a shadow4 optical element object of type optical surface.

        Returns
        -------
        instance of shadow4.s4_optical_surfaces.s4_conic.S4Conic
        """
        return S4Conic.initialize_as_plane()

class S4CurvedOpticalElementDecorator(S4OpticalElementDecorator):
    """
    Creates the shadow4 optical element decorator for a curved surface.

    Parameters
    ----------
    surface_calculation : int, optional
        flag:
            0 = SurfaceCalculation.INTERNAL,
            1 = SurfaceCalculation.EXTERNAL.
    is_cylinder : int, optional
        flag:
            0=No (there is revolution symmetry along Y)
            1=Yes (flat surface along X or Y).
    curved_surface_shape : syned shape instance, optional

    """
    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 curved_surface_shape=None
                 ):

        S4OpticalElementDecorator.__init__(self)

        self._surface_calculation  = surface_calculation
        self._is_cylinder          = is_cylinder
        self._curved_surface_shape = curved_surface_shape
    
    def get_surface_shape_instance(self):
        """
        Returns the SYNED instance of type Shape.

        Returns
        -------
        instance of shape class.
        """
        return self._curved_surface_shape

##################################################
# CURVED OPTICAL ELEMENTS
##################################################

class S4SphereOpticalElementDecorator(S4CurvedOpticalElementDecorator):
    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL, # external=1, internal=0
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 radius=0.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):
        """
        Creates the shadow4 optical element decorator for a Sphere surface.

        Parameters
        ----------
        surface_calculation : int, optional
            flag:
                0 = SurfaceCalculation.INTERNAL,
                1 = SurfaceCalculation.EXTERNAL.
        is_cylinder : int, optional
            flag:
                0=No (there is revolution symmetry along Y)
                1=Yes (flat surface along X or Y).
        cylinder_direction : int (as defined by Direction), optional
            TANGENTIAL = 0, SAGITTAL = 1.
        convexity : int (as defined by Convexity), optional
            NONE = -1, UPWARD = 0, DOWNWARD = 1.
        radius : float, optional
            The radius of the sphere.
        p_focus : float, optional
            For surface_calculation=1, the p distance (from source plane to center of the optical element).
        q_focus : float, optional
            For surface_calculation=1, the q distance (from the center of the optical element to the image plane).
        grazing_angle : float, optional
            For surface_calculation=1, the grazing angle in rad.
        """
        if surface_calculation == SurfaceCalculation.EXTERNAL:
            if is_cylinder: curved_surface_shape = SphericalCylinder.create_spherical_cylinder_from_radius(radius, convexity, cylinder_direction)
            else:           curved_surface_shape = Sphere.create_sphere_from_radius(radius, convexity)
        else:
            if is_cylinder: curved_surface_shape = SphericalCylinder.create_spherical_cylinder_from_p_q(p_focus, q_focus, grazing_angle, convexity, cylinder_direction)
            else:           curved_surface_shape = Sphere.create_sphere_from_p_q(p_focus, q_focus, grazing_angle, convexity)

        S4CurvedOpticalElementDecorator.__init__(self,
                                                 surface_calculation=surface_calculation,
                                                 is_cylinder=is_cylinder,
                                                 curved_surface_shape=curved_surface_shape)

    def get_optical_surface_instance(self):
        """
        Returns a shadow4 optical element object of type optical surface.

        Returns
        -------
        instance of shadow4.s4_optical_surfaces.s4_conic.S4Conic
        """
        surface_shape = self.get_surface_shape_instance()

        switch_convexity = 0 if surface_shape.get_convexity() == Convexity.DOWNWARD else 1

        radius = surface_shape.get_radius()
        if isinstance(surface_shape, SphericalCylinder):
            cylindrical = 1
            cylangle = 0.0 if surface_shape.get_cylinder_direction() == Direction.TANGENTIAL else (0.5 * numpy.pi)
        elif isinstance(surface_shape, Sphere):
            cylindrical = 0
            cylangle    = 0.0

        print(">>>>> S4SphereOpticalElement.get_optical_surface_instance(): R, cyl, cyl_angle, optical element, ", radius, cylindrical, cylangle, surface_shape)
        out = S4Conic.initialize_as_sphere_from_curvature_radius(radius, cylindrical=cylindrical, cylangle=cylangle, switch_convexity=switch_convexity)
        print(">>>>> Sphere ccc", out.ccc)
        return out

class S4EllipsoidOpticalElementDecorator(S4CurvedOpticalElementDecorator):
    """
    Creates the shadow4 optical element decorator for an Ellipsoid surface.

    Parameters
    ----------
    surface_calculation : int, optional
        flag:
            0 = SurfaceCalculation.INTERNAL,
            1 = SurfaceCalculation.EXTERNAL.
    is_cylinder : int, optional
        flag:
            0=No (there is revolution symmetry along Y)
            1=Yes (flat surface along X or Y).
    cylinder_direction : int (as defined by Direction), optional
        TANGENTIAL = 0, SAGITTAL = 1.
    convexity : int (as defined by Convexity), optional
        NONE = -1, UPWARD = 0, DOWNWARD = 1.
    min_axis : float, optional
        For surface_calculation=0, The minor axis of the ellipsoid (2a).
    maj_axis : float, optional
        For surface_calculation=0, The major axis of the ellipsoid (2b)
    p_focus : float, optional
        For surface_calculation=1,2, the p distance (from source plane to center of the optical element).
    q_focus : float, optional
        For surface_calculation=1, the q distance (from the center of the optical element to the image plane).
    grazing_angle : float, optional
        For surface_calculation=1, the grazing angle in rad.
    """
    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 min_axis=0.0,
                 maj_axis=0.0,
                 pole_to_focus=0.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):
        if surface_calculation == SurfaceCalculation.EXTERNAL:
            if is_cylinder: curved_surface_shape = EllipticalCylinder.create_elliptical_cylinder_from_axes(min_axis, maj_axis, pole_to_focus, convexity, cylinder_direction)
            else:           curved_surface_shape = Ellipsoid.create_ellipsoid_from_axes(min_axis, maj_axis, pole_to_focus, convexity)
        else:
            if is_cylinder: curved_surface_shape = EllipticalCylinder.create_elliptical_cylinder_from_p_q(p_focus, q_focus, grazing_angle, convexity, cylinder_direction)
            else:           curved_surface_shape = Ellipsoid.create_ellipsoid_from_p_q(p_focus, q_focus, grazing_angle, convexity)

        S4CurvedOpticalElementDecorator.__init__(self, surface_calculation, is_cylinder, curved_surface_shape)

    def get_optical_surface_instance(self): # todo: update this one like hyperboloid
        """
        Returns a shadow4 optical element object of type optical surface.

        Returns
        -------
        instance of shadow4.s4_optical_surfaces.s4_conic.S4Conic
        """
        surface_shape = self.get_surface_shape_instance()

        switch_convexity = 0 if surface_shape.get_convexity() == Convexity.DOWNWARD else 1

        if isinstance(surface_shape, EllipticalCylinder):
            print(">>>>> EllipticalCylinder optical element", surface_shape)
            cylindrical = 1
            cylangle = 0.0 if surface_shape.get_cylinder_direction() == Direction.TANGENTIAL else (0.5 * numpy.pi)
        elif isinstance(surface_shape, Ellipsoid):
            print(">>>>> Ellipsoid optical element", surface_shape)
            cylindrical = 0
            cylangle    = 0.0

        out = S4Conic.initialize_as_ellipsoid_from_focal_distances(surface_shape.get_p_focus(),
                                                                    surface_shape.get_q_focus(),
                                                                    surface_shape.get_grazing_angle(),
                                                                    cylindrical=cylindrical,
                                                                    cylangle=cylangle,
                                                                    switch_convexity=switch_convexity)
        print(">>>>> Ellipsoid ccc", out.ccc)
        return out

class S4HyperboloidOpticalElementDecorator(S4CurvedOpticalElementDecorator):
    """
    Creates the shadow4 optical element decorator for an Hyperboloid surface.

    Parameters
    ----------
    surface_calculation : int, optional
        flag:
            0 = SurfaceCalculation.INTERNAL,
            1 = SurfaceCalculation.EXTERNAL.
    is_cylinder : int, optional
        flag:
            0=No (there is revolution symmetry along Y)
            1=Yes (flat surface along X or Y).
    cylinder_direction : int (as defined by Direction), optional
        NONE = -1, UPWARD = 0, DOWNWARD = 1.
    convexity : int (as defined by Convexity), optional
        NONE = -1, UPWARD = 0, DOWNWARD = 1.
    min_axis : float, optional
        For surface_calculation=0, The minor axis of the ellipsoid (2a).
    maj_axis : float, optional
        For surface_calculation=0, The major axis of the ellipsoid (2b)
    p_focus : float, optional
        For surface_calculation=1,2, the p distance (from source plane to center of the optical element).
    q_focus : float, optional
        For surface_calculation=1, the q distance (from the center of the optical element to the image plane).
    grazing_angle : float, optional
        For surface_calculation=1, the grazing angle in rad.
    """
    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 min_axis=0.0,
                 maj_axis=0.0,
                 pole_to_focus=0.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):
        if surface_calculation == SurfaceCalculation.EXTERNAL:
            if is_cylinder: curved_surface_shape = HyperbolicCylinder.create_hyperbolic_cylinder_from_axes(min_axis, maj_axis, pole_to_focus, convexity, cylinder_direction)
            else:           curved_surface_shape = Hyperboloid.create_hyperboloid_from_axes(min_axis, maj_axis, pole_to_focus, convexity)
        else:
            if is_cylinder: curved_surface_shape = HyperbolicCylinder.create_hyperbolic_cylinder_from_p_q(p_focus, q_focus, grazing_angle, convexity, cylinder_direction)
            else:           curved_surface_shape = Hyperboloid.create_hyperboloid_from_p_q(p_focus, q_focus, grazing_angle, convexity)

        S4CurvedOpticalElementDecorator.__init__(self, surface_calculation, is_cylinder, curved_surface_shape)

    def get_optical_surface_instance(self):
        """
        Returns a shadow4 optical element object of type optical surface.

        Returns
        -------
        instance of shadow4.s4_optical_surfaces.s4_conic.S4Conic
        """
        surface_shape = self.get_surface_shape_instance()

        switch_convexity = 0 if surface_shape.get_convexity() == Convexity.DOWNWARD else 1

        if isinstance(surface_shape, HyperbolicCylinder):
            print(">>>>> HyperbolicCylinder optical element", surface_shape)
            cylindrical = 1
            cylangle = 0.0 if surface_shape.get_cylinder_direction() == Direction.TANGENTIAL else (0.5 * numpy.pi)
        elif isinstance(surface_shape, Hyperboloid):
            print(">>>>> Hyperboloid optical element", surface_shape)
            cylindrical = 0
            cylangle    = 0.0

        out = S4Conic.initialize_as_hyperboloid_from_focal_distances(surface_shape.get_p_focus(),
                                                                      surface_shape.get_q_focus(),
                                                                      surface_shape.get_grazing_angle(),
                                                                      cylindrical=cylindrical,
                                                                      cylangle=cylangle,
                                                                      switch_convexity=switch_convexity)
        print(">>>>> Hyperboloid ccc", out.ccc)
        return out

class S4ToroidOpticalElementDecorator(S4CurvedOpticalElementDecorator):
    """
    Creates the shadow4 optical element decorator for an Toroid surface.

    Parameters
    ----------
    surface_calculation : int, optional
        flag:
            0 = SurfaceCalculation.INTERNAL,
            1 = SurfaceCalculation.EXTERNAL.
    min_radius : float, optional
        The minor axis of the toroid.
    maj_radius : float, optional
        The major axis of the toroid.
    p_focus : float, optional
        For surface_calculation=1, the p distance (from source plane to center of the optical element).
    q_focus : float, optional
        For surface_calculation=1, the q distance (from the center of the optical element to the image plane).
    grazing_angle : float, optional
        For surface_calculation=1, the grazing angle in rad.
    """
    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 min_radius=0.0,
                 maj_radius=0.0,
                 f_torus=0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):
        if surface_calculation == SurfaceCalculation.EXTERNAL:
            curved_surface_shape = Toroid.create_toroid_from_radii(min_radius, maj_radius)
        else:
            curved_surface_shape = Toroid.create_toroid_from_p_q(p_focus, q_focus, grazing_angle)

        S4CurvedOpticalElementDecorator.__init__(self, surface_calculation, False, curved_surface_shape)

        self._f_torus = f_torus

    def get_optical_surface_instance(self):
        """
        Returns a shadow4 optical element object of type optical surface.

        Returns
        -------
        instance of shadow4.s4_optical_surfaces.s4_conic.S4Toroid
        """
        surface_shape = self.get_surface_shape_instance()

        print(">>>>> Toroidal optical element (syned stored) %g   %g" % (surface_shape.get_maj_radius(), surface_shape.get_min_radius()) )
        print(">>>>>       setting S4Toroid r_maj=%f, r_min=%g, f_torus=%d" % (
            surface_shape.get_maj_radius() - surface_shape.get_min_radius(),
            surface_shape.get_min_radius(),
            self._f_torus))


        return S4Toroid(
            r_maj=surface_shape.get_maj_radius() - surface_shape.get_min_radius(),
            r_min=surface_shape.get_min_radius(),
            f_torus=self._f_torus)

class S4ParaboloidOpticalElementDecorator(S4CurvedOpticalElementDecorator):
    """
    Creates the shadow4 optical element decorator for an Paraboloid surface.

    Parameters
    ----------
    surface_calculation : int, optional
        flag:
            0 = SurfaceCalculation.INTERNAL,
            1 = SurfaceCalculation.EXTERNAL.
    is_cylinder : int, optional
        flag:
            0=No (there is revolution symmetry along Y)
            1=Yes (flat surface along X or Y).
    cylinder_direction : int (as defined by Direction), optional
        NONE = -1, UPWARD = 0, DOWNWARD = 1.
    convexity : int (as defined by Convexity), optional
        NONE = -1, UPWARD = 0, DOWNWARD = 1.
    parabola_parameter : float, optional
        For surface_calculation=0, The parabola parameter PARAM (y^2 = 2 PARAM z)
    at_infinity : int, optional
        For surface_calculation=0, flag to indicate that
            0: the source is at infinity (focusing paraboloid) = Side.SOURCE.
            1: the image is at infinity (collimating paraboloid) = Side.IMAGE.
    pole_to_focus : float, optional
        For surface_calculation=0, the p or q distance (from focus to center of the optical element).
    p_focus : float, optional
        For surface_calculation=1, the q distance (from source plane to the center of the optical element).
    q_focus : float, optional
        For surface_calculation=1, the q distance (from the center of the optical element to the image plane).
    grazing_angle : float, optional
        For surface_calculation=1, the grazing angle in rad.
    """
    def __init__(self,
                 surface_calculation=SurfaceCalculation.INTERNAL,
                 is_cylinder=False,
                 cylinder_direction=Direction.TANGENTIAL,
                 convexity=Convexity.UPWARD,
                 parabola_parameter=0.0,
                 at_infinity=Side.SOURCE,
                 pole_to_focus=0.0,
                 p_focus=0.0,
                 q_focus=0.0,
                 grazing_angle=0.0,
                 ):
        if surface_calculation == SurfaceCalculation.EXTERNAL:
            if is_cylinder:
                curved_surface_shape = ParabolicCylinder.create_parabolic_cylinder_from_parabola_parameter(
                    parabola_parameter, at_infinity, pole_to_focus, convexity, cylinder_direction)
            else:
                curved_surface_shape = Paraboloid.create_paraboloid_from_parabola_parameter(
                    parabola_parameter, at_infinity, pole_to_focus, convexity)
        else:
            if is_cylinder:
                curved_surface_shape = ParabolicCylinder.create_parabolic_cylinder_from_p_q(
                    p_focus, q_focus, grazing_angle, at_infinity, convexity, cylinder_direction)
            else:
                curved_surface_shape = Paraboloid.create_paraboloid_from_p_q(
                    p_focus, q_focus, grazing_angle, at_infinity, convexity)

        S4CurvedOpticalElementDecorator.__init__(self, surface_calculation, is_cylinder, curved_surface_shape)

    def get_optical_surface_instance(self):
        """
        Returns a shadow4 optical element object of type optical surface.

        Returns
        -------
        instance of shadow4.s4_optical_surfaces.s4_conic.S4Conic
        """
        surface_shape = self.get_surface_shape_instance()

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

        out = S4Conic.initialize_as_paraboloid_from_focal_distances(p, q, surface_shape.get_grazing_angle(), cylindrical=cylindrical, cylangle=cylangle, switch_convexity=switch_convexity)

        print(">>>>> Paraboloid ccc", out.ccc)
        return out

class S4ConicOpticalElementDecorator(S4CurvedOpticalElementDecorator):

    def __init__(self, conic_coefficients=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0]):
        """
        Creates the shadow4 optical element decorator for a Conic surface.

        Parameters
        ----------
        conic_coefficients : list, ndarray, optional
            The list of the 10 conic coefficients.
        """
        S4CurvedOpticalElementDecorator.__init__(
            self,
            surface_calculation=SurfaceCalculation.INTERNAL,
            is_cylinder=False,
            curved_surface_shape=None if conic_coefficients is None else Conic(conic_coefficients=conic_coefficients))

    def get_optical_surface_instance(self):
        """
        Returns a shadow4 optical element object of type optical surface.

        Returns
        -------
        instance of shadow4.s4_optical_surfaces.s4_conic.S4Conic
        """
        surface_shape = self.get_surface_shape_instance()
        out = S4Conic.initialize_from_coefficients(surface_shape.get_conic_coefficients())
        print(">>>>> Conic optical element")
        print(">>>>> Conic ccc", out.ccc)
        return out


class S4NumericalMeshOpticalElementDecorator(S4CurvedOpticalElementDecorator):
    """
    Creates the shadow4 optical element decorator for a Mesh surface.

    Parameters
    ----------
    xx : ndarray, optional
        the 1D array with the X points.
    yy : ndarray, optional
        the 1D array with the Y points.
    zz : ndarray, optional
        the 2D [shape Nx,Ny] array with the Z points.
    surface_data_file : str, optional
        the name of the h5 file with the mesh.
    """
    def __init__(self, xx=None, yy=None, zz=None, surface_data_file=None):
        S4CurvedOpticalElementDecorator.__init__(self,
                                                 surface_calculation=SurfaceCalculation.INTERNAL,
                                                 is_cylinder=False,
                                                 curved_surface_shape = NumericalMesh(xx, yy, zz, surface_data_file))

    def get_optical_surface_instance(self):
        """
        Returns a shadow4 optical element object of type optical surface.

        Returns
        -------
        instance of shadow4.s4_optical_surfaces.s4_conic.S4Mesh
        """
        surface_shape = self.get_surface_shape_instance()

        print(">>>>> SurfaceData optical element")
        numerical_mesh = S4Mesh()

        if surface_shape.has_surface_data(): numerical_mesh.load_surface_data(surface_shape)
        elif surface_shape.has_surface_data_file():
            filename, file_extension = os.path.splitext(surface_shape._surface_data_file)

            if file_extension.lower() in [".h5", ".hdf", ".hdf5"]: numerical_mesh.load_h5file(surface_shape._surface_data_file)
            else:                                                  numerical_mesh.load_file(surface_shape._surface_data_file) # 3 columns ASCII

        return numerical_mesh

##################################################
# LENS OPTICAL ELEMENTS
##################################################

class S4RefractiveLensOpticalElementDecorator(S4CurvedOpticalElementDecorator):
    """
    Creates the shadow4 optical element decorator for a Refractive lens.

    Parameters
    ----------
    surface_shape : int, optional
        Flag for the single lens surface shape: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
    convex_to_the_beam : int, optional
        Flag for convexity of the first interface exposed to the beam
        For surface_shape (1,2): 0=No, 1=Yes, or in other words 0=Concave, 1=Convex
    cylinder_angle : int, optional
        Flag for cylindrical symmetry.
        For surface_shape (1,2): 0=not cylindricaL, 1=meridional 2=sagittal
    ri_calculation_mode : int, optional
        Flag for the used source of refraction indices and absorption coefficients:
           0=User
           1=prerefl file
           2=direct calculation using xraylib
           3=direct calculation using dabax
    prerefl_file : str, optional
        For ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
    refraction_index : float, optional
        For ri_calculation_mode=1: the real part of the refraction coefficient Real[n].
    attenuation_coefficient : float, optional
        For ri_calculation_mode=1: the linear absorption coefficient, in cm^-1 (real)..=
    radius : float, optional
        For surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid).
    conic_coefficients : list or ndarray, optional
        For surface_shape = 3: the 10 conic coefficients.
    """
    def __init__(self,
                 surface_shape=1,      # now: 0=plane, 1=sphere, 2=parabola, 3=conic coefficients
                                       # (in shadow3: 1=sphere 4=paraboloid, 5=plane)
                 convex_to_the_beam=1, # for surface_shape (1,2): convexity of the first interface exposed to the beam 0=No, 1=Yes
                                       # the second interface has opposite convexity
                 cylinder_angle=0,     # for surface_shape (1,2): 0=not cylindricaL, 1=meridional 2=sagittal
                 ri_calculation_mode=0,       # source of refraction indices and absorption coefficients
                                 # 0=User
                                 # 1=prerefl file
                                 # 2=direct calculation using xraylib
                                 # 3=direct calculation using dabax
                 prerefl_file=None,    # for ri_calculation_mode=0: file name (from prerefl) to get the refraction index.
                 refraction_index=1.0, # for ri_calculation_mode=1: n (real)
                 attenuation_coefficient=0.0, # for ri_calculation_mode=1: mu in cm^-1 (real)
                 radius=500e-6,        # for surface_shape=(1,2): lens radius [m] (for spherical, or radius at the tip for paraboloid)
                 conic_coefficients=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0],   # for surface_shape = 3: the conic coefficients
                 ):


        conic_coefficients = self._get_conic_coefficients(surface_shape,
                                                          radius,
                                                          cylinder_angle,
                                                          convex_to_the_beam,
                                                          conic_coefficients)

        curved_surface_shape = [Conic(conic_coefficients=conic_coefficients[0]),
                                Conic(conic_coefficients=conic_coefficients[1])]

        S4CurvedOpticalElementDecorator.__init__(self,
                                                 surface_calculation=SurfaceCalculation.EXTERNAL,
                                                 curved_surface_shape=curved_surface_shape)

        self._ri_calculation_mode     = ri_calculation_mode
        self._prerefl_file            = prerefl_file
        self._refraction_index        = refraction_index
        self._attenuation_coefficient = attenuation_coefficient

        print(">>>>> conic_coefficients: ", conic_coefficients)

    def _get_conic_coefficients(self, surface_shape, radius, cylinder_angle, convex_to_the_beam, conic_coefficients):
        if surface_shape == 0:   conic_coefficients_1 = [0, 0, 0, 0, 0, 0, 0, 0, -1, 0]
        elif surface_shape == 1: conic_coefficients_1 = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0 * radius, 0.0]
        elif surface_shape == 2: conic_coefficients_1 = [1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0 * radius, 0.0]
        elif surface_shape == 3: conic_coefficients_1 = conic_coefficients.copy()
        else: return None

        if   cylinder_angle == 1: conic_coefficients_1[0] = conic_coefficients_1[3] = conic_coefficients_1[5] = conic_coefficients_1[6] = 0
        elif cylinder_angle == 2: conic_coefficients_1[1] = conic_coefficients_1[3] = conic_coefficients_1[4] = conic_coefficients_1[7] = 0

        if convex_to_the_beam == 1: conic_coefficients_1[8] *= -1

        conic_coefficients_2 = conic_coefficients_1.copy()
        conic_coefficients_2[8] *= -1

        return [conic_coefficients_1, conic_coefficients_2]

    def get_optical_surface_instance(self):
        """
        Returns a shadow4 optical element object of type optical surface.

        Returns
        -------
        list
            instances of shadow4.s4_optical_surfaces.s4_conic.S4Conic
        """
        surface_shapes = self.get_surface_shape_instance()

        return [S4Conic.initialize_from_coefficients(numpy.array(surface_shapes[0].get_conic_coefficients())),
                S4Conic.initialize_from_coefficients(numpy.array(surface_shapes[1].get_conic_coefficients()))]
