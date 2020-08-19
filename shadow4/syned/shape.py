
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

class Conic(SurfaceShape):
    def __init__(self, 
                 conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
                 convexity=Convexity.UPWARD):
        SurfaceShape.__init__(self, convexity)

        self._conic_coefficients = conic_coefficients

class Plane(SurfaceShape):
    def __init__(self):
        SurfaceShape.__init__(self, convexity=Convexity.NONE)

class Sphere(SurfaceShape):
    def __init__(self, radius=1.0, convexity=Convexity.UPWARD):
        SurfaceShape.__init__(self, convexity=Convexity.NONE)
        self._radius = radius

    def get_radius(self):
        return self._radius

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        self.self._radius = Sphere.get_radius_from_p_q(p, q, grazing_angle)

    @classmethod
    def get_radius_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003):
        # 1/p + 1/q = 2/(R cos(pi/2 - gr.a.))
        return (2*p*q/(p+q))/numpy.sin(grazing_angle)

class SphericalCylinder(Sphere, Cylinder):
    def __init__(self, 
                 radius=1.0, 
                 convexity=Convexity.UPWARD, 
                 cylinder_direction=Direction.TANGENTIAL):
        Sphere.__init__(self, radius, convexity)
        Cylinder.__init__(self, cylinder_direction)

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        if self._cylinder_direction == Direction.TANGENTIAL:
            self._conic_coefficients[9] = -Sphere.get_radius_from_p_q(p, q, grazing_angle)**2
        elif self._cylinder_direction == Direction.SAGITTAL:
            self._conic_coefficients[9] = -SphericalCylinder.get_radius_from_p_q_sagittal(p, q, grazing_angle)**2

    @classmethod
    def get_radius_from_p_q_sagittal(cls, p=2.0, q=1.0, grazing_angle=0.003):
        # 1/p + 1/q = 2 cos(pi/2 - gr.a.)/r
        return (2*p*q/(p+q))*numpy.sin(grazing_angle)

class Ellipsoid(SurfaceShape):
    def __init__(self, min_axis=0.0, maj_axis=0.0, convexity=Convexity.UPWARD):
        SurfaceShape.__init__(self, convexity)

        self._min_axis = min_axis
        self._maj_axis = maj_axis

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        self._min_axis, self._maj_axis = Ellipsoid.get_axis_from_p_q(p, q, grazing_angle)

    def get_p_q(self, grazing_angle=0.003):
        return Ellipsoid.get_p_q_from_axis(self._min_axis, self._maj_axis, grazing_angle)

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
                 convexity=Convexity.UPWARD, 
                 cylinder_direction=Direction.TANGENTIAL):
        Ellipsoid.__init__(self, min_axis, maj_axis, convexity)
        Cylinder.__init__(self, cylinder_direction)

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        if self._cylinder_direction == Direction.SAGITTAL:
            raise NotImplementedError("Operation not possible for SAGITTAL direction")

        return super().initialize_from_p_q(p, q, grazing_angle)

    def get_p_q(self, grazing_angle=0.003):
        if self._cylinder_direction == Direction.SAGITTAL:
            raise NotImplementedError("Operation not possible for SAGITTAL direction")

        return super().get_p_q(grazing_angle)

class Paraboloid(SurfaceShape):
    def __init__(self, 
                 parabola_parameter=0.0, 
                 convexity=Convexity.UPWARD):
        SurfaceShape.__init__(self, convexity)

        self._parabola_parameter = parabola_parameter

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003, at_infinity=Side.SOURCE):
        self._parabola_parameter = Paraboloid.get_parabola_parameter_from_p_q(p, q, grazing_angle, at_infinity)

    @classmethod
    def get_parabola_parameter_from_p_q(cls, p=2.0, q=1.0, grazing_angle=0.003, at_infinity=Side.SOURCE):
        if at_infinity == Side.IMAGE:
            return 2*p*(numpy.sin(grazing_angle))**2
        elif at_infinity == Side.SOURCE:
            return 2*q*(numpy.sin(grazing_angle))**2


class ParabolicCylinder(Paraboloid, Cylinder):
    def __init__(self, 
                 parabola_parameter=0.0, 
                 convexity=Convexity.UPWARD, 
                 cylinder_direction=Direction.TANGENTIAL):
        Paraboloid.__init__(self, parabola_parameter, convexity)
        Cylinder.__init__(self, cylinder_direction)

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003, at_infinity=Side.SOURCE):
        if self._cylinder_direction == Direction.SAGITTAL:
            raise NotImplementedError("Operation not possible for SAGITTAL direction")

        return super().initialize_from_p_q(p, q, grazing_angle, at_infinity)

class Hyperboloid(SurfaceShape):
    def __init__(self, min_axis=0.0, maj_axis=0.0, convexity=Convexity.UPWARD):
        SurfaceShape.__init__(self, convexity)

        self._min_axis = min_axis
        self._maj_axis = maj_axis

    def initialize_from_p_q(self, p=2.0, q=1.0, grazing_angle=0.003):
        raise NotImplementedError("TBD")

    def get_p_q(self, grazing_angle=0.003):
        raise NotImplementedError("TBD")

class HyperbolicCylinder(Hyperboloid, Cylinder):
    def __init__(self, 
                 min_axis=0.0, 
                 maj_axis=0.0, 
                 convexity=Convexity.UPWARD, 
                 cylinder_direction=Direction.TANGENTIAL):
        Hyperboloid.__init__(self, min_axis, maj_axis, convexity)
        Cylinder.__init__(self, cylinder_direction)

class Toroidal(SurfaceShape):
    def __init__(self, min_radius=0.0, maj_radius=0.0):
        SurfaceShape.__init__(self, convexity=Convexity.NONE)
        
        self._min_radius = min_radius
        self._maj_radius = maj_radius

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("min_radius"         , "Minor radius r   ", "m" ),
                    ("maj_radius"         , "Major radius R (optical=R+r)", "m" ),
            ] )
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

class NumericalMesh(SurfaceShape):
    def __init__(self, h5file=""):
        SurfaceShape.__init__(self, convexity=Convexity.NONE)
        self._h5file = h5file

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

    # ell = Ellipsoid()
    # ell.initialize_from_p_q(20, 10, 0.2618)
    #
    # print (ell._min_axis/2, ell._maj_axis/2)
    #
    # ell = Ellipsoid(min_axis=ell._min_axis, maj_axis=ell._maj_axis)
    #
    # print(ell.get_p_q(0.2618))


    # circle = Circle(3.0)
    #
    # print(circle.get_radius(),circle.get_center())
    # print(circle.get_boundaries())




    # patches = MultiplePatch()
    #
    # patches.append_rectangle(-0.02,-0.01,-0.001,0.001)
    # patches.append_rectangle(0.01,0.02,-0.001,0.001)
    #
    # print(patches.get_number_of_patches(),patches.get_boundaries())
    #
    # for patch in patches.get_patches():
    #     print(patch.info())
    #
    # print("Patch 0 is: ",patches.get_name_of_patch(0))
    # print("Patch 1 is: ",patches.get_name_of_patch(1))
    # print(patches.get_boundaries())


    double_rectangle = DoubleRectangle()
    double_rectangle.set_boundaries(-0.02,-0.01,-0.001,0.001,0.01,0.02,-0.001,0.001)
    print("Rectangle 0 is: ",double_rectangle.get_name_of_patch(0))
    print("Rectangle 1 is: ",double_rectangle.get_name_of_patch(1))
    print(double_rectangle.get_boundaries())
