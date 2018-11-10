
import numpy
from syned.syned_object import SynedObject


class DistributionGeneric(SynedObject):
    def __init__(self):
        SynedObject.__init__(self)

    def get_dimension(self):
        raise Exception("To be implemented in the subclasses")

    def get_sampled_points(self):
        raise Exception("To be implemented in the subclasses")

#
# main distribution subclasses:
#      Distribution1D
#      Distribution2D
class Distribution2D(DistributionGeneric):
    def __init__(self):
        DistributionGeneric.__init__(self)

    def get_dimension(self):
        return 2


class Distribution1D(DistributionGeneric):
    def __init__(self):
        pass
        
    def get_dimension(self):
        return 1


#
# Subclasses for Distribution2D used for spatial type sampling
#  ["Point","Rectangle","Ellipse","Gaussian"]

class Point2D(Distribution2D):
    def __init__(self, h_center=0.0,v_center=0.0):
        self._h_center = h_center
        self._v_center = v_center

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("h_center"         , "h (center) ", "" ),
                    ("v_center"         , "v (center) ", "" ),
            ] )

    def get_sampled_points(self,N):
        return numpy.zeros(N)+self._h_center,numpy.zeros(N)+self._v_center


class Rectangle2D(Distribution2D):
    def __init__(self, h_min,h_max,v_min,v_max):
        self._h_min = h_min
        self._h_max = h_max
        self._v_min = v_min
        self._v_max = v_max
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("h_min"         , "h (width) minimum (signed)   ", "" ),
                    ("h_max"         , "h (width) maximum (signed)   ", "" ),
                    ("v_min"         , "v (length) minimum (signed)  ", "" ),
                    ("v_max"         , "v (length) maximum (signed)  ", "" ),
            ] )

        # return ["Point","Rectangle","Ellipse","Gaussian"]
        # return ["Flat","Uniform","Gaussian","Conical"]

    def get_sampled_points(self,N):
        return Uniform1D.sample(N,self._h_min,self._h_max),Uniform1D.sample(N,self._v_min,self._v_max)

    @classmethod
    def sample(cls,N,h_min,h_max,v_min,v_max):
        return Rectangle2D(h_min,h_max,v_min,v_max).get_sampled_points(N)


class Ellipse2D(Distribution2D):
    def __init__(self, h_min,h_max,v_min,v_max):
        self._h_min = h_min
        self._h_max = h_max
        self._v_min = v_min
        self._v_max = v_max
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("h_min"         , "h (width) minimum (signed)   ", "" ),
                    ("h_max"         , "h (width) maximum (signed)   ", "" ),
                    ("v_min"         , "v (length) minimum (signed)  ", "" ),
                    ("v_max"         , "v (length) maximum (signed)  ", "" ),
            ] )

        # return ["Point","Rectangle","Ellipse","Gaussian"]
        # return ["Flat","Uniform","Gaussian","Conical"]

    def get_sampled_points(self,N):
        # ! C
        # ! C Elliptical source **
        # ! C Uses a transformation algorithm to generate a uniform variate distribution
        # ! C
        phi = numpy.pi * 2 * numpy.random.random(N)
        radius = numpy.sqrt(numpy.random.random(N))
        x = 0.5 *(self._h_max+self._h_min) + 0.5 * (self._h_max-self._h_min) * radius * numpy.cos(phi)
        y = 0.5 *(self._v_max+self._v_min) + 0.5 * (self._v_max-self._v_min) * radius * numpy.sin(phi)
        return x,y

    @classmethod
    def sample(cls,N,h_min,h_max,v_min,v_max):
        return Ellipse2D(h_min,h_max,v_min,v_max).get_sampled_points(N)

class Gaussian2D(Distribution2D):
    def __init__(self, sigma_h,sigma_v):
        self._sigma_h = sigma_h
        self._sigma_v = sigma_v

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("sigma_h"         , "h (width) sigma", "" ),
                    ("sigma_v"         , "v (width) sigma", "" ),
            ] )

        # return ["Point","Rectangle","Ellipse","Gaussian"]
        # return ["Flat","Uniform","Gaussian","Conical"]

    def get_sampled_points(self,N):
        return Gaussian1D.sample(N,self._sigma_h),Gaussian1D.sample(N,self._sigma_v)

    @classmethod
    def sample(cls,N,sigma_h,sigma_v):
        return Gaussian2D(sigma_h,sigma_v).get_sampled_points(N)


#
# Subclasses for Distribution2D used for angle emission sampling
#  ["Flat","Uniform","Gaussian","Conical"]
#


class Flat2D(Rectangle2D):
    def __init__(self,h_min,h_max,v_min,v_max):
        super().__init__(h_min,h_max,v_min,v_max)


class Uniform2D(Distribution2D):
    def __init__(self, h_min,h_max,v_min,v_max):
        self._h_min = h_min
        self._h_max = h_max
        self._v_min = v_min
        self._v_max = v_max
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("h_min"         , "h (width) minimum (signed)   ", "" ),
                    ("h_max"         , "h (width) maximum (signed)   ", "" ),
                    ("v_min"         , "v (length) minimum (signed)  ", "" ),
                    ("v_max"         , "v (length) maximum (signed)  ", "" ),
            ] )

        # return ["Point","Rectangle","Ellipse","Gaussian"]
        # return ["Flat","Uniform","Gaussian","Conical"]

    def get_sampled_points(self,N):
        # ! C
        # ! C   Uniform distribution ( Isotrope emitter )
        # ! C
        XMAX1 =   numpy.tan(self._h_min)
        XMAX2 =   numpy.tan(self._h_max)
        ZMAX1 =   numpy.tan(self._v_min)
        ZMAX2 =   numpy.tan(self._v_max)
        XRAND = numpy.random.random(N) * (XMAX1 - XMAX2) + XMAX2
        ZRAND = numpy.random.random(N) * (ZMAX1 - ZMAX2) + ZMAX2
        THETAR  = numpy.arctan(numpy.sqrt(XRAND**2+ZRAND**2))
        PHIR = numpy.arctan2(ZRAND,XRAND)
        DIREC1  = numpy.cos(PHIR) * numpy.sin(THETAR)
        # DIREC2  = numpy.cos(THETAR)
        DIREC3  = numpy.sin(PHIR) * numpy.sin(THETAR)
        return DIREC1,DIREC3

    @classmethod
    def sample(cls,N,h_min,h_max,v_min,v_max):
        return Uniform2D(h_min,h_max,v_min,v_max).get_sampled_points(N)




class Conical2D(Distribution2D):
    def __init__(self, cone_max=10e-6,cone_min=0.0):
        self._cone_max = cone_max
        self._cone_min = cone_min
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("cone_max"         , "max angle for cone semiaperture  ", "" ),
                    ("cone_min"         , "max angle for cone semiaperture  ", "" ),
            ] )

    def get_sampled_points(self,N):
        # ! C   Now generates a set of rays along a cone centered about the normal,
        # ! C   plus a ray along the normal itself.
        # ! C
        # IF (FGRID.EQ.1.OR.FGRID.EQ.3) THEN
        #   ANGLE =   TWOPI*GRID(4,ITIK)*(IDO_VX-1)/IDO_VX
        # ELSE
        # ANGLE =   TWOPI*GRID(4,ITIK)
        # END IF
        # ! C temp fix -- 16 Jan 1987
        # ! C      ANG_CONE =   CONE_MIN +
        # ! C     $ (CONE_MAX - CONE_MIN)*GRID(6,ITIK)
        # ANG_CONE =   COS(CONE_MIN) - GRID(6,ITIK)*(COS(CONE_MIN)-COS(CONE_MAX))
        # ANG_CONE =  ACOS(ANG_CONE)
        # DIREC(1) =   SIN(ANG_CONE)*COS(ANGLE)
        # DIREC(2) =   COS(ANG_CONE)
        # DIREC(3) =   SIN(ANG_CONE)*SIN(ANGLE)

        ANGLE = 2 * numpy.pi * numpy.random.random(N)
        ANG_CONE = numpy.cos(self._cone_min) - numpy.random.random(N) * \
                                               (numpy.cos(self._cone_min)-numpy.cos(self._cone_max))
        ANG_CONE = numpy.arccos(ANG_CONE)
        DIREC1 = numpy.sin(ANG_CONE) * numpy.cos(ANGLE)
        # DIREC2 = numpy.cos(ANG_CONE)
        DIREC3 = numpy.sin(ANG_CONE) * numpy.sin(ANGLE)

        return DIREC1,DIREC3

    @classmethod
    def sample(cls,N,cone_max=10e-6,cone_min=0.0):
        return Conical2D(cone_max=cone_max,cone_min=cone_min).get_sampled_points(N)







# class NumericalMesh2D(Distribution2D):
#     def __init__(self,function=None):
#         self.function = function


#
# subclasses for Distribution1D
#


class Uniform1D(Distribution1D):
    def __init__(self, x_min=-0.010, x_max=0.010):

        self._x_min  = x_min
        self._x_max  = x_max


        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("x_min"         , "minimum (signed)", "" ),
                    ("x_max"         , "maximum (signed)", "" ),
            ] )

    def get_sampled_points(self,N):
        return numpy.random.random(N) * (self._x_max-self._x_min) + self._x_min

    @classmethod
    def sample(cls,N=1000,x_min=-0.010, x_max=0.010):
        return Uniform1D(x_min=x_min, x_max=x_max).get_sampled_points(N)

class Gaussian1D(Distribution1D):
    def __init__(self, sigma=1e-3, center=0.0):
        super().__init__()

        self._sigma  = sigma
        self._center  = center


        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("sigma"         , "sigma", "" ),
                    ("center"        , "center", "" ),
            ] )

    def get_sampled_points(self,N):
        return numpy.random.normal(loc=self._center, scale=self._sigma, size=N)

    @classmethod
    def sample(cls,N=1000,sigma=0.25,center=0.0):
        return Gaussian1D(sigma=sigma,center=center).get_sampled_points(N)

if __name__=="__main__":

    from srxraylib.plot.gol import plot,plot_scatter

    do_plot = False

    #
    # constructors
    #
    u = Uniform1D(-10,10)
    print(u.info())

    u2 = Uniform2D(-10,10,-5,5)
    print(u2.info())

    #
    # Gaussian1D
    #
    N = 1000
    g = Gaussian1D(sigma=0.25)
    sampled_gaussian = g.get_sampled_points(N)
    if do_plot:
        plot(numpy.arange(N),sampled_gaussian)
    assert ( (sampled_gaussian.std() - 0.25) < 0.05)

    if do_plot:
        plot(numpy.arange(N),Gaussian1D.sample(sigma=0.25,N=N),title="Gaussian")

    #
    # Uniform1D
    #
    N = 1000
    g = Uniform1D(10,20)
    sampled = g.get_sampled_points(N)
    if do_plot:
        plot(numpy.arange(N),sampled_gaussian,title="Uniform")

    print(">>>>>",numpy.abs(sampled.mean() - 15 ))
    assert ( numpy.abs(sampled.mean() - 15 ) < 0.2)


    if do_plot:
        plot(numpy.arange(N),Uniform1D.sample(N,10,20))

    #
    # Elipse
    #
    N = 1000
    g = Ellipse2D(-4,0,-3,-1)
    x,y = g.get_sampled_points(N)
    # plot_scatter(x,y,xrange=[-6,6],yrange=[-6,6])

    #
    # Rectangle2D
    #

    dx,dy = Rectangle2D.sample(5000,-2.5e-6,2.5e-6,-0.5e-6,0.5e-6,)
    print(dx)
    if do_plot:
        plot_scatter(1e6*dx,1e6*dy,title="Rectangle")

    #
    # Uniform2D
    #

    dx,dy = Uniform2D.sample(5000,-2.5e-6,2.5e-6,-0.5e-6,0.5e-6,)
    print(dx)
    if do_plot:
        plot_scatter(1e6*dx,1e6*dy,title="Uniform")

    #
    # Flat2D
    #

    dx,dy = Flat2D.sample(5000,-2.5e-6,2.5e-6,-0.5e-6,0.5e-6,)

    if do_plot:
        plot_scatter(1e6*dx,1e6*dy,title="Flat2D")

    #
    # Gaussian2D
    #

    dx,dy = Gaussian2D.sample(5000,2.5e-6,0.5e-6,)
    if do_plot:
        plot_scatter(1e6*dx,1e6*dy,title="Gaussian2D")

    #
    # Conical2D
    #

    dx,dy = Conical2D.sample(5000,2.5e-6,0.5e-6,)
    if do_plot:
        plot_scatter(1e6*dx,1e6*dy,title="Gaussian2D")