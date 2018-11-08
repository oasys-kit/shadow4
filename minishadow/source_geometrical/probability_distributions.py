
import numpy
from syned.syned_object import SynedObject
from collections import OrderedDict




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
        DistributionGeneric.__init__(self)
        
    def get_dimension(self):
        return 1

#
# Subclasses for Distribution2D
#

class RectangleUniform2D(Distribution2D):
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


class NumericalMesh2D(Distribution2D):
    def __init__(self,function=None):
        self.function = function


#
# subclasses for Distribution1D
#


class Uniform1D(Distribution1D):
    def __init__(self, x_min=-0.010, x_max=0.010):
        super().__init__()

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

    u2 = RectangleUniform2D(-10,10,-5,5)
    print(u2.info())

    #
    # Gaussian
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
    # Uniform
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