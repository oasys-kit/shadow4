"""
Samplers for typical mathematical probability distributions (1D and 2D).
"""
import numpy
from syned.syned_object import SynedObject

class DistributionGeneric(SynedObject):
    """
    Base class for a mathematical distribution.

    Note:
        It inherits from SynedObject to use the syned mechanism of displaying parameters.

    """
    def __init__(self):
        SynedObject.__init__(self)

    def get_dimension(self):
        """
        Returns the dimension. To be defined in the derived class.

        Raises
        ------
            NotImplementedError
        """
        raise NotImplementedError("To be implemented in a derived class.")

    def get_sampled_points(self):
        """
        Returns the number of sampling points. To be defined in the derived class.

        Raises
        ------
            NotImplementedError
        """
        raise NotImplementedError("To be implemented in the subclasses")

#
# main distribution subclasses:
#      Distribution1D
#      Distribution2D
#
class Distribution2D(DistributionGeneric):
    """
    Defines a generic 2D mathematical distribution class.
    """
    def __init__(self):
        DistributionGeneric.__init__(self)

    def get_dimension(self):
        """
        Returns the dimension of the distribution.

        Returns
        -------
        int
            returns 2.
        """
        return 2


class Distribution1D(DistributionGeneric):
    """
    Defines a generic 1D mathematical distribution class.
    """
    def __init__(self):
        pass
        
    def get_dimension(self):
        """
        Returns the dimension of the distribution.

        Returns
        -------
        int
            returns 1.
        """
        return 1


#
# Subclasses for Distribution2D used for spatial type sampling
#  ["Point","Rectangle","Ellipse","Gaussian"]
#
class Point2D(Distribution2D):
    def __init__(self, h_center=0.0, v_center=0.0):
        """
        Defines a point distribution in 2D.

        Parameters
        ----------
        h_center : float, optional
            The position of the point (abscissa).
        v_center : float, optional
            The position of the point (ordinate).
        """
        self._h_center = h_center
        self._v_center = v_center

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("h_center"         , "h (center) ", "" ),
                    ("v_center"         , "v (center) ", "" ),
            ] )

    def get_sampled_points(self,N):
        """
        Returns the sampled points.

        Notes:
        -----
        The multiple resulting points are always identical, as it is a point.

        Parameters
        ----------
        N : int
            The number of points to be sampled.

        Returns
        -------
        tuple
            (H,V) The arrays for the H and V.
        """
        return numpy.zeros(N) + self._h_center, numpy.zeros(N) + self._v_center


class Rectangle2D(Distribution2D):
    def __init__(self, h_min, h_max, v_min, v_max):
        """
        Defines a rectangular 2D mathematical distribution.

        Parameters
        ----------
        h_min : float
            The minimum coordinate of the rectangle in the horizontal direction.
        h_max : float
            The maximum coordinate of the rectangle in the horizontal direction.
        v_min : float
            The minimum coordinate of the rectangle in the vertical direction.
        v_max : float
            The maximum coordinate of the rectangle in the vertical direction.
        """
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

    def get_sampled_points(self, N):
        """
        Returns the sampled points.

        Parameters
        ----------
        N : int
            The number of points to be sampled.

        Returns
        -------
        tuple
            (H,V) The arrays for the H and V.
        """
        return Uniform1D.sample(N, self._h_min, self._h_max), Uniform1D.sample(N, self._v_min, self._v_max)

    @classmethod
    def sample(cls, N, h_min, h_max, v_min, v_max):
        """
        Returns sampled points for a 2D rectangular distribution.

        Parameters
        ----------
        N : int
            The number of points to be sampled.
        h_min : float
            The minimum coordinate of the rectangle in the horizontal direction
        h_max : float
            The maximum coordinate of the rectangle in the horizontal direction
        v_min : float
            The minimum coordinate of the rectangle in the vertical direction
        v_max : float
            The maximum coordinate of the rectangle in the vertical direction

        Returns
        -------
        tuple
            (H,V) The arrays for the H and V.
        """
        return Rectangle2D(h_min, h_max, v_min, v_max).get_sampled_points(N)


class Ellipse2D(Distribution2D):
    """
    Defines an ellipse 2D mathematical distribution.

    Parameters
    ----------
    h_min : float
        The minimum coordinate of the ellipse in the horizontal direction.
    h_max : float
        The maximum coordinate of the ellipse in the horizontal direction.
    v_min : float
        The minimum coordinate of the ellipse in the vertical direction.
    v_max : float
        The maximum coordinate of the ellipse in the vertical direction.
    """
    def __init__(self, h_min, h_max, v_min, v_max):
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
        # return ["Flat","Uniform","Gaussian","Cone"]

    def get_sampled_points(self, N):
        """
        Returns the sampled points.

        Parameters
        ----------
        N : int
            The number of points to be sampled.

        Returns
        -------
        tuple
            (H,V) The arrays for the H and V.
        """
        # ! C Elliptical source **
        # ! C Uses a transformation algorithm to generate a uniform variate distribution
        phi = numpy.pi * 2 * numpy.random.random(N)
        radius = numpy.sqrt(numpy.random.random(N))
        x = 0.5 *(self._h_max+self._h_min) + 0.5 * (self._h_max-self._h_min) * radius * numpy.cos(phi)
        y = 0.5 *(self._v_max+self._v_min) + 0.5 * (self._v_max-self._v_min) * radius * numpy.sin(phi)
        return x,y

    @classmethod
    def sample(cls, N, h_min, h_max, v_min, v_max):
        """
        Returns sampled points for a 2D ellipse distribution.

        Parameters
        ----------
        N : int
            The number of points to be sampled.
        h_min : float
            The minimum coordinate of the ellipse in the horizontal direction
        h_max : float
            The maximum coordinate of the ellipse in the horizontal direction
        v_min : float
            The minimum coordinate of the ellipse in the vertical direction
        v_max : float
            The maximum coordinate of the ellipse in the vertical direction

        Returns
        -------
        tuple
            (H,V) The arrays for the H and V.
        """
        return Ellipse2D(h_min,h_max,v_min,v_max).get_sampled_points(N)

class Gaussian2D(Distribution2D):
    """
    Defines a Gaussian 2D mathematical distribution.

    Parameters
    ----------
    sigma_h : float
        The sigma in the horizontal direction.
    sigma_v : float
        The sigma in the vertical direction.
    """
    def __init__(self, sigma_h, sigma_v):
        self._sigma_h = sigma_h
        self._sigma_v = sigma_v
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("sigma_h"         , "h (width) sigma", "" ),
                    ("sigma_v"         , "v (width) sigma", "" ),
            ] )


    def get_sampled_points(self, N):
        """
        Returns the sampled points.

        Parameters
        ----------
        N : int
            The number of points to be sampled.

        Returns
        -------
        tuple
            (H,V) The arrays for the H and V.
        """
        return Gaussian1D.sample(N, self._sigma_h), Gaussian1D.sample(N, self._sigma_v)

    @classmethod
    def sample(cls, N, sigma_h, sigma_v):
        """
        Returns sampled points for a 2D Gaussian distribution.

        Parameters
        ----------
        N : int
            The number of points to be sampled.
        sigma_h : float
            The Gaussian sigma in the horizontal direction.
        sigma_v : float
            The Gaussian sigma in the vertical direction.

        Returns
        -------
        tuple
            (H,V) The arrays for the H and V.
        """
        return Gaussian2D(sigma_h,sigma_v).get_sampled_points(N)


#
# Subclasses for Distribution2D used for angle emission sampling
#  ["Flat","Uniform","Gaussian","Cone"]
#


class Flat2D(Rectangle2D):
    """
    Defines a flat 2D mathematical distribution (the same as Rectangle2D).

    "Flat" means that the ray divergence distribution is constant versus the angles with the YZ and XY planes.
    Strictly speaking, the "angles" with the YZ and XY planes are indeed the direction cosines with the X and Z axis.
    In the small angle approximation, theta=sin(theta).

    Parameters
    ----------
    h_min : float
        The minimum coordinate of the rectangle in the horizontal direction.
    h_max : float
        The maximum coordinate of the rectangle in the horizontal direction.
    v_min : float
        The minimum coordinate of the rectangle in the vertical direction.
    v_max : float
        The maximum coordinate of the rectangle in the vertical direction.
    """
    def __init__(self, h_min, h_max, v_min, v_max):
        super().__init__(h_min, h_max, v_min, v_max)


class Uniform2D(Distribution2D):
    """
    Defines a uniform 2D mathematical distribution.

    "Uniform" means that the rays will illuminate homogeneously a screen at a given distance of the source.
    This corresponds to the isotropic emitter.

    Parameters
    ----------
    h_min : float
        The minimum angular coordinate in the horizontal direction.
    h_max : float
        The maximum angular coordinate in the horizontal direction.
    v_min : float
        The minimum angular coordinate in the vertical direction.
    v_max : float
        The maximum angular coordinate in the vertical direction.
    """
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

    def get_sampled_points(self, N):
        """
        Returns the sampled points.

        Parameters
        ----------
        N : int
            The number of points to be sampled.

        Returns
        -------
        tuple
            (H,V) The arrays for the H and V.
        """
        # ! C   Uniform distribution ( Isotrope emitter )
        XMAX1 =   numpy.tan(self._h_min)
        XMAX2 =   numpy.tan(self._h_max)
        ZMAX1 =   numpy.tan(self._v_min)
        ZMAX2 =   numpy.tan(self._v_max)
        XRAND = numpy.random.random(N) * (XMAX1 - XMAX2) + XMAX2
        ZRAND = numpy.random.random(N) * (ZMAX1 - ZMAX2) + ZMAX2
        THETAR  = numpy.arctan(numpy.sqrt(XRAND**2 + ZRAND**2))
        PHIR = numpy.arctan2(ZRAND, XRAND)
        DIREC1  = numpy.cos(PHIR) * numpy.sin(THETAR)
        DIREC3  = numpy.sin(PHIR) * numpy.sin(THETAR)
        return DIREC1, DIREC3

    @classmethod
    def sample(cls, N, h_min, h_max, v_min, v_max):
        """
        Returns sampled points for a 2D Uniform distribution.

        Parameters
        ----------
        N : int
            The number of points to be sampled.
        h_min : float
            The minimum angular coordinate in the horizontal direction.
        h_max : float
            The maximum angular coordinate in the horizontal direction.
        v_min : float
            The minimum angular coordinate in the vertical direction.
        v_max : float
            The maximum angular coordinate in the vertical direction.

        Returns
        -------
        tuple
            (H,V) The arrays for the H and V.
        """
        return Uniform2D(h_min,h_max,v_min,v_max).get_sampled_points(N)


class Cone2D(Distribution2D):
    """
    Defines the 2D Cone mathematical distribution.

    Parameters
    ----------
    cone_max : float, optional
        The maximum aperture of the cone in rad.
    cone_min : float, optional
        The minimum aperture of the cone in rad.
    """
    def __init__(self, cone_max=10e-6, cone_min=0.0):
        self._cone_max = cone_max
        self._cone_min = cone_min
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("cone_max"         , "max angle for cone semiaperture  ", "" ),
                    ("cone_min"         , "max angle for cone semiaperture  ", "" ),
            ] )

    def get_sampled_points(self, N):
        """
        Returns the sampled points.

        Parameters
        ----------
        N : int
            The number of points to be sampled.

        Returns
        -------
        tuple
            (H,V) The arrays for the H and V.
        """
        # ! C   Now generates a set of rays along a cone centered about the normal, plus a ray along the normal itself.
        ANGLE = 2 * numpy.pi * numpy.random.random(N)
        ANG_CONE = numpy.cos(self._cone_min) - numpy.random.random(N) * \
                                               (numpy.cos(self._cone_min)-numpy.cos(self._cone_max))
        ANG_CONE = numpy.arccos(ANG_CONE)
        DIREC1 = numpy.sin(ANG_CONE) * numpy.cos(ANGLE)
        DIREC3 = numpy.sin(ANG_CONE) * numpy.sin(ANGLE)
        return DIREC1,DIREC3

    @classmethod
    def sample(cls,N,cone_max=10e-6,cone_min=0.0):
        def sample(cls, N, h_min, h_max, v_min, v_max):
            """
            Returns sampled points for a 2D Cone distribution.

            Parameters
            ----------
            N : int
                The number of points to be sampled.
            cone_max : float, optional
                The maximum aperture of the cone in rad.
            cone_min : float, optional
                The minimum aperture of the cone in rad.

            Returns
            -------
            tuple
                (H,V) The arrays for the H and V.
            """
        return Cone2D(cone_max=cone_max,cone_min=cone_min).get_sampled_points(N)


#
# subclasses for Distribution1D
#
class Uniform1D(Distribution1D):
    def __init__(self, x_min=-0.010, x_max=0.010):
        """
        Defines a 1D uniform (flat) distribution.

        Parameters
        ----------
        x_min : float, optional
            The minimum coordinate.
        x_max : float, optional
            The maximum coordinate.
        """
        self._x_min  = x_min
        self._x_max  = x_max
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("x_min"         , "minimum (signed)", "" ),
                    ("x_max"         , "maximum (signed)", "" ),
            ] )

    def get_sampled_points(self, N):
        """
        Returns the sampled points.

        Parameters
        ----------
        N : int
            The number of points to be sampled.

        Returns
        -------
        numpy array
            The arrays with the N sampled points.
        """
        return numpy.random.random(N) * (self._x_max-self._x_min) + self._x_min

    @classmethod
    def sample(cls, N=1000, x_min=-0.010, x_max=0.010):
        """
        Returns sampled points for a 1D Uniform distribution.

        Parameters
        ----------
        N : int
            The number of points to be sampled.
        x_min : float, optional
            The minimum coordinate.
        x_max : float, optional
            The maximum coordinate.

        Returns
        -------
        numpy array
            The arrays with the N sampled points.
        """
        return Uniform1D(x_min=x_min, x_max=x_max).get_sampled_points(N)

class Gaussian1D(Distribution1D):
    """
    Defines a 1D Gaussian distribution.

    Parameters
    ----------
    sigma : float, optional
        The sigma of the Gaussian.
    center : float, optional
        The center of the Gaussian.
    """
    def __init__(self, sigma=1e-3, center=0.0):
        super().__init__()
        self._sigma  = sigma
        self._center  = center
        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
                    ("sigma"         , "sigma", "" ),
                    ("center"        , "center", "" ),
            ] )

    def get_sampled_points(self, N):
        """
        Returns the sampled points.

        Parameters
        ----------
        N : int
            The number of points to be sampled.

        Returns
        -------
        numpy array
            The arrays with the N sampled points.
        """
        return numpy.random.normal(loc=self._center, scale=self._sigma, size=N)

    @classmethod
    def sample(cls, N=1000, sigma=0.25, center=0.0):
        """
        Returns sampled points for a 1D Uniform distribution.

        Parameters
        ----------
        N : int
            The number of points to be sampled.
        sigma : float, optional
            The sigma of the Gaussian.
        center : float, optional
            The center of the Gaussian.

        Returns
        -------
        numpy array
            The arrays with the N sampled points.
        """
        return Gaussian1D(sigma=sigma, center=center).get_sampled_points(N)

if __name__=="__main__":

    from srxraylib.plot.gol import plot,plot_scatter

    do_plot = True

    #
    # constructors
    #
    u = Uniform1D(-10, 10)
    print(u.info())
    u2 = Uniform2D(-10, 10, -5, 5)
    print(u2.info())

    #
    # Gaussian1D
    #
    N = 1000
    g = Gaussian1D(sigma=0.25)
    print(g.info())
    sampled_gaussian = g.get_sampled_points(N)
    assert ( (sampled_gaussian.std() - 0.25) < 0.05)
    if do_plot: plot(numpy.arange(N), Gaussian1D.sample(sigma=0.25, N=N), title="Gaussian")

    #
    # Uniform1D
    #
    N = 1000
    g = Uniform1D(10, 20)
    print(g.info())
    sampled = g.get_sampled_points(N)
    if do_plot: plot(numpy.arange(N),sampled_gaussian,title="Uniform")
    print(">>>>>",numpy.abs(sampled.mean() - 15 ))
    assert ( numpy.abs(sampled.mean() - 15 ) < 0.3)
    if do_plot: plot(numpy.arange(N), Uniform1D.sample(N, 10, 20))

    #
    # Ellipse
    #
    N = 1000
    g = Ellipse2D(-4, 0, -3, -1)
    print(g.info())
    x,y = g.get_sampled_points(N)
    if do_plot: plot_scatter(1e6*x, 1e6*y, title="Ellipse")

    #
    # Rectangle2D
    #
    dx, dy = Rectangle2D.sample(5000, -2.5e-6, 2.5e-6, -0.5e-6, 0.5e-6,)
    if do_plot: plot_scatter(1e6*dx, 1e6*dy, title="Rectangle")

    #
    # Uniform2D
    #

    dx, dy = Uniform2D.sample(5000, -2.5e-6, 2.5e-6, -0.5e-6, 0.5e-6,)
    if do_plot: plot_scatter(1e6*dx, 1e6*dy, title="Uniform")

    #
    # Flat2D
    #

    dx, dy = Flat2D.sample(5000, -2.5e-6, 2.5e-6, -0.5e-6, 0.5e-6,)
    if do_plot: plot_scatter(1e6*dx, 1e6*dy, title="Flat2D")

    #
    # Gaussian2D
    #

    dx, dy = Gaussian2D.sample(5000, 2.5e-6, 0.5e-6,)
    if do_plot: plot_scatter(1e6*dx, 1e6*dy, title="Gaussian2D")

    #
    # Cone2D
    #
    dx, dy = Cone2D.sample(5000, 2.5e-6, 0.5e-6,)
    if do_plot: plot_scatter(1e6*dx, 1e6*dy, title="Cone2D")