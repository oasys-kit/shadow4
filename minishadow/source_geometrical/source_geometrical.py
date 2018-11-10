import numpy
from probability_distributions import Rectangle2D, Ellipse2D, Gaussian2D
from probability_distributions import Flat2D, Uniform2D, Conical2D


class SourceGeometrical(object):
    def __init__(self,
                    spatial_type="Point",
                    angular_distribution = "Flat",
                    energy_distribution = "Single line",
                 ):

        self.set_spatial_type_by_name(spatial_type)
        self.set_angular_distribution_by_name(angular_distribution)
        # self.energy_distribution = energy_distribution



    @classmethod
    def spatial_type_list(cls):
        # fsour	= 0 - spatial source type/shape in X-Z plane.  Options are:
        # 			  point (0), rectangle (1), ellipse (2), gaussian (3).
        return ["Point","Rectangle","Ellipse","Gaussian"]

    def set_spatial_type_by_name(self,name):
        if name == "Point":
            self.set_spatial_type_point()
        elif name == "Rectangle":
            self.set_spatial_type_rectangle()
        elif name == "Ellipse":
            self.set_spatial_type_ellipse()
        elif name == "Gaussian":
            self.set_spatial_type_gaussian()
        else:
            raise Exception("Wrong spatial type: %s"%name)

    def set_spatial_type_point(self):
        self.spatial_type = "Point"

    def set_spatial_type_rectangle(self,width=2.0,height=1.0):
        # wxsou	=  0.0000000000000000E+00 - for fsour=1,2; source width (X).
        # wzsou	=  0.0000000000000000E+00 - for fsour=1,2; source height (Z).
        self.spatial_type = "Rectangle"
        self.__wxsou = width
        self.__wzsou = height

    def set_spatial_type_ellipse(self,width=2.0,height=1.0):
        # wxsou	=  0.0000000000000000E+00 - for fsour=1,2; source width (X).
        # wzsou	=  0.0000000000000000E+00 - for fsour=1,2; source height (Z).
        self.spatial_type = "Ellipse"
        self.__wxsou = width
        self.__wzsou = height

    def set_spatial_type_gaussian(self,sigma_h=2.0,sigma_v=1.0):
        # sigmax	=  0.0000000000000000E+00 - for fsour=3; sigma in X
        # sigmaz	=  0.0000000000000000E+00 - for fsour=3; sigma in Z
        self.spatial_type = "Gaussian"
        self.__sigmax = sigma_h
        self.__sigmaz = sigma_v

    @classmethod
    def angular_distribution_list(cls):
         # fdistr	= 2 - defines source angle distribution types:
         # 		  Available options are: flat(1),uniform(2),
         # 		  gaussian(3), synchrotron(4), conical(5), exact
         # 		  synchrotron(6).
        return ["Flat","Uniform","Gaussian","Conical"]

        # cone_max	=  0.0000000000000000E+00 - for fdistr=5; maximum half
        # 					    divergence.
        # cone_min	=  0.0000000000000000E+00 - for fdistr=5; minimum half
        # 					    divergence.
        #
        # sigdix	=  0.0000000000000000E+00 - for fdistr=3; sigma (radians)for horizontal
        # 				    divergence (gaussian angle distribution).
        # sigdiz	=  0.0000000000000000E+00 - for fdistr=3; sigma (radians) for vertical
        # 				    divergence (gaussian angle distribution).
        # hdiv1	=  0.0000000000000000E+00 - horizontal divergence in +X (radians).
        # hdiv2	=  0.0000000000000000E+00 - horizontal divergence in -X (radians).
        # vdiv1	=  6.0000000000000002E-05 - vertical divergence in +Z (radians).
        # vdiv2	=  6.0000000000000002E-05 - vertical divergence in -Z (radians).

    def set_angular_distribution_by_name(self,name):
        if name == "Flat":
            self.set_angular_distribution_flat()
        elif name == "Uniform":
            self.set_angular_distribution_uniform()
        elif name == "Gaussian":
            self.set_angular_distribution_gaussian()
        elif name == "Conical":
            self.set_angular_distribution_conical()
        else:
            raise Exception("Wrong angular distribution: %s"%name)



    # WARNING: in shadow4 limits are signed!!!!
    def set_angular_distribution_flat(self,hdiv1=-5e-6,hdiv2=5e-6,vdiv1=-0.5e-6,vdiv2=0.5e-6):
        self.angular_distribution = "Flat"
        self.__hdiv1 = hdiv1
        self.__hdiv2 = hdiv2
        self.__vdiv1 = vdiv1
        self.__vdiv2 = vdiv2

    # WARNING: in shadow4 limits are signed!!!!
    def set_angular_distribution_uniform(self,hdiv1=-5e-6,hdiv2=5e-6,vdiv1=-0.5e-6,vdiv2=0.5e-6):
        self.angular_distribution = "Uniform"
        self.__hdiv1 = hdiv1
        self.__hdiv2 = hdiv2
        self.__vdiv1 = vdiv1
        self.__vdiv2 = vdiv2

    def set_angular_distribution_gaussian(self,sigdix=1e-6,sigdiz=1e-6):
        self.angular_distribution = "Gaussian"
        self.__sigdix = sigdix
        self.__sigdiz = sigdiz

    def set_angular_distribution_conical(self,cone_max=10e-6,cone_min=0.0):
        self.angular_distribution = "Conical"
        self.__cone_max = cone_max
        self.__cone_min = cone_min

    @classmethod
    def energy_distribution_list(cls):
        return ["Single line","Several lines","Uniform","Relative intensities","Gaussian","User defined"]

    @classmethod
    def _sample_rays_default(cls,N=5000):
        rays = numpy.zeros((N,18))
        rays[:,4] = 1.0  # vy
        rays[:,6] = 1.0  # Ex
        rays[:,9] = 1.0  # flag
        rays[:,10] = 2 * numpy.pi / (1e-10 * 100) # wavenumber
        rays[:,11] = numpy.arange(N) + 1          # index
        return rays

    def calculate_rays(self,N=5000):

        rays = self._sample_rays_default(N)

        print(">> Spatial type: %s"%(self.spatial_type))

        #
        # spatial type
        #
        if self.spatial_type == "Point":
            pass
        elif self.spatial_type == "Rectangle":
            rays[:,0],rays[:,2] = Rectangle2D.sample(N,
                                    -0.5*self.__wxsou,
                                    +0.5*self.__wxsou,
                                    -0.5*self.__wzsou,
                                    +0.5*self.__wzsou)
        elif self.spatial_type == "Ellipse":
            rays[:,0],rays[:,2] = Ellipse2D.sample(N,
                                    -0.5*self.__wxsou,
                                    +0.5*self.__wxsou,
                                    -0.5*self.__wzsou,
                                    +0.5*self.__wzsou)
        elif self.spatial_type == "Gaussian":
            rays[:,0],rays[:,2] = Gaussian2D.sample(N,
                                    self.__sigmax,
                                    self.__sigmaz)
        else:
            raise Exception("Bad value of spatial_type")

        #
        # angular distribution
        #
        print(">> Angular distribution: %s"%(self.angular_distribution))

        if self.angular_distribution == "Flat":
            rays[:,3],rays[:,5] = Flat2D.sample(N,
                                    self.__hdiv1,
                                    self.__hdiv2,
                                    self.__vdiv1,
                                    self.__vdiv2)
            rays[:,4] = numpy.sqrt(-rays[:,3]**2 - rays[:,5]**2 + 1.0)
        elif self.angular_distribution == "Uniform":
            rays[:,3],rays[:,5] = Uniform2D.sample(N,
                                    self.__hdiv1,
                                    self.__hdiv2,
                                    self.__vdiv1,
                                    self.__vdiv2)
            rays[:,4] = numpy.sqrt(-rays[:,3]**2 - rays[:,5]**2 + 1.0)
        elif self.angular_distribution == "Gaussian":
            rays[:,3],rays[:,5] = Gaussian2D.sample(N,
                                    self.__sigdix,
                                    self.__sigdiz)
            rays[:,4] = numpy.sqrt(-rays[:,3]**2 - rays[:,5]**2 + 1.0)
        elif self.angular_distribution == "Conical":
            rays[:,3],rays[:,5] = Conical2D.sample(N,
                                    self.__cone_max,
                                    self.__cone_min)
            rays[:,4] = numpy.sqrt(-rays[:,3]**2 - rays[:,5]**2 + 1.0)
        else:
            raise Exception("Bad value of angular_distribution")

        return rays

if __name__ == "__main__":

    from srxraylib.plot.gol import plot,plot_scatter

    # rectangle
    gs = SourceGeometrical()
    gs.set_spatial_type_rectangle(4,2)
    rays = gs.calculate_rays(1000)

    for i in range(6):
        print("Limits for column %d : %g,%g"%(1+i,rays[:,i].min(),rays[:,i].max()))
    # plot_scatter(rays[:,0],rays[:,2],xrange=[-3,3],yrange=[-3,3])


    # Gaussian
    gs = SourceGeometrical()
    gs.set_spatial_type_gaussian(2e-6,1e-6)
    rays = gs.calculate_rays(1000)

    for i in range(6):
        print("Std for column %d : %g"%(i+1,rays[:,i].std()))

    # Ellipse
    gs = SourceGeometrical()
    gs.set_spatial_type_ellipse(2,1)
    rays = gs.calculate_rays(1000)
    # plot_scatter(rays[:,0],rays[:,2],xrange=[-3,3],yrange=[-3,3])

    # Uniform divergence
    gs = SourceGeometrical()
    print(">>>>",rays[:,3].shape)
    gs.set_angular_distribution_uniform(-10e-6,5e-6,-3e-6,6e-6)
    rays = gs.calculate_rays(5000)
    # plot_scatter(1e6*rays[:,3],1e6*rays[:,5],xrange=[-10,10],yrange=[-10,10])

    # Gaussian divergence
    gs = SourceGeometrical()
    print(">>>>",rays[:,3].shape)
    gs.set_angular_distribution_gaussian(2e-5,1e-5)

    rays = gs.calculate_rays(5000)
    # plot_scatter(1e6*rays[:,3],1e6*rays[:,5],title="Gaussian div")

    # Conical divergence
    gs = SourceGeometrical()
    print(">>>>",rays[:,3].shape)
    gs.set_angular_distribution_conical(2e-5,1e-5)

    rays = gs.calculate_rays(5000)
    plot_scatter(1e6*rays[:,3],1e6*rays[:,5],title="Conical div")