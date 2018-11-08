import numpy
from probability_distributions import Uniform1D


class SourceGeometrical(object):
    def __init__(self,
                    spatial_type="Point",
                    angular_distribution = "Flat",
                    energy_distribution = "Single line",
                 ):

        self.spatial_type = spatial_type
        self.angular_distribution = angular_distribution
        self.energy_distribution = energy_distribution

        self.rays = None

    @classmethod
    def spatial_type_list(cls):
        return ["Point","Rectangle","Ellipse","Gaussian"]

    @classmethod
    def angular_distribution_list(cls):
        return ["Flat","Uniform","Gaussian","Conical"]

    @classmethod
    def energy_distribution_list(cls):
        return ["Single line","Several lines","Uniform","Relative intensities","Gaussian","User defined"]

    def _sample_rays_default(self,N):
        self.rays = numpy.zeros((N,18))
        self.rays[:,4] = 1.0  # vy
        self.rays[:,6] = 1.0  # Ex
        self.rays[:,9] = 1.0  # flag
        self.rays[:,10] = 2 * numpy.pi / (1e-10 * 100) # wavenumber
        self.rays[:,11] = numpy.arange(N) + 1          # index

    def calculate_rays(self,N=1000):

        self._sample_rays_default(N)

        print(">> Spatial type: %s"%(self.spatial_type))

        if self.spatial_type == "Point":
            pass
        elif self.spatial_type == "Rectangle":
            self.rays[:,0] = Uniform1D.sample(N,-1.0,1.0)
            self.rays[:,2] = Uniform1D.sample(N,-2.0,2.0)
        elif self.spatial_type == "Ellipse":
            raise Exception(NotImplemented)
        elif self.spatial_type == "Gaussian":
            raise Exception(NotImplemented)
        else:
            raise Exception("Bad value of spatial_type")


if __name__ == "__main__":


    gs = SourceGeometrical(spatial_type="Rectangle")
    gs.calculate_rays(1000)

    for i in range(6):
        print("Limits for column %d : %f,%f"%(1+i,gs.rays[:,i].min(),gs.rays[:,i].max()))

