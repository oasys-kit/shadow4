"""
This is a Gaussian source in both space and divergence coordinates.
It is redundant with source_geometrical
"""
import numpy
from shadow4.beam.s4_beam import S4Beam
from shadow4.sources.s4_light_source_base import S4LightSourceBase

class SourceGaussian(S4LightSourceBase):

    def __init__(self,
                 sigmaX,
                 sigmaY,
                 sigmaZ,
                 sigmaXprime,
                 sigmaZprime,
                 real_space_center,
                 direction_space_center,
                 name="Undefined",
                 nrays=5000,
                 seed=123456,
                 ):

        super().__init__(name=name, nrays=nrays, seed=seed)
        # self._number_of_rays = number_of_rays
        self._sigmaX = sigmaX
        self._sigmaY = sigmaY
        self._sigmaZ = sigmaZ
        self._sigmaXprime = sigmaXprime
        self._sigmaZprime = sigmaZprime
        self._real_space_center = real_space_center
        self._direction_space_center = direction_space_center
        # self._seed = seed

        if seed != 0:
            numpy.random.seed(seed)


    @classmethod
    def initialize_from_keywords(cls,
                 sigmaX=1.0e-6,
                 sigmaY=1.0e-6,
                 sigmaZ=0.0,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0],
                 nrays=10000,
                 seed=12345,
                 ):
        return SourceGaussian(
                 sigmaX,
                 sigmaY,
                 sigmaZ,
                 sigmaXprime,
                 sigmaZprime,
                 real_space_center,
                 direction_space_center,
                 nrays=nrays,
                 seed=seed)

    @classmethod
    def initialize_point_source(cls,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0],
                 nrays=10000,
                 seed=12345,
                 ):

        return SourceGaussian(
                 0.0,
                 0.0,
                 0.0,
                 sigmaXprime,
                 sigmaZprime,
                 real_space_center,
                 direction_space_center,
                 nrays=nrays,
                 seed=seed)


    @classmethod
    def initialize_collimated_source(cls,
                 sigmaX=1.0,
                 sigmaY=0.0,
                 sigmaZ=1.0,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0],
                 nrays = 10000,
                 seed=12345,
                 ):

        return SourceGaussian(
                 sigmaX,
                 sigmaY,
                 sigmaZ,
                 0.0,
                 0.0,
                 real_space_center,
                 direction_space_center,
                 nrays=nrays,
                 seed=seed)


    #
    # getters
    #
    def get_number_of_points(self):
        return self.get_nrays()

    def get_sigmas_real_space(self):
        return self._sigmaX,self._sigmaY,self._sigmaZ

    def get_sigmas_direction_space(self):
        return self._sigmaXprime,self._sigmaZprime

    def get_arrays_real_space(self):

        if self._sigmaX > 0.0:
            x = numpy.random.normal(self._real_space_center[0],self._sigmaX,self.get_number_of_points())
        else:
            x = numpy.zeros(self.get_number_of_points())

        if self._sigmaY > 0.0:
            y = numpy.random.normal(self._real_space_center[1],self._sigmaY,self.get_number_of_points())
        else:
            y = numpy.zeros(self.get_number_of_points())

        if self._sigmaZ > 0.0:
            z = numpy.random.normal(self._real_space_center[2],self._sigmaZ,self.get_number_of_points())
        else:
            z = numpy.zeros(self.get_number_of_points())

        return x,y,z

    def get_arrays_direction_space(self):
        if self._sigmaXprime > 0:
            x = numpy.random.normal(self._direction_space_center[0],self._sigmaXprime,self.get_number_of_points())
        else:
            x = numpy.zeros(self.get_number_of_points())

        if self._sigmaZprime > 0:
            z = numpy.random.normal(self._direction_space_center[1],self._sigmaZprime,self.get_number_of_points())
        else:
            z = numpy.zeros(self.get_number_of_points())
        return x,z



    def get_volume_divergences(self):
        """
        Returns an array (3,npoints) with xp,yp,zp (first index 0,1,2, respectively) with the
        direction vectors
        :return: xpypzp array
        """
        XP,ZP = self.get_arrays_direction_space()
        YP = numpy.sqrt(1 - XP**2 - ZP**2 )
        tmp = numpy.vstack((XP.flatten(),YP.flatten(),ZP.flatten()))
        return tmp

    def get_volume_real_space(self):
        """
        Returns an array (3,npoints) with x,y,z (first index 0,1,2, respectively) with the
        spatial coordinates
        :return: xyz
        """
        X,Y,Z = self.get_arrays_real_space()
        return numpy.vstack((X.flatten(),Y.flatten(),Z.flatten()))

    def get_volume(self):
        """
        Returns an array (6,npoints) with x,y,z,xp,yp,zp (first index 0,1,2,3,4,5 respectively) with the
        spatial and direction coordinates
        :return: xyzxpypzp array
        """

        if self.get_seed() != 0:
            numpy.random.seed(self.get_seed())

        v1 = self.get_volume_real_space()
        v2 = self.get_volume_divergences()

        V1x = v1[0,:].copy().flatten()
        V1y = v1[1,:].copy().flatten()
        V1z = v1[2,:].copy().flatten()
        V2x = v2[0,:].copy().flatten()
        V2y = v2[1,:].copy().flatten()
        V2z = v2[2,:].copy().flatten()

        return numpy.vstack((V1x,V1y,V1z,V2x,V2y,V2z)).T


    #
    # info
    #
    def info(self):
        """
        Returns an array of strings with info.
        :return:
        """
        f2dot35 = 2*numpy.sqrt(2*numpy.log(2))
        txt = ""

        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        txt += "Number of rays: %d \n"%self.get_number_of_points()
        txt += "Gausian source: \n"
        txt += "   Real space: \n"
        txt += "         Horizontal: sigmaX = %5.3f um; FWHM X = %5.3f um\n"%(1e6*self._sigmaX,1e6*f2dot35*self._sigmaX)
        txt += "         Vertical:   sigmaZ = %5.3f um; FWHM Z = %5.3f um\n"%(1e6*self._sigmaZ,1e6*f2dot35*self._sigmaZ)
        txt += "         Depth:      sigmaY = %5.3f um; FWHM Y = %5.3f um\n"%(1e6*self._sigmaY,1e6*f2dot35*self._sigmaY)
        txt += "   Direction space (divergences): \n"
        txt += "         Horizontal: sigmaXprime = %5.3f urad; FWHM X = %5.3f urad\n"%(1e6*self._sigmaXprime,1e6*f2dot35*self._sigmaXprime)
        txt += "         Vertical:   sigmaZprime = %5.3f urad; FWHM Z = %5.3f urad\n"%(1e6*self._sigmaZprime,1e6*f2dot35*self._sigmaZprime)
        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'n"

        return txt


    #
    # interfaces
    #

    def get_beam(self,wavelength=1e-10):
        """
        Returns a Beam
        :param wavelength: the photon wavelength in m
        :return:
        """

        rays = numpy.zeros((self.get_number_of_points(),18))
        rays[:,0:6] = self.get_volume()
        rays[:,6] = 1.0 # Es
        rays[:,9] = 1   # flag
        rays[:,10] = 2 * numpy.pi / (wavelength * 1e2) # wavenumber
        rays[:,11] = numpy.arange(self.get_number_of_points(),dtype=float) # index
        return S4Beam.initialize_from_array(rays)


if __name__ == "__main__":
    pass
