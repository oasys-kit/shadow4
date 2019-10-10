
import numpy
from shadow4.beam.beam import Beam


class SourceGaussian(object):

    def __init__(self,
                 number_of_rays,
                 sigmaX,
                 sigmaY,
                 sigmaZ,
                 sigmaXprime,
                 sigmaZprime,
                 real_space_center,
                 direction_space_center
                 ):

        self._number_of_rays = number_of_rays
        self._sigmaX = sigmaX
        self._sigmaY = sigmaY
        self._sigmaZ = sigmaZ
        self._sigmaXprime = sigmaXprime
        self._sigmaZprime = sigmaZprime
        self._real_space_center = real_space_center
        self._direction_space_center = direction_space_center


    @classmethod
    def initialize_from_keywords(cls,
                 number_of_rays=10000,
                 sigmaX=1.0e-6,
                 sigmaY=1.0e-6,
                 sigmaZ=0.0,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0]
                                 ):
        return SourceGaussian(
                 number_of_rays,
                 sigmaX,
                 sigmaY,
                 sigmaZ,
                 sigmaXprime,
                 sigmaZprime,
                 real_space_center,
                 direction_space_center)

    @classmethod
    def initialize_point_source(cls,
                 number_of_rays=10000,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0] ):

        return SourceGaussian(
                 number_of_rays,
                 0.0,
                 0.0,
                 0.0,
                 sigmaXprime,
                 sigmaZprime,
                 real_space_center,
                 direction_space_center)


    @classmethod
    def initialize_collimated_source(cls,
                 number_of_rays=10000,
                 sigmaX=1.0,
                 sigmaY=0.0,
                 sigmaZ=1.0,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0]
                                 ):

        return SourceGaussian(
                 number_of_rays,
                 sigmaX,
                 sigmaY,
                 sigmaZ,
                 0.0,
                 0.0,
                 real_space_center,
                 direction_space_center)


    #
    # getters
    #
    def get_number_of_points(self):
        return self._number_of_rays

    def get_sigmas_real_space(self):
        return self._sigmaX,self._sigmaY,self._sigmaZ

    def get_sigmas_direction_space(self):
        return self._sigmaXprime,self._sigmaZprime

    def get_arrays_real_space(self):

        if self._sigmaX > 0.0:
            x = numpy.random.normal(self._real_space_center[0],self._sigmaX,self._number_of_rays)
        else:
            x = numpy.zeros(self._number_of_rays)

        if self._sigmaY > 0.0:
            y = numpy.random.normal(self._real_space_center[1],self._sigmaY,self._number_of_rays)
        else:
            y = numpy.zeros(self._number_of_rays)

        if self._sigmaZ > 0.0:
            z = numpy.random.normal(self._real_space_center[2],self._sigmaZ,self._number_of_rays)
        else:
            z = numpy.zeros(self._number_of_rays)

        return x,y,z

    def get_arrays_direction_space(self):
        if self._sigmaXprime > 0:
            x = numpy.random.normal(self._direction_space_center[0],self._sigmaXprime,self._number_of_rays)
        else:
            x = numpy.zeros(self._number_of_rays)

        if self._sigmaZprime > 0:
            z = numpy.random.normal(self._direction_space_center[1],self._sigmaZprime,self._number_of_rays)
        else:
            z = numpy.zeros(self._number_of_rays)
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
        txt += "Number of rays: %d \n"%self._number_of_rays
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
    def get_beam_shadow3(self,wavelength=1e-10):
        """
        Returns a shadow3 beam
        :param wavelength: the photon wavelength in m
        :return:
        """
        import Shadow
        beam = Shadow.Beam(self.get_number_of_points())
        beam.rays[:,0:6] = self.get_volume().T
        beam.rays[:,6:18] = 0.0
        beam.rays[:,6] = 1.0 # Es
        beam.rays[:,9] = 1   # flag
        beam.rays[:,10] = 2 * numpy.pi / (wavelength * 1e2) # wavenumber
        beam.rays[:,11] = numpy.arange(self.get_number_of_points(),dtype=float) # index
        return beam

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
        return Beam.initialize_from_array(rays)


if __name__ == "__main__":

    # a = SourceGaussian.initialize_point_source(number_of_rays=10000,
    #              sigmaXprime=3e-6,
    #              sigmaZprime=1e-6,
    #              real_space_center=[0.0,0.0,0.0],
    #              direction_space_center=[0.0,0.0])
    # print(a.info())
    #
    # x,y,z = a.get_arrays_real_space()
    # print("x.shape:",x.shape)
    # print("y.shape:",y.shape)
    # print("z.shape:",z.shape)
    #
    # xp,zp = a.get_arrays_direction_space()
    # print("xp.shape:",xp.shape)
    # print("zp.shape:",zp.shape)
    #
    #
    # VP = a.get_volume_divergences()
    # print("VP",VP.shape,VP.size)
    #
    #
    # Vx = a.get_volume_real_space()
    # print("Vx: ",Vx.shape)
    #
    # V = a.get_volume()
    # print("V: ",V.shape)
    #
    #
    # beam_shadow3 = a.get_beam_shadow3()
    # beam_shadow3.write("begin.dat")
    #
    # import Shadow
    # Shadow.ShadowTools.plotxy(beam_shadow3,4,6)
    #
    #
    #


    # #
    # #
    # a = SourceGaussian.initialize_collimated_source(number_of_rays=10000,
    #              sigmaX=1.0,
    #              sigmaY=0.0,
    #              sigmaZ=1.0,
    #              real_space_center=[0.0,0.0,0.0],
    #              direction_space_center=[0.0,0.0] )
    #
    # print(a.info())
    # beam_shadow3 = a.get_beam_shadow3()
    # import Shadow
    # Shadow.ShadowTools.plotxy(beam_shadow3,1,3)
    #
    # beam = a.get_beam()
    # from srxraylib.plot.gol import plot_scatter
    # plot_scatter(beam.get_column(1),beam.get_column(3))
    # print(beam.info())

    #
    #
    a = SourceGaussian.initialize_from_keywords(number_of_rays=10000,
                 sigmaX=0.0,
                 sigmaY=1.0e-3,
                 sigmaZ=0.0,
                 sigmaXprime=0.0,
                 sigmaZprime=0.0,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=[0.0,0.0] )

    print(a.info())

    beam = a.get_beam()
    from srxraylib.plot.gol import plot_scatter, set_qt
    set_qt()
    plot_scatter(beam.get_column(2),beam.get_column(1))
    print(beam.info())

    print(isinstance(a,SourceGaussian ))