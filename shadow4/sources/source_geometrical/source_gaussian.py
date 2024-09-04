"""
This is a Gaussian source in both space and divergence coordinates.

It may be redundant with SourceGeometrical.

Notes
-----
Not used in OASYS. Not creating python scripts.

"""
import numpy
from shadow4.beam.s4_beam import S4Beam
from shadow4.sources.s4_light_source_base import S4LightSourceBase
from shadow4.tools.arrayofvectors import vector_default_efields

class SourceGaussian(S4LightSourceBase):
    """
    Defines a Gaussian source in 3D.

    Parameters
    ----------
    sigmaX : float, optional
        The sigma in X direction (width).
    sigmaY : float, optional
        The sigma in Y direction (depth).
    sigmaZ : float, optional
        The sigma in Z direction (height).
    sigmaXprime : float, optional
        The divergence in X direction in rad.
    sigmaZprime : float, optional
        The divergence in Z direction in rad.
    real_space_center : list, tuple or numpy array, optional
        The 3 coordinates of the center in real space.
    direction_space_center : list, tuple or numpy array, optional
        The 2 coordinates of the center in divergence space (X,Z).
    name : str, optional
        A name.
    nrays : int, optional
        Number of rays generated using SourceGaussian.get_beam()
    seed : int, optional
        Seed for the Monte Carlo generator.
    """
    def __init__(self,
                 name="Undefined",
                 nrays=5000,
                 seed=123456,
                 sigmaX=1.0e-6,
                 sigmaY=1.0e-6,
                 sigmaZ=0.0,
                 sigmaXprime=1e-6,
                 sigmaZprime=1e-6,
                 real_space_center=None,
                 direction_space_center=None,
                 ):
        if real_space_center is None:
            real_space_center = [0.0, 0.0, 0.0]

        if direction_space_center is None:
            direction_space_center = [0.0, 0.0]

        super().__init__(name=name, nrays=nrays, seed=seed)
        self._sigmaX = sigmaX
        self._sigmaY = sigmaY
        self._sigmaZ = sigmaZ
        self._sigmaXprime = sigmaXprime
        self._sigmaZprime = sigmaZprime
        self._real_space_center      = numpy.array(real_space_center)       # must be defined as numpy array to allow syned file i/o
        self._direction_space_center = numpy.array(direction_space_center)  # must be defined as numpy array to allow syned file i/o

        if seed != 0:
            numpy.random.seed(seed)

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
            ("sigmaX","The sigma in X direction (width)",""),
            ("sigmaY","The sigma in Y direction (depth)",""),
            ("sigmaZ","The sigma in Z direction (height)",""),
            ("sigmaXprime","The divergence in X direction in rad",""),
            ("sigmaZprime","The divergence in Z direction in rad",""),
            ("real_space_center","The 3 coordinates of the center in real space",""),
            ("direction_space_center","The 2 coordinates of the center in divergence space (X,Z)",""),
            ("name","A name",""),
            ("nrays"," Number of rays generated using SourceGaussian.get_beam(",""),
            ("seed","Seed for the Monte Carlo generator",""),
            ] )

    @classmethod
    def initialize_from_keywords(cls,
                                name="Undefined",
                                nrays=10000,
                                seed=12345,
                                sigmaX=1.0e-6,
                                sigmaY=1.0e-6,
                                sigmaZ=0.0,
                                sigmaXprime=1e-6,
                                sigmaZprime=1e-6,
                                real_space_center=None,
                                direction_space_center=None,
                                ):
        """
        Creates a Gaussian source in 3D.

        Parameters
        ----------
        sigmaX : float, optional
            The sigma in X direction (width).
        sigmaY : float, optional
            The sigma in Y direction (depth).
        sigmaZ : float, optional
            The sigma in Z direction (height).
        sigmaXprime : float, optional
            The divergence in X direction in rad.
        sigmaZprime : float, optional
            The divergence in Z direction in rad.
        real_space_center : list, optional
            The 3 coordinates of the center in real space.
        direction_space_center : list, optional
            The 2 coordinates of the center in divergence space (X,Z).
        name : str, optional
            A name.
        nrays : int, optional
            Number of rays generated using SourceGaussian.get_beam()
        seed : int, optional
            Seed for the Monte Carlo generator.

        Returns
        -------
        instance of SourceGaussian.
        """
        return SourceGaussian(
                                name                   = name,
                                nrays                  = nrays,
                                seed                   = seed,
                                sigmaX                 = sigmaX,
                                sigmaY                 = sigmaY,
                                sigmaZ                 = sigmaZ,
                                sigmaXprime            = sigmaXprime,
                                sigmaZprime            = sigmaZprime,
                                real_space_center      = None,
                                direction_space_center = None,
                                )

    @classmethod
    def initialize_point_source(cls,
                                name="Undefined",
                                nrays=10000,
                                seed=12345,
                                sigmaXprime=1e-6,
                                sigmaZprime=1e-6,
                                real_space_center=None,
                                direction_space_center=None,
                                ):
        """
        Creates a Gaussian source with zero dimension in real space.

        Parameters
        ----------
        real_space_center : list, optional
            The 3 coordinates of the center in real space.
        direction_space_center : list, optional
            The 2 coordinates of the center in divergence space (X,Z).
        name : str, optional
            A name.
        nrays : int, optional
            Number of rays generated using SourceGaussian.get_beam()
        seed : int, optional
            Seed for the Monte Carlo generator.

        Returns
        -------
        instance of SourceGaussian.
        """
        return SourceGaussian(
                 name                   = name,
                 nrays                  = nrays,
                 seed                   = seed,
                 sigmaX                 = 0.0,
                 sigmaY                 = 0.0,
                 sigmaZ                 = 0.0,
                 sigmaXprime            = sigmaXprime,
                 sigmaZprime            = sigmaZprime,
                 real_space_center      = real_space_center,
                 direction_space_center = direction_space_center,
        )


    @classmethod
    def initialize_collimated_source(cls,
                 name="Undefined",
                 nrays=10000,
                 seed=12345,
                 sigmaX=1.0,
                 sigmaY=0.0,
                 sigmaZ=1.0,
                 real_space_center=None,
                 direction_space_center=None,
                 ):
        """
        Creates a Gaussian source with zero dimension in divergence space space.

        Parameters
        ----------
        sigmaX : float, optional
            The sigma in X direction (width).
        sigmaY : float, optional
            The sigma in Y direction (depth).
        sigmaZ : float, optional
            The sigma in Z direction (height).
        name : str, optional
            A name.
        nrays : int, optional
            Number of rays generated using SourceGaussian.get_beam()
        seed : int, optional
            Seed for the Monte Carlo generator.

        Returns
        -------
        instance of SourceGaussian.
        """
        return SourceGaussian(
                                name                   = name,
                                nrays                  = nrays,
                                seed                   = seed,
                                sigmaX                 = sigmaX,
                                sigmaY                 = sigmaY,
                                sigmaZ                 = sigmaZ,
                                sigmaXprime            = 0.0,
                                sigmaZprime            = 0.0,
                                real_space_center      = real_space_center,
                                direction_space_center = direction_space_center,
                                )

    #
    # getters
    #
    def get_number_of_points(self):
        """
        Returns the number of rays.

        Returns
        -------
        int
        """
        return self.get_nrays()

    def get_sigmas_real_space(self):
        """
        Returns the sigmas in real space.

        Returns
        -------
        tuple
            (sigmaX, sigmaY, sigmaZ)
        """
        return self._sigmaX, self._sigmaY, self._sigmaZ

    def get_sigmas_direction_space(self):
        """
        Returns the sigmas in divergence space.

        Returns
        -------
        tuple
            (sigmaX', sigmaZ')
        """
        return self._sigmaXprime, self._sigmaZprime

    def _get_arrays_real_space(self):
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

    def _get_arrays_direction_space(self):
        if self._sigmaXprime > 0:
            x = numpy.random.normal(self._direction_space_center[0],self._sigmaXprime,self.get_number_of_points())
        else:
            x = numpy.zeros(self.get_number_of_points())

        if self._sigmaZprime > 0:
            z = numpy.random.normal(self._direction_space_center[1],self._sigmaZprime,self.get_number_of_points())
        else:
            z = numpy.zeros(self.get_number_of_points())
        return x,z



    def _get_volume_divergences(self):
        # Returns an array (3,npoints) with xp,yp,zp (first index 0,1,2, respectively) with the direction vectors
        XP,ZP = self._get_arrays_direction_space()
        YP = numpy.sqrt(1 - XP**2 - ZP**2 )
        tmp = numpy.vstack((XP.flatten(),YP.flatten(),ZP.flatten()))
        return tmp

    def _get_volume_real_space(self):
        # Returns an array (3,npoints) with x,y,z (first index 0,1,2, respectively) with the spatial coordinates
        X,Y,Z = self._get_arrays_real_space()
        return numpy.vstack((X.flatten(),Y.flatten(),Z.flatten()))

    def _get_volume(self):
        # Returns an array (6,npoints) with x,y,z,xp,yp,zp (first index 0,1,2,3,4,5 respectively) with the
        # spatial and direction coordinates
        if self.get_seed() != 0:
            numpy.random.seed(self.get_seed())

        v1 = self._get_volume_real_space()
        v2 = self._get_volume_divergences()

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
    def get_info(self):
        """
        Returns an array of strings with info.

        Returns
        -------
        str
        """

        f2dot35 = 2*numpy.sqrt(2*numpy.log(2))
        txt = ""

        txt += "SourceGaussian: \n"
        txt += "   Real space: \n"
        txt += "         Horizontal: sigmaX = %5.3f um; FWHM X = %5.3f um\n"%(1e6*self._sigmaX,1e6*f2dot35*self._sigmaX)
        txt += "         Vertical:   sigmaZ = %5.3f um; FWHM Z = %5.3f um\n"%(1e6*self._sigmaZ,1e6*f2dot35*self._sigmaZ)
        txt += "         Depth:      sigmaY = %5.3f um; FWHM Y = %5.3f um\n"%(1e6*self._sigmaY,1e6*f2dot35*self._sigmaY)
        txt += "   Direction space (divergences): \n"
        txt += "         Horizontal: sigmaXprime = %5.3f urad; FWHM X = %5.3f urad\n"%(1e6*self._sigmaXprime,1e6*f2dot35*self._sigmaXprime)
        txt += "         Vertical:   sigmaZprime = %5.3f urad; FWHM Z = %5.3f urad\n"%(1e6*self._sigmaZprime,1e6*f2dot35*self._sigmaZprime)

        return txt


    #
    # interfaces
    #

    def get_beam(self, wavelength=1e-10):
        """
        Returns an instance of S4Beam with the sampled rays.

        Parameters
        ----------
        wavelength : float, optional
            The photon wavelength in Angstroms.

        Returns
        -------
        instance of S4Beam

        """
        rays = numpy.zeros((self.get_number_of_points(),18))
        rays[:,0:6] = self._get_volume()
        # rays[:,6] = 1.0 # Es
        rays[:,9] = 1   # flag
        rays[:,10] = 2 * numpy.pi / (wavelength * 1e2) # wavenumber
        rays[:,11] = numpy.arange(self.get_number_of_points(),dtype=float) # index

        DIREC = rays[:,3:6]
        A_VEC, AP_VEC = vector_default_efields(DIREC)
        rays[:, 6:9] = A_VEC
        rays[:, 15:18] = AP_VEC

        return S4Beam.initialize_from_array(rays)


if __name__ == "__main__":
    src = SourceGaussian.initialize_from_keywords(
                 nrays=10000,
                 sigmaX=1.0e-6,
                 sigmaY=0.0,
                 sigmaZ=1.0e-6,
                 sigmaXprime=0.02,
                 sigmaZprime=0.02,
                 real_space_center=numpy.array([0.0,0.0,0.0]),
                 direction_space_center=[0.0,0.0]
                                 )
    print(src.info())

    beam = S4Beam()
    beam.generate_source(src)
    beam.set_photon_energy_eV(1000.0)
    print(beam.info())

    print("check orthogonality", beam.efields_orthogonal())