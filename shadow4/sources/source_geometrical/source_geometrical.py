import numpy
import scipy.constants as codata

from shadow4.sources.source_geometrical.probability_distributions import Rectangle2D, Ellipse2D, Gaussian2D
from shadow4.sources.source_geometrical.probability_distributions import Flat2D, Uniform2D, Conical2D

from srxraylib.util.inverse_method_sampler import Sampler1D


class SourceGeometrical(object):
    def __init__(self,
                    spatial_type="Point",
                    angular_distribution = "Flat",
                    energy_distribution = "Single line",
                 ):

        self.set_spatial_type_by_name(spatial_type)
        self.set_angular_distribution_by_name(angular_distribution)
        self.set_energy_distribution_by_name(energy_distribution)


    #
    # spatial type
    #
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
        self.__fsour = 0

    def set_spatial_type_rectangle(self,width=2.0,height=1.0):
        # wxsou	=  0.0000000000000000E+00 - for fsour=1,2; source width (X).
        # wzsou	=  0.0000000000000000E+00 - for fsour=1,2; source height (Z).
        self.spatial_type = "Rectangle"
        self.__wxsou = width
        self.__wzsou = height
        self.__fsour = 1

    def set_spatial_type_ellipse(self,width=2.0,height=1.0):
        # wxsou	=  0.0000000000000000E+00 - for fsour=1,2; source width (X).
        # wzsou	=  0.0000000000000000E+00 - for fsour=1,2; source height (Z).
        self.spatial_type = "Ellipse"
        self.__wxsou = width
        self.__wzsou = height
        self.__fsour = 2

    def set_spatial_type_gaussian(self,sigma_h=2.0,sigma_v=1.0):
        # sigmax	=  0.0000000000000000E+00 - for fsour=3; sigma in X
        # sigmaz	=  0.0000000000000000E+00 - for fsour=3; sigma in Z
        self.spatial_type = "Gaussian"
        self.__sigmax = sigma_h
        self.__sigmaz = sigma_v
        self.__fsour = 3

    #
    # angular distribution
    #
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
        self.__fdist = 1

    # WARNING: in shadow4 limits are signed!!!!
    def set_angular_distribution_uniform(self,hdiv1=-5e-6,hdiv2=5e-6,vdiv1=-0.5e-6,vdiv2=0.5e-6):
        self.angular_distribution = "Uniform"
        self.__hdiv1 = hdiv1
        self.__hdiv2 = hdiv2
        self.__vdiv1 = vdiv1
        self.__vdiv2 = vdiv2
        self.__fdist = 2

    def set_angular_distribution_gaussian(self,sigdix=1e-6,sigdiz=1e-6):
        self.angular_distribution = "Gaussian"
        self.__sigdix = sigdix
        self.__sigdiz = sigdiz
        self.__fdist = 3

    def set_angular_distribution_conical(self,cone_max=10e-6,cone_min=0.0):
        self.angular_distribution = "Conical"
        self.__cone_max = cone_max
        self.__cone_min = cone_min
        self.__fdist = 5

    #
    # energy distribution
    #


    #
    # f_color = 1 - photon energy distribution type.  Options are:
    #   single energy (1),
    #   multiple discrete energies, up to 10 energies (2),
    #   uniform energy distribution (3),
    #
    @classmethod
    def energy_distribution_list(cls):
        return ["Single line","Several lines","Uniform","Relative intensities","Gaussian","User defined"]

    def set_energy_distribution_by_name(self,name):
        if name == "Single line":
            self.set_energy_distribution_singleline()
        elif name == "Several lines":
            self.set_energy_distribution_severallines()
        elif name == "Uniform":
            self.set_energy_distribution_uniform()
        elif name == "Relative intensities":
            self.set_energy_distribution_relativeintensities()
        elif name == "Gaussian":
            self.set_energy_distribution_gaussian()
        elif name == "User defined":
            self.set_energy_distribution_userdefined()
        else:
            raise Exception("Wrong angular distribution: %s"%name)

    def _set_energy_distribution_unit(self,name='eV'):
        if name == 'eV':
            self.__f_phot = 0
        elif name == 'A':
            self.__f_phot = 1
        else:
            raise Exception("Bar name for energy units, valid names are: eV, A")

    def set_energy_distribution_singleline(self,value=1000,unit='eV'):
        # f_phot	=             0 - defines whether the photon energy will be
        # 			  specified in eV (0) or Angstroms (1).
        self.energy_distribution = "Single line"
        self._set_energy_distribution_unit(unit)
        #WARNING: ph1, ph2, etc. become ph (array)
        self.__ph = [value]
        self.__f_color = 1


    def set_energy_distribution_severallines(self,values=[1000.0,2000.0],unit='eV'):
        self.energy_distribution = "Several lines"
        self._set_energy_distribution_unit(unit)
        #WARNING: ph1, ph2, etc. become ph (array)
        self.__ph = values
        self.__f_color = 2

    def set_energy_distribution_uniform(self,value_min=1000.0,value_max=2000.0,unit='eV'):
        self.energy_distribution = "Uniform"
        self._set_energy_distribution_unit(unit)
        #WARNING: ph1, ph2, etc. become ph (array)
        self.__ph = [value_min,value_max]
        self.__f_color = 3

    def set_energy_distribution_relativeintensities(self,values=[1000.0,2000.0],weights=[1.0,2.0],unit='eV'):
        self.energy_distribution = "Relative intensities"
        self._set_energy_distribution_unit(unit)
        #WARNING: ph1, ph2, etc. become ph (array)
        self.__ph = values
        self.__rl = weights
        self.__f_color = 4 # not in source.nml

    # WARNING: limits suppressed
    def set_energy_distribution_gaussian(self,center=1000.0,sigma=10.0,unit='eV'):
        self.energy_distribution = "Gaussian"
        self._set_energy_distribution_unit(unit)
        #WARNING: ph1, ph2, etc. become ph (array)
        self.__ph = [center,sigma]
        self.__f_color = 5  # new

    def set_energy_distribution_userdefined(self,spectrum_abscissas,spectrum_ordinates,unit='eV'):
        self.energy_distribution = "User defined"
        self._set_energy_distribution_unit(unit)
        self.__ph_spectrum_abscissas = spectrum_abscissas
        self.__ph_spectrum_ordinates = spectrum_ordinates
        self.__f_color = 6  # new

    # conversors: wavenumber is in cm^(-1), wavenumber in m, energy in eV
    @classmethod
    def _energy_to_wavelength(cls, photon_energy):
        wavelength = codata.h * codata.c / codata.e / photon_energy
        return wavelength

    @classmethod
    def _wavelength_to_wavenumber(cls, wavelength):
        return 2 * numpy.pi / (wavelength * 1e2)

    @classmethod
    def _energy_to_wavenumber(cls, photon_energy):
        return cls._wavelength_to_wavenumber(cls._energy_to_wavelength(photon_energy))

    # the other way around
    @classmethod
    def _wavenumber_to_wavelength(cls, wavenumber):
        return 2 * numpy.pi / (wavenumber * 1e2)

    @classmethod
    def _wavelength_to_energy(cls, wavelength):
        return codata.h * codata.c / codata.e / wavelength

    @classmethod
    def _wavenumber_to_energy(cls,wavenumber):
        return cls._wavelength_to_energy(cls._wavenumber_to_wavelength(wavenumber))



    #
    # sampler
    #
    @classmethod
    def _sample_rays_default(cls,N=5000):
        rays = numpy.zeros((N,18))
        rays[:,4] = 1.0  # vy
        rays[:,6] = 1.0  # Ex
        rays[:,9] = 1.0  # flag
        rays[:,10] = 2 * numpy.pi / (1e-10 * 100) # wavenumber
        rays[:,11] = numpy.arange(N) + 1          # index
        return rays

    def calculate_rays(self,N=5000,POL_DEG=1.0,POL_ANGLE=0.0,F_COHER=False):

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

        #
        # energy distribution
        # ["Single line","Several lines","Uniform","Relative intensities","Gaussian","User defined"]
        #
        print(">> energy distribution: ",self.energy_distribution)

        if self.energy_distribution == "Single line":
            if self.__f_phot == 0:
                rays[:,10] = self._energy_to_wavenumber(self.__ph[0])
            else:
                rays[:,10] = self._wavelength_to_wavenumber(self.__ph[0])
        elif self.energy_distribution == "Several lines":
            values = numpy.array(self.__ph)
            n_test =   (numpy.random.random(N) * values.size).astype(int)
            sampled_values = values[n_test]
            if self.__f_phot == 0:
                rays[:,10] = self._energy_to_wavenumber(sampled_values)
            else:
                rays[:,10] = self._wavelength_to_wavenumber(sampled_values)
        elif self.energy_distribution == "Relative intensities":
            values = numpy.array(self.__ph)
            relative_intensities = numpy.array(self.__rl)
            # ! C
            # ! C Normalize so that each energy has a probability and so that the sum
            # ! C of the probabilities of all the energies is 1.
            # ! C
            relative_intensities /= relative_intensities.sum()
            # ! C
            # ! C Arrange the probabilities so that they comprise the (0,1) interval,
            # ! C e.g. (energy1,0.3), (energy2, 0.1), (energy3, 0.6) is translated to
            # ! C 0.0, 0.3, 0.4, 1.0. Then a random number falling in an interval
            # ! C assigned to a certain energy results in the ray being assigned that
            # ! C photon energy.
            # ! C

            TMP_B = 0
            for i in range(relative_intensities.size):
                TMP_B += relative_intensities[i]
                relative_intensities[i] = TMP_B

            sampled_values = numpy.zeros(N)
            for i in range(N):
                DPS_RAN3 = numpy.random.random()
                if (DPS_RAN3 > 0. and DPS_RAN3 <= relative_intensities[0]):
                    sampled_values[i] = values[0]

                for j in range(1,values.size):
                    if (DPS_RAN3 > relative_intensities[j-1] and DPS_RAN3 <= relative_intensities[j]):
                        sampled_values[i] = values[j]

            if self.__f_phot == 0:
                rays[:,10] = self._energy_to_wavenumber(sampled_values)
            else:
                rays[:,10] = self._wavelength_to_wavenumber(sampled_values)

        elif self.energy_distribution == "Gaussian":
            sampled_values = numpy.random.normal(loc=self.__ph[0], scale=self.__ph[1], size=N)
            if self.__f_phot == 0:
                rays[:,10] = self._energy_to_wavenumber(sampled_values)
            else:
                rays[:,10] = self._wavelength_to_wavenumber(sampled_values)
        elif self.energy_distribution == "User defined":
            sampler = Sampler1D(self.__ph_spectrum_ordinates,self.__ph_spectrum_abscissas)
            sampled_values = sampler.get_n_sampled_points(N)
            # sampled_values, hy, hx = sampler.get_n_sampled_points_and_histogram(N)
            # plot(hx,hy)

            if self.__f_phot == 0:
                rays[:,10] = self._energy_to_wavenumber(sampled_values)
            else:
                rays[:,10] = self._wavelength_to_wavenumber(sampled_values)
        else:
            raise Exception("Bad value of energy_distribution")

        #
        # polarization
        #
        # ! C
        # ! C  ---------------------------------------------------------------------
        # ! C                 POLARIZATION
        # ! C
        # ! C   Generates the polarization of the ray. This is defined on the
        # ! C   source plane, so that A_VEC is along the X-axis and AP_VEC is along Z-axis.
        # ! C   Then care must be taken so that A will be perpendicular to the ray
        # ! C   direction.
        # ! C
        # ! C
        # A_VEC(1) = 1.0D0
        # A_VEC(2) = 0.0D0
        # A_VEC(3) = 0.0D0

        DIREC = rays[:, 3:6].copy()
        A_VEC = numpy.zeros_like(DIREC)
        A_VEC[:, 0] = 1.0

        # ! C
        # ! C   Rotate A_VEC so that it will be perpendicular to DIREC and with the
        # ! C   right components on the plane.
        # ! C
        # CALL CROSS (A_VEC,DIREC,A_TEMP)
        A_TEMP = self._cross(A_VEC, DIREC)
        # CALL CROSS (DIREC,A_TEMP,A_VEC)
        A_VEC = self._cross(DIREC, A_TEMP)
        # CALL NORM (A_VEC,A_VEC)
        A_VEC = self._norm(A_VEC)
        # CALL CROSS (A_VEC,DIREC,AP_VEC)
        AP_VEC = self._cross(A_VEC, DIREC)
        # CALL NORM (AP_VEC,AP_VEC)
        AP_VEC = self._norm(AP_VEC)

        #
        # obtain polarization for each ray (interpolation)
        #

        DENOM = numpy.sqrt(1.0 - 2.0 * POL_DEG + 2.0 * POL_DEG ** 2)
        AX = POL_DEG / DENOM
        for i in range(3):
            A_VEC[:, i] *= AX

        AZ = (1.0 - POL_DEG) / DENOM
        for i in range(3):
            AP_VEC[:, i] *= AZ

        rays[:, 6:9] = A_VEC
        rays[:, 15:18] = AP_VEC

        #
        # ! C
        # ! C Now the phases of A_VEC and AP_VEC.
        # ! C

        #

        if F_COHER == 1:
            PHASEX = 0.0
        else:
            PHASEX = numpy.random.random(N) * 2 * numpy.pi

        PHASEZ = PHASEX + POL_ANGLE

        rays[:, 13] = PHASEX
        rays[:, 14] = PHASEZ

        # set flag (col 10)
        rays[:, 9] = 1.0

        print("<><>",dir(self))
        return rays

    def _cross(self,u,v):
        # w = u X v
        # u = array (npoints,vector_index)

        w = numpy.zeros_like(u)
        w[:,0] = u[:,1] * v[:,2] - u[:,2] * v[:,1]
        w[:,1] = u[:,2] * v[:,0] - u[:,0] * v[:,2]
        w[:,2] = u[:,0] * v[:,1] - u[:,1] * v[:,0]

        return w

    def _norm(self,u):
        # w = u / |u|
        # u = array (npoints,vector_index)
        u_norm = numpy.zeros_like(u)
        uu = numpy.sqrt( u[:,0]**2 + u[:,1]**2 + u[:,2]**2)
        for i in range(3):
            u_norm[:,i] = uu
        return u / u_norm

if __name__ == "__main__":

    from srxraylib.plot.gol import plot,plot_scatter, set_qt
    set_qt()

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
    # plot_scatter(1e6*rays[:,3],1e6*rays[:,5],title="Conical div")


    # energy monochromatic
    gs = SourceGeometrical()
    # gs.set_energy_distribution_singleline(2000.0,unit='eV')
    # gs.set_energy_distribution_severallines([1000.0,2000.0,3000.,8000],unit='eV')
    # gs.set_energy_distribution_relativeintensities([1000.0,2000.0,3000.,8000],[1.,2.,3,4],unit='eV')
    gs.set_energy_distribution_gaussian(10000.0,2000.0,unit='eV')

    rays = gs.calculate_rays(5000)
    ev = gs._wavenumber_to_energy(rays[:,10])
    # plot_scatter(rays[:,11],ev,title="Energy: xxx")


    # energy external spectrum
    gs = SourceGeometrical()
    x = numpy.linspace(1000.0,100000,2000)
    y = numpy.exp(- (x-50000)**2 / 2 / 10000**2 )
    # plot(x,y)
    gs.set_energy_distribution_userdefined(x,y,unit='eV')
    rays = gs.calculate_rays(5000)
    ev = gs._wavenumber_to_energy(rays[:,10])
    plot_scatter(rays[:,11],ev,title="Energy: sampled from numerical spectrum")

    print(dir(gs))

