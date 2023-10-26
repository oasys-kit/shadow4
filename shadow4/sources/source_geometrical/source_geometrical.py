import numpy
import scipy.constants as codata

from shadow4.beam.s4_beam import S4Beam

from shadow4.sources.source_geometrical.probability_distributions import Rectangle2D, Ellipse2D, Gaussian2D
from shadow4.sources.source_geometrical.probability_distributions import Flat2D, Uniform2D, Cone2D

from srxraylib.util.inverse_method_sampler import Sampler1D
from shadow4.sources.s4_light_source_base import S4LightSourceBase

class SourceGeometrical(S4LightSourceBase):
    def __init__(self,
                    name="Undefined",
                    spatial_type="Point",
                    angular_distribution = "Flat",
                    energy_distribution = "Single line",
                    depth_distribution = "Off",
                    nrays=5000,
                    seed=1234567,
                 ):
        super().__init__(name=name, nrays=nrays, seed=seed)
        self.set_spatial_type_by_name(spatial_type)  # see SourceGeometrical.spatial_type_list()
        self.set_angular_distribution_by_name(angular_distribution) # see SourceGeometrical.angular_distribution_list()
        self.set_energy_distribution_by_name(energy_distribution) # see SourceGeometrical.energy_distribution_list()
        self.set_depth_distribution_by_name(depth_distribution)
        self.set_polarization(polarization_degree=1.0, phase_diff=0.0, coherent_beam=0)


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
         # 		  gaussian(3), synchrotron(4), cone(5), exact
         # 		  synchrotron(6).
        return ["Flat","Uniform","Gaussian","Cone","Collimated"]

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
        elif name == "Cone":
            self.set_angular_distribution_cone()
        elif name == "Collimated":
            self.set_angular_distribution_collimated()
        else:
            raise Exception("Wrong angular distribution: %s"%name)


    def set_angular_distribution_collimated(self,):
        self.angular_distribution = "Flat"
        self.__hdiv1 = 0.0
        self.__hdiv2 = 0.0
        self.__vdiv1 = 0.0
        self.__vdiv2 = 0.0
        self.__fdist = 1

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

    def set_angular_distribution_cone(self,cone_max=10e-6,cone_min=0.0):
        self.angular_distribution = "Cone"
        self.__cone_max = cone_max
        self.__cone_min = cone_min
        self.__fdist = 4

    def set_polarization(self,
                         polarization_degree=1, # cos_s / (cos_s + sin_s)
                         phase_diff=0.000000,   # in rad, 0=linear, pi/2=elliptical/right
                         coherent_beam=0,       # 0=No (random phase s: col 14)  1=Yes (zero s: col 14)
                         ):
        self.__pol_deg = polarization_degree
        self.__pol_angle = phase_diff
        self.__f_foher = coherent_beam
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

    def set_depth_distribution_by_name(self, name, value=0.0):
        if name == "Off":
            self.set_depth_distribution(depth=0, source_depth_y=value)
        elif name == "Uniform":
            self.set_depth_distribution(depth=1, source_depth_y=value)
        elif name == "Gaussian":
            self.set_depth_distribution(depth=2, source_depth_y=value)
        else:
            raise Exception("Wrong depth distribution: %s"%name)

    def set_depth_distribution(self,
                               depth=0,             #0=off, 1=uniform, 2=Gaussian
                               source_depth_y=0.0,  #width if depth=1, sigma if depth=2
                               ):

        self.__fsource_depth = depth
        self.__wysou = source_depth_y

        if self.__fsource_depth == 0:
            self.depth_distribution = "Off"
        elif self.__fsource_depth == 1:
            self.depth_distribution = "Uniform"
        elif self.__fsource_depth == 2:
            self.depth_distribution = "Gaussian"
        else:
            raise Exception("Not implemented depth type (index=%d)" % self.__fsour)

    def set_depth_distribution_off(self):
        self.set_depth_distribution(0)

    def set_depth_distribution_uniform(self, value):
        self.set_depth_distribution(1, value)

    def set_depth_distribution_gaussian(self, value):
        self.set_depth_distribution(2, value)

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

    def calculate_beam(self):
        rays = self.calculate_rays()
        return S4Beam.initialize_from_array(rays)

    def calculate_rays(self):

        if self.get_seed() != 0:
            numpy.random.seed(self.get_seed())

        N = self.get_nrays()

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
        # depth distribution
        #
        print(">> Depth distribution: %s"%(self.depth_distribution))

        if self.depth_distribution == "Off":
            pass
        elif self.depth_distribution == "Uniform":
            rays[:,1] = (numpy.random.rand(N) - 0.5) * self.__wysou
        elif self.depth_distribution == "Gaussian":
            rays[:,1] = numpy.random.normal(loc=0.0, scale=self.__wysou, size=N)
        else:
            raise Exception("Bad value of depth_distribution")


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
        elif self.angular_distribution == "Cone":
            rays[:,3],rays[:,5] = Cone2D.sample(N,
                                    self.__cone_max,
                                    self.__cone_min)
            rays[:,4] = numpy.sqrt(-rays[:,3]**2 - rays[:,5]**2 + 1.0)
        else:
            raise Exception("Bad value of angular_distribution")



        #
        # energy distribution
        # ["Single line","Several lines","Uniform","Relative intensities","Gaussian","User defined"]
        #
        print(">> Energy distribution: ",self.energy_distribution)

        if self.energy_distribution == "Single line":
            if self.__f_phot == 0:
                rays[:,10] = self._energy_to_wavenumber(self.__ph[0])
            else:
                rays[:,10] = self._wavelength_to_wavenumber(self.__ph[0] * 1e-10)
        elif self.energy_distribution == "Several lines":
            values = numpy.array(self.__ph)
            n_test =   (numpy.random.random(N) * values.size).astype(int)
            sampled_values = values[n_test]
            if self.__f_phot == 0:
                rays[:,10] = self._energy_to_wavenumber(sampled_values)
            else:
                rays[:,10] = self._wavelength_to_wavenumber(sampled_values * 1e-10)
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
                rays[:,10] = self._wavelength_to_wavenumber(sampled_values * 1e-10)
        elif self.energy_distribution == "Uniform":
            sampled_values = self.__ph[0] + (self.__ph[1]-self.__ph[0]) * numpy.random.rand(N)
            if self.__f_phot == 0:
                rays[:,10] = self._energy_to_wavenumber(sampled_values)
            else:
                rays[:,10] = self._wavelength_to_wavenumber(sampled_values * 1e-10)
        elif self.energy_distribution == "Gaussian":
            sampled_values = numpy.random.normal(loc=self.__ph[0], scale=self.__ph[1], size=N)
            if self.__f_phot == 0:
                rays[:,10] = self._energy_to_wavenumber(sampled_values)
            else:
                rays[:,10] = self._wavelength_to_wavenumber(sampled_values * 1e-10)
        elif self.energy_distribution == "User defined":
            sampler = Sampler1D(self.__ph_spectrum_ordinates,self.__ph_spectrum_abscissas)
            sampled_values = sampler.get_n_sampled_points(N)
            # sampled_values, hy, hx = sampler.get_n_sampled_points_and_histogram(N)
            # plot(hx,hy)

            if self.__f_phot == 0:
                rays[:,10] = self._energy_to_wavenumber(sampled_values)
            else:
                rays[:,10] = self._wavelength_to_wavenumber(sampled_values * 1e-10)
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

        DENOM = numpy.sqrt(1.0 - 2.0 * self.__pol_deg + 2.0 * self.__pol_deg ** 2)
        AX = self.__pol_deg / DENOM
        for i in range(3):
            A_VEC[:, i] *= AX

        AZ = (1.0 - self.__pol_deg) / DENOM
        for i in range(3):
            AP_VEC[:, i] *= AZ


        rays[:, 6:9] = A_VEC
        rays[:, 15:18] = AP_VEC
        #
        # ! C
        # ! C Now the phases of A_VEC and AP_VEC.
        # ! C

        #

        if self.__f_foher == 1:
            PHASEX = 0.0
        else:
            PHASEX = numpy.random.random(N) * 2 * numpy.pi

        PHASEZ = PHASEX + self.__pol_angle

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

    def get_beam(self):
        rays =  self.calculate_rays()
        return S4Beam.initialize_from_array(rays)

    def to_python_code(self):

        txt = ""

        txt += "\n#\n#\n#"

        txt += "\nfrom shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical"

        txt += "\nlight_source = SourceGeometrical(name='%s', nrays=%d, seed=%d)" % \
               (self.get_name(), self.get_nrays(), self.get_seed())

        # spatial type
        if self.__fsour == 0:  # point
            txt += "\nlight_source.set_spatial_type_point()"
        elif self.__fsour == 1:  # rectangle
            txt += "\nlight_source.set_spatial_type_rectangle(width=%f,height=%f)" % \
                   (self.__wxsou, self.__wzsou)
        elif self.__fsour == 2:  # ellipse
            txt += "\nlight_source.set_spatial_type_ellipse(width=%f,height=%f)" % \
                   (self.__wxsou, self.__wzsou)
        elif self.__fsour == 3:  # Gaussian
            txt += "\nlight_source.set_spatial_type_gaussian(sigma_h=%s,sigma_v=%f)" % \
                   (self.__sigmax, self.__sigmaz)

        # depth
        if self.__fsource_depth == 0:  # off
            txt += "\nlight_source.set_depth_distribution_off()"
        elif self.__fsource_depth == 1:  # uniform
            txt += "\nlight_source.set_depth_distribution_uniform(%f)" % (self.__wysou)
        elif self.__fsource_depth == 2:  # gaussian
            txt += "\nlight_source.set_depth_distribution_gaussian(%f)" % (self.__wysou)

        # divergence
        if self.__fdist == 1:  # flat
            txt += "\nlight_source.set_angular_distribution_flat(hdiv1=%f,hdiv2=%f,vdiv1=%f,vdiv2=%f)" % \
                   (self.__hdiv1, self.__hdiv2, self.__vdiv1, self.__vdiv2)
        elif self.__fdist == 2:  # Uniform
            txt += "\nlight_source.set_angular_distribution_uniform(hdiv1=%f,hdiv2=%f,vdiv1=%f,vdiv2=%f)" % \
                   (self.__hdiv1, self.__hdiv2, self.__vdiv1, self.__vdiv2)

        elif self.__fdist == 3:  # Gaussian
            txt += "\nlight_source.set_angular_distribution_gaussian(sigdix=%f,sigdiz=%f)" % \
                   (self.__sigdix, self.__sigdiz)
        elif self.__fdist == 4:  # cone
            txt += "\nlight_source.set_angular_distribution_cone(cone_max=%f,cone_min=%f)" % \
                   (self.__cone_max, self.__cone_min)
        elif self.__fdist == 5:  # Zero (collimated) - New in shadow4
            txt += "\nlight_source.set_angular_distribution_collimated()"


        # energy
        unit = ['eV', 'A'][self.__f_phot]

        if self.__f_color == 1:  # "Single line":
            txt += "\nlight_source.set_energy_distribution_singleline(%f, unit='%s')" % \
                   (self.__ph[0], unit)
        elif self.__f_color == 2:  # "Several lines":
            nlines = (numpy.array(self.__ph)).size
            ff = "["
            for i in range(nlines):
                ff += "%f," % self.__ph[i]
            ff += "]"
            txt += "\nlight_source.set_energy_distribution_severallines(values=%s, unit='%s')" % (ff, unit)
        elif self.__f_color == 3:  # "Uniform":
            txt += "\nlight_source.set_energy_distribution_uniform(value_min=%f, value_max=%f, unit='%s')" % \
                   (self.__ph[0], self.__ph[1], unit)
        elif self.__f_color == 4:  # "Relative intensities":
            nlines = (numpy.array(self.__ph)).size
            ff = "["
            ww = "["
            for i in range(nlines):
                ff += "%f," % self.__ph[i]
                ww += "%f," % self.__rl[i]
            ff += "]"
            ww += "]"
            txt += "\nlight_source.set_energy_distribution_relativeintensities(values=%s, weights=%s, unit='%s')" % \
                   (ff, ww, unit)
        elif self.__f_color == 5:  # "Gaussian":
            txt += "\nlight_source.set_energy_distribution_gaussian(center=%f, sigma=%f, unit='%s')" % \
                (self.__ph[0], self.__ph[1], unit)
        elif self.__f_color == 6:  # "User defined":
            # a = numpy.loadtxt(self.user_defined_file)
            txt += "\nlight_source.set_energy_distribution_userdefined() #### TODO: COMPLETE" 

        #polarization/coherence
        txt += "\nlight_source.set_polarization(polarization_degree=%f, phase_diff=%f, coherent_beam=%s)" % \
               (self.__pol_deg, self.__pol_angle, self.__f_foher)

        txt += "\nbeam = light_source.get_beam()"

        return txt
if __name__ == "__main__":
    a = SourceGeometrical(depth_distribution="Gaussian")
    print(a.to_python_code())
