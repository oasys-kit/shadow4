"""
python version of the mirror reflectivity code and refractive index calculations in shadow.
"""
import numpy
import scipy.constants as codata
from shadow4.tools.logger import is_verbose, is_debug, set_verbose

tocm = codata.h * codata.c / codata.e * 1e2 # 12398.419739640718e-8

class PreRefl(object):

    def __init__(self):
        """
        Constructor.

        Get reflectivities using one of these options:
        * Using preprocessor file: ref = PreRefl(); ref.read_preprocessor_file('myprerefl.dat') ; ref.get_refraction_index() ; ref;get_attenuation_coefficient().
        * PreRefl.get_refraction_index_external_xraylib(<kwds>) PreRefl.get_attenuation_coefficient_external_xraylib(<kwds>)
        * PreRefl.get_refraction_index_external_dabax(<kwds>) PreRefl.get_attenuation_coefficient_external_dabax(<kwds>)

        """
        self.prerefl_dict = None

    def read_preprocessor_file(self, filename):
        """
        Reads a preprocessor (prerefl) file. The same as in shadow3.

        Parameters
        ----------
        filename : str
            The name of the filename.
        """
        fp = open(filename) # Open file on read mode
        lines = fp.read().split("\n") # Create a list containing all lines
        fp.close() # Close file

        index_pointer = 0
        mylist = lines[index_pointer].split("    ")
        QMIN = float(mylist[0])
        QMAX = float(mylist[1])
        QSTEP = float(mylist[2])
        DEPTH0 = float(mylist[3])


        index_pointer += 1
        NREFL = int(lines[index_pointer])

        ZF1 = numpy.zeros(NREFL)
        ZF2 = numpy.zeros(NREFL)
        for i in range(NREFL):
            index_pointer += 1
            ZF1[i] = float(lines[index_pointer])
        for i in range(NREFL):
            index_pointer += 1
            ZF2[i] = float(lines[index_pointer])

        self.prerefl_dict = {"QMIN":QMIN,"QMAX":QMAX,"QSTEP":QSTEP,"DEPTH0":DEPTH0,"NREFL":NREFL,"ZF1":ZF1,"ZF2":ZF2}

    def preprocessor_info(self):
        """
        Prints some information.
        """
        print("\n========================================")
        print("         preprocesor PreRefl info         ")
        for k in self.prerefl_dict.keys():
            try:
                print(k,self.prerefl_dict[k][0])
            except:
                print(k,self.prerefl_dict[k])

        print("QMIN: %f  EMIN: %f "%(self.prerefl_dict["QMIN"], self.prerefl_dict["QMIN"] * tocm / (2*numpy.pi)))
        print("QMAX: %f  EMAX: %f "%(self.prerefl_dict["QMAX"], self.prerefl_dict["QMAX"] * tocm / (2*numpy.pi)))
        print("========================================")

    def info(self):
        """
        Prints some information.
        """
        return self.preprocessor_info()

    def get_refraction_index(self, energy1):
        """
        Returns the complex refraction index.

        Parameters
        ----------
        energy1 : float or numpy array
            The photon energy or array of energies in eV.

        Returns
        -------
        complex or numpy array
            The complex refraction index.
        """

        wnum = 2 * numpy.pi * energy1 / tocm

        QSTEP = self.prerefl_dict["QSTEP"]
        QMIN = self.prerefl_dict["QMIN"]

        index1 = (wnum - QMIN) / QSTEP
        index1 = numpy.array(index1).astype(int)

        if (index1.max() + 1) > (self.prerefl_dict["NREFL"] - 1):
            raise Exception("Error: Photon energy above (or equal) tabulated upper limit.")

        if index1.min() < 0:
            raise Exception("Error: Photon energy below tabulated lower limit.")

        WNUM0 = QSTEP * index1 + QMIN
        DEL_X = wnum - WNUM0
        DEL_X = DEL_X / QSTEP

        ALFA = self.prerefl_dict["ZF1"][index1] + (self.prerefl_dict["ZF1"][index1 + 1] - self.prerefl_dict["ZF1"][index1]) * DEL_X
        GAMMA = self.prerefl_dict["ZF2"][index1] + (self.prerefl_dict["ZF2"][index1 + 1] - self.prerefl_dict["ZF2"][index1]) * DEL_X

        refraction_index = (1.0 - ALFA / 2) + (GAMMA / 2)*1j

        if is_verbose():
            try:
                print("   prerefl_test: calculates refraction index for a given energy")
                print("                 using a file created by the prerefl preprocessor.")
                print("   \n")

                print("------------------------------------------------------------------------")
                print("Inputs: ")
                print("   energy [eV]:                       ", wnum[0] * tocm / (2 * numpy.pi))
                print("   wavelength [A]:                    ", (1.0 / wnum[0]) * 2 * numpy.pi * 1e8)
                print("   wavenumber (2 pi/lambda) [cm^-1]:  ", wnum[0])
                # wnum = 2*numpy.pi * energy1 / tocm
                print("Outputs: ")
                print("   refraction index = (1-delta) + i*beta : ")
                print("   delta:                          ", 1.0 - refraction_index[0].real)
                print("   beta:                           ", refraction_index[0].imag)
                print("   real(n):                        ", refraction_index[0].real)
                print("   attenuation coef [cm^-1]:       ", 2 * refraction_index[0].imag * wnum[0])
                print("The given results correspond to ray index 0")

                print("------------------------------------------------------------------------")
            except:
                print("   prerefl_test: calculates refraction index for a given energy")
                print("                 using a file created by the prerefl preprocessor.")
                print("   \n")

                print("------------------------------------------------------------------------" )
                print("Inputs: " )
                print("   energy [eV]:                       ",wnum * tocm / (2 * numpy.pi))
                print("   wavelength [A]:                    ",(1.0 / wnum) * 2 * numpy.pi * 1e8 )
                print("   wavenumber (2 pi/lambda) [cm^-1]:  ",wnum )
                # wnum = 2*numpy.pi * energy1 / tocm
                print("Outputs: " )
                print("   refraction index = (1-delta) + i*beta : " )
                print("   delta:                          ", 1.0 - refraction_index.real )
                print("   beta:                           ", refraction_index.imag)
                print("   real(n):                        ", refraction_index.real )
                print("   attenuation coef [cm^-1]:       ", 2*refraction_index.imag*wnum )
                print("------------------------------------------------------------------------" )


        return refraction_index

    def get_attenuation_coefficient(self, energy1):
        """
        Returns the attenuation coefficient.

        Parameters
        ----------
        energy1 : float or numpy array
            The photon energy or array of energies in eV.

        Returns
        -------
        float or numpy array
            The attenuation coefficient in cm^-1.
        """
        refraction_index = self.get_refraction_index(energy1)
        wnum = 2 * numpy.pi * energy1 / tocm

        return 2 * refraction_index.imag * wnum

    @classmethod
    def get_attenuation_coefficient_external_xraylib(self,
                                                     photon_energy_ev=10000.0,
                                                     material="SiC",
                                                     density=3.217,):

        try:    import xraylib
        except: raise ImportError("xraylib not available")

        photon_energy_ev_array = numpy.array(photon_energy_ev)
        attenuation_coefficient = numpy.zeros_like(photon_energy_ev_array, dtype=float)

        for i,photon_energy_ev in enumerate(photon_energy_ev_array):
            attenuation_coefficient[i] = xraylib.CS_Total_CP(material, photon_energy_ev*1e-3) * density
        return attenuation_coefficient

    @classmethod
    def get_attenuation_coefficient_external_dabax(self,
                                                     photon_energy_ev=10000.0,
                                                     material="SiC",
                                                     density=3.217,
                                                     dabax=None,
                                                     ):
        """
        Standalone method to return the attenuation coefficient.

        Parameters
        ----------
        photon_energy_ev : float or numpy array
            The photon energy or array of energies in eV.
        material : str, optional
            The symbol/formula of the material.
        density : float, optional
            The material density in g/cm3.
        dabax : None or instance of DabaxXraylib
            A pointer to the dabax library. Use None for default.
        Returns
        -------
        float or numpy array
            The array with attenuation coefficient in cm^-1.
        """
        from dabax.dabax_xraylib import DabaxXraylib
        if isinstance(dabax, DabaxXraylib):
            dx = dabax
        else:
            dx = DabaxXraylib()

        return dx.CS_Total_CP(material, photon_energy_ev*1e-3) * density

    @classmethod
    def get_refraction_index_external_xraylib(self,
                                              photon_energy_ev=10000.0,
                                              material="SiC",
                                              density=3.217,):
        """
        Standalone method to return the complex refraction index using xraylib.

        Parameters
        ----------
        photon_energy_ev : float or numpy array
            The photon energy or array of energies in eV.
        material : str, optional
            The symbol/formula of the material.
        density : float, optional
            The material density in g/cm3.

        Returns
        -------
        float or numpy array
            The array with attenuation coefficient in cm^-1.
        """
        try:    import xraylib
        except: raise ImportError("xraylib not available")

        if isinstance(photon_energy_ev, (float, int)):
            return xraylib.Refractive_Index(material, photon_energy_ev * 1e-3, density)
        else:
            photon_energy_ev_array = numpy.array(photon_energy_ev)
            refraction_index = numpy.zeros_like(photon_energy_ev_array, dtype=complex)

            for i,photon_energy_ev in enumerate(photon_energy_ev_array):
                refraction_index[i] = xraylib.Refractive_Index(material, photon_energy_ev * 1e-3, density)
            return refraction_index

    @classmethod
    def get_refraction_index_real_external_xraylib(cls,
                                                     photon_energy_ev=10000.0,
                                                     material="SiC",
                                                     density=3.217,):
        """
        Standalone method to return the real part of the refraction index using xraylib.

        Parameters
        ----------
        photon_energy_ev : float or numpy array
            The photon energy or array of energies in eV.
        material : str, optional
            The symbol/formula of the material.
        density : float, optional
            The material density in g/cm3.

        Returns
        -------
        float or numpy array
            The array with attenuation coefficient in cm^-1.
        """

        return cls.get_refraction_index_external_xraylib(
                                                     photon_energy_ev=photon_energy_ev,
                                                     material=material,
                                                     density=density).real

    @classmethod
    def get_refraction_index_external_dabax(self,
                                            photon_energy_ev=10000.0,
                                            material="SiC",
                                            density=3.217,
                                            dabax=None):

        from dabax.dabax_xraylib import DabaxXraylib
        if isinstance(dabax, DabaxXraylib): dx = dabax
        else:                               dx = DabaxXraylib()

        return dx.Refractive_Index_Re(material, photon_energy_ev * 1e-3, density) + \
               1j * dx.Refractive_Index_Im(material, photon_energy_ev * 1e-3, density)

    @classmethod
    def get_refraction_index_real_external_dabax(self,
                                                 photon_energy_ev=10000.0,
                                                 material="SiC",
                                                 density=3.217,
                                                 dabax=None):
        """
        Standalone method to return the complex refraction index using dabax.

        Parameters
        ----------
        photon_energy_ev : float or numpy array
            The photon energy or array of energies in eV.
        material : str, optional
            The symbol/formula of the material.
        density : float, oprional
            The material density in g/cm3.
        dabax : None or instance of DabaxXraylib
            A pointer to the dabax library. Use None for default.

        Returns
        -------
        float or numpy array
            The array with attenuation coefficient in cm^-1.
        """
        from dabax.dabax_xraylib import DabaxXraylib
        if isinstance(dabax, DabaxXraylib): dx = dabax
        else:                               dx = DabaxXraylib()

        return dx.Refractive_Index_Re(material, photon_energy_ev * 1e-3, density)

    def reflectivity_fresnel(self,
                             photon_energy_ev=10000.0,
                             grazing_angle_mrad=3.0,
                             roughness_rms_A=0.0,
                             method=2,
                             ):
        """
        Calculates the reflectivity (intensity) of an interface using Fresnel formulas.

        Parameters
        ----------
        photon_energy_ev : float or numpy array
            The photon energy in eV.
        grazing_angle_mrad : float or numpy array
            The grazing incident angle in mrad.
        roughness_rms_A : float or numpy array
            The rouughness RMS in Angstroms,
        method : int, optional
            0=Born&Wolf, 1=Parratt, 2=shadow3  (avoid using 0 or 1, experimental!!)

        Returns
        -------
        tuple
            (rs, rp, runp) the s-polarized, p-pol and unpolarized reflectivities
        """
        rs, rp = self.reflectivity_amplitudes_fresnel(photon_energy_ev=photon_energy_ev,
                                                      grazing_angle_mrad=grazing_angle_mrad,
                                                      roughness_rms_A=roughness_rms_A,
                                                      method=method)

        return numpy.abs(rs)**2, numpy.abs(rp)**2, numpy.abs(0.5 * (rs + rp))**2,

    def reflectivity_amplitudes_fresnel(self,
                                        photon_energy_ev=10000.0,
                                        grazing_angle_mrad=3.0,
                                        roughness_rms_A=0.0,
                                        method=0 # 0=born & wolf, 1=parratt, 2=shadow3 (avoid using 0 or 1, experimental!!)
                                        ):
        """
        Calculates the reflectivity (amplitude) of an interface using Fresnel formulas.

        Parameters
        ----------
        photon_energy_ev : float or numpy array
            The photon energy in eV.
        grazing_angle_mrad : float or numpy array
            The grazing incident angle in mrad.
        roughness_rms_A : float or numpy array
            The rouughness RMS in Angstroms,
        method : int, optional
            0=Born&Wolf, 1=Parratt, 2=shadow3  (avoid using 0 or 1, experimental!!)

        Returns
        -------
        tuple
            (rs, rp, runp) the s-polarized, p-pol and unpolarized reflectivities
        """
        refraction_index_2 = self.get_refraction_index(photon_energy_ev)
        refraction_index_1 = numpy.ones_like(refraction_index_2)

        return PreRefl.reflectivity_amplitudes_fresnel_external(
                                                 photon_energy_ev=photon_energy_ev,
                                                 refraction_index_1=refraction_index_1,
                                                 refraction_index_2=refraction_index_2,
                                                 grazing_angle_mrad=grazing_angle_mrad,
                                                 roughness_rms_A=roughness_rms_A,
                                                 method=method)

    # to be used externally
    @classmethod
    def reflectivity_amplitudes_fresnel_external_xraylib(cls,
                                                         photon_energy_ev=10000.0,
                                                         coating_material="SiC",
                                                         coating_density=3.217,
                                                         grazing_angle_mrad=3.0,
                                                         roughness_rms_A=0.0,
                                                         method=0,
                                                         # 0=born & wolf, 1=parratt, 2=shadow3 (avoid using 0 or 1, experimental!!) ):
                                                         ):
        """
        Standalone method to calculate the reflectivity (amplitude) of an interface using Fresnel formulas using xraylib.

        Parameters
        ----------
        photon_energy_ev : float or numpy array
            The photon energy in eV.
        grazing_angle_mrad : float or numpy array
            The grazing incident angle in mrad.
        roughness_rms_A : float or numpy array
            The rouughness RMS in Angstroms,
        method : int, optional
            0=Born&Wolf, 1=Parratt, 2=shadow3  (avoid using 0 or 1, experimental!!)
        coating_material : int, optional
            The symbol/formula of the coating material.
        coating_density : int, optional
            The density in g/cm3 of the coating material.
        Returns
        -------
        tuple
            (rs, rp, runp) the s-polarized, p-pol and unpolarized reflectivities
        """
        try:    import xraylib
        except: raise ImportError("xraylib not available")

        photon_energy_ev_array = numpy.array(photon_energy_ev)
        refraction_index_2 = numpy.zeros_like(photon_energy_ev_array, dtype=complex)

        for i,photon_energy_ev in enumerate(photon_energy_ev_array):
            refraction_index_2[i] = xraylib.Refractive_Index_Re(coating_material, photon_energy_ev*1e-3, coating_density) + \
                                1j * xraylib.Refractive_Index_Im(coating_material, photon_energy_ev * 1e-3, coating_density)

        refraction_index_1 = numpy.ones_like(refraction_index_2)

        return PreRefl.reflectivity_amplitudes_fresnel_external(
                                                 photon_energy_ev=photon_energy_ev,
                                                 refraction_index_1=refraction_index_1,
                                                 refraction_index_2=refraction_index_2,
                                                 grazing_angle_mrad=grazing_angle_mrad,
                                                 roughness_rms_A=roughness_rms_A,
                                                 method=method)

    @classmethod
    def reflectivity_amplitudes_fresnel_external_dabax(cls,
                                                         photon_energy_ev=10000.0,
                                                         coating_material="SiC",
                                                         coating_density=3.217,
                                                         grazing_angle_mrad=3.0,
                                                         roughness_rms_A=0.0,
                                                         method=0, # 0=born & wolf, 1=parratt, 2=shadow3 (avoid using 0 or 1, experimental!!) ):
                                                         dabax=None,
                                                       ):
        """
        Standalone method to calculate the reflectivity (amplitude) of an interface using Fresnel formulas using dabax.

        Parameters
        ----------
        photon_energy_ev : float or numpy array
            The photon energy in eV.
        grazing_angle_mrad : float or numpy array
            The grazing incident angle in mrad.
        roughness_rms_A : float or numpy array
            The rouughness RMS in Angstroms,
        method : int, optional
            0=Born&Wolf, 1=Parratt, 2=shadow3  (avoid using 0 or 1, experimental!!)
        dabax : None or instance of DabaxXraylib
            A pointer to the dabax library. Use None for default.
        coating_material : int, optional
            The symbol/formula of the coating material.
        coating_density : int, optional
            The density in g/cm3 of the coating material.
        Returns
        -------
        tuple
            (rs, rp, runp) the s-polarized, p-pol and unpolarized reflectivities
        """
        from dabax.dabax_xraylib import DabaxXraylib
        if isinstance(dabax, DabaxXraylib):
            dx = dabax
        else:
            dx = DabaxXraylib()

        refraction_index_2 = dx.Refractive_Index_Re(coating_material, photon_energy_ev * 1e-3, coating_density) + \
                                1j * dx.Refractive_Index_Im(coating_material, photon_energy_ev * 1e-3, coating_density)

        refraction_index_1 = numpy.ones_like(refraction_index_2)

        return PreRefl.reflectivity_amplitudes_fresnel_external(
                                                 photon_energy_ev=photon_energy_ev,
                                                 refraction_index_1=refraction_index_1,
                                                 refraction_index_2=refraction_index_2,
                                                 grazing_angle_mrad=grazing_angle_mrad,
                                                 roughness_rms_A=roughness_rms_A,
                                                 method=method)


    @classmethod
    def reflectivity_amplitudes_fresnel_external(cls,
                                                 photon_energy_ev=10000.0,
                                                 refraction_index_1=1.0,
                                                 refraction_index_2=1.0,
                                                 grazing_angle_mrad=3.0,
                                                 roughness_rms_A=0.0,
                                                 method=0, # 0=born & wolf, 1=parratt, 2=shadow3
                                                 ):
        """
        Standalone method to calculate the reflectivity (amplitude) of an interface using Fresnel formulas using optical
        constants entered by user.

        This is the basic routine for the Fresnel equations.


        Parameters
        ----------
        photon_energy_ev : float or numpy array
            The photon energy in eV. This is only used for the Debye-Waller factor, therefore it is not used if
            roughness_rms_A=0.
        grazing_angle_mrad : float or numpy array
            The grazing incident angle in mrad.
        roughness_rms_A : float or numpy array
            The rouughness RMS in Angstroms,
        method : int, optional
            The equations used for the calculation:
            0=Born&Wolf [1]_ , 1=Parratt [2]_, 2=shadow3 [3]_
            Note that Parratt and shadow3 are only good for X-rays. Parratt gives bad p-pol.
            shadow3 does not calculate correct phase. Therefore, prefer Born&Wolf.
        coating_material : int, optional
            The symbol/formula of the coating material.
        coating_density : int, optional
            The density in g/cm3 of the coating material.
        refraction_index_1 : float or numpy array
            the refraction index (complex) for medium 1 (object).
        refraction_index_2 : float or numpy array
            the refraction index (complex) for medium 2 (image).

        Returns
        -------
        tuple
            (rs, rp) the s-polarized and p-pol and amplitude reflectivities

        References
    ----------
        .. [1] Born and Wolf "Principles of Optics" 6th edition, pag 40, eqs 21.
        .. [2] Parratt LG (1954) Surface studies of solids by total reflection of X-rays. Phys Rev 95(2):359–369.
        .. [3] Equations implemented in Shadow3.
        """

        theta1g = grazing_angle_mrad * 1e-3     # in rad
        theta1 = numpy.pi / 2 - theta1g
        rough1 = roughness_rms_A

        sin_theta1 = numpy.sin(theta1)
        cos_theta1 = numpy.sqrt(1 - sin_theta1**2, dtype=complex)

        sin_theta2 = refraction_index_1 * sin_theta1 / refraction_index_2
        cos_theta2 = numpy.sqrt(1 - sin_theta2**2, dtype=complex)


        if method == 0: # Born & Wolf 6th edition, pag 40
            if is_verbose(): print("Fresnel equations from B&W")
            rs1 = refraction_index_1 * cos_theta1 - refraction_index_2 * cos_theta2
            rs2 = refraction_index_1 * cos_theta1 + refraction_index_2 * cos_theta2
            rs = rs1 / rs2

            rp1 = refraction_index_2 * cos_theta1 - refraction_index_1 * cos_theta2
            rp2 = refraction_index_2 * cos_theta1 + refraction_index_1 * cos_theta2

            rp = rp1 / rp2

        elif method == 1: # parratt
            if is_verbose(): print("Fresnel equations from PARRATT")
            phi = theta1g
            delta1 = 1.0 - numpy.real(refraction_index_1)
            beta1 = numpy.imag(refraction_index_1)
            delta2 = 1.0 - numpy.real(refraction_index_2)
            beta2 = numpy.imag(refraction_index_2)
            f1 = numpy.sqrt(phi ** 2 - 2 * delta1 - 2 * 1j * beta1, dtype=complex)
            f2 = numpy.sqrt(phi ** 2 - 2 * delta2 - 2 * 1j * beta2, dtype=complex)
            rs = (f1 - f2) / (f1 + f2)
            rp = rs # TODO complete!

        elif method == 2:
            if is_verbose(): print("Fresnel equations from Shadow3")
            alpha = 2 * (1.0 - refraction_index_2.real)
            gamma = 2 * refraction_index_2.imag

            rho = (numpy.sin(theta1g)) ** 2 - alpha
            rho += numpy.sqrt((numpy.sin(theta1g) ** 2 - alpha) ** 2 + gamma ** 2)
            rho = numpy.sqrt(rho / 2)

            rs1 = 4 * (rho ** 2) * (numpy.sin(theta1g) - rho) ** 2 + gamma ** 2
            rs2 = 4 * (rho ** 2) * (numpy.sin(theta1g) + rho) ** 2 + gamma ** 2

            rs = rs1 / rs2

            ratio1 = 4 * rho ** 2 * (rho * numpy.sin(theta1g) - numpy.cos(theta1g) ** 2) ** 2 + gamma ** 2 * numpy.sin(
                theta1g) ** 2
            ratio2 = 4 * rho ** 2 * (rho * numpy.sin(theta1g) + numpy.cos(theta1g) ** 2) ** 2 + gamma ** 2 * numpy.sin(
                theta1g) ** 2
            ratio = ratio1 / ratio2

            rp = rs * ratio
            rs = numpy.sqrt(rs, dtype = complex)
            rp = numpy.sqrt(rp, dtype = complex)

        if roughness_rms_A != 0.0:
            wavelength_m = codata.h * codata.c / codata.e / photon_energy_ev
            debyewaller = numpy.exp( -(4.0 * numpy.pi * numpy.sin(theta1g) * rough1 / (wavelength_m * 1e10))**2 )
        else:
            debyewaller = 1.0

        # runp = 0.5 * (rs + rp)
        # return numpy.abs(rs)**2 * debyewaller, numpy.abs(rp)**2 * debyewaller, numpy.abs(runp)**2 * debyewaller

        return rs * numpy.sqrt(debyewaller), rp * numpy.sqrt(debyewaller)


    #
    # this is copied (and sligtly cleaned) from shadow3 python preprocessors
    #
    @classmethod
    def prerefl(cls,
                interactive=True,
                SYMBOL="SiC",
                DENSITY=3.217,
                FILE="prerefl.dat",
                E_MIN=100.0,
                E_MAX=20000.0,
                E_STEP=100.0,
                materials_library=None
                ):
        """
        Creates an instance of PreRefl with parameters initialized from the keyword parameters and the
        prerefl preprocessor file. It uses xraylib for accessing the optical constants.

        Parameters
        ----------
        interactive : bool, optional
            Set True for running interactively in the terminal and answer the questions (like in shadow2).
        SYMBOL : str, optional
            The material symbol/formula.
        DENSITY : float, optional
            The density of the material in g/cm3.
        FILE : str, optional
            The file name (output) for the prerefl preprocessor file.
        E_MIN : float, optional
            The minimum photon energy in eV (for creating the tabulated refraction index).
        E_MAX : float, optional
            The maximum photon energy in eV (for creating the tabulated refraction index).
        E_STEP : float optional
            The photon energy step in eV.

        Returns
        -------
        instance of PreRefl
        """

        if materials_library is None:
            try:    import xraylib as materials_library
            except:
                from dabax.dabax_xraylib import DabaxXraylib
                materials_library = DabaxXraylib()

        if interactive:
            # input section
            print("prerefl: Preprocessor for mirrors - python+xraylib version")
            iMaterial = input("Enter material expression (symbol,formula): ")
            density = input("Density [ g/cm3 ] ?")
            density = float(density)

            estart = input("Enter starting photon energy: ")
            estart = float(estart)

            efinal = input("Enter end photon energy: ")
            efinal = float(efinal)

            estep = input("Enter step photon energy:")
            estep = float(estep)

            out_file = input("Output file : ")
        else:
            iMaterial = SYMBOL
            density = DENSITY
            estart = E_MIN
            efinal = E_MAX
            estep = E_STEP
            out_file = FILE

        twopi = numpy.pi * 2
        npoint = int((efinal - estart) / estep + 1)
        depth0 = density / 2.0
        qmin = estart / tocm * twopi
        qmax = efinal / tocm * twopi
        qstep = estep / tocm * twopi

        f = open(out_file, 'wt')
        f.write(("%20.11e " * 4 + "\n") % tuple([qmin, qmax, qstep, depth0]))
        f.write("%i \n" % int(npoint))
        for i in range(npoint):
            energy = (estart + estep * i) * 1e-3
            tmp = 2e0 * (1e0 - materials_library.Refractive_Index_Re(iMaterial, energy, density))
            f.write("%e \n" % tmp)
        for i in range(npoint):
            energy = (estart + estep * i) * 1e-3
            tmp2 = 2e0 * (materials_library.Refractive_Index_Im(iMaterial, energy, density))
            f.write("%e \n" % tmp2)
        print("File written to disk: %s" % out_file)
        f.close()

        # test (not needed)
        itest = 0
        if itest:
            cdtest = materials_library.CompoundParser(iMaterial)
            print("    ", iMaterial, " contains %i atoms and %i elements" % (cdtest['nAtomsAll'], cdtest['nElements']))
            for i in range(cdtest['nElements']):
                print("    Element %i: %lf %%" % (cdtest['Elements'][i], cdtest['massFractions'][i] * 100.0))
            print(qmin, qmax, qstep, depth0)
            print(npoint)
            for i in range(npoint):
                energy = (estart + estep * i) * 1e-3
                qq = qmin + qstep * i
                print(energy, qq, \
                      2e0 * (1e0 - materials_library.Refractive_Index_Re(iMaterial, energy, density)), \
                      2e0 * (materials_library.Refractive_Index_Im(iMaterial, energy, density)))

        return None

    #
    # this will create prerefl file from a file prepared at https://henke.lbl.gov/tmp/xray8378.dat
    #
    @classmethod
    def prerefl_cxro(cls, input_file="https://henke.lbl.gov/tmp/xray8378.dat", output_file="prerefl.dat",
                     E_MIN=None, E_MAX=None, NPOINTS=1000):
        """
        Creates an instance of PreRefl with parameters initialized from the keyword parameters and a file
        with energy, alpha and beta created by the CXRO web site https://henke.lbl.gov/optical_constants/getdb2.html

        When running the cxro site, request a "text file" and copy the link below "If your file does not appear, click here."



        Parameters
        ----------
        input_file : str, optional
            The URL from where the data is downloaded.
        output_file : str, optional
            The file name (output) for the prerefl preprocessor file.
        E_MIN : float, optional
            The minimum photon energy in eV (for creating the tabulated refraction index).
        E_MAX : float, optional
            The maximum photon energy in eV (for creating the tabulated refraction index).
        NPOINTS : int optional
            The number of points for the photon energy array.

        Returns
        -------
        instance of PreRefl
        """
        a = numpy.loadtxt(input_file, skiprows=2)
        if a.shape[1] != 3:
            raise Exception("Bar input. Take output URL from https://henke.lbl.gov/optical_constants/getdb2.html ")

        energy0 = a[:,0]
        delta0  = a[:,1]
        beta0   = a[:,2]

        if E_MIN is None:
            E_MIN = energy0.min()

        if E_MAX is None:
            E_MAX = energy0.max()

        if energy0[0] > E_MIN:
            raise Exception("File min(energy) = %g greater than limit %g, cannot interpolate" % (energy0[0], E_MIN))

        if energy0[-1] < E_MAX:
            raise Exception("File max(energy) = %g smaller than limit %g, cannot interpolate" % (energy0[-1], E_MAX))

        # read density from header
        if "http" in input_file:
            import urllib.request as urllib
            fp = urllib.urlopen(input_file)
            lines = fp.readlines()
            mylist = lines[0].decode('utf-8').split("=")
            density = float(mylist[1])
        else:
            fp = open(input_file) # Open file on read mode
            lines = fp.read().split("\n") # Create a list containing all lines
            fp.close() # Close file
            mylist = lines[0].split("=")
            density = float(mylist[1])

        energy = numpy.linspace(E_MIN, E_MAX, int(NPOINTS))
        delta = numpy.interp(energy, energy0, delta0)
        beta  = numpy.interp(energy, energy0, beta0)

        twopi = numpy.pi * 2

        npoint = energy.size
        depth0 = density / 2.0

        qmin = energy[0] / tocm * twopi
        qmax = energy[-1] / tocm * twopi
        qstep = (energy[1] - energy[0]) / tocm * twopi

        f = open(output_file, 'wt')
        f.write(("%20.11e " * 4 + "\n") % tuple([qmin, qmax, qstep, depth0]))
        f.write("%i \n" % int(npoint))
        for i in range(npoint):
            tmp = 2e0 * delta[i]
            f.write("%e \n" % tmp)
        for i in range(npoint):
            tmp2 = 2e0 * beta[i]
            f.write("%e \n" % tmp2)
        print("File written to disk: %s" % output_file)
        f.close()

        return None


    #
    # this will create prerefl file from a file from https://refractiveindex.info
    # (it needs that the python package refractiveindex ( https://pypi.org/project/refractiveindex ) is installed)
    #
    @classmethod
    def prerefl_refractiveindexinfo(cls,
                    shelf='main', book='SiO2', page='Franta',
                    output_file="prerefl.dat",
                    WAVELENGTH_MIN=100e-9, WAVELENGTH_MAX=500e-9, NPOINTS=1000, density=2.2,
                    use_absorption=1, invert_refraction_index=0):
        """
        Creates an instance of PreRefl with parameters initialized from the keyword parameters and data
        from https://refractiveindex.info/ .

        Parameters
        ----------
        shelf : str
            The shelf name (see [1]_ and [2]_).
        book : str
            The book name (see [1]_ and [2]_).
        page : str
            The page name (see [1]_ and [2]_).
        output_file : str, optional
            The file name (output) for the prerefl preprocessor file.
        WAVELENGTH_MIN : float, optional
            The minimum value for the wavelength in m.
        WAVELENGTH_MAX: float, optional
            The maximum value for the wavelength in m.
        NPOINTS : int, optional
            The number of points.
        density : float, optional
            The density in g/cm^3.
        use_absorption : int, optional
            Set to zero for using zero absorption (zero imaginary part of refraction index).
        invert_refraction_index : int, optional
            Set to 1 to invert the refraction index (to be used for example, for total internal reflection)

        Returns
        -------
        instance of PreRefl

        References
        ----------
        .. [1] https://github.com/toftul/refractiveindex
        .. [2] https://refractiveindex.info/

        """

        try: from refractiveindex import RefractiveIndexMaterial
        except: raise ImportError("refractiveindex not available")

        depth0 = density / 2.0

        mat = RefractiveIndexMaterial(shelf=shelf, book=book, page=page)

        qmax = 2 * numpy.pi / (WAVELENGTH_MIN * 100)  # in cm^-1
        qmin = 2 * numpy.pi / (WAVELENGTH_MAX * 100)  # in cm^-1
        QARRAY = numpy.linspace(qmin, qmax, NPOINTS)

        N = []
        K = []
        for i in range(NPOINTS):
            wavelength_nm = 2 * numpy.pi / QARRAY[i] * 1e7
            if invert_refraction_index:
                N.append(1.0 / mat.get_refractive_index(wavelength_nm))
            else:
                N.append(mat.get_refractive_index(wavelength_nm))
            if use_absorption:
                K.append(mat.get_extinction_coefficient(wavelength_nm))
            else:
                K.append(0.0)

        delta = 1 - numpy.array(N)
        beta  = numpy.array(K)
        qstep = numpy.abs(QARRAY[1] - QARRAY[0])


        f = open(output_file, 'wt')
        f.write(("%20.11e " * 4 + "\n") % tuple([qmin, qmax, qstep, depth0]))
        f.write("%i \n" % int(NPOINTS))
        for i in range(NPOINTS):
            tmp = 2e0 * delta[i]
            f.write("%e \n" % tmp)
        for i in range(NPOINTS):
            tmp2 = 2e0 * beta[i]
            f.write("%e \n" % tmp2)
        print("File written to disk: %s" % output_file)
        f.close()

        return None

if __name__ == "__main__":
    from srxraylib.plot.gol import plot

    # set_verbose()

    if False:
        from srxraylib.plot.gol import plot

        #
        # refractive index
        #


        prerefl_file = "Be5_30.dat"
        PreRefl.prerefl(interactive=False, SYMBOL="Be", DENSITY=1.848, FILE=prerefl_file,
                        E_MIN=5000.0, E_MAX=30000.0, E_STEP=100.0)



        a = PreRefl()


        a.read_preprocessor_file(prerefl_file)

        refraction_index = a.get_refraction_index(10000.0)

        a.preprocessor_info()

    if False: # energy scan
        from srxraylib.plot.gol import plot

        #
        # mirror reflectivity
        #

        a = PreRefl()

        prerefl_file = "Rh5_50.dat"
        PreRefl.prerefl(interactive=False, SYMBOL="Rh", DENSITY=12.41, FILE=prerefl_file,
                        E_MIN=5000.0, E_MAX=50000.0, E_STEP=100.0)

        a.read_preprocessor_file(prerefl_file)

        Energy = numpy.linspace(5000.0,40000.0,100)

        RS0 = numpy.zeros_like(Energy)
        RS5 = numpy.zeros_like(Energy)

        a.preprocessor_info()

        # array inputs

        Grazing_angle = numpy.ones_like(Energy) * 3.0

        rs0_m0, rp0_m0 = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=Grazing_angle, photon_energy_ev=Energy,
                                                     roughness_rms_A=0.0, method=0)
        rs1_m0, rp1_m0 = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=Grazing_angle, photon_energy_ev=Energy,
                                                     roughness_rms_A=5.0, method=0)

        rs0_m1, rp0_m1 = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=Grazing_angle, photon_energy_ev=Energy,
                                                     roughness_rms_A=0.0, method=1)
        rs1_m1, rp1_m1 = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=Grazing_angle, photon_energy_ev=Energy,
                                                     roughness_rms_A=5.0, method=1)

        rs0_m2, rp0_m2 = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=Grazing_angle, photon_energy_ev=Energy,
                                                     roughness_rms_A=0.0, method=2)
        rs1_m2, rp1_m2 = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=Grazing_angle, photon_energy_ev=Energy,
                                                     roughness_rms_A=5.0, method=2)

        plot(Energy, numpy.abs(rs0_m0) ** 2, Energy, numpy.abs(rs1_m0) ** 2,
             Energy, numpy.abs(rs0_m1) ** 2, Energy, numpy.abs(rs1_m1) ** 2,
             Energy, numpy.abs(rs0_m2) ** 2, Energy, numpy.abs(rs1_m2) ** 2,
             ylog=True,
             legend=["method 0 - no roughness", "method 0 - 5A RMS roughness",
                     "method 1 - no roughness", "method 1 - 5A RMS roughness",
                     "method 2 - no roughness", "method 2 - 5A RMS roughness"], title="Rh mirror 3 mrad Energy Scan")


    if False: # theta scan
        from srxraylib.plot.gol import plot

        #
        # mirror reflectivity vs angle
        #


        angle_mrad = numpy.linspace(0.0, 10.0, 100)

        RS0 = numpy.zeros_like(angle_mrad)
        RS5 = numpy.zeros_like(angle_mrad)


        rs0_m0, rp0_m0 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=0.0, method=0,
                                refraction_index_1=1.0, refraction_index_2=0.999995+9.88072e-08j)
        rs1_m0, rp1_m0 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=5.0, method=0,
                                refraction_index_1=1.0,
                                refraction_index_2=0.999995 + 9.88072e-08j)

        rs0_m1, rp0_m1 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=0.0, method=1,
                                refraction_index_1=1.0, refraction_index_2=0.999995+9.88072e-08j)
        rs1_m1, rp1_m1 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=5.0, method=1,
                                refraction_index_1=1.0,
                                refraction_index_2=0.999995 + 9.88072e-08j)

        rs0_m2, rp0_m2 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=0.0, method=2,
                                refraction_index_1=1.0, refraction_index_2=0.999995+9.88072e-08j)
        rs1_m2, rp1_m2 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=5.0, method=2,
                                refraction_index_1=1.0,
                                refraction_index_2=0.999995 + 9.88072e-08j)

        plot(angle_mrad, numpy.abs(rs0_m0) ** 2, angle_mrad, numpy.abs(rs1_m0) ** 2,
             angle_mrad, numpy.abs(rs0_m1) ** 2, angle_mrad, numpy.abs(rs1_m1) ** 2,
             angle_mrad, numpy.abs(rs0_m2) ** 2, angle_mrad, numpy.abs(rs1_m2) ** 2,
             ylog=True, xtitle="Grazing angle [mrad]",
             legend=["method 0 - no roughness", "method 0 - 5A RMS roughness",
                     "method 1 - no roughness", "method 1 - 5A RMS roughness",
                     "method 2 - no roughness", "method 2 - 5A RMS roughness"],
             marker=["+", "+", None, None, None, None],
             linestyle=["", "", None, None, None, None],
             title="Rh mirror @ 20 keV angle scan; external refraction index",
             show=0)

        plot(
             angle_mrad, numpy.abs(rs0_m0) ** 2, angle_mrad, numpy.abs(rp0_m0) ** 2,
             angle_mrad, numpy.abs(rs0_m1) ** 2, angle_mrad, numpy.abs(rp0_m1) ** 2,
             angle_mrad, numpy.abs(rs0_m2) ** 2, angle_mrad, numpy.abs(rp0_m2) ** 2,
             ylog=0, xtitle="Incidence angle (grazing) [mrad]",
             legend=["method 0 - s", "method 0 - p",
                     "method 1 - s", "method 1 - p !!!",
                     "method 2 - s", "method 2 - p"],
             marker=["+","+",None,None,None,None],
             linestyle=["","",None,None,None,None],
             title="Rh mirror @ 20 keV angle scan; external refraction index", yrange=[0,1],
             show=0)

        # TODO: there are many differences in phase
        plot(
             angle_mrad, numpy.angle(rs0_m0), angle_mrad, numpy.angle(rp0_m0),
             angle_mrad, numpy.angle(rs0_m1), angle_mrad, numpy.angle(rp0_m1),
             # angle_mrad, numpy.angle(rs0_m2), angle_mrad, numpy.angle(rp0_m2),
             ylog=0, xtitle="Incidence angle (grazing) [mrad]", ytitle="Phase angle [rad]",
             legend=["method 0 - s", "method 0 - p",
                     "method 1 - s", "method 1 - p !!!",
                     "method 2 - s", "method 2 - p"],
             marker=["x","+",None,None,None,None],
             linestyle=["","",None,None,None,None],
             title="Rh mirror @ 20 keV angle scan; external refraction index")

    if False: # theta scan GLASS https://www.rp-photonics.com/fresnel_equations.html , see also fig 1.12 (pag 44) B&W 6th ed.
        from srxraylib.plot.gol import plot

        #
        # mirror reflectivity vs angle
        #


        angle_deg = numpy.linspace(0.0, 90.0, 100)
        angle_mrad = 1e3 * numpy.radians(90.0 - angle_deg)

        RS0 = numpy.zeros_like(angle_deg)
        RS5 = numpy.zeros_like(angle_deg)


        rs0_m0, rp0_m0 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=0.0, method=0,
                                refraction_index_1=1.0, refraction_index_2=1.52)
        rs1_m0, rp1_m0 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=5.0, method=0,
                                refraction_index_1=1.0, refraction_index_2=1.52)

        rs0_m1, rp0_m1 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=0.0, method=1,
                                refraction_index_1=1.0, refraction_index_2=1.52)
        rs1_m1, rp1_m1 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=5.0, method=1,
                                refraction_index_1=1.0, refraction_index_2=1.52)

        rs0_m2, rp0_m2 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=0.0, method=2,
                                refraction_index_1=1.0, refraction_index_2=1.52)
        rs1_m2, rp1_m2 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=5.0, method=2,
                                refraction_index_1=1.0, refraction_index_2=1.52)

        plot(
             angle_deg, numpy.abs(rs0_m0) ** 2, angle_deg, numpy.abs(rp0_m0) ** 2,
             angle_deg, numpy.abs(rs0_m1) ** 2, angle_deg, numpy.abs(rp0_m1) ** 2,
             angle_deg, numpy.abs(rs0_m2) ** 2, angle_deg, numpy.abs(rp0_m2) ** 2,
             ylog=0, xtitle="Incidence angle (to normal) [deg]",
             legend=["method 0 - s", "method 0 - p",
                     "method 1 - s", "method 1 - p !!!",
                     "method 2 - s", "method 2 - p"],
             marker=["+","+",None,None,None,None],
             linestyle=["","",None,None,None,None],
             title="GLASS n=1.52", xrange=[0,90], yrange=[0,1],
             show=0)

        # degree of polarization (op.cit. pag 44 eq 42)
        P_m0 = numpy.abs((numpy.abs(rp0_m0) ** 2 - numpy.abs(rs0_m0) ** 2) / (numpy.abs(rp0_m0) ** 2 + numpy.abs(rs0_m0) ** 2))
        P_m1 = numpy.abs((numpy.abs(rp0_m1) ** 2 - numpy.abs(rs0_m1) ** 2) / (numpy.abs(rp0_m1) ** 2 + numpy.abs(rs0_m1) ** 2))
        P_m2 = numpy.abs((numpy.abs(rp0_m2) ** 2 - numpy.abs(rs0_m2) ** 2) / (numpy.abs(rp0_m2) ** 2 + numpy.abs(rs0_m2) ** 2))

        print(">>>>", P_m0)
        plot(
             angle_deg, P_m0,
             angle_deg, P_m1,
             angle_deg, P_m2,
             ylog=0, xtitle="Incidence angle (to normal) [deg]",
             legend=["method 0",
                     "method 1 !!!",
                     "method 2"],
             marker=["+",None,None],
             linestyle=["",None,None],
             title="GLASS n=1.52 POLARIZATION DEGREE", xrange=[0,90],
             show=0)

        # phase diff (op.cit. pag 50)
        delta_m0 = numpy.angle(rs0_m0) - numpy.angle(rp0_m0)
        delta_m1 = numpy.angle(rs0_m1) - numpy.angle(rp0_m1)
        delta_m2 = numpy.angle(rs0_m2) - numpy.angle(rp0_m2)

        plot(
             angle_deg, delta_m0,
             angle_deg, delta_m1,
             angle_deg, delta_m2,
             ylog=0, xtitle="Incidence angle (to normal) [deg]",
             legend=["method 0",
                     "method 1",
                     "method 2"],
             marker=["+",None,None],
             linestyle=["",None,None],
             title="GLASS n=1.52 delta (phase s - phase p)", xrange=[0,90],
             show=1)

    if False: # Fresnel Rhomb  GLASS op cit  fig 1.16 (pag 450-51).
        from srxraylib.plot.gol import plot

        #
        # mirror reflectivity vs angle
        #


        angle_deg = numpy.linspace(0.0, 90.0, 1000)
        angle_mrad = 1e3 * numpy.radians(90.0 - angle_deg)

        refraction_index_1 = 1.51
        refraction_index_2 = 1.0

        RS0 = numpy.zeros_like(angle_deg)
        RS5 = numpy.zeros_like(angle_deg)

        n12 = refraction_index_2 / refraction_index_1

        rs0_m0, rp0_m0 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=0.0, method=0,
                                refraction_index_1=refraction_index_1, refraction_index_2=refraction_index_2)
        rs1_m0, rp1_m0 = PreRefl.reflectivity_amplitudes_fresnel_external(grazing_angle_mrad=angle_mrad, photon_energy_ev=20000,
                                roughness_rms_A=5.0, method=0,
                                refraction_index_1=refraction_index_1, refraction_index_2=refraction_index_2)


        plot(
             angle_deg, numpy.abs(rs0_m0) ** 2, angle_deg, numpy.abs(rp0_m0) ** 2,
             ylog=0, xtitle="Incidence angle (to normal) [deg]",
             marker=["+","+",None,None,None,None],
             title="FRESNEL RHOMB - REFLECTANCE", xrange=[0,90], yrange=[0,1],
             show=0)

        # degree of polarization (op.cit. pag 44 eq 42)
        P_m0 = numpy.abs((numpy.abs(rp0_m0) ** 2 - numpy.abs(rs0_m0) ** 2) / (numpy.abs(rp0_m0) ** 2 + numpy.abs(rs0_m0) ** 2))

        plot(
             angle_deg, P_m0,
             ylog=0, xtitle="Incidence angle (to normal) [deg]",
             title="FRESNEL RHOMB - POLARIZATION DEGREE", xrange=[0,90],
             show=0)

        # phase diff (op.cit. pag 50)
        delta_m0 = numpy.angle(rs0_m0) - numpy.angle(rp0_m0)

        angle = numpy.radians(angle_deg)
        arg_eq61 = numpy.cos(angle) * numpy.sqrt(numpy.sin(angle) ** 2 - n12 ** 2) / numpy.sin(angle) ** 2
        ii = numpy.nanargmax(arg_eq61)
        print("Max of eq 61 found at [deg]: ", angle_deg[ii], arg_eq61[ii])
        delta_eq61 = 2 * numpy.arctan(arg_eq61)


        plot(
             angle_deg, numpy.degrees(delta_m0),
             angle_deg, numpy.degrees(delta_eq61),
             angle_deg, angle_deg * 0 + 45,
             ylog=0, xtitle="Incidence angle (to normal) [deg]", ytitle="Phase diff delta [deg]",
             title="FRESNEL RHOMB - delta (phase s - phase p)", xrange=[0,90],
             legend=["fresnel", "eq61 TOTAL REFLECTION", "45"],
             show=1)


    if False:
        from srxraylib.plot.gol import plot

        #
        # mirror reflectivity
        #

        a = PreRefl()

        prerefl_file = "Rh5_50.dat"
        PreRefl.prerefl(interactive=False, SYMBOL="Be", DENSITY=1.848, FILE=prerefl_file,
                        E_MIN=5000.0, E_MAX=50000.0, E_STEP=100.0)

        a.read_preprocessor_file(prerefl_file)

        Energy = numpy.linspace(5000.0,40000.0,100)

        RS0 = numpy.zeros_like(Energy)
        RS5 = numpy.zeros_like(Energy)

        a.preprocessor_info()

        #
        # scalar inputs
        #
        process_phase = True
        method = 0

        for ii,ee in enumerate(Energy):
            if process_phase:
                ss, pp = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=3.0, photon_energy_ev=ee,
                                                           roughness_rms_A=0.0, method=method)
                RS0[ii] = numpy.abs(ss)**2
                ss, pp = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=3.0, photon_energy_ev=ee,
                                                           roughness_rms_A=5.0, method=method)
                RS5[ii] = numpy.abs(ss)**2
            else:
                rs, rp, rav = a.reflectivity_fresnel(grazing_angle_mrad=3.0,photon_energy_ev=ee,
                                                     roughness_rms_A=0.0, method=method)
                RS0[ii] = rs
                rs, rp, rav = a.reflectivity_fresnel(grazing_angle_mrad=3.0,photon_energy_ev=ee,
                                                     roughness_rms_A=5.0, method=method)
                RS5[ii] = rs

        plot(Energy,RS0,Energy,RS5,ylog=True,legend=["no roughness","5A RMS roughness"],title="Rh5_50.dat scalar inputs")



        # array inputs

        Grazing_angle = numpy.ones_like(Energy) * 3.0

        if process_phase:
            rs0, rp0 = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=Grazing_angle, photon_energy_ev=Energy,
                                                         roughness_rms_A=0.0)
            rs1, rp1 = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=Grazing_angle, photon_energy_ev=Energy,
                                                         roughness_rms_A=5.0)

            plot(Energy, numpy.abs(rs0)**2, Energy, numpy.abs(rs1)**2, ylog=True, legend=["no roughness", "5A RMS roughness"], title="array inputs")

        else:
            rs0, rp0, rav0 = a.reflectivity_fresnel(grazing_angle_mrad=Grazing_angle, photon_energy_ev=Energy, roughness_rms_A=0.0)
            rs1, rp1, rav1 = a.reflectivity_fresnel(grazing_angle_mrad=Grazing_angle, photon_energy_ev=Energy, roughness_rms_A=5.0)

            plot(Energy, rs0, Energy, rs1, ylog=True, legend=["no roughness", "5A RMS roughness"],title="Rh5_50.dat array inputs")


    if False: # compare with cxro i) Au at low energies (differences observed)

        prerefl_file = "reflec1.dat"
        prerefl_file_cxro = "reflec_cxro.dat"


        PreRefl.prerefl(interactive=False, SYMBOL="Au", DENSITY=19.3, FILE=prerefl_file,
                        E_MIN=110.0, E_MAX=501.0, E_STEP=1.0)

        PreRefl.prerefl_cxro(input_file="https://henke.lbl.gov/tmp/xray8383.dat", # https://henke.lbl.gov/tmp/xray8354.dat
                             output_file=prerefl_file_cxro,
                             E_MIN=110, E_MAX=501,  NPOINTS=500)

        Energy = numpy.linspace(110.0,500.0,1000)
        grazing_angle_mrad = 175.5 # 10.055409304545947 degrees

        a = PreRefl()
        a.read_preprocessor_file(filename=prerefl_file)
        rs, rp, rav = a.reflectivity_fresnel(grazing_angle_mrad=grazing_angle_mrad, photon_energy_ev=Energy,
                                             roughness_rms_A=0.0, method=2)

        a2 = PreRefl()
        a2.read_preprocessor_file(filename=prerefl_file_cxro)
        rs2, rp2, rav2 = a2.reflectivity_fresnel(grazing_angle_mrad=grazing_angle_mrad, photon_energy_ev=Energy,
                                             roughness_rms_A=0.0, method=2)

        plot(Energy, rs,
             Energy, rs2, legend=["xraylib","cxro"],
             xtitle="Photon energy [eV]", ytitle="Reflectivity", title="COMPARISON WITH CXRO - LOW ENERGY - Au@%g mrad" % grazing_angle_mrad)

    if False: # refractiveindex.info site

        prerefl_file_cxro = "reflec_refractiveindexinfo.dat"



        PreRefl.prerefl_refractiveindexinfo(shelf='main', book='SiO2', page='Franta',
                                output_file=prerefl_file_cxro,
                    WAVELENGTH_MIN=5000e-9, WAVELENGTH_MAX=15000e-9, NPOINTS=1000, density=2.2)
        #           WAVELENGTH_MIN = 0.025e-6, WAVELENGTH_MAX = 125e-6, NPOINTS = 100000, density = 2.2)

        a2 = PreRefl()
        a2.read_preprocessor_file(filename=prerefl_file_cxro)

        wavelength_nm = numpy.linspace(6000, 14000, 1500)
        energy_eV = codata.h * codata.c / codata.e / (wavelength_nm * 1e-9)

        print(">>>> eV: ", energy_eV)
        nn = a2.get_refraction_index(energy_eV)
        plot(wavelength_nm, nn.real,
             wavelength_nm, nn.imag,
             legend=["real(n)","imag(n)"],
             xtitle="Wavelength [nm]", title="data from database at: https://refractiveindex.info")



    if False: # compare with cxro ii) Si at high energies (similar)

        prerefl_file = "reflec1.dat"
        prerefl_file_cxro = "reflec_cxro.dat"


        PreRefl.prerefl(interactive=False, SYMBOL="Si", DENSITY=2.33, FILE=prerefl_file,
                        E_MIN=1000.0, E_MAX=25000.0, E_STEP=25.0)

        PreRefl.prerefl_cxro(input_file="https://henke.lbl.gov/tmp/xray8434.dat", # https://henke.lbl.gov/tmp/xray8354.dat
                             output_file=prerefl_file_cxro,
                             E_MIN=1000, E_MAX=25000,  NPOINTS=500)

        Energy = numpy.linspace(5000.0,20000.0,1000)
        grazing_angle_mrad = 5.0 # 10.055409304545947 degrees

        a = PreRefl()
        a.read_preprocessor_file(filename=prerefl_file)
        rs, rp, rav = a.reflectivity_fresnel(grazing_angle_mrad=grazing_angle_mrad, photon_energy_ev=Energy,
                                             roughness_rms_A=0.0, method=2)

        a2 = PreRefl()
        a2.read_preprocessor_file(filename=prerefl_file_cxro)
        rs2, rp2, rav2 = a2.reflectivity_fresnel(grazing_angle_mrad=grazing_angle_mrad, photon_energy_ev=Energy,
                                             roughness_rms_A=0.0, method=2)

        plot(Energy, rs,
             Energy, rs2, legend=["xraylib","cxro"],
             xtitle="Photon energy [eV]", ytitle="Reflectivity", title="COMPARISON WITH CXRO - HIGH ENERGY - Si@%g mrad" % grazing_angle_mrad)

    if False: # comparing with f1f2_calc

        # """

        prerefl_file = "reflec1.dat"

        PreRefl.prerefl(interactive=False, SYMBOL="Au", DENSITY=19.3, FILE=prerefl_file,
                        E_MIN=110.0, E_MAX=501.0, E_STEP=1.0)
        #
        Energy = numpy.linspace(110.0, 500.0, 1000)
        grazing_angle_mrad = 175.5
        #
        # PreRefl.prerefl(interactive=False, SYMBOL="Au", DENSITY=19.3, FILE=prerefl_file,
        #                 E_MIN=110.0, E_MAX=501.0, E_STEP=1.0)

        a = PreRefl()
        a.read_preprocessor_file(prerefl_file)
        a.get_refraction_index(110.0)

        RS0 = numpy.zeros_like(Energy)

        #
        # scalar inputs
        #
        process_phase = True
        method = 0

        for ii,ee in enumerate(Energy):
            if process_phase:
                ss, pp = a.reflectivity_amplitudes_fresnel(grazing_angle_mrad=grazing_angle_mrad, photon_energy_ev=ee,
                                                           roughness_rms_A=0.0, method=method)
                RS0[ii] = numpy.abs(ss)**2
            else:
                rs, rp, rav = a.reflectivity_fresnel(grazing_angle_mrad=grazing_angle_mrad,photon_energy_ev=ee,
                                                     roughness_rms_A=0.0, method=method)
                RS0[ii] = rs


        try:
            from xoppylib.scattering_functions.f1f2_calc import f1f2_calc
            try: import xraylib
            except: "xraylib not available"
        except:
            print("Please install xoppylib (to be avoided)")

        aa = f1f2_calc("Au", Energy, theta=grazing_angle_mrad*1e-3, F=8, density=None, rough=0.0,
                       material_constants_library=xraylib)
        print(aa.shape)
        plot(Energy,RS0,
             Energy,aa,
             legend=["xraylib",'f1f2_calc'],
             xtitle="Photon energy [eV]", ytitle="Reflectivity",
             title="COMPARISON WITH f1f2_calc - LOW ENERGY - Si@%g mrad" % grazing_angle_mrad)


    if False:
        prerefl_file = "reflec1.dat"

        PreRefl.prerefl(interactive=False, SYMBOL="Au", DENSITY=19.3, FILE=prerefl_file,
                        E_MIN=1000.0, E_MAX=5000.0, E_STEP=100.0)


        a = PreRefl()
        a.read_preprocessor_file(prerefl_file)


        energies = (2000, 3000, 4000)
        print("Energies: ", energies)
        method = 0

        for energy in energies:
            print("+++++++++++++++++++++++++++++ Energy: ", energy)
            print(
                a.reflectivity_amplitudes_fresnel(photon_energy_ev=energy,
                                                         grazing_angle_mrad=3.0,
                                                         roughness_rms_A=0.0,
                                                         method=method # 0=born & wolf, 1=parratt, 2=shadow3
                                                        )
            )

        energies = numpy.array((2000,3000,4000,4000))
        print("+++++++++++++++++++++++++++++ Energy Scan: ", energies)
        print(
                        a.reflectivity_amplitudes_fresnel(photon_energy_ev=energies,
                                                                 grazing_angle_mrad=numpy.array((3.0,3.0,3.0,1.0)),
                                                                 roughness_rms_A=0.0,
                                                                 method=method # 0=born & wolf, 1=parratt, 2=shadow3
                                                                )
                    )

        tmp_xrl = PreRefl.reflectivity_amplitudes_fresnel_external_xraylib(
                                          photon_energy_ev=energies,
                                          coating_material="Au",
                                          coating_density=19.3,
                                          grazing_angle_mrad=numpy.array((3.0, 3.0, 3.0, 1.0)),
                                          roughness_rms_A=0.0,
                                          method=method  # 0=born & wolf, 1=parratt, 2=shadow3
                                          )

        tmp_dx = PreRefl.reflectivity_amplitudes_fresnel_external_dabax(
                                          photon_energy_ev=energies,
                                          coating_material="Au",
                                          coating_density=19.3,
                                          grazing_angle_mrad=numpy.array((3.0, 3.0, 3.0, 1.0)),
                                          roughness_rms_A=0.0,
                                          method=method  # 0=born & wolf, 1=parratt, 2=shadow3
                                          )

        print(">>>> tmp_xrl", tmp_xrl)
        print(">>>> tmp_dx", tmp_dx)

    if False :
        prerefl_file = "reflec1.dat"

        PreRefl.prerefl(interactive=False, SYMBOL="Au", DENSITY=19.3, FILE=prerefl_file,
                        E_MIN=1000.0, E_MAX=5000.0, E_STEP=100.0)


        a = PreRefl()
        a.read_preprocessor_file(prerefl_file)


        energies = (2000,3000,4000)
        print("Energies: ", energies)
        for energy in energies:
            print(
                a.get_attenuation_coefficient(energy))

        tmp_xrl = PreRefl.get_attenuation_coefficient_external_xraylib(
                                          photon_energy_ev=numpy.array((2000, 3000, 4000, 4000)),
                                          material="Au",
                                          density=19.3,
                                            )

        tmp_dx = PreRefl.get_attenuation_coefficient_external_dabax(
                                          photon_energy_ev=numpy.array((2000, 3000, 4000, 4000)),
                                          material="Au",
                                          density=19.3,
                                            )

        print(">>>> tmp_xrl", tmp_xrl)
        print(">>>> tmp_dx", tmp_dx)


        if True:
            print("Refraction index Au (prerefl file): ", a.get_refraction_index(numpy.array(energies)))
            n_xrl = PreRefl.get_refraction_index_external_xraylib(
                                                photon_energy_ev=numpy.array((2000, 3000, 4000, 4000)),
                                                material="Au",
                                                density=19.3,
                                                       )

            n_dx = PreRefl.get_refraction_index_external_dabax(
                                                photon_energy_ev=numpy.array((2000, 3000, 4000, 4000)),
                                                material="Au",
                                                density=19.3,
                                                dabax=None,
                                                       )
            print("Refraction index Au (xraylib): ", n_xrl)
            print("Refraction index Au (dabax): ", n_dx)

        if True:
            prerefl_file = "reflec1.dat"
            PreRefl.prerefl(interactive=False, SYMBOL="Al", DENSITY=2.6989, FILE=prerefl_file,
                            E_MIN=13000.0, E_MAX=15000.0, E_STEP=100.0)

            a = PreRefl()
            a.read_preprocessor_file(prerefl_file)

            energies = (14000, 14001, 14002)
            n_xrl = PreRefl.get_refraction_index_external_xraylib(
                                                photon_energy_ev=numpy.array(energies),
                                                material="Al",
                                                density=2.6989,
                                                       )

            n_dx = PreRefl.get_refraction_index_external_dabax(
                                                photon_energy_ev=numpy.array(energies),
                                                material="Al",
                                                density=2.6989,
                                                dabax=None,
                                                       )

            a_xrl = PreRefl.get_attenuation_coefficient_external_xraylib(
                                                photon_energy_ev=numpy.array(energies),
                                                material="Al",
                                                density=2.6989,
                                                       )

            a_dx = PreRefl.get_attenuation_coefficient_external_dabax(
                                                photon_energy_ev=numpy.array(energies),
                                                material="Al",
                                                density=2.6989,
                                                dabax=None,
                                                       )
            print("Attenuation Al (prerefl file): ", a.get_attenuation_coefficient(numpy.array(energies)))
            print("Attenuation Al (xraylib): ", a_xrl)
            print("Attenuation Al (dabax): ", a_dx)