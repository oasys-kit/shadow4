"""
python version of the mirror reflectivity code and refractive index calculations in shadow.
"""
# TODO: vectorization
import numpy
import scipy.constants as codata
from shadow4.tools.logger import is_verbose, is_debug

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

    def get_refraction_index(self,energy1):
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

        wnum = 2*numpy.pi * energy1 / tocm

        QSTEP = self.prerefl_dict["QSTEP"]
        QMIN = self.prerefl_dict["QMIN"]

        index1 = (wnum - QMIN) / QSTEP
        index1 = numpy.array(index1).astype(int)

        if index1.max() > self.prerefl_dict["NREFL"]:
            raise Exception("Error: Photon energy above tabulated upper limit.")

        if index1.min() < 0:
            raise Exception("Error: Photon energy below tabulated lower limit.")

        # WNUM0 = QSTEP * int(index1) + QMIN
        WNUM0 = QSTEP * index1 + QMIN
        DEL_X = wnum - WNUM0
        DEL_X = DEL_X / QSTEP

        # index1 = int(index1)

        ALFA = self.prerefl_dict["ZF1"][index1] + (self.prerefl_dict["ZF1"][index1 + 1] - self.prerefl_dict["ZF1"][index1]) * DEL_X
        GAMMA = self.prerefl_dict["ZF2"][index1] + (self.prerefl_dict["ZF2"][index1 + 1] - self.prerefl_dict["ZF2"][index1]) * DEL_X

        refraction_index = (1.0 - ALFA / 2) + (GAMMA / 2)*1j

        if is_verbose():
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

        import xraylib
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
        import xraylib
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
                                                     dabax=None,
                                                     ):

        from dabax.dabax_xraylib import DabaxXraylib
        if isinstance(dabax, DabaxXraylib):
            dx = dabax
        else:
            dx = DabaxXraylib()

        return dx.Refractive_Index_Re(material, photon_energy_ev * 1e-3, density) + \
               1j * dx.Refractive_Index_Im(material, photon_energy_ev * 1e-3, density)

    @classmethod
    def get_refraction_index_real_external_dabax(self,
                                                     photon_energy_ev=10000.0,
                                                     material="SiC",
                                                     density=3.217,
                                                     dabax=None,
                                                     ):
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
        if isinstance(dabax, DabaxXraylib):
            dx = dabax
        else:
            dx = DabaxXraylib()

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
                                        method=2 # 0=born & wolf, 1=parratt, 2=shadow3 (avoid using 0 or 1, experimental!!)
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
                                                         method=2,
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
        import xraylib
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
                                                         method=2, # 0=born & wolf, 1=parratt, 2=shadow3 (avoid using 0 or 1, experimental!!) ):
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
                                                 method=2, # 0=born & wolf, 1=parratt, 2=shadow3 (avoid using 0 or 1, experimental!!) ):
                                                 ):
        """
        Standalone method to calculate the reflectivity (amplitude) of an interface using Fresnel formulas using optical
        constants entered by user.

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
        refraction_index_1 : float or numpy array
            the refraction index (complex) for medium 1 (object).
        refraction_index_2 : float or numpy array
            the refraction index (complex) for medium 2 (image).

        Returns
        -------
        tuple
            (rs, rp) the s-polarized and p-pol and amplitude reflectivities
        """

        theta1g = grazing_angle_mrad * 1e-3     # in rad
        theta1 = numpy.pi / 2 - theta1g
        rough1 = roughness_rms_A

        sin_theta1 = numpy.sin(theta1)
        cos_theta1 = numpy.sqrt(1 - sin_theta1**2, dtype=complex)

        sin_theta2 = refraction_index_1 * sin_theta1 / refraction_index_2
        cos_theta2 = numpy.sqrt(1 - sin_theta2**2, dtype=complex)


        if method == 0: # born & wolf
            rs1 = refraction_index_2 * cos_theta1 - refraction_index_1 * cos_theta2
            rs2 = refraction_index_2 * cos_theta1 + refraction_index_1 * cos_theta2
            rs = rs1 / rs2

            rp1 = refraction_index_1 * cos_theta1 - refraction_index_2 * cos_theta2
            rp2 = refraction_index_1 * cos_theta1 + refraction_index_2 * cos_theta2
            rp = rp1 / rp2

        elif method == 1: # parratt
            phi = theta1g
            delta1 = 1.0 - numpy.real(refraction_index_1)
            beta1 = numpy.imag(refraction_index_1)
            delta2 = 1.0 - numpy.real(refraction_index_2)
            beta2 = numpy.imag(refraction_index_2)
            f1 = numpy.sqrt(phi ** 2 - 2 * delta1 - 2 * 1j * beta1, dtype=complex)
            f2 = numpy.sqrt(phi ** 2 - 2 * delta2 - 2 * 1j * beta2, dtype=complex)
            rs = (f1 - f2) / (f1 + f2)
            rp = rs

        elif method == 2:
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

        wavelength_m = codata.h * codata.c / codata.e / photon_energy_ev

        debyewaller = numpy.exp( -(4.0 * numpy.pi * numpy.sin(theta1g) * rough1 / (wavelength_m * 1e10))**2 )

        # runp = 0.5 * (rs + rp)
        # return numpy.abs(rs)**2 * debyewaller, numpy.abs(rp)**2 * debyewaller, numpy.abs(runp)**2 * debyewaller

        return rs * numpy.sqrt(debyewaller), rp * numpy.sqrt(debyewaller)


    #
    # this is copied (ans sligtly ceaned) from shadow3 python preprocessors
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
        import xraylib

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
            tmp = 2e0 * (1e0 - xraylib.Refractive_Index_Re(iMaterial, energy, density))
            f.write("%e \n" % tmp)
        for i in range(npoint):
            energy = (estart + estep * i) * 1e-3
            tmp2 = 2e0 * (xraylib.Refractive_Index_Im(iMaterial, energy, density))
            f.write("%e \n" % tmp2)
        print("File written to disk: %s" % out_file)
        f.close()

        # test (not needed)
        itest = 0
        if itest:
            cdtest = xraylib.CompoundParser(iMaterial)
            print("    ", iMaterial, " contains %i atoms and %i elements" % (cdtest['nAtomsAll'], cdtest['nElements']))
            for i in range(cdtest['nElements']):
                print("    Element %i: %lf %%" % (cdtest['Elements'][i], cdtest['massFractions'][i] * 100.0))
            print(qmin, qmax, qstep, depth0)
            print(npoint)
            for i in range(npoint):
                energy = (estart + estep * i) * 1e-3
                qq = qmin + qstep * i
                print(energy, qq, \
                      2e0 * (1e0 - xraylib.Refractive_Index_Re(iMaterial, energy, density)), \
                      2e0 * (xraylib.Refractive_Index_Im(iMaterial, energy, density)))

        return None

    #
    # this is copied (ans sligtly ceaned) from shadow3 python preprocessors
    #
    @classmethod
    def prerefl_cxro(cls, input_file="https://henke.lbl.gov/tmp/xray8378.dat", output_file="prerefl.dat",
                     E_MIN=None, E_MAX=None, NPOINTS=1000):
        """
        Creates an instance of PreRefl with parameters initialized from the keyword parameters and the
        prerefl preprocessor file. It uses the CXRO web site for accessing the optical constants.

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

if __name__ == "__main__":
    from srxraylib.plot.gol import plot
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
        process_phase = False
        method = 2

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

        plot(Energy,RS0,Energy,RS5,ylog=True,legend=["no roughness","5A RMS roughness"],title="scalar inputs")



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

            plot(Energy, rs0, Energy, rs1, ylog=True, legend=["no roughness", "5A RMS roughness"],title="array inputs")


    if False: # compare with cxro

        prerefl_file = "reflec1.dat"
        prerefl_file_cxro = "reflec_cxro.dat"


        PreRefl.prerefl(interactive=False, SYMBOL="Au", DENSITY=19.3, FILE=prerefl_file,
                        E_MIN=110.0, E_MAX=501.0, E_STEP=1.0)

        PreRefl.prerefl_cxro(input_file="https://henke.lbl.gov/tmp/xray8378.dat",
                             output_file=prerefl_file_cxro,
                             E_MIN=110, E_MAX=501,  NPOINTS=500)

        Energy = numpy.linspace(110.0,500.0,1000)
        grazing_angle_mrad = 175.5

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
             xtitle="Photon energy [eV]", ytitle="Reflectivity", title="Au@%g mrad" % grazing_angle_mrad)

    if False: # comparing with f1f2_calc

        # """

        prerefl_file = "reflec1.dat"

        PreRefl.prerefl(interactive=False, SYMBOL="Au", DENSITY=19.3, FILE=prerefl_file,
                        E_MIN=110.0, E_MAX=501.0, E_STEP=1.0)

        Energy = numpy.linspace(110.0, 500.0, 1000)
        grazing_angle_mrad = 175.5

        PreRefl.prerefl(interactive=False, SYMBOL="Au", DENSITY=19.3, FILE=prerefl_file,
                        E_MIN=110.0, E_MAX=501.0, E_STEP=1.0)

        a = PreRefl()
        a.read_preprocessor_file(prerefl_file)
        a.get_refraction_index(110.0)

        RS0 = numpy.zeros_like(Energy)

        #
        # scalar inputs
        #
        process_phase = False
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
            import xraylib
        except:
            print("Please install xoppylib (to be avoided)")

        aa = f1f2_calc("Au", Energy, theta=grazing_angle_mrad*1e-3, F=8, density=None, rough=0.0,
                       material_constants_library=xraylib)
        print(aa.shape)
        plot(Energy,RS0,
             Energy,aa,
             legend=["RS0",'f1f2_calc'])


    if False:
        prerefl_file = "reflec1.dat"

        PreRefl.prerefl(interactive=False, SYMBOL="Au", DENSITY=19.3, FILE=prerefl_file,
                        E_MIN=1000.0, E_MAX=5000.0, E_STEP=100.0)


        a = PreRefl()
        a.read_preprocessor_file(prerefl_file)


        energies = (2000,3000,4000)
        print("Energies: ", energies)
        for energy in energies:
            print(
                a.reflectivity_amplitudes_fresnel(photon_energy_ev=energy,
                                                         grazing_angle_mrad=3.0,
                                                         roughness_rms_A=0.0,
                                                         method=2 # 0=born & wolf, 1=parratt, 2=shadow3
                                                        )
            )

        print(
                        a.reflectivity_amplitudes_fresnel(photon_energy_ev=numpy.array((2000,3000,4000,4000)),
                                                                 grazing_angle_mrad=numpy.array((3.0,3.0,3.0,1.0)),
                                                                 roughness_rms_A=0.0,
                                                                 method=2 # 0=born & wolf, 1=parratt, 2=shadow3
                                                                )
                    )

        tmp_xrl = PreRefl.reflectivity_amplitudes_fresnel_external_xraylib(
                                          photon_energy_ev=numpy.array((2000, 3000, 4000, 4000)),
                                          coating_material="Au",
                                          coating_density=19.3,
                                          grazing_angle_mrad=numpy.array((3.0, 3.0, 3.0, 1.0)),
                                          roughness_rms_A=0.0,
                                          method=2  # 0=born & wolf, 1=parratt, 2=shadow3
                                          )

        tmp_dx = PreRefl.reflectivity_amplitudes_fresnel_external_dabax(
                                          photon_energy_ev=numpy.array((2000, 3000, 4000, 4000)),
                                          coating_material="Au",
                                          coating_density=19.3,
                                          grazing_angle_mrad=numpy.array((3.0, 3.0, 3.0, 1.0)),
                                          roughness_rms_A=0.0,
                                          method=2  # 0=born & wolf, 1=parratt, 2=shadow3
                                          )

        print(">>>> tmp_xrl", tmp_xrl)
        print(">>>> tmp_dx", tmp_dx)

    if True:
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


        if False:
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
            # print("Refraction index Al (prerefl file): ", a.get_refraction_index(numpy.array(energies)))
            # print("Refraction index Al (xraylib): ", n_xrl)
            # print("Refraction index Al (dabax): ", n_dx)

            print("Attenuation index Al (prerefl file): ", a.get_attenuation_coefficient(numpy.array(energies)))
            print("Attenuation index Al (xraylib): ", a_xrl)
            print("Attenuation index Al (dabax): ", a_dx)