# todo: DELETE? this class is not used... use crystalpy instead...
#
# Bragg() a class to manipulate the shadow bragg preprocessor file v1
#
#  uses xoppylib that contains bragg preprocessor readers
#
import numpy
import scipy.constants as codata
from xoppylib.crystals.bragg_preprocessor_file_io import bragg_preprocessor_file_v1_read

class Bragg(object):
    def __init__(self, preprocessor_file=None, preprocessor_dictionary=None):
        self._preprocessor_file = preprocessor_file
        self._preprocessor_dictionary = preprocessor_dictionary

    #
    # setters and getters
    #
    def set_preprocessor_file(self, filename):
        self._preprocessor_file = filename

    def get_preprocessor_file(self):
        return self._preprocessor_file

    def get_preprocessor_dictionary(self):
        return self._preprocessor_dictionary

    #
    # loader
    #
    @classmethod
    def create_from_preprocessor_file(cls, filename):
        out = Bragg(preprocessor_file=filename)
        out.load_preprocessor_file()
        return out

    def load_preprocessor_file(self):
        out_dict = bragg_preprocessor_file_v1_read(self._preprocessor_file)
        self._preprocessor_dictionary = out_dict

    #
    # extract values
    #


    def __F_elastic(self, ratio):
        CA = self._preprocessor_dictionary["fit_a"]
        CB = self._preprocessor_dictionary["fit_b"]
        FOA = CA[2] * ratio ** 2 + CA[1] * ratio + CA[0]
        FOB = CB[2] * ratio ** 2 + CB[1] * ratio + CB[0]
        return FOA, FOB

    def structure_factor(self, energy, ratio=None):

        if ratio is None:
            ratio = self.ratio(energy, theta=None)

        F1A, F2A, F1B, F2B = self.__interpolate(energy)
        FOA, FOB = self.__F_elastic(ratio)
        FA = FOA + F1A + 1j * F2A
        FB = FOB + F1B + 1j * F2B

        I_LATT = self._preprocessor_dictionary["i_latt"]

        GA = self._preprocessor_dictionary["ga.real"] + 1j * self._preprocessor_dictionary["ga.imag"]
        GB = self._preprocessor_dictionary["gb.real"] + 1j * self._preprocessor_dictionary["gb.imag"]
        GA_BAR = self._preprocessor_dictionary["ga_bar.real"] + 1j * self._preprocessor_dictionary["ga_bar.imag"]
        GB_BAR = self._preprocessor_dictionary["gb_bar.real"] + 1j * self._preprocessor_dictionary["gb_bar.imag"]
        ATNUM_A = self._preprocessor_dictionary["zeta_a"]
        ATNUM_B = self._preprocessor_dictionary["zeta_b"]
        TEMPER = self._preprocessor_dictionary["temper"]


        # RN = self._preprocessor_dictionary["one_over_volume_times_electron_radius_in_cm"]

        if (I_LATT == 0):
            # ABSORP = 2.0 * RN * R_LAM0 * (4.0*(DIMAG(FA)+DIMAG(FB)))
            F_0 = 4*((F1A + ATNUM_A + F1B + ATNUM_B) + 1j*(F2A + F2B))
        elif (I_LATT == 1):
            # ABSORP = 2.0D0*RN*R_LAM0*(4.0D0*(DIMAG(FA)+DIMAG(FB)))
            F_0 = 4*((F1A + ATNUM_A + F1B + ATNUM_B) + 1j*(F2A + F2B))
        elif (I_LATT == 2):
            FB	 = 0.0 + 0.0j
            # ABSORP = 2.0D0*RN*R_LAM0*(4.0D0*DIMAG(FA))
            F_0 = 4*(F1A + ATNUM_A + 1j*F2A)
        elif (I_LATT == 3):
            # ABSORP = 2.0D0*RN*R_LAM0*(DIMAG(FA)+DIMAG(FB))
            F_0 = (F1A + ATNUM_A + F1B + ATNUM_B) + CI*(F2A + F2B)
        elif (I_LATT == 4):
            FB = 0.0 + 0.0j
            # ABSORP = 2.0D0*RN*R_LAM0*(2.0D0*(DIMAG(FA)))
            F_0 = 2*(F1A+ 1j*F2A)
        elif (I_LATT == 5):
            FB = 0.0 + 0.0j
            # ABSORP = 2.0D0*RN*R_LAM0*(4.0D0*(DIMAG(FA)))
            F_0 = 4*(F1A + 1j*F2A )

        # ! C
        # ! C FH and FH_BAR are the structure factors for (h,k,l) and (-h,-k,-l).
        # ! C
        # ! C srio, Added TEMPER here (95/01/19)
        FH 	= ( (GA * FA) + (GB * FB) ) * TEMPER
        FH_BAR	= ( (GA_BAR * FA) + (GB_BAR * FB) ) * TEMPER
        STRUCT = numpy.sqrt(FH * FH_BAR)
        return F_0, FH, FH_BAR, STRUCT, FA, FB

    def info(self, energy):

        F_0, FH, FH_BAR, STRUCT, FA, FB = self.structure_factor(energy)
        ratio = self.ratio(energy, theta=None)
        darwin_s, darwin_p = self.darwin_halfwidth(energy)

        info_text = ""
        if isinstance(energy, float) or isinstance(energy, int):
            info_text += "\nPhoton energy: %g eV" % energy
            info_text += "\nsin(theta)/wavelength = %g A**-1" % ratio
            info_text += "\n    For atom A, fo + f' + if'' = (%g + %g i) " %  (FA.real, FA.imag)
            info_text += "\n    For atom B, fo + f' + if'' = (%g + %g i) " %  (FB.real, FB.imag)
            info_text += "\n    Structure factor F0: (%g + %g i)" %  (F_0.real, F_0.imag)
            info_text += "\n    Structure factor FH: (%g + %g i)" %  (FH.real, FH.imag)
            info_text += "\n    Structure factor FH_BAR (%g + %g i): " % (FH_BAR.real, FH_BAR.imag)
            info_text += "\n    Structure factor STRUCT: (%g + %g i): " % (STRUCT.real, STRUCT.imag)
            info_text += "\n    dSpacing = %g m" % self.d_spacing()
            info_text += "\n    unit cell volume = %g m^3" % (self.unit_cell_volume())
            info_text += "\n    wavelength = %g m" % self.wavelength(energy)
            info_text += "\n    bragg angle = %g rad = %g deg " % (self.bragg_angle(energy), self.bragg_angle(energy) * 180 / numpy.pi)
            info_text += "\n    darwin_halfwidth_s = %f urad " % (1e6 * darwin_s)
            info_text += "\n    darwin_halfwidth_p = %f urad " % (1e6 * darwin_p)
        else:
            for i, ener in enumerate(energy):
                info_text += "\nPhoton energy: %g eV" % ener
                info_text += "\nsin(theta)/wavelength = %g A**-1" % ratio[i]
                info_text += "\n    For atom A, fo + f' + if'' = (%g + i %g) " %  (FA[i].real, FA[i].imag)
                info_text += "\n    For atom B, fo + f' + if'' = (%g + i %g) " %  (FB[i].real, FB[i].imag)
                info_text += "\n    Structure factor F0: (%g + %g i)" %  (F_0[i].real, F_0[i].imag)
                info_text += "\n    Structure factor FH: (%g + %g i)" %  (FH[i].real, FH[i].imag)
                info_text += "\n    Structure factor FH_BAR (%g + %g i): " % (FH_BAR[i].real, FH_BAR[i].imag)
                info_text += "\n    Structure factor STRUCT: (%g + %g i): " % (STRUCT[i].real, STRUCT[i].imag)
                info_text += "\n    dSpacing = %g m" % self.d_spacing()
                info_text += "\n    unit cell volume = %g m^3" % (self.unit_cell_volume())
                info_text += "\n    wavelength = %g m" % self.wavelength(ener)
                info_text += "\n    bragg angle = %g rad = %g deg " % (self.bragg_angle(ener), self.bragg_angle(ener) * 180 / numpy.pi)
                info_text += "\n    darwin_halfwidth_s = %f urad " % (1e6 * darwin_s[i])
                info_text += "\n    darwin_halfwidth_p = %f urad " % (1e6 * darwin_p[i])
        return info_text

    def wavelength(self, energy):
        return codata.h * codata.c / codata.e / energy

    def wavelength_in_A(self, energy):
        return 1e10 * self.wavelength(energy)

    def ratio(self, energy, theta=None):
        """
        sin(theta) / lambda in A**-1

        Parameters
        ----------
        energy

        Returns
        -------

        """

        if theta is None:
            theta = self.bragg_angle(energy)

        return numpy.sin(theta) / ( self.wavelength_in_A(energy))

    def darwin_halfwidth_s(self, energy):
        return self.darwin_halfwidth(energy)[0]

    def darwin_halfwidth_p(self, energy):
        return self.darwin_halfwidth(energy)[1]

    def darwin_halfwidth(self, energy):
        if isinstance(energy, int): energy = float(energy)
        RN = self.get_preprocessor_dictionary()["one_over_volume_times_electron_radius_in_cm"]
        R_LAM0 = self.wavelength(energy) * 1e2
        F_0, FH, FH_BAR, STRUCT, FA, FB = self.structure_factor(energy, ratio=None)
        STRUCT = numpy.sqrt( FH * FH_BAR)
        TEMPER = self.get_preprocessor_dictionary()["temper"]
        GRAZE = self.bragg_angle(energy)
        SSVAR	= RN*(R_LAM0**2)*STRUCT*TEMPER/numpy.pi/numpy.sin(2.0*GRAZE)
        SPVAR = SSVAR * numpy.abs(numpy.cos(2.0 * GRAZE))
        return SSVAR.real, SPVAR.real

    def __energy_index(self, energy1):
        Energy = self._preprocessor_dictionary["Energy"]
        energy = numpy.array(energy1)
        if energy.size == 1:
            if (energy < Energy.min()) or (energy > Energy.max()):
                return -100
            ll = numpy.where(Energy > energy)[0][0]
        else:
            ll = numpy.zeros(energy.size, dtype=int)
            for i, ener in enumerate(energy):
                if (ener < Energy.min()) or (ener > Energy.max()):
                    ll[i] = -100
                else:
                    ll[i] = numpy.where( Energy > ener)[0][0]

        NENER = numpy.array(ll)

        if (NENER < 0).any():
            raise Exception("Cannot interpolate: energy outside limits")

        return NENER

    def __interpolate(self, PHOT):
        ENERGY = self._preprocessor_dictionary["Energy"]
        NENER = self.__energy_index(PHOT) - 1
        FP_A  = self._preprocessor_dictionary["F1a"]
        FP_B  = self._preprocessor_dictionary["F1b"]
        FPP_A = self._preprocessor_dictionary["F2a"]
        FPP_B = self._preprocessor_dictionary["F2b"]
        F1A	=  FP_A[NENER] +  (FP_A[NENER+1] -  FP_A[NENER]) *  (PHOT - ENERGY[NENER]) / (ENERGY[NENER+1] - ENERGY[NENER])
        F2A	= FPP_A[NENER] + (FPP_A[NENER+1] - FPP_A[NENER]) *  (PHOT - ENERGY[NENER]) / (ENERGY[NENER+1] - ENERGY[NENER])
        F1B	=  FP_B[NENER] +  (FP_B[NENER+1] -  FP_B[NENER]) *  (PHOT - ENERGY[NENER]) / (ENERGY[NENER+1] - ENERGY[NENER])
        F2B	= FPP_B[NENER] + (FPP_B[NENER+1] - FPP_B[NENER]) *  (PHOT - ENERGY[NENER]) / (ENERGY[NENER+1] - ENERGY[NENER])
        return F1A, F2A, F1B, F2B


    def d_spacing(self):
        """
        Returns the lattice spacing d in m
        :return: Lattice spacing.
        """

        return self._preprocessor_dictionary["dspacing_in_cm"] * 1e-2

    def bragg_angle(self, energy):
        """
        Returns the Bragg angle for a given energy.
        :param energy: Energy to calculate the Bragg angle for.
        :return: Bragg angle.
        """

        return numpy.arcsin(self.wavelength(energy) / 2 / self.d_spacing())

    def unit_cell_volume(self):
        """
        Returns the unit cell volume.

        :return: Unit cell volume
        """
        # Retrieve unit cell volume from xraylib.
        one_over_volume_times_electron_radius_in_cm = self.get_preprocessor_dictionary()["one_over_volume_times_electron_radius_in_cm"]
        # codata_e2_mc2 = 2.81794032e-15 = Classical electron radius in S.I.
        codata_e2_mc2 = codata.hbar * codata.alpha / codata.m_e / codata.c
        one_over_volume_in_cm = one_over_volume_times_electron_radius_in_cm / (codata_e2_mc2 * 1e2)
        unit_cell_volume = 1.0 / one_over_volume_in_cm * (1e-2)**3
        return unit_cell_volume


if __name__ == "__main__":
    import numpy
    import xraylib


    from xoppylib.crystals.create_bragg_preprocessor_file_v1 import create_bragg_preprocessor_file_v1

    SHADOW_FILE = "bragg_v1_xraylib.dat"
    tmp = create_bragg_preprocessor_file_v1(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1,
                                                           K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                                                           TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0,
                                                           E_STEP=100.0, SHADOW_FILE=SHADOW_FILE,
                                                           material_constants_library=xraylib)

    b = Bragg.create_from_preprocessor_file(SHADOW_FILE)
    tmp = b.get_preprocessor_dictionary()

    for key in tmp.keys():
        print("---------------", key)

    ener = numpy.linspace(8000.0,9000,2)
    ratio = 0.1595
    F_0, FH, FH_BAR, STRUCT, FA, FB = b.structure_factor(ener, ratio=None)
    info_text = b.info(ener)
    print(info_text)
