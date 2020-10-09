
import numpy
import xraylib
import scipy.constants as codata

class Bragg(object):
    def __init__(self, preprocessor_file=None, preprocessor_dictionary=None):
        self._preprocessor_file = preprocessor_file
        self._preprocessor_dictionary = preprocessor_dictionary

    def set_preprocessor_file(self, filename):
        self._preprocessor_file = filename

    def get_preprocessor_file(self):
        return self._preprocessor_file

    def get_preprocessor_dictionary(self):
        return self._preprocessor_dictionary

    def run_preprocessor(self, interactive=True, DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
              TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0, SHADOW_FILE="bragg.dat"):

        out_dict = Bragg.bragg(interactive=interactive, DESCRIPTOR=DESCRIPTOR,
                               H_MILLER_INDEX=H_MILLER_INDEX, K_MILLER_INDEX=K_MILLER_INDEX, L_MILLER_INDEX=L_MILLER_INDEX,
                               TEMPERATURE_FACTOR=TEMPERATURE_FACTOR, E_MIN=E_MIN, E_MAX=E_MAX, E_STEP=E_STEP,
                               SHADOW_FILE=SHADOW_FILE)

        return Bragg.__init__(self, preprocessor_file=SHADOW_FILE, preprocessor_dictionary=out_dict)

    def __parse_line(self, line, remove=[]):
        if len(remove) > 0:
            for str1 in remove:
                line = line.replace(str1, " ")
        line = " ".join(line.split())
        variables = line.split(" ")
        return variables

    def load_preprocessor_file(self, filename=None):
        if filename is None:
            filename = self._preprocessor_file
        else:
            self._preprocessor_file = filename

        f = open(filename, 'r')
        lines = f.read().splitlines()
        f.close()

        out_dict = {}

        line_index = 0
        variables = self.__parse_line(lines[line_index])
        # print(">>>>>>>>>> variables", variables)
        out_dict["i_latt"] = int(variables[0])
        out_dict["one_over_volume_times_electron_radius_in_cm"] = float(variables[1])
        out_dict["dspacing_in_cm"] = float(variables[2])

        line_index += 1
        variables = self.__parse_line(lines[line_index])
        out_dict["zeta_a"] = int(variables[0])
        out_dict["zeta_b"] = int(variables[1])
        out_dict["temper"] = float(variables[2])

        # line_index += 1

        line_index += 1
        variables = self.__parse_line(lines[line_index], remove=["(",")",","])
        # print(">>>>>>>>>> variables", variables)
        out_dict["ga.real"] = float(variables[0])
        out_dict["ga.imag"] = float(variables[1])

        line_index += 1
        # print(">>>>>>>>>> variables", variables)
        variables = self.__parse_line(lines[line_index], remove=["(", ")", ","])
        out_dict["ga_bar.real"] = float(variables[0])
        out_dict["ga_bar.imag"] = float(variables[1])

        line_index += 1
        # print(">>>>>>>>>> variables", variables)
        variables = self.__parse_line(lines[line_index], remove=["(", ")", ","])
        out_dict["gb.real"] = float(variables[0])
        out_dict["gb.imag"] = float(variables[1])

        line_index += 1
        line = lines[line_index]
        line = line.replace("(", "")
        line = line.replace(")", "")
        line = line.replace(" ", "")
        variables = line.split(",")
        # print(">>>>>>>>>> variables", variables)
        out_dict["gb_bar.real"] = float(variables[0])
        out_dict["gb_bar.imag"] = float(variables[1])

        line_index += 1
        variables = self.__parse_line(lines[line_index])
        # print(">>>>>>>>>> variables", variables)
        out_dict["fit_a"] = []
        for variable in variables:
            out_dict["fit_a"].append(float(variable))

        line_index += 1
        variables = self.__parse_line(lines[line_index])
        # print(">>>>>>>>>> variables", variables)
        out_dict["fit_b"] = []
        for variable in variables:
            out_dict["fit_b"].append(float(variable))

        line_index += 1
        variables = self.__parse_line(lines[line_index])
        npoint = int(variables[0])
        out_dict["npoint"] = npoint

        line_index += 1
        variables = self.__parse_line(" ".join(lines[line_index:]))


        Energy = numpy.zeros(npoint)
        F1a = numpy.zeros(npoint)
        F2a = numpy.zeros(npoint)
        F1b = numpy.zeros(npoint)
        F2b = numpy.zeros(npoint)
        iacc = -1
        for i in range(npoint):
            iacc += 1
            Energy[i] = variables[iacc]
            iacc += 1
            F1a[i] = variables[iacc]
            iacc += 1
            F2a[i] = variables[iacc]
            iacc += 1
            F1b[i] = variables[iacc]
            iacc += 1
            F2b[i] = variables[iacc]

        out_dict["Energy"] = Energy
        out_dict["F1a"] = F1a
        out_dict["F2a"] = F2a
        out_dict["F1b"] = F1b
        out_dict["F2b"] = F2b

        self._preprocessor_dictionary = out_dict

    def write_preprocessor_file(self, filename=None):
        if filename is not None: self.set_preprocessor_file(filename)
        self.dump_preprocessor_file(self.get_preprocessor_dictionary(), fileout=self.get_preprocessor_file())

    @classmethod
    def dump_preprocessor_file(cls, out_dict, fileout=""):
        if fileout != "":
            f = open(fileout, 'wt')

            f.write("%i " % out_dict["i_latt"])  # flag ZincBlende
            # f.write("%e " % ((1e0 / volume_in_cm3) * (codata_e2_mc2 * 1e2))) # 1/V*electronRadius
            f.write("%e " % (out_dict["one_over_volume_times_electron_radius_in_cm"]))
            f.write("%e " % out_dict["dspacing_in_cm"])
            f.write("\n")
            f.write("%i " % out_dict["zeta_a"])
            f.write("%i " % out_dict["zeta_b"])
            f.write("%e " % out_dict["temper"])  # temperature parameter
            f.write("\n")
            f.write("(%20.11e,%20.11e ) \n" % (out_dict["ga.real"], out_dict["ga.imag"]))
            f.write("(%20.11e,%20.11e ) \n" % (out_dict["ga_bar.real"], out_dict["ga_bar.imag"]))
            f.write("(%20.11e,%20.11e ) \n" % (out_dict["gb.real"], out_dict["gb.imag"]))
            f.write("(%20.11e,%20.11e ) \n" % (out_dict["gb_bar.real"], out_dict["gb_bar.imag"]))
            f.write("%e %e %e  \n" % (out_dict["fit_a"][0], out_dict["fit_a"][1], out_dict["fit_a"][2]  ))
            f.write("%e %e %e  \n" % (out_dict["fit_b"][0], out_dict["fit_b"][1], out_dict["fit_b"][2]  ))
            f.write(("%i \n") % out_dict["npoint"])
            for i in range(out_dict["npoint"]):
                f.write(("%20.11e %20.11e %20.11e \n %20.11e %20.11e \n") % ( \
                    out_dict["Energy"][i],
                    out_dict["F1a"][i],
                    out_dict["F2a"][i],
                    out_dict["F1b"][i],
                    out_dict["F2b"][i]
                ))
            f.close()
            print("File written to disk: %s" % fileout)



    @classmethod
    def create_from_preprocessor_file(cls, filename):
        out = Bragg(preprocessor_file=filename)
        out.load_preprocessor_file()
        return out

    @classmethod
    def bragg(cls, interactive=True, DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
              TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0, SHADOW_FILE="bragg.dat"):
        """
         SHADOW preprocessor for crystals - python+xraylib version

         -"""

        # codata_e2_mc2 = 2.81794032e-15 = Classical electron radius in S.I.
        codata_e2_mc2 = codata.hbar * codata.alpha / codata.m_e / codata.c

        if interactive:
            print("bragg: SHADOW preprocessor for crystals - python+xraylib version")
            fileout = input("Name of output file : ")

            print(" bragg (python) only works now for ZincBlende Cubic structures. ")
            print(" Valid descriptor are: ")
            print("     Si (alternatively Si_NIST, Si2) ")
            print("     Ge")
            print("     Diamond")
            print("     GaAs, GaSb, GaP")
            print("     InAs, InP, InSb")
            print("     SiC")

            descriptor = input("Name of crystal descriptor : ")

            print("Miller indices of crystal plane of reeflection.")
            miller = input("H K L: ")
            miller = miller.split()
            hh = int(miller[0])
            kk = int(miller[1])
            ll = int(miller[2])

            temper = input("Temperature (Debye-Waller) factor (set 1 for default): ")
            temper = float(temper)

            emin = input("minimum photon energy (eV): ")
            emin = float(emin)
            emax = input("maximum photon energy (eV): ")
            emax = float(emax)
            estep = input("energy step (eV): ")
            estep = float(estep)

        else:
            fileout    = SHADOW_FILE
            descriptor = DESCRIPTOR
            hh         = int(H_MILLER_INDEX)
            kk         = int(K_MILLER_INDEX)
            ll         = int(L_MILLER_INDEX)
            temper     = float(TEMPERATURE_FACTOR)
            emin       = float(E_MIN)
            emax       = float(E_MAX)
            estep      = float(E_STEP)


        #
        # end input section, start calculations
        #
        cryst = xraylib.Crystal_GetCrystal(descriptor)
        volume = cryst['volume']
        volume_in_cm3 = volume * 1e-8 * 1e-8 * 1e-8  # in cm^3
        one_over_volume_times_electron_radius_in_cm = (1e0 / volume_in_cm3) * (codata_e2_mc2 * 1e2)

        # test crystal data - not needed

        if (cryst == None):
            raise Exception("Undefined crystal")
        print("  Unit cell dimensions are %f %f %f" % (cryst['a'], cryst['b'], cryst['c']))
        print("  Unit cell angles are %f %f %f" % (cryst['alpha'], cryst['beta'], cryst['gamma']))
        print("  Unit cell volume is %f A^3" % volume)
        print("  Atoms at:")
        print("     Z  fraction    X        Y        Z")
        for i in range(cryst['n_atom']):
            atom = cryst['atom'][i]
            print("    %3i %f %f %f %f" % (atom['Zatom'], atom['fraction'], atom['x'], atom['y'], atom['z']))
        print("  ")

        # dspacing
        dspacing = xraylib.Crystal_dSpacing(cryst, hh, kk, ll)
        dspacing_in_cm = dspacing * 1e-8
        atom = cryst['atom'] # Z's

        ga = (1e0 + 0j) + numpy.exp(1j * numpy.pi * (hh + kk)) \
             + numpy.exp(1j * numpy.pi * (hh + ll)) \
             + numpy.exp(1j * numpy.pi * (kk + ll))
        gb = ga * numpy.exp(1j * numpy.pi * 0.5 * (hh + kk + ll))
        ga_bar = ga.conjugate()
        gb_bar = gb.conjugate()
        zetas = numpy.array([atom[0]["Zatom"], atom[7]["Zatom"]])
        zeta_a = zetas[0]
        zeta_b = zetas[1]
        npoint = int((emax - emin) / estep + 1)
        for i,zeta in enumerate(zetas):
            xx01 = 1e0 / 2e0 / dspacing
            xx00 = xx01 - 0.1
            xx02 = xx01 + 0.1
            yy00 = xraylib.FF_Rayl(int(zeta), xx00)
            yy01 = xraylib.FF_Rayl(int(zeta), xx01)
            yy02 = xraylib.FF_Rayl(int(zeta), xx02)
            xx = numpy.array([xx00, xx01, xx02])
            yy = numpy.array([yy00, yy01, yy02])
            fit = numpy.polyfit(xx, yy, 2)
            if i == 0:
                fit_a = fit[::-1] # (tuple(fit[::-1].tolist()))
            elif i == 1:
                fit_b = fit[::-1] # (tuple(fit[::-1].tolist()))
            else:
                raise Exception("Unknown crystal structure")

        Energy = numpy.zeros(npoint)
        F1a = numpy.zeros(npoint)
        F1b = numpy.zeros(npoint)
        F2a = numpy.zeros(npoint)
        F2b = numpy.zeros(npoint)


        for i in range(npoint):
            energy = (emin + estep * i)
            Energy[i] = energy
            F1a[i] = xraylib.Fi(int(zetas[0]), energy * 1e-3)
            F2a[i] = abs(xraylib.Fii(int(zetas[0]), energy * 1e-3))
            F1b[i] = xraylib.Fi(int(zetas[1]), energy * 1e-3)
            F2b[i] = abs(xraylib.Fii(int(zetas[1]), energy * 1e-3))

        out_dict = {}
        # inputs
        out_dict["fileout"]    = fileout
        out_dict["descriptor"] = descriptor
        out_dict["hh"]         = hh
        out_dict["kk"]         = kk
        out_dict["ll"]         = ll
        out_dict["temper"]     = temper
        out_dict["emin"]       = emin
        out_dict["emax"]       = emax
        out_dict["estep"]      = estep
        # outputs
        out_dict["i_latt"] = 0
        out_dict["one_over_volume_times_electron_radius_in_cm"] = one_over_volume_times_electron_radius_in_cm
        out_dict["dspacing_in_cm"] = dspacing_in_cm
        out_dict["zeta_a"] = zeta_a
        out_dict["zeta_b"] = zeta_b
        out_dict["temper"] = temper
        out_dict["ga.real"] = ga.real
        out_dict["ga.imag"] = ga.imag
        out_dict["ga_bar.real"] = ga_bar.real
        out_dict["ga_bar.imag"] = ga_bar.imag

        out_dict["gb.real"] = gb.real
        out_dict["gb.imag"] = gb.imag
        out_dict["gb_bar.real"] = gb_bar.real
        out_dict["gb_bar.imag"] = gb_bar.imag

        out_dict["fit_a"] = fit_a
        out_dict["fit_b"] = fit_b
        out_dict["npoint"] = npoint
        out_dict["Energy"] = Energy
        out_dict["F1a"] = F1a
        out_dict["F2a"] = F2a
        out_dict["F1b"] = F1b
        out_dict["F2b"] = F2b

        Bragg.dump_preprocessor_file(out_dict, fileout=fileout)

        return out_dict


    def __energy_index(self, energy):
        Energy = self._preprocessor_dictionary["Energy"]
        if isinstance(energy, int):
            energy = float(energy)
        if isinstance(energy, float):
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
        # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>> NENER", NENER, type(NENER), NENER.any())

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

    def F_elastic(self, ratio):
        CA = self._preprocessor_dictionary["fit_a"]
        CB = self._preprocessor_dictionary["fit_b"]
        FOA = CA[2] * ratio ** 2 + CA[1] * ratio + CA[0]
        FOB = CB[2] * ratio ** 2 + CB[1] * ratio + CB[0]
        return FOA, FOB

    def structure_factor(self, energy, ratio=None):

        if ratio is None:
            ratio = self.ratio(energy, theta=None)

        F1A, F2A, F1B, F2B = self.__interpolate(energy)
        FOA, FOB = self.F_elastic(ratio)
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
            info_text += "\n    dSpacing = %g m" % self.dSpacing()
            info_text += "\n    unit cell volume = %g A^3" % (self.unitcellVolume())
            info_text += "\n    wavelength = %g m" % self.wavelength(energy)
            info_text += "\n    bragg angle = %g rad = %g deg " % ( self.angleBragg(energy), self.angleBragg(energy) * 180 / numpy.pi)
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
                info_text += "\n    dSpacing = %g m" % self.dSpacing()
                info_text += "\n    unit cell volume = %g A^3" % (self.unitcellVolume())
                info_text += "\n    wavelength = %g m" % self.wavelength(ener)
                info_text += "\n    bragg angle = %g rad = %g deg " % ( self.angleBragg(ener), self.angleBragg(ener) * 180 / numpy.pi)
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
            theta = self.angleBragg(energy)

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
        GRAZE = self.angleBragg(energy)
        SSVAR	= RN*(R_LAM0**2)*STRUCT*TEMPER/numpy.pi/numpy.sin(2.0*GRAZE)
        SPVAR = SSVAR * numpy.abs(numpy.cos(2.0 * GRAZE))
        return SSVAR.real, SPVAR.real

    # #
    # # methods as in crystalpy
    # #
    # def F0(self, energy, ratio=None):
    #     """
    #     Calculate F0 from Zachariasen.
    #     :param energy: photon energy in eV.
    #     :return: F0
    #     """
    #     F_0, FH, FH_BAR, STRUCT, FA, FB = self.structure_factor(energy, ratio)
    #     return F_0
    #
    # def FH(self, energy, ratio=None):
    #     """
    #     Calculate FH from Zachariasen.
    #     :param energy: photon energy in eV.
    #     :return: FH
    #     """
    #     F_0, FH, FH_BAR, STRUCT, FA, FB = self.structure_factor(energy, ratio)
    #     return FH
    #
    # def FH_bar(self, energy, ratio=None):
    #     """
    #     Calculate FH_bar from Zachariasen.
    #     :param energy: photon energy in eV.
    #     :return: FH_bar
    #     """
    #     F_0, FH, FH_BAR, STRUCT, FA, FB = self.structure_factor(energy, ratio)
    #     return FH_BAR
    #
    # def dSpacing(self):
    #     """
    #     Returns the lattice spacing d in A
    #     :return: Lattice spacing.
    #     """
    #
    #     return self._preprocessor_dictionary["dspacing_in_cm"] * 1e8
    #
    # def angleBragg(self, energy):
    #     """
    #     Returns the Bragg angle for a given energy.
    #     :param energy: Energy to calculate the Bragg angle for.
    #     :return: Bragg angle.
    #     """
    #
    #     return numpy.arcsin( self.wavelength_in_A(energy) / 2 / self.dSpacing())
    #
    # def unitcellVolume(self):
    #     """
    #     Returns the unit cell volume.
    #
    #     :return: Unit cell volume
    #     """
    #     # Retrieve unit cell volume from xraylib.
    #     one_over_volume_times_electron_radius_in_cm = self.get_preprocessor_dictionary()["one_over_volume_times_electron_radius_in_cm"]
    #     # codata_e2_mc2 = 2.81794032e-15 = Classical electron radius in S.I.
    #     codata_e2_mc2 = codata.hbar * codata.alpha / codata.m_e / codata.c
    #     one_over_volume_in_cm = one_over_volume_times_electron_radius_in_cm / (codata_e2_mc2 * 1e2)
    #     unit_cell_volume = 1.0 / one_over_volume_in_cm * (1e-2)**3
    #     return unit_cell_volume * (1e10)**3

if __name__ == "__main__":
    import numpy
    out = Bragg.bragg(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
          TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0, SHADOW_FILE="bragg.dat")

    #==========================================================================================
    a = Bragg()
    a.run_preprocessor(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
          TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0, SHADOW_FILE="bragg.dat")
    out = a.get_preprocessor_dictionary()


    # ==========================================================================================
    b = Bragg.create_from_preprocessor_file("bragg.dat")
    tmp = b.get_preprocessor_dictionary()

    for key in tmp.keys():
        print("---------------", key)
        print(out[key] - tmp[key])

    print(a.F0(10000.0), a.FH(10000.0), a.FH_bar(10000.0))

    # ==========================================================================================
    b = Bragg.create_from_preprocessor_file("bragg_xop.dat")
    ener = numpy.linspace(8000.0,9000,2)
    ratio = 0.1595
    F_0, FH, FH_BAR, STRUCT, FA, FB = b.structure_factor(ener, ratio=None)
    info_text = b.info(ener)
    print(info_text)
