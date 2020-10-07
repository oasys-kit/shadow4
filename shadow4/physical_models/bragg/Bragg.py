
import numpy
import xraylib
import scipy.constants as codata

class Bragg(object):
    def __init__(self, preprocessor_file=None, preprocessor_dictionary=None):
        self._preprocessor_file = preprocessor_file
        self._preprocessor_dictionary = preprocessor_dictionary

    def run_preprocessor(self, interactive=True, DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
              TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0, SHADOW_FILE="bragg.dat"):

        out_dict = Bragg.bragg(interactive=interactive, DESCRIPTOR=DESCRIPTOR,
                               H_MILLER_INDEX=H_MILLER_INDEX, K_MILLER_INDEX=K_MILLER_INDEX, L_MILLER_INDEX=L_MILLER_INDEX,
                               TEMPERATURE_FACTOR=TEMPERATURE_FACTOR, E_MIN=E_MIN, E_MAX=E_MAX, E_STEP=E_STEP,
                               SHADOW_FILE=SHADOW_FILE)

        return Bragg.__init__(self, preprocessor_file=SHADOW_FILE, preprocessor_dictionary=out_dict)

    def set_preprocessor_file(self, filename):
        self._preprocessor_file = filename

    def get_preprocessor_file(self):
        return self._preprocessor_file

    def get_preprocessor_dictionary(self):
        return self._preprocessor_dictionary

    def load_preprocessor_file(self, filename=None):
        if filename is None:
            filename = self._preprocessor_file
        else:
            self._preprocessor_file = filename

        f = open(filename, 'r')
        # lines = f.readlines()
        lines = f.read().splitlines()
        f.close()

        print(lines)
        out_dict = {}

        line_index = 0
        line = lines[line_index]
        line = " ".join(line.split())
        variables = line.split(" ")
        print(">>>>>>>>>> variables", variables)
        out_dict["i_latt"] = int(variables[0])
        out_dict["one_over_volume_times_electron_radius_in_cm"] = float(variables[1])
        out_dict["dspacing_in_A"] = float(variables[2])

        line_index += 1
        line = lines[line_index]
        line = " ".join(line.split())
        variables = line.split(" ")
        out_dict["zeta_a"] = int(variables[0])
        out_dict["zeta_b"] = int(variables[1])
        out_dict["temper"] = float(variables[2])

        # line_index += 1

        line_index += 1
        line = lines[line_index]
        line = line.replace("(", "")
        line = line.replace(")", "")
        line = line.replace(" ", "")
        variables = line.split(",")
        # print(">>>>>>>>>> variables", variables)
        out_dict["ga.real"] = float(variables[0])
        out_dict["ga.imag"] = float(variables[1])

        line_index += 1
        line = lines[line_index]
        line = line.replace("(", "")
        line = line.replace(")", "")
        line = line.replace(" ", "")
        variables = line.split(",")
        # print(">>>>>>>>>> variables", variables)
        out_dict["ga_bar.real"] = float(variables[0])
        out_dict["ga_bar.imag"] = float(variables[1])

        line_index += 1
        line = lines[line_index]
        line = line.replace("(", "")
        line = line.replace(")", "")
        line = line.replace(" ", "")
        variables = line.split(",")
        # print(">>>>>>>>>> variables", variables)
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
        line = lines[line_index]
        line = " ".join(line.split())
        variables = line.split(" ")
        # print(">>>>>>>>>> variables", variables)
        out_dict["fit_a"] = []
        for variable in variables:
            out_dict["fit_a"].append(float(variable))

        line_index += 1
        line = lines[line_index]
        line = " ".join(line.split())
        variables = line.split(" ")
        # print(">>>>>>>>>> variables", variables)
        out_dict["fit_b"] = []
        for variable in variables:
            out_dict["fit_b"].append(float(variable))

        line_index += 1
        line = lines[line_index]
        line = " ".join(line.split())
        variables = line.split(" ")
        npoint = int(variables[0])
        out_dict["npoint"] = npoint

        line_index += 1
        text = " ".join(lines[line_index:])
        variables = text.split()
        # print(">>>> text: ", text)
        # print(">>>> text: ", variables)


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
        # retrieve physical constants needed
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
        dspacing_in_A = dspacing * 1e-8
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
            # print "zeta: ",zeta
            # print "z,xx,YY: ",zeta,xx,yy
            # print "fit: ",fit[::-1] # reversed coeffs
            # print "fit-tuple: ",(tuple(fit[::-1].tolist())) # reversed coeffs
            # print("fit-tuple: %e %e %e  \n" % (tuple(fit[::-1].tolist())) ) # reversed coeffs
            # f.write("%e %e %e  \n" % (tuple(fit[::-1].tolist())))  # reversed coeffs

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
        out_dict["dspacing_in_A"] = dspacing_in_A
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



        if fileout != "":
            f = open(fileout, 'wt')

            f.write("%i " % out_dict["i_latt"])  # flag ZincBlende
            # f.write("%e " % ((1e0 / volume_in_cm3) * (codata_e2_mc2 * 1e2))) # 1/V*electronRadius
            f.write("%e " % (out_dict["one_over_volume_times_electron_radius_in_cm"]))
            f.write("%e " % out_dict["dspacing_in_A"])
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
            for i in range(npoint):
                f.write(("%20.11e %20.11e %20.11e \n %20.11e %20.11e \n") % ( \
                    out_dict["Energy"][i],
                    out_dict["F1a"][i],
                    out_dict["F2a"][i],
                    out_dict["F1b"][i],
                    out_dict["F2b"][i]
                ))
            f.close()
            print("File written to disk: %s" % fileout)


        return out_dict

    def F0(self, energy):
        """
        Calculate F0 from Zachariasen.
        :param energy: photon energy in eV.
        :return: F0
        """
        F_0 = 0
        return F_0

    def FH(self, energy):
        """
        Calculate FH from Zachariasen.
        :param energy: photon energy in eV.
        :return: FH
        """
        F_H = 0
        return F_H

    def FH_bar(self, energy):
        """
        Calculate FH_bar from Zachariasen.
        :param energy: photon energy in eV.
        :return: FH_bar
        """
        F_H_bar = 0

        return F_H_bar

    def energy_index(self, energy):
        Energy = self._preprocessor_dictionary["Energy"]
        if (energy < Energy.min()) or (energy > Energy.max()):
            return -100
            # raise Exception("Energy %f outside the defined interval [%f, %f]" % (energy, Energy.min(), Energy.max() ))

        ll = numpy.where( Energy > energy)
        # print(">>> ll: ", ll, Energy.size, Energy.min(), Energy.max())
        return ll[0][0]

    def interpolate(self, PHOT):
        ENERGY = self._preprocessor_dictionary["Energy"]
        NENER = self.energy_index(PHOT) - 1
        FP_A = self._preprocessor_dictionary["F1a"]
        FP_B = self._preprocessor_dictionary["F1b"]
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
        # print(">>> CA", CA, ratio)
        FOA = CA[2] * ratio ** 2 + CA[1] * ratio + CA[0]
        FOB = CB[2] * ratio ** 2 + CB[1] * ratio + CB[0]
        # print("FOA, FOB", FOA, FOB)
        return FOA, FOB

    def structure_factor(self, energy, ratio):

        F1A, F2A, F1B, F2B = self.interpolate(energy)
        FOA, FOB = self.F_elastic(ratio)
        FA = FOA + F1A + 1j * F2A
        FB = FOB + F1B + 1j * F2B
        print(">> For atom A, fo + f' + if'' = ", FA)
        print(">> For atom B, fo + f' + if'' = ", FB)

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
        # ELSE IF (I_LATT.EQ.1) THEN
        # ABSORP = 2.0D0*RN*R_LAM0*(4.0D0*(DIMAG(FA)+DIMAG(FB)))
        # F_0 = 4*((F1A + ATNUM_A + F1B + ATNUM_B) + CI*(F2A + F2B))
        # ELSE IF (I_LATT.EQ.2) THEN
        # FB	 = (0.0D0,0.0D0)
        # ABSORP = 2.0D0*RN*R_LAM0*(4.0D0*DIMAG(FA))
        # F_0 = 4*(F1A + ATNUM_A + CI*F2A)
        # ELSE IF (I_LATT.EQ.3) THEN
        # ABSORP = 2.0D0*RN*R_LAM0*(DIMAG(FA)+DIMAG(FB))
        # F_0 = (F1A + ATNUM_A + F1B + ATNUM_B) + CI*(F2A + F2B)
        # ELSE IF (I_LATT.EQ.4) THEN
        # FB     = (0.0D0,0.0D0)
        # ABSORP = 2.0D0*RN*R_LAM0*(2.0D0*(DIMAG(FA)))
        # F_0 = 2*(F1A+ CI*F2A )
        # ELSE IF (I_LATT.EQ.5) THEN
        # FB     = (0.0D0,0.0D0)
        # ABSORP = 2.0D0*RN*R_LAM0*(4.0D0*(DIMAG(FA)))
        # F_0 = 4*(F1A + CI*F2A )
        # END IF


        # ! C
        # ! C FH and FH_BAR are the structure factors for (h,k,l) and (-h,-k,-l).
        # ! C
        # ! C srio, Added TEMPER here (95/01/19)
        FH 	= ( (GA * FA) + (GB * FB) ) * TEMPER
        FH_BAR	= ( (GA_BAR * FA) + (GB_BAR * FB) ) * TEMPER

        return F_0, FH, FH_BAR

if __name__ == "__main__":
    import numpy
    # out = Bragg.bragg(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
    #       TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0, SHADOW_FILE="bragg.dat")

    #==========================================================================================
    # a = Bragg()
    # a.run_preprocessor(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
    #       TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0, SHADOW_FILE="bragg.dat")
    # out = a.get_preprocessor_dictionary()
    #
    #
    # b = Bragg.create_from_preprocessor_file("bragg.dat")
    # # b.load_preprocessor_file("bragg.dat")
    # tmp = b.get_preprocessor_dictionary()
    #
    # for key in tmp.keys():
    #     print("---------------", key)
    #     print(out[key] - tmp[key])
    #
    # print(a.F0(10000.0), a.FH(10000.0), a.FH_bar(10000.0))

    # ==========================================================================================
    b = Bragg.create_from_preprocessor_file("bragg_xop.dat")
    # for energy in numpy.linspace(4000, 16000, 100):
    #     print("Energy index for %f : %f" % (energy, b.energy_index(energy)))
    # ratio = 0.1595 # sin_theta_over_lambda
    ener = numpy.linspace(8000.0,9000,10)
    ratio = 0.1595
    F0, FH, FH_BAR = b.structure_factor(ener, ratio)
    print("Structure factor F0: (%g, %g)" % (F0.real, F0.imag))
    print("Structure factor FH: (%g, %g)" % (FH.real, FH.imag))
    print("Structure factor FH_BAR: (%g, %g)" % (FH_BAR.real, FH_BAR.imag))
    struct = numpy.sqrt( FH * FH_BAR)
    print("Structure factor STRUCT: (%g, %g)" % (struct.real, struct.imag))