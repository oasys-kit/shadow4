"""
python version of the multilayer code in shadow.

The stack is as follows: Vacuum+[Odd,Even]xn+Substrate

>>> #                  vacuum   ")
>>> #       |------------------------------|  \
>>> #       |          odd (n)             |  |
>>> #       |------------------------------|  | BILAYER # n
>>> #       |          even (n)            |  |
>>> #       |------------------------------|  /
>>> #       |          .                   |
>>> #       |          .                   |
>>> #       |          .                   |
>>> #       |------------------------------|  \
>>> #       |          odd (1)             |  |
>>> #       |------------------------------|  | BILAYER # 1
>>> #       |          even (1)            |  |
>>> #       |------------------------------|  /
>>> #       |                              |
>>> #       |///////// substrate //////////|
>>> #       |                              |

"""

# TODO: graded multilayers

import numpy
import scipy.constants as codata
from srxraylib.util.h5_simple_writer import H5SimpleWriter

tocm = codata.h * codata.c / codata.e*1e2 # 12398.419739640718e-8

class MLayer(object):
    """
    Constructor.

    Several ways of using MLayer:
    * ml = MLayer.pre_mlayer(<keywords>)
    * ml = MLayer(); ml.read_preprocessor_file('mymlayer.dat')
    * ml = MLayer.initialize_from_bilayer_stack(<keywords>)

    """
    def __init__(self):
        self.using_pre_mlayer = False # False=use craylib, True=use preprocessor
        self.pre_mlayer_dict = None

    def read_preprocessor_file(self, filename):
        """
        Reads a preprocessor (pre_mlayer) file. The same as in shadow3.

        Parameters
        ----------
        filename : str
            The name of the filename.
        """
        out_dict = {}

        fp = open(filename) # Open file on read mode
        lines = fp.read().split("\n") # Create a list containing all lines
        fp.close() # Close file

        index_pointer = 0
        np = int(lines[index_pointer])
        out_dict['np'] = np
        energy = numpy.zeros(np)

        index_pointer += 1
        mylist = lines[index_pointer].split()
        for i in range(np):
                energy[i] = float(mylist[i])

        out_dict["energy"] = energy

        delta_s = numpy.zeros_like(energy)
        beta_s = numpy.zeros_like(energy)

        delta_e = numpy.zeros_like(energy)
        beta_e = numpy.zeros_like(energy)

        delta_o = numpy.zeros_like(energy)
        beta_o = numpy.zeros_like(energy)

        for i in range(np):
            index_pointer += 1
            mylist = lines[index_pointer].strip().split()
            delta_s[i] = float(mylist[0])
            beta_s[i] = float(mylist[1])

        for i in range(np):
            index_pointer += 1
            mylist = lines[index_pointer].strip().split()
            delta_e[i] = float(mylist[0])
            beta_e[i] = float(mylist[1])

        for i in range(np):
            index_pointer += 1
            mylist = lines[index_pointer].strip().split()
            delta_o[i] = float(mylist[0])
            beta_o[i] = float(mylist[1])


        out_dict["delta_s"] = delta_s
        out_dict["beta_s"] = beta_s

        out_dict["delta_e"] = delta_e
        out_dict["beta_e"] = beta_e

        out_dict["delta_o"] = delta_o
        out_dict["beta_o"] = beta_o

        # #! srio@esrf.eu 2012-06-07 Nevot-Croce ML roughness model implemented.
        # #! By convention, starting from the version that includes ML roughness
        # #! we set NPAR negative, in order to assure compatibility with old
        # #! versions. If NPAR<0, roughness data are read, if NPAR>0 no roughness.

        index_pointer += 1
        npair = int(lines[index_pointer])

        out_dict["npair"] = npair

        thick = numpy.zeros(numpy.abs(npair))
        gamma1 = numpy.zeros_like(thick)
        mlroughness1 = numpy.zeros_like(thick)
        mlroughness2 = numpy.zeros_like(thick)

        for i in range(numpy.abs(npair)):
            index_pointer += 1
            mylist = lines[index_pointer].strip().split()
            thick[i] = float(mylist[0])
            gamma1[i] = float(mylist[1])
            mlroughness1[i] = float(mylist[2])
            mlroughness2[i] = float(mylist[3])

        out_dict["thick"] = thick
        out_dict["gamma1"] = gamma1
        out_dict["mlroughness1"] = mlroughness1
        out_dict["mlroughness2"] = mlroughness2

        index_pointer += 1
        igrade = int(lines[index_pointer])

        out_dict["igrade"] = igrade

        if igrade == 1:
            index_pointer += 1
            fgrade = int(lines[index_pointer])

            out_dict["fgrade"] = fgrade

        elif igrade == 2:  # igrade=2,

            index_pointer += 1
            mylist = lines[index_pointer].strip().split("   ")
            a0 = float(mylist[0])
            a1 = float(mylist[1])
            a2 = float(mylist[2])
            a3 = float(mylist[3])

            out_dict["a0"] = a0
            out_dict["a1"] = a1
            out_dict["a2"] = a2
            out_dict["a3"] = a3

        self.pre_mlayer_dict = out_dict
        self.using_pre_mlayer = True

    # this is copied from shadow3 python preprocessors
    @classmethod
    def pre_mlayer(cls,
                   interactive=False,
                   FILE="pre_mlayer.dat",
                   E_MIN=5000.0, E_MAX=20000.0,
                   S_DENSITY=2.33, S_MATERIAL="Si",
                   E_DENSITY=2.40, E_MATERIAL="B4C",
                   O_DENSITY=9.40, O_MATERIAL="Ru",
                   N_PAIRS=70,
                   THICKNESS=33.1,
                   GAMMA=0.483,
                   ROUGHNESS_EVEN=3.3, ROUGHNESS_ODD=3.1,
                   GRADE_DEPTH=0,
                   GRADE_SURFACE=0,
                   FILE_DEPTH="myfile_depth.dat",
                   FILE_SHADOW="mlayer1.sha",
                   FILE_THICKNESS="mythick.dat",
                   FILE_GAMMA="mygamma.dat",
                   AA0=1.0, AA1=0.0, AA2=0.0, AA3=0.0):
        """
        Creates an instance of MLayer with parameters initialized from the keyword parameters and the
        pre_mlayer preprocessor file. It uses xraylib for accessing the optical constants.

        Parameters
        ----------
        interactive : bool, optional
            Set True for running interactively in the terminal and answer the questions (like in shadow2).
        FILE : str, optional
            File name to be created (preprocessor pre_mlayer file).
        E_MIN : float, optional
            The minimum photon energy in eV (for creating the tabulated refraction index).
        E_MAX : float, optional
            The maximum photon energy in eV (for creating the tabulated refraction index).
        S_DENSITY : float, optional
            Density (in g/cm3) for the SUBSTRATE material.
        S_MATERIAL : float, optional
            Material formula for the SUBSTRATE material.
        E_DENSITY : float, optional
            Density (in g/cm3) for the EVEN material.
        E_MATERIAL : float, optional
            Material formula for the EVEN material.
        O_DENSITY : float, optional
            Density (in g/cm3) for the ODD material.
        O_MATERIAL : float, optional
            Material formula for the ODD material.
        N_PAIRS : int, optional
            The number of bilayers of the multilayer.
        THICKNESS : float, optional
            The thickness of a bilayer in Angstroms.
        GAMMA : float, optional
            The gamma factor thickness(even) / (thickness(odd) + thickness(even)) for the multilayer.
        ROUGHNESS_EVEN : float, optional
            The rounghness of the EVEN layers in A.
        ROUGHNESS_ODD : float, optional
            The rounghness of the ODD layers in A.
        GRADE_DEPTH : int, optional
            Flag to indicated that the multilayes has a graded depth.
        GRADE_SURFACE : int, optional
            Flag for multilayer graded over the surface:
             * 0: No.
             * 1: t and/or gamma graded over the surface (input spline files with t and gamma gradient).
             * 2: t graded over the surface (input quadratic fit to t gradient).

        FILE_DEPTH : str, optional
            ????
        FILE_SHADOW : str, optional
            ????
        FILE_THICKNESS : str, optional
            ????
        FILE_GAMMA : str, optional
            ????
        AA0 : float, optional
            Coefficient when GRADED_SURFACE=1
        AA1 : float, optional
            Coefficient when GRADED_SURFACE=1
        AA2 : float, optional
            Coefficient when GRADED_SURFACE=1
        AA3 : float, optional
            Coefficient when GRADED_SURFACE=1

        Returns
        -------
        instance of MLayer

        Notes
        -----
        This method requires xraylib to access to the optical constants.
        """
        import xraylib

        # input section
        if interactive:
            print("pre_mlayer: SHADOW preprocessor for multilayers - python+xraylib version")
            fileout = input("Name of output file : ")
            estart = input("Photon energy (eV) from : ")
            estart = float(estart)
            efinal = input("                     to : ")
            efinal = float(efinal)

            print("  ")
            print("The stack is as follows: ")
            print("      ")
            print("                 vacuum   ")
            print("      |------------------------------|  \   ")
            print("      |          odd (n)             |  |   ")
            print("      |------------------------------|  | BILAYER # n   ")
            print("      |          even (n)            |  |   ")
            print("      |------------------------------|  /   ")
            print("      |          .                   |   ")
            print("      |          .                   |   ")
            print("      |          .                   |   ")
            print("      |------------------------------|  \   ")
            print("      |          odd (1)             |  |   ")
            print("      |------------------------------|  | BILAYER # 1   ")
            print("      |          even (1)            |  |   ")
            print("      |------------------------------|  /   ")
            print("      |                              |   ")
            print("      |///////// substrate //////////|   ")
            print("      |                              |   ")
            print("      ")
            print(" ")

            # substrate
            matSubstrate = input("Specify the substrate material : ")
            denSubstrate = input("Specify the substrate density [g/cm^3] : ")
            denSubstrate = float(denSubstrate)

            print("Right above the substrate is the even layer material")
            matEven = input("Specify the even layer material : ")
            denEven = input("Specify the even layer density [g/cm^3] : ")
            denEven = float(denEven)

            print("Odd layer material is on top of the even layer.")
            matOdd = input("Specify the odd layer material : ")
            denOdd = input("Specify the odd layer density [g/cm^3] : ")
            denOdd = float(denOdd)

            #! By convention, starting from the version that includes ML roughness
            #! we set NPAR negative, in order to assure compatibility with old
            #! versions. If NPAR<0, roughness data are read, if NPAR>0 no roughness.
            npair = input("No. of layer pairs : ")
            npair = int(npair)

            print(" ")
            print("Starting from the substrate surface, specify the thickness t :")
            print("      t = t(odd) + t(even)        in Angstroms,")
            print("and the gamma ratio :")
            print("      t(even) / (t(odd) + t(even))")
            print("for EACH bilayer.")
            print(" ")
            print("Type two -1 whenever you want the remaining layers ")
            print("to assume the thickness, gamma ratio and roughnesses of the previous one.")
            print(" ")

            #define variables
            thick=[0e0]*npair
            gamma1=[0e0]*npair
            mlroughness1=[0e0]*npair
            mlroughness2=[0e0]*npair

            for i in range(npair):
                tmps = ("thickness [A], gamma ratio, roughness even [A] and roughness odd [A] of bilayer %i: \n"% (i+1) )
                tmp = input(tmps)
                tmp1 = tmp.split()
                if ((i != 0) and (int(float(tmp1[0])) == -1)):
                    thick[i:(npair-1)] = [thick[i-1]] * (npair-i)
                    gamma1[i:(npair-1)] = [gamma1[i-1]] * (npair-i)
                    mlroughness1[i:(npair-1)] = [mlroughness1[i-1]] * (npair-i)
                    mlroughness2[i:(npair-1)] = [mlroughness2[i-1]] * (npair-i)
                    break
                else:
                    thick[i] = float(tmp1[0])
                    gamma1[i] = float(tmp1[1])
                    mlroughness1[i] = float(tmp1[2])
                    mlroughness2[i] = float(tmp1[3])

            print("***************************************************")
            print("  Is the multilayer graded over the surface? ")
            print("      0: No ")
            print("      1: t and/or gamma graded over the surface ")
            print("         (input spline files with t and gamma gradient")
            print("      2: t graded over the surface ")
            print("         (input quadratic fit to t gradient)")
            print("      ")

            igrade = input("Is t and/or gamma graded over the surface [0=No/1=Yes] ? ")
            igrade = int(igrade)
            if igrade == 1:
                print("Generation of the spline coefficients for the t and gamma factors")
                print("over the surface.")
                print("Then GRADE_MLAYER should be run to generate the spline ")
                print("coefficients for the t and gamma factors over the surface.")
                print("Here just type in the file name that WILL be used to store")
                print("the spline coefficients :")
                fgrade = input("File name (output from grade_mlayer: ")
            elif igrade == 2:  # igrade=2, coefficients
                print("A second degree polynomial fit of the thickness grading")
                print("must be available:")
                print("t(y) = BILATER_THICHNESS(y)/BILAYER_THICKNESS(y=0)")
                print("t(y) = a0 + a1*y + a2*(y^2) + a3*(y^3)  ")
                print("a0 (constant term) ")
                print("a1 (slope term) ")
                print("a2 (quadratic term) ")
                print("a3 (cubic term) ")
                tmp = input("Enter a0, a1, a2, a3: ")
                tmp = tmp.split()
                a0 = float(tmp[0])
                a1 = float(tmp[1])
                a2 = float(tmp[2])
                a3 = float(tmp[3])
        else:
            #--- From input keywords...
            fileout = FILE
            estart = float(E_MIN)
            efinal = float(E_MAX)

            # substrate
            matSubstrate = S_MATERIAL
            denSubstrate = float(S_DENSITY)

            matEven = E_MATERIAL
            denEven = float(E_DENSITY)

            matOdd = O_MATERIAL
            denOdd = float(O_DENSITY)

            npair = int(N_PAIRS)

            #define variables
            thick=[0e0]*npair
            gamma1=[0e0]*npair
            mlroughness1=[0e0]*npair
            mlroughness2=[0e0]*npair

            for i in range(npair):
                thick[i] = float(THICKNESS)
                gamma1[i] = float(GAMMA)
                mlroughness1[i] = float(ROUGHNESS_EVEN)
                mlroughness2[i] = float(ROUGHNESS_ODD)

            igrade = int(GRADE_SURFACE)

            #TODO: check if needed file_gamma
            fgrade = FILE_THICKNESS # raw_input("File name (output from grade_mlayer: ")

            a0 = float(AA0)
            a1 = float(AA1)
            a2 = float(AA2)
            a3 = float(AA3)

        elfactor = numpy.log10(1.0e4/30.0)/300.0
        istart = int(numpy.log10(estart/30.0e0)/elfactor + 1)
        ifinal = int(numpy.log10(efinal/30.0e0)/elfactor + 2)
        np = int(ifinal - istart) + 1


        f = open(fileout, 'wt')


        pre_mlayer_dict = {}

        f.write("%i \n" % np)
        pre_mlayer_dict["np"] = np


        ENERGY = numpy.zeros(np)
        for i in range(np):
            energy = 30e0*numpy.power(10,elfactor*(istart+i-1))
            f.write("%e " % energy)
            ENERGY[i] = energy
        f.write( "\n")
        pre_mlayer_dict["energy"] = ENERGY

        DELTA = numpy.zeros(np)
        BETA = numpy.zeros(np)
        for i in range(np):  #substrate
            energy = 30e0*numpy.power(10,elfactor*(istart+i-1)) *1e-3 # in keV!!
            delta = 1e0-xraylib.Refractive_Index_Re(matSubstrate,energy,denSubstrate)
            beta = xraylib.Refractive_Index_Im(matSubstrate,energy,denSubstrate)
            DELTA[i] = delta
            BETA[i] = beta
            f.write( ("%26.17e "*2+"\n") % tuple([delta,beta]) )
        pre_mlayer_dict["delta_s"] = DELTA
        pre_mlayer_dict["beta_s"] = BETA

        DELTA = numpy.zeros(np)
        BETA = numpy.zeros(np)
        for i in range(np): #even
            energy = 30e0*numpy.power(10,elfactor*(istart+i-1)) *1e-3 # in keV!!
            delta = 1e0-xraylib.Refractive_Index_Re(matEven,energy,denEven)
            beta = xraylib.Refractive_Index_Im(matEven,energy,denEven)
            DELTA[i] = delta
            BETA[i] = beta
            f.write( ("%26.17e  "*2+"\n") % tuple([delta,beta]) )
        pre_mlayer_dict["delta_e"] = DELTA
        pre_mlayer_dict["beta_e"] = BETA

        DELTA = numpy.zeros(np)
        BETA = numpy.zeros(np)
        for i in range(np): #odd
            energy = 30e0*numpy.power(10,elfactor*(istart+i-1)) *1e-3 # in keV!!
            delta = 1e0-xraylib.Refractive_Index_Re(matOdd,energy,denOdd)
            beta = xraylib.Refractive_Index_Im(matOdd,energy,denOdd)
            DELTA[i] = delta
            BETA[i] = beta
            f.write( ("%26.17e "*2+"\n") % tuple([delta,beta]) )
        pre_mlayer_dict["delta_o"] = DELTA
        pre_mlayer_dict["beta_o"] = BETA


        #! srio@esrf.eu 2012-06-07 Nevot-Croce ML roughness model implemented.
        #! By convention, starting from the version that includes ML roughness
        #! we set NPAR negative, in order to assure compatibility with old
        #! versions. If NPAR<0, roughness data are read, if NPAR>0 no roughness.
        f.write("%i \n" % -npair)
        pre_mlayer_dict["npair"] = -npair

        for i in range(npair):
            f.write( ("%26.17e "*4+"\n") % tuple([thick[i],gamma1[i],mlroughness1[i],mlroughness2[i]]) )
        pre_mlayer_dict["thick"]        = numpy.array(thick)
        pre_mlayer_dict["gamma1"]       = numpy.array(gamma1)
        pre_mlayer_dict["mlroughness1"] = numpy.array(mlroughness1)
        pre_mlayer_dict["mlroughness2"] = numpy.array(mlroughness2)

        f.write("%i \n" % igrade)
        pre_mlayer_dict["igrade"] = igrade
        if igrade == 1:
            f.write("%s \n" % fgrade)
            pre_mlayer_dict["fgrade"] = fgrade
        elif igrade == 2:  # igrade=2, coefficients
            f.write("%f  %f  %f  %f\n"%(a0,a1,a2,a3))
            pre_mlayer_dict["a0"] = a0
            pre_mlayer_dict["a1"] = a1
            pre_mlayer_dict["a2"] = a2
            pre_mlayer_dict["a3"] = a3

        f.close()
        print("File written to disk: %s" % fileout)


        out = MLayer()
        out.pre_mlayer_dict = pre_mlayer_dict
        out.using_pre_mlayer = True
        return out


    @classmethod
    def initialize_from_bilayer_stack(cls,
            material_S="Si", density_S=None, roughness_S=0.0,
            material_E="B4C",density_E=None, roughness_E=0.0,
            material_O="Ru", density_O=None, roughness_O=0.0,
            bilayer_pairs=70,
            bilayer_thickness=33.1,
            bilayer_gamma=0.483,
            ):
        """
        Creates an instance of MLayer with the main parameters of a multilater. In this case, the calculations
        are done without using the preprocessor data, and the optical constants are supplied by xraylib at run time.

        Parameters
        ----------
        material_S : str, optional
            The symbol of the SUBSTRATE material.
        density_S : float, optional
            The density in g/cm3 for the SUBSTRATE layers.
        roughness_S : float, optional
            The roughness in A for the SUBSTRATE layers.
        material_E : str, optional
            The symbol of the EVEN layer material.
        density_E : float, optional
            The density in g/cm3 for the EVEN layers.
        roughness_E : float, optional
            The roughness in A for the EVEN layers.
        material_O : str, optional
            The symbol of the ODD layer material.
        density_O : float, optional
            The density in g/cm3 for the ODD layers.
        roughness_O : float, optional
            The roughness in A for the ODD layers.
        bilayer_pairs : int, optional
            The number of bilayers.
        bilayer_thickness : float, optional
            The thickness in A of the bilayers.
        bilayer_gamma : float, optional
            The gamma factor of the bilayer thickness(even) / (thichness(odd) + thichness(even)).

        Returns
        -------
        instance of MLayer
        """

        import xraylib

        npair = int(bilayer_pairs)

        #define variables
        thick        = []
        gamma1       = []
        mlroughness1 = []
        mlroughness2 = []

        for i in range(npair):
            thick.append(bilayer_thickness)
            gamma1.append(bilayer_gamma)
            mlroughness1.append(roughness_E)
            mlroughness2.append(roughness_O)

        pre_mlayer_dict = {}

        pre_mlayer_dict["np"] = bilayer_pairs


        #! srio@esrf.eu 2012-06-07 Nevot-Croce ML roughness model implemented.
        #! By convention, starting from the version that includes ML roughness
        #! we set NPAR negative, in order to assure compatibility with old
        #! versions. If NPAR<0, roughness data are read, if NPAR>0 no roughness.
        pre_mlayer_dict["npair"] = -npair

        pre_mlayer_dict["thick"]        = numpy.array(thick)
        pre_mlayer_dict["gamma1"]       = numpy.array(gamma1)
        pre_mlayer_dict["mlroughness1"] = numpy.array(mlroughness1)
        pre_mlayer_dict["mlroughness2"] = numpy.array(mlroughness2)

        #These keys are not in the original pre_mlayer_dict
        pre_mlayer_dict["material1"] = material_E
        pre_mlayer_dict["material2"] = material_O
        pre_mlayer_dict["materialS"] = material_S

        pre_mlayer_dict["roughnessS"] = roughness_S

        pre_mlayer_dict["density1"] = density_E
        pre_mlayer_dict["density2"] = density_O
        pre_mlayer_dict["densityS"] = density_S

        if pre_mlayer_dict["densityS"] is None:
            try:
                pre_mlayer_dict["densityS"] = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(pre_mlayer_dict["materialS"]))
                print("Using density for substrate (%s): %f"%(pre_mlayer_dict["materialS"],
                      pre_mlayer_dict["densityS"]))
            except:
                raise Exception("Failed to load density for material: %s"%(pre_mlayer_dict["material1"]))
        if pre_mlayer_dict["density1"] is None:
            try:
                pre_mlayer_dict["density1"] = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(pre_mlayer_dict["material1"]))
                print("Using density for layer 1 (even) (%s): %f" % (pre_mlayer_dict["material1"],
                      pre_mlayer_dict["density1"]))
            except:
                raise Exception("Failed to load density for material: %s"%(pre_mlayer_dict["material1"]))
        if pre_mlayer_dict["density2"] is None:
            try:
                pre_mlayer_dict["density2"] = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(pre_mlayer_dict["material2"]))
                print("Using density for layer 2 (odd) (%s): %f" % (pre_mlayer_dict["material2"],
                      pre_mlayer_dict["density2"]))
            except:
                raise Exception("Failed to load density for material: %s"%(pre_mlayer_dict["material2"]))

        if isinstance(pre_mlayer_dict["densityS"],str):
            pre_mlayer_dict["densityS"] = float(pre_mlayer_dict["densityS"])
        if isinstance(pre_mlayer_dict["density1"],str):
            pre_mlayer_dict["density1"] = float(pre_mlayer_dict["density1"])
        if isinstance(pre_mlayer_dict["density2"], str):
            pre_mlayer_dict["density2"] = float(pre_mlayer_dict["density2"])


        # fill unused keys
        pre_mlayer_dict["energy"] = None
        pre_mlayer_dict["delta_s"] = None
        pre_mlayer_dict["beta_s"] = None
        pre_mlayer_dict["delta_e"] = None
        pre_mlayer_dict["beta_e"] = None
        pre_mlayer_dict["delta_o"] = None
        pre_mlayer_dict["beta_o"] = None
        pre_mlayer_dict["igrade"] = None
        if pre_mlayer_dict["igrade"] == 1:
            pre_mlayer_dict["fgrade"] = None
        elif pre_mlayer_dict["igrade"] == 2:  # igrade=2, coefficients
            pre_mlayer_dict["a0"] = None
            pre_mlayer_dict["a1"] = None
            pre_mlayer_dict["a2"] = None
            pre_mlayer_dict["a3"] = None

        # return
        out = MLayer()
        out.pre_mlayer_dict = pre_mlayer_dict
        out.using_pre_mlayer = False
        return out

    @classmethod
    def initialize_from_bilayer_stack_in_compressed_format(cls,
            structure='[Pd/B4C]x150+Si',
            density_O=None, roughness_O=0.0,
            density_E=None, roughness_E=0.0,
            density_S=None, roughness_S=0.0,
            bilayer_thickness=33.1,
            bilayer_gamma=0.483,
            ):
        """
        Creates an instance of MLayer with the main parameters of a multilater. In this case, the calculations
        are done without using the preprocessor data, and the optical constants are supplied by xraylib at run time.

        Parameters
        ----------
        material_S : str, optional
            The symbol of the SUBSTRATE material.
        density_S : float, optional
            The density in g/cm3 for the SUBSTRATE layers.
        roughness_S : float, optional
            The roughness in A for the SUBSTRATE layers.
        material_E : str, optional
            The symbol of the EVEN layer material.
        density_E : float, optional
            The density in g/cm3 for the EVEN layers.
        roughness_E : float, optional
            The roughness in A for the EVEN layers.
        material_O : str, optional
            The symbol of the ODD layer material.
        density_O : float, optional
            The density in g/cm3 for the ODD layers.
        roughness_O : float, optional
            The roughness in A for the ODD layers.
        bilayer_pairs : int, optional
            The number of bilayers.
        bilayer_thickness : float, optional
            The thickness in A of the bilayers.
        bilayer_gamma : float, optional
            The gamma factor of the bilayer thickness(even) / (thichness(odd) + thichness(even)).

        Returns
        -------
        instance of MLayer
        """

        import xraylib

        #
        # parse
        #

        i0 = structure.find('[')
        i1 = structure.find('/')
        i2 = structure.find(']')
        i3 = structure.find('x')
        i4 = structure.find('+')
        if i0 * i2 * i3 * i4 < 0:
            raise Exception("Bad format for multilayer structure %s (must be: %s)" % (structure, ['[W/B]x50+Si']))
        material_O = structure[(i0+1):(i1)]
        material_E = structure[(i1+1):(i2)]
        npair =  int(structure[(i3+1):(i4)])
        material_S = structure[(i4+1)::]

        print("[%s,%s]x%d+%s" % (material_O, material_E, npair, material_S))

        #define variables
        thick        = []
        gamma1       = []
        mlroughness1 = []
        mlroughness2 = []

        for i in range(npair):
            thick.append(bilayer_thickness)
            gamma1.append(bilayer_gamma)
            mlroughness1.append(roughness_E)
            mlroughness2.append(roughness_O)

        pre_mlayer_dict = {}

        pre_mlayer_dict["np"] = npair


        #! srio@esrf.eu 2012-06-07 Nevot-Croce ML roughness model implemented.
        #! By convention, starting from the version that includes ML roughness
        #! we set NPAR negative, in order to assure compatibility with old
        #! versions. If NPAR<0, roughness data are read, if NPAR>0 no roughness.
        pre_mlayer_dict["npair"] = -npair

        pre_mlayer_dict["thick"]        = numpy.array(thick)
        pre_mlayer_dict["gamma1"]       = numpy.array(gamma1)
        pre_mlayer_dict["mlroughness1"] = numpy.array(mlroughness1)
        pre_mlayer_dict["mlroughness2"] = numpy.array(mlroughness2)

        #These keys are not in the original pre_mlayer_dict
        pre_mlayer_dict["material1"] = material_E
        pre_mlayer_dict["material2"] = material_O
        pre_mlayer_dict["materialS"] = material_S

        pre_mlayer_dict["roughnessS"] = roughness_S

        pre_mlayer_dict["density1"] = density_E
        pre_mlayer_dict["density2"] = density_O
        pre_mlayer_dict["densityS"] = density_S

        if pre_mlayer_dict["densityS"] is None:
            try:
                pre_mlayer_dict["densityS"] = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(pre_mlayer_dict["materialS"]))
                print("Using density for substrate (%s): %f"%(pre_mlayer_dict["materialS"],
                      pre_mlayer_dict["densityS"]))
            except:
                raise Exception("Failed to load density for material: %s"%(pre_mlayer_dict["material1"]))
        if pre_mlayer_dict["density1"] is None:
            try:
                pre_mlayer_dict["density1"] = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(pre_mlayer_dict["material1"]))
                print("Using density for layer 1 (even) (%s): %f" % (pre_mlayer_dict["material1"],
                      pre_mlayer_dict["density1"]))
            except:
                raise Exception("Failed to load density for material: %s"%(pre_mlayer_dict["material1"]))
        if pre_mlayer_dict["density2"] is None:
            try:
                pre_mlayer_dict["density2"] = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(pre_mlayer_dict["material2"]))
                print("Using density for layer 2 (odd) (%s): %f" % (pre_mlayer_dict["material2"],
                      pre_mlayer_dict["density2"]))
            except:
                raise Exception("Failed to load density for material: %s"%(pre_mlayer_dict["material2"]))

        if isinstance(pre_mlayer_dict["densityS"],str):
            pre_mlayer_dict["densityS"] = float(pre_mlayer_dict["densityS"])
        if isinstance(pre_mlayer_dict["density1"],str):
            pre_mlayer_dict["density1"] = float(pre_mlayer_dict["density1"])
        if isinstance(pre_mlayer_dict["density2"], str):
            pre_mlayer_dict["density2"] = float(pre_mlayer_dict["density2"])


        # fill unused keys
        pre_mlayer_dict["energy"] = None
        pre_mlayer_dict["delta_s"] = None
        pre_mlayer_dict["beta_s"] = None
        pre_mlayer_dict["delta_e"] = None
        pre_mlayer_dict["beta_e"] = None
        pre_mlayer_dict["delta_o"] = None
        pre_mlayer_dict["beta_o"] = None
        pre_mlayer_dict["igrade"] = None
        if pre_mlayer_dict["igrade"] == 1:
            pre_mlayer_dict["fgrade"] = None
        elif pre_mlayer_dict["igrade"] == 2:  # igrade=2, coefficients
            pre_mlayer_dict["a0"] = None
            pre_mlayer_dict["a1"] = None
            pre_mlayer_dict["a2"] = None
            pre_mlayer_dict["a3"] = None

        # return
        out = MLayer()
        out.pre_mlayer_dict = pre_mlayer_dict
        out.using_pre_mlayer = False
        return out

    def scan(self, h5file="",
            energyN=51, energy1=5000.0, energy2=20000.0,
            thetaN=1, theta1=0.75, theta2=0.75, verbose=0):
        """
        this method computes the multilayer reflectivity vs photon energy and/or angle.
        The reflectivity values are for amplitude.

        Parameters
        ----------
        energyN : int, optional
            Number of points in photon energy.
        energy1 : float, optional
            The minimum photon energy in eV.
        energy2 : float, optional
            The maximum photon energy in eV.
        thetaN : int, optional
            Number of points in incident grazing angle.
        theta1 : float, optional
            The minimum value of the incident grazing angle in deg.
        theta2 : float, optional
            The minimum value of the incident grazing angle in deg.
        h5file : str, optional
            Name of the h5 file to dump calculated data. Set to "" to avoid creating file.
        verbose : int, optional
            if verbose=1, and thetaN=angleN=1 a text with results is printed .

        Returns
        -------
        tuple
            (R_S_array, R_P_array, energy_array, theta_array) numpy arrays for the sigma reflectivity (energyN, thetaN),
            pi reflectivity (energyN, thetaN), photon energy (energyN) and incident grazing angle (thetaN).

        """

        if self.pre_mlayer_dict is None:
            raise Exception("load preprocessor file before!")

        energy_array = numpy.linspace(energy1, energy2, energyN)
        theta_array = numpy.linspace(theta1, theta2, thetaN)

        ENERGY = numpy.outer(energy_array, numpy.ones_like(theta_array))
        THETA = numpy.outer(numpy.ones_like(energy_array), theta_array)

        ENERGY_flatten = ENERGY.flatten()
        THETA_flatten = THETA.flatten()

        sin_ref = numpy.sin(THETA_flatten * numpy.pi / 180)
        COS_POLE = 1.0
        k_what = 1
        R_S, R_P, _, _ = self._reflec(ENERGY_flatten, sin_ref, COS_POLE, k_what)

        R_S_array =   R_S.reshape((energyN, thetaN))
        R_P_array =   R_P.reshape((energyN, thetaN))

        if verbose:
            if ((thetaN == 1) and (energyN == 1)):
                print("------------------------------------------------------------------------")
                print("Inputs: ")
                print("   for E=", energy1, "eV: ")
                print("   energy [eV]:                       ", energy_array[0])
                print("   grazing angle [deg]:               ", theta_array[0])
                print("   wavelength [A]:                    ", (tocm / theta_array[0] * 1e8))
                print("Outputs: ")
                print("   R_S:                          ", R_P_array)
                print("   R_P:                          ", R_P_array)
                print("------------------------------------------------------------------------")

        if h5file != "":
            h5_initialize = True
            try:
                if h5_initialize:
                    h5w = H5SimpleWriter.initialize_file(h5file, creator="xoppy_multilayer.py")
                else:
                    h5w = H5SimpleWriter(h5file, None)
                h5_entry_name = "MLayer"
                h5w.create_entry(h5_entry_name, nx_default="reflectivity-s")
                if energyN == 1:
                    h5w.add_dataset(theta_array, R_S_array[0], dataset_name="reflectivity-s", entry_name=h5_entry_name,
                                    title_x="Grazing angle [deg]", title_y="Reflectivity-s")
                elif thetaN == 1:
                    h5w.add_dataset(energy_array, R_S_array[:,0], dataset_name="reflectivity-s", entry_name=h5_entry_name,
                                    title_x="Photon energy [eV]", title_y="Reflectivity-s")
                else:
                    h5w.add_image(R_S_array, energy_array, theta_array, image_name="EnergyAngleScan",
                                  entry_name=h5_entry_name,
                                  title_x="Photon Energy [eV]",
                                  title_y="Grazing Angle [deg]")

                h5w.create_entry("parameters", root_entry=h5_entry_name, nx_default=None)
                for key in self.pre_mlayer_dict.keys():
                    try:
                        h5w.add_key(key, self.pre_mlayer_dict[key], entry_name=h5_entry_name + "/parameters")
                    except:
                        pass

                print("File written to disk: %s" % h5file)
            except:
                raise Exception("ERROR writing h5 file")


        return R_S_array, R_P_array, energy_array, theta_array

    def _reflec(self, PHOT_ENER, SIN_REF, COS_POLE=1.0, K_WHAT=1.0):
        # ! C+++
        # ! C	SUBROUTINE	REFLEC
        # ! C
        # ! C	PURPOSE		To compute the local reflectivity of a mirror or
        # ! C                     multilayer. Also compute filter transmittivity.
        # ! C
        # ! C
        # ! C	ARGUMENTS	[ I ] PIN	: (x,y,z) of the intercept
        # ! C			[ I ] wnum 	: wavenumber (cm-1) => changed to PHOT_ENER array in eV
        # ! C			[ I ] sin_ref	: sine of angle from surface
        # ! C			[ I ] cos_pole  : cosine of angle of normal from pole
        # ! C			[ O ] R_P 	: p-pol reflection coefficient
        # ! C			[ O ] R_S 	: s-pol    "  "
        # ! C
        # ! C---

        nn = PHOT_ENER.size

        DELS = numpy.zeros(nn)
        BETS = numpy.zeros(nn)
        DELE = numpy.zeros(nn)
        BETE = numpy.zeros(nn)
        DELO = numpy.zeros(nn)
        BETO = numpy.zeros(nn)

        NIN = self.pre_mlayer_dict["np"]
        NPAIR = numpy.abs(self.pre_mlayer_dict["npair"])

        XLAM = tocm / PHOT_ENER * 1.0e8  # Angstrom

        gamma1 = self.pre_mlayer_dict["gamma1"]
        t_oe = self.pre_mlayer_dict["thick"]

        # gamma1 = ratio t(even)/(t(odd)+t(even))  of each layer pair
        t_e = gamma1 * t_oe
        t_o = (1.0 - gamma1) * t_oe

        mlroughness1 = self.pre_mlayer_dict["mlroughness1"]
        mlroughness2 = self.pre_mlayer_dict["mlroughness2"]

        is_monochromatic = numpy.all(PHOT_ENER == PHOT_ENER[0])

        if self.using_pre_mlayer:
            ENER = self.pre_mlayer_dict["energy"]
            wnum = 2 * numpy.pi * ENER / tocm
            QMIN = wnum[0]
            QSTEP = wnum[1] - wnum[0]

            DELTA_S = self.pre_mlayer_dict["delta_s"]
            DELTA_E = self.pre_mlayer_dict["delta_e"]
            DELTA_O = self.pre_mlayer_dict["delta_o"]
            BETA_S = self.pre_mlayer_dict["beta_s"]
            BETA_E = self.pre_mlayer_dict["beta_e"]
            BETA_O = self.pre_mlayer_dict["beta_o"]

            i_grade = self.pre_mlayer_dict["igrade"]

            # TODO graded ml
            #         ! C
            #         ! C Is the multilayer thickness graded ?
            #         ! C
            #         read    (iunit,*)   i_grade
            #         ! 0=None
            #         ! 1=spline files
            #         ! 2=quadic coefficients
            #
            #         ! spline
            #         if (i_grade.eq.1) then
            #           read  (iunit,'(a)') file_grade
            #           OPEN  (45, FILE=adjustl(FILE_GRADE), STATUS='OLD', &
            #                 FORM='UNFORMATTED', IOSTAT=iErr)
            #           ! srio added test
            #           if (iErr /= 0 ) then
            #             print *,"REFLEC: File not found: "//trim(adjustl(file_grade))
            #             print *,'Error: REFLEC: File not found. Aborted.'
            #             ! stop 'File not found. Aborted.'
            #           end if
            #
            #           READ  (45) NTX, NTY
            #           READ  (45) TX,TY
            #           !DO 205 I = 1, NTX
            #           !DO 205 J = 1, NTY
            #           DO I = 1, NTX
            #             DO J = 1, NTY
            #               READ  (45) TSPL(1,I,1,J),TSPL(1,I,2,J),    & ! spline for t
            #                          TSPL(2,I,1,J),TSPL(2,I,2,J)
            #             END DO
            #           END DO
            #
            #           READ (45) NGX, NGY
            #           READ (45) GX,GY
            #           DO I = 1, NGX
            #             DO J = 1, NGY
            #               READ (45) GSPL(1,I,1,J),GSPL(1,I,2,J),    & ! spline for gamma
            #                         GSPL(2,I,1,J),GSPL(2,I,2,J)
            #             END DO
            #           END DO
            #
            #           CLOSE (45)
            #         end if
            #
            #         if (i_grade.eq.2) then  ! quadric coefficients
            #           !
            #           ! laterally gradded multilayer
            #           !
            #
            #           ! srio@esrf.eu added cubic term (requested B Meyer, LNLS)
            #           read(iunit,*,IOSTAT=iErr) lateral_grade_constant,lateral_grade_slope, &
            #                         lateral_grade_quadratic,lateral_grade_cubic
            #
            #         end if
            #
            #         close(unit=iunit)
            #         tfilm = absor
            #         RETURN
            #     END IF
            # END IF

            ELFACTOR = numpy.log10(1.0e4 / 30.0e0) / 300.0e0

            if is_monochromatic: # just avoid the loop
                i = 0
                index1  = int(numpy.log10(PHOT_ENER[i]/ENER[0])/ELFACTOR)
                DELS[:]  = DELTA_S[index1] + (DELTA_S[index1 + 1] - DELTA_S[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
                BETS[:]  =  BETA_S[index1] + ( BETA_S[index1 + 1] -  BETA_S[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
                DELE[:]  = DELTA_E[index1] + (DELTA_E[index1 + 1] - DELTA_E[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
                BETE[:]  =  BETA_E[index1] + ( BETA_E[index1 + 1] -  BETA_E[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
                DELO[:]  = DELTA_O[index1] + (DELTA_O[index1 + 1] - DELTA_O[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
                BETO[:]  =  BETA_O[index1] + ( BETA_O[index1 + 1] -  BETA_O[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
            else:
                for i in range(nn):
                    index1  = int(numpy.log10(PHOT_ENER[i]/ENER[0])/ELFACTOR)
                    DELS[i]  = DELTA_S[index1] + (DELTA_S[index1+1] - DELTA_S[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
                    BETS[i]  =  BETA_S[index1] + ( BETA_S[index1+1] -  BETA_S[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
                    DELE[i]  = DELTA_E[index1] + (DELTA_E[index1+1] - DELTA_E[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
                    BETE[i]  =  BETA_E[index1] + ( BETA_E[index1+1] -  BETA_E[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
                    DELO[i]  = DELTA_O[index1] + (DELTA_O[index1+1] - DELTA_O[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
                    BETO[i]  =  BETA_O[index1] + ( BETA_O[index1+1] -  BETA_O[index1]) * (PHOT_ENER[i] - ENER[index1]) / (ENER[index1 + 1] - ENER[index1])
        else: # not using preprocessor, using xraylib
            import xraylib
            if is_monochromatic:
                DELS[:]  = 1.0 - xraylib.Refractive_Index_Re(self.pre_mlayer_dict["materialS"], 1e-3 * PHOT_ENER[0], self.pre_mlayer_dict["densityS"])
                BETS[:]  =       xraylib.Refractive_Index_Im(self.pre_mlayer_dict["materialS"], 1e-3 * PHOT_ENER[0], self.pre_mlayer_dict["densityS"])
                DELE[:]  = 1.0 - xraylib.Refractive_Index_Re(self.pre_mlayer_dict["material1"], 1e-3 * PHOT_ENER[0], self.pre_mlayer_dict["density1"])
                BETE[:]  =       xraylib.Refractive_Index_Im(self.pre_mlayer_dict["material1"], 1e-3 * PHOT_ENER[0], self.pre_mlayer_dict["density1"])
                DELO[:]  = 1.0 - xraylib.Refractive_Index_Re(self.pre_mlayer_dict["material2"], 1e-3 * PHOT_ENER[0], self.pre_mlayer_dict["density2"])
                BETO[:]  =       xraylib.Refractive_Index_Im(self.pre_mlayer_dict["material2"], 1e-3 * PHOT_ENER[0], self.pre_mlayer_dict["density2"])
            else:
                for i in range(nn):
                    DELS[i]  = 1.0 - xraylib.Refractive_Index_Re(self.pre_mlayer_dict["materialS"], 1e-3 * PHOT_ENER[i], self.pre_mlayer_dict["densityS"])
                    BETS[i]  =       xraylib.Refractive_Index_Im(self.pre_mlayer_dict["materialS"], 1e-3 * PHOT_ENER[i], self.pre_mlayer_dict["densityS"])
                    DELE[i]  = 1.0 - xraylib.Refractive_Index_Re(self.pre_mlayer_dict["material1"], 1e-3 * PHOT_ENER[i], self.pre_mlayer_dict["density1"])
                    BETE[i]  =       xraylib.Refractive_Index_Im(self.pre_mlayer_dict["material1"], 1e-3 * PHOT_ENER[i], self.pre_mlayer_dict["density1"])
                    DELO[i]  = 1.0 - xraylib.Refractive_Index_Re(self.pre_mlayer_dict["material2"], 1e-3 * PHOT_ENER[i], self.pre_mlayer_dict["density2"])
                    BETO[i]  =       xraylib.Refractive_Index_Im(self.pre_mlayer_dict["material2"], 1e-3 * PHOT_ENER[i], self.pre_mlayer_dict["density2"])



        TFACT = 1.0
        GFACT = 1.0

        #TODO graded multilayers
        #         IF (I_GRADE.EQ.1) THEN
        #             XIN = PIN(1)
        #             YIN = PIN(2)
        #             CALL DBCEVL (TX,NTX,TY,NTY,TSPL,i101,XIN,YIN,PDS,IER)
        #             IF (IER.NE.0) THEN
        #               CALL      MSSG ('REFLEC','Spline error # ',IER)
        #               RETURN
        #             END IF
        #             TFACT = PDS(1)
        #             ! C
        #             CALL DBCEVL (GX,NGX,GY,NGY,GSPL,i101,XIN,YIN,PDS,IER)
        #             IF (IER.NE.0) THEN
        #               CALL MSSG ('REFLEC','Spline error # ',IER)
        #               RETURN
        #             END IF
        #             GFACT = PDS(1)
        #         ELSE IF (I_GRADE.EQ.2) THEN
        #             TFACT = lateral_grade_constant+ &
        #                     lateral_grade_slope*pin(2) + &
        #                     lateral_grade_quadratic*pin(2)*pin(2) + &
        #                     lateral_grade_cubic*pin(2)*pin(2)*pin(2)
        #         ELSE
        #         END IF

        R_S, R_P, PHASES, PHASEP = self._fresnel(TFACT, GFACT, NPAIR, SIN_REF, COS_POLE, XLAM,
                                             DELO, DELE, DELS, BETO, BETE, BETS, t_o, t_e, mlroughness1, mlroughness2)

        return R_S, R_P, PHASES, PHASEP


    def reflectivity(self, grazing_angle_deg, photon_energy_ev):
        """
        Computes the reflectivity (in amplitude) as a function of the incident angle and photon energy.

        Parameters
        ----------
        grazing_angle_deg : numpy array
            The array with the incident gracing angle in deg.
        photon_energy_ev : numpy array
            The array with the photon energy in eV.

        Returns
        -------
        tuple
            (R_S_array, R_P_array, phase_S, phase_P) arrays with reflectivities (in amplitude) and phases.
        """
        if self.pre_mlayer_dict is None:
            raise Exception("load preprocessor file before!")

        R_S_array = numpy.zeros_like(grazing_angle_deg)
        R_P_array = numpy.zeros_like(grazing_angle_deg)
        k_what = 1

        vectorized = 1
        if vectorized == 0:
            for i in range(grazing_angle_deg.size):
                theta = grazing_angle_deg[i]
                energy = photon_energy_ev[i]
                sin_ref = numpy.sin(theta * numpy.pi / 180)
                wnum = 2 * numpy.pi * energy / tocm
                COS_POLE = 1.0
                R_S, R_P, tmp, phases, phasep = self._reflec(wnum, sin_ref, COS_POLE, k_what)
                R_S_array[i] = R_S
                R_P_array[i] = R_P
        else:
            sin_ref = numpy.sin(grazing_angle_deg * numpy.pi / 180)
            COS_POLE = 1.0
            R_S_array, R_P_array, phase_S, phase_P = self._reflec(photon_energy_ev, sin_ref, COS_POLE, k_what)

        return R_S_array, R_P_array, phase_S, phase_P

    @classmethod
    def _fresnel(cls, TFACT, GFACT, NPAIR, SIN_REF, COS_POLE, XLAM,
                delo, dele, dels, beto, bete, bets, t_o, t_e, mlroughness1, mlroughness2):

        # !C------------------------------------------------------------------------------
        # !C  subroutine FRESNEL
        # !C------------------------------------------------------------------------------
        # !c  compute x-ray/u.v. reflection efficiency of multilayers
        # !c
        # !c  inputs:
        # !c         tfact   : used for ML with graded thickness (thickness coeff)
        # !c         gfact   : used for ML with graded thickness (gamma coeff)
        # !c         n       : number of bilayers
        # !c         sin_ref :  sin of angle of incidence (grazing??)
        # !c         cos_pole:  cos of angle of between normal and pole??
        # !c
        # !c        delo,dele,dels = parameter delta odd, even, substrate respectively
        # !c        belo,bele,bels = parametro beta odd, even, substrate respectively
        # !c            1.0 - delo - i*beto = complex refractive index (odd)
        # !c            1.0 - dele - i*bete = complex refractive index (even)
        # !c        t_o = thickness of odd layers (a)
        # !c        t_e = thickness of even layers (a)
        # !c outputs:
        # !c         ans = S polarization  reflectivity
        # !c         anp = P polarization  reflectivity
        # !c         phaseS = change of phase S
        # !c         phaseP = change di phase P

        # !c
        # !c----------------------------------------------------------------------------
        # !C
        # !C
        # !C               vacuum
        # !C    |------------------------------|  \
        # !C    |          odd (n)             |  |
        # !C    |------------------------------|  | BILAYER # n
        # !C    |          even (n)            |  |
        # !C    |------------------------------|  /
        # !C    |          .                   |
        # !C    |          .                   |
        # !C    |          .                   |
        # !C    |------------------------------|  \
        # !C    |          odd (1)             |  |
        # !C    |------------------------------|  | BILAYER # 1
        # !C    |          even (1)            |  |
        # !C    |------------------------------|  /
        # !C    |                              |
        # !C    |///////// substrate //////////|
        # !C    |                              |
        # !C
        # !c----------------------------------------------------------------------------
        # !c----------------------------------------------------------------------------

        ci = 0.+1.0j

        # ! (refraction index "odd,even,substrate")**2

        ro2 = (1.0 - delo - ci * beto)**2
        re2 = (1.0 - dele - ci * bete)**2
        rs2 = (1.0 - dels - ci * bets)**2

        # ! angles
        SIN_REF2 = SIN_REF**2
        COS_REF2 = 1.0 - SIN_REF2

        fo = ro2 - COS_REF2
        fe = re2 - COS_REF2
        refv = SIN_REF2
        xmfv = 0.0

        fv = refv + ci * xmfv
        fs = rs2 - COS_REF2

        fo = numpy.sqrt(fo) # complex!!
        fe = numpy.sqrt(fe) # complex!!
        fv = numpy.sqrt(fv) # complex!!
        fs = numpy.sqrt(fs) # complex!!


        # ! Fresnel formula "S" (in function of incidence angle and critical angle)
        ffe = (fe-fo)/(fe+fo)
        ffo = -ffe
        ffv = (fv-fo)/(fv+fo)
        ffs = (fe-fs)/(fe+fs)
        # ! Fresnel formula "P" (in function of incidence angle and critical angle)
        ffep = (fe/re2-fo/ro2)/(fe/re2+fo/ro2)
        ffop = -ffep
        ffvp = (fv-fo/ro2)/(fv+fo/ro2)
        ffsp = (fe/re2-fs/rs2)/(fe/re2+fs/rs2)

        if NPAIR == 0: # now there is only substrate and vacuum
            fe = fv
            fo = fs
            # ! Fresnel formula "S" (in function of incidence angle and critical angle)
            ffe = (fe - fo) / (fe + fo)
            ffo = -ffe
            ffv = (fv - fo) / (fv + fo)
            ffs = (fe - fs) / (fe + fs)
            # ! Fresnel formula "P" (in function of incidence angle and critical angle)
            ffep = (fe / re2 - fo / ro2) / (fe / re2 + fo / ro2)
            ffop = -ffep
            ffvp = (fv - fo / ro2) / (fv + fo / ro2)
            ffsp = (fe / re2 - fs / rs2) / (fe / re2 + fs / rs2)





        # ! another way
        # ! ro=(1.0D0-delo-ci*beto)
        # ! re=(1.0D0-dele-ci*bete)
        # ! rs=(1.0D0-dels-ci*bets)
        # !
        # ! !
        # ! cos_ref = sqrt(1.0D0-sin_ref**2)   ! in vacuum
        # ! !!! snell (top to bottom propagation)
        # ! cos_o = (1.0D0/ro)*cos_ref  ! in odd medium
        # ! cos_e = (ro/re)*cos_o       ! in even medium
        # ! cos_s = (re/rs)*cos_e       ! in substrate medium
        # !
        # ! sin_o = mysqrt(1.0d0 - cos_o**2)
        # ! sin_e = mysqrt(1.0d0 - cos_e**2)
        # ! sin_s = mysqrt(1.0d0 - cos_s**2)
        # !
        # ! ! even->odd interface
        # ! ffe =  (re*sin_e - ro*sin_o)/ &  ! e->o
        # !        (re*sin_e + ro*sin_o)
        # ! ffo=-ffe
        # ! ffv = (sin_ref - ro*sin_o )/ &   ! v->o
        # !       (sin_ref + ro*sin_o )
        # ! ! even->substrate interface
        # ! ffs = (re*sin_e - rs*sin_s )/ &  ! e->s
        # !       (re*sin_e + rs*sin_s)
        # !
        # ! ! p-polarization
        # ! ffep = (re*sin_o - ro*sin_e)/ &  ! e->o
        # !        (re*sin_o + ro*sin_e)
        # !
        # ! ffop=-ffep
        # ! ffvp = (sin_o - ro*sin_ref )/&   !v->o
        # !        (sin_o + ro*sin_ref )
        # !
        # ! ffsp = (re*sin_s - rs*sin_e )/ &  ! e->s
        # !        (re*sin_s + rs*sin_e)

        r = 0.0 + 0.0j
        rp = 0.0 + 0.0j
        prefact = (8.*(numpy.pi**2.)) / (XLAM**2)


        # !c Nevot-Croce roughness
        # !c DO NOT include refraction index in the roughness formula

        sigma_s2 = 0.0 # ! sigma_s**2.0 !roughn. substrate
        sigma_v2 = 0.0 # ! sigma_v**2.0!roughn. vacuum


        # ! loop over the bilayers
        # ! remember that "even" is the bottom sublayer
        for j in range(NPAIR): #   =1,n   ! n is the number of bilayers
            # ! C
            # ! C compute the thickness for the odd and even material :
            # ! C
            ao = -ci * (numpy.pi * fo * t_o[j] * COS_POLE / XLAM)
            ae = -ci * (numpy.pi * fe * t_e[j] * COS_POLE / XLAM)
            ao = numpy.exp(ao)
            ae = numpy.exp(ae)

            if j != 0:
                sigma_e2 = mlroughness1[j]**2.0 #!roughn. even layer
                arg_e = fo * fe * sigma_e2 / (numpy.sqrt(ro2) * numpy.sqrt(re2))
                fnevot_e = numpy.exp(-prefact * arg_e)
                r = (ae**4) * (r + ffe * fnevot_e) / (r * ffe * fnevot_e + 1.0)
                rp = (ae**4) * (rp + ffep * fnevot_e) / (rp * ffep * fnevot_e + 1.0)
            else:
                # ! layer on top of substrate
                arg_s = fe * fs * sigma_s2 / (numpy.sqrt(re2) * numpy.sqrt(rs2))
                fnevot_s = numpy.exp(-prefact * arg_s)
                r = (ae**4.0) * (r + ffs * fnevot_s) / (r * ffs * fnevot_s + 1.0)
                rp = (ae**4.0) * (rp + ffsp * fnevot_s) / (rp * ffsp * fnevot_s + 1.0)

            # ! odd layer (top sublayer)
            sigma_o2 = mlroughness2[j]**2.0 #!roughn. odd layer
            arg_o = fo * fe * sigma_o2 / (numpy.sqrt(ro2) * numpy.sqrt(re2))
            fnevot_o = numpy.exp(-prefact * arg_o)
            r = (ao**4.0) * (r + ffo * fnevot_o) / (r * ffo * fnevot_o + 1.0)
            rp = (ao**4.0) * (rp + ffop * fnevot_o) / (rp * ffop * fnevot_o + 1.0)

        # !
        # ! vacuum interface
        # !
        arg_v = fo * fv * sigma_v2 / numpy.sqrt(ro2)
        fnevot_v = numpy.exp(-prefact * arg_v)
        r = (r + ffv * fnevot_v) / (r * ffv * fnevot_v + 1.0)
        rp = (rp + ffvp * fnevot_v) / (rp * ffvp * fnevot_v + 1.0)
        #
        # !
        # ! calculate phases
        # !
        PHASES = numpy.arctan2(r.imag,r.real)
        ans = numpy.abs(r)

        PHASEP = numpy.arctan2(rp.imag,rp.real)
        anp = numpy.abs(rp)

        # print("\>>>>>> fresnel inputs\n TFACT, GFACT, NPAIR, SIN_REF, COS_POLE, XLAM, delo, dele, dels, beto, bete, bets, t_o, t_e, mlroughness1, mlroughness2 = ",
        #       TFACT, GFACT, NPAIR, SIN_REF, COS_POLE, XLAM,
        #       delo, dele, dels, beto, bete, bets, t_o[0], t_e[0], mlroughness1[0], mlroughness2[0])
        #
        # print(">>>>>> fresnel output\n ans, anp, PHASES, PHASEP = ", ans, anp, PHASES, PHASEP)
        return ans, anp, PHASES, PHASEP

if __name__ == "__main__":

    from srxraylib.plot.gol import plot

    if 0:
        a = MLayer.pre_mlayer(
            interactive=False,
            FILE="pre_mlayer.dat",
            E_MIN=110.0, E_MAX=500.0,
            O_DENSITY=7.19, O_MATERIAL="Cr",  # odd: closer to vacuum
            E_DENSITY=3.00, E_MATERIAL="Sc",  # even: closer to substrate
            S_DENSITY=2.33, S_MATERIAL="Si",  # substrate
            GRADE_DEPTH=0,
            N_PAIRS=50,
            THICKNESS=22.0,
            GAMMA=10.0/22.0,  #  gamma ratio  =  t(even) / (t(odd) + t(even))")
            ROUGHNESS_EVEN=0.0,
            ROUGHNESS_ODD=0.0,
            FILE_DEPTH="myfile_depth.dat",
            GRADE_SURFACE=0,
            FILE_SHADOW="mlayer1.sha",
            FILE_THICKNESS="mythick.dat",
            FILE_GAMMA="mygamma.dat",
            AA0=1.0,AA1=0.0,AA2=0.0,AA3=0.0)

        b = MLayer()
        b.read_preprocessor_file("pre_mlayer.dat")


        #
        # energy scan
        #
        rs, rp, e, t = a.scan(h5file="",
                energyN=100,energy1=300.0,energy2=500.0,
                thetaN=1,theta1=45.0,theta2=45.0)

        print(rs.shape,rp.shape,e.shape,t.shape)

        plot(e,rs[:,0],xtitle="Photon energy [eV]",ytitle="Reflectivity")

        #
        # theta scan
        #
        rs, rp, e, t = a.scan(h5file="",
                energyN=1,energy1=400.0,energy2=401.0,
                thetaN=1000,theta1=40.0,theta2=50.0)

        print(rs.shape,rp.shape,e.shape,t.shape)

        plot(t,rs[0],xtitle="angle [deg]",ytitle="Reflectivity",ylog=False)

        #
        # single point
        #
        a.scan(h5file="",
                energyN=1,energy1=398.0,thetaN=1,theta1=45.0)


    if 0: # prepate optical element: MLayer.reflectivity()
        if 1:
            b = MLayer.pre_mlayer(
                interactive=False,
                FILE="mlayer1.dat",
                E_MIN=30000.0, E_MAX=100000.0,
                O_DENSITY=19.3, O_MATERIAL="W",  # odd: closer to vacuum
                E_DENSITY=2.37, E_MATERIAL="B",  # even: closer to substrate
                S_DENSITY=2.33, S_MATERIAL="Si",  # substrate
                GRADE_DEPTH=0,
                N_PAIRS=350,
                THICKNESS=30.19,
                GAMMA=0.666,  #  gamma ratio  =  t(even) / (t(odd) + t(even))")
                ROUGHNESS_EVEN=0.0,
                ROUGHNESS_ODD=0.0,
                FILE_DEPTH="myfile_depth.dat",
                GRADE_SURFACE=0,
                FILE_SHADOW="tmp.sha",
                FILE_THICKNESS="tmp.thi",
                FILE_GAMMA="tmp.gam",
                AA0=1.0, AA1=0.0, AA2=0.0, AA3=0.0)

            # b = MLayer()
            # # b.read_preprocessor_file("mlayer1.dat") # "/home/srio/Oasys/mlayer.dat")
            # b.read_preprocessor_file("/home/srio/Oasys/mlayer.dat")

        else:

            b = MLayer.initialize_from_bilayer_stack(
                                              material_S="Si", density_S=None, roughness_S=0.0,
                                              material_E="B", density_E=None, roughness_E=0.0,
                                              material_O="W", density_O=None, roughness_O=0.0,
                                              bilayer_pairs=350,
                                              bilayer_thickness=30.19,
                                              bilayer_gamma=0.666,
                                              )


        print(">>>> using preprocessor data: ", b.using_pre_mlayer)

        #
        # energy scan
        #
        rs, rp, e, t = b.scan(h5file="",
                              energyN=501, energy1=41000, energy2=49000,
                              thetaN=1, theta1=0.26356, theta2=0.26356)
        plot(e, rs[:, 0], xtitle="Photon energy [eV]", ytitle="Reflectivity", title="ENERGY scan", )
        energy = e.copy()
        #
        # theta scan
        #
        rs, rp, e, t = b.scan(h5file="",
                              energyN=1, energy1=45000, energy2=45000,
                              thetaN=1000, theta1=0.001, theta2=1.100)

        print(rs.shape, rp.shape, e.shape, t.shape)

        plot(t, rs[0], xtitle="angle [deg]", ytitle="Reflectivity", title="THETA scan", ylog=False)

        #
        # single point
        #
        b.scan(h5file="",
               energyN=1, energy1=45000, thetaN=1, theta1=0.26356)

        #
        # MLayer.reflectivity() scan (for shadow4) scan
        #
        rs1, rp1, phase_s1, phase_p1= b.reflectivity(numpy.ones_like(energy) * 0.26356, energy)
        plot(energy, rs1, xtitle="Photon energy [eV]", ytitle="Reflectivity", title="MIXED scan", )

    if 0: # compressed format

        b = MLayer.initialize_from_bilayer_stack_in_compressed_format(
        structure='[Pd/B4C]x150+Si',
        density_O=None,  roughness_O=0.0,
        density_E=2.52,  roughness_E=0.0,
        density_S=None,  roughness_S=0.0,
        bilayer_thickness=25.14,
        bilayer_gamma=0.5,
        )

        # b = MLayer()
        # # b.read_preprocessor_file("mlayer1.dat") # "/home/srio/Oasys/mlayer.dat")
        # b.read_preprocessor_file("/home/srio/Oasys/mlayer.dat")

        print(">>>> using preprocessor data: ", b.using_pre_mlayer)

        #
        # theta scan
        #

        npoints = 1000

        # compare with IMD
        a = numpy.loadtxt('/home/srio/Oasys/Pd_B4C.txt', skiprows=17)

        rs, rp, e, t = b.scan(h5file="",
                              energyN=1, energy1=8047.80, energy2=8047.8,
                              thetaN=npoints, theta1=0.001, theta2=8.0)
        plot(
             t, rs[0] ** 2,
             a[:, 0], a[:, 1],
             xtitle="angle [deg]", ytitle="Reflectivity", title="THETA scan", ylog=1, yrange=[1e-8, 1],
             legend=['S4 s-pol','IMD s-pol'])


        rs1, rp1, phase_s1, phase_p1 = b.reflectivity(numpy.linspace(0,8,npoints), numpy.zeros(npoints)+8047.80)
        plot(
             numpy.linspace(0,8,npoints), rs1 ** 2,
             a[:, 0], a[:, 1],
             xtitle="angle [deg]", ytitle="Reflectivity", title="THETA scan", ylog=1, yrange=[1e-8, 1],
             legend=['S4 reflectivity s-pol','IMD s-pol'])

        # Note that the phases in IMD are reversed!
        plot(
             numpy.linspace(0,8,npoints), phase_s1,
             a[:, 0], numpy.radians(-a[:, 3]),
             xtitle="angle [deg]", ytitle="Phase S [rad]", title="THETA scan", ylog=0,
             legend=['S4 reflectivity s-pol','IMD s-pol'])

        plot(
             numpy.linspace(0,8,npoints), phase_p1,
             a[:, 0], numpy.radians(-a[:, 4]),
             xtitle="angle [deg]", ytitle="Phase P [rad]", title="THETA scan", ylog=0,
             legend=['S4 reflectivity s-pol','IMD s-pol'])



    if 1: #
        b = MLayer()
        b.read_preprocessor_file("/home/srio/Oasys/mlayer.dat")
        print(">>>> using preprocessor data: ", b.using_pre_mlayer)

        #
        # energy scan
        #
        npoints = 1000
        energy = numpy.linspace(41000, 49000, npoints)
        theta = numpy.zeros(npoints) + 0.26356
        rs1, rp1, phase_s1, phase_p1 = b.reflectivity(theta, energy)

        # rs, rp, e, t = b.scan(h5file="",
        #                       energyN=501, energy1=41000, energy2=49000,
        #                       thetaN=1, theta1=0.26356, theta2=0.26356)
        # plot(e, rs[:, 0], xtitle="Photon energy [eV]", ytitle="Reflectivity", title="ENERGY scan", )

        plot(energy, rs1**2, xtitle="Photon energy [eV]", ytitle="Reflectivity", title="ENERGY scan", )
