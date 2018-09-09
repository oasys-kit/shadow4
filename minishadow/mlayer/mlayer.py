__author__ = 'srio'

#
import numpy
import scipy.constants as codata
# import math
#



class MLayer(object):

    def __init__(self):

        self.pre_mlayer_dict = None



    def read_preprocessor_file(self,filename):

        out_dict = {}


        fp = open(filename) # Open file on read mode
        lines = fp.read().split("\n") # Create a list containing all lines
        fp.close() # Close file


        index_pointer = 0
        np = int(lines[index_pointer])
        out_dict['np'] = np


        # print(">>>>",np)
        # f.write("%i \n" % np)
        energy = numpy.zeros(np)

        index_pointer += 1
        mylist = lines[index_pointer].split(" ")
        for i in range(np):
            energy[i] = float(mylist[i])

        out_dict["energy"] = energy

        # print(energy,energy.shape)
        #
        delta_s = numpy.zeros_like(energy)
        beta_s = numpy.zeros_like(energy)

        delta_e = numpy.zeros_like(energy)
        beta_e = numpy.zeros_like(energy)

        delta_o = numpy.zeros_like(energy)
        beta_o = numpy.zeros_like(energy)

        for i in range(np):
            index_pointer += 1
            mylist = lines[index_pointer].strip().split(" ")
            delta_s[i] = float(mylist[0])
            beta_s[i] = float(mylist[0])

        for i in range(np):
            index_pointer += 1
            mylist = lines[index_pointer].strip().split(" ")
            delta_e[i] = float(mylist[0])
            beta_e[i] = float(mylist[0])

        for i in range(np):
            index_pointer += 1
            mylist = lines[index_pointer].strip().split(" ")
            delta_o[i] = float(mylist[0])
            beta_o[i] = float(mylist[0])


        out_dict["delta_s"] = delta_s
        out_dict["beta_s"] = beta_s

        out_dict["delta_e"] = delta_e
        out_dict["beta_e"] = beta_e

        out_dict["delta_o"] = delta_o
        out_dict["beta_o"] = beta_o

        # print(delta_s,delta_s.shape)
        # print(beta_s,beta_s.shape)
        #
        # for i in range(np):
        #     energy = 30e0*math.pow(10,elfactor*(istart+i-1)) *1e-3 # in keV!!
        #     delta = 1e0-xraylib.Refractive_Index_Re(matEven,energy,denEven)
        #     beta = xraylib.Refractive_Index_Im(matEven,energy,denEven)
        #     f.write( ("%26.17e  "*2+"\n") % tuple([delta,beta]) )
        #
        # for i in range(np):
        #     energy = 30e0*math.pow(10,elfactor*(istart+i-1)) *1e-3 # in keV!!
        #     delta = 1e0-xraylib.Refractive_Index_Re(matOdd,energy,denOdd)
        #     beta = xraylib.Refractive_Index_Im(matOdd,energy,denOdd)
        #     f.write( ("%26.17e "*2+"\n") % tuple([delta,beta]) )
        #
        #
        # #! srio@esrf.eu 2012-06-07 Nevot-Croce ML roughness model implemented.
        # #! By convention, starting from the version that includes ML roughness
        # #! we set NPAR negative, in order to assure compatibility with old
        # #! versions. If NPAR<0, roughness data are read, if NPAR>0 no roughness.
        # f.write("%i \n" % -npair)

        index_pointer += 1
        npair = int(lines[index_pointer])

        out_dict["npair"] = npair

        # print(">>> npair",npair)

        thick = numpy.zeros(numpy.abs(npair))
        gamma1 = numpy.zeros_like(thick)
        mlroughness1 = numpy.zeros_like(thick)
        mlroughness2 = numpy.zeros_like(thick)

        for i in range(-npair):
            index_pointer += 1
            mylist = lines[index_pointer].strip().split("   ")
            thick[i] = float(mylist[0])
            gamma1[i] = float(mylist[1])
            mlroughness1[i] = float(mylist[2])
            mlroughness2[i] = float(mylist[3])

        out_dict["thick"] = thick
        out_dict["gamma1"] = gamma1
        out_dict["mlroughness1"] = mlroughness1
        out_dict["mlroughness2"] = mlroughness2

        #
        #
        # for i in range(npair):
        #     f.write( ("%26.17e "*4+"\n") % tuple([thick[i],gamma1[i],mlroughness1[i],mlroughness2[i]]) )
        #

        index_pointer += 1
        igrade = int(lines[index_pointer])

        # print(">>> igrade",igrade)

        out_dict["igrade"] = igrade


        # f.write("%i \n" % igrade)
        if igrade == 1:
            # f.write("%s \n" % fgrade)
            index_pointer += 1
            fgrade = int(lines[index_pointer])

            # print(">>> fgrade",fgrade)
            out_dict["fgrade"] = fgrade

        elif igrade == 2:  # igrade=2,
        # coefficients
        #     f.write("%f  %f  %f  %f\n"%(a0,a1,a2,a3))
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



    @classmethod
    def pre_mlayer(self, interactive=False, FILE="pre_mlayer.dat",E_MIN=5000.0,E_MAX=20000.0,
                   S_DENSITY=2.33,S_MATERIAL="Si",
                   E_DENSITY=2.40,E_MATERIAL="B4C",
                   O_DENSITY=9.40,O_MATERIAL="Ru",
                   GRADE_DEPTH=0,N_PAIRS=70,THICKNESS=33.1,GAMMA=0.483,
                   ROUGHNESS_EVEN=3.3,ROUGHNESS_ODD=3.1,
                   FILE_DEPTH="myfile_depth.dat",GRADE_SURFACE=0,FILE_SHADOW="mlayer1.sha",
                   FILE_THICKNESS="mythick.dat",FILE_GAMMA="mygamma.dat",AA0=1.0,AA1=0.0,AA2=0.0,AA3=0.0):

        """
         SHADOW preprocessor for multilayers - python+xraylib version
         -"""

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
            #--

        ###--------------------------------------------------------------------------------------

        elfactor = math.log10(1.0e4/30.0)/300.0
        istart = int(math.log10(estart/30.0e0)/elfactor + 1)
        ifinal = int(math.log10(efinal/30.0e0)/elfactor + 2)
        np = int(ifinal - istart) + 1

        f = open(fileout, 'wt')
        f.write("%i \n" % np)
        for i in range(np):
            energy = 30e0*math.pow(10,elfactor*(istart+i-1))
            f.write("%e " % energy)
        f.write( "\n")

        for i in range(np):
            energy = 30e0*math.pow(10,elfactor*(istart+i-1)) *1e-3 # in keV!!
            delta = 1e0-xraylib.Refractive_Index_Re(matSubstrate,energy,denSubstrate)
            beta = xraylib.Refractive_Index_Im(matSubstrate,energy,denSubstrate)
            f.write( ("%26.17e "*2+"\n") % tuple([delta,beta]) )

        for i in range(np):
            energy = 30e0*math.pow(10,elfactor*(istart+i-1)) *1e-3 # in keV!!
            delta = 1e0-xraylib.Refractive_Index_Re(matEven,energy,denEven)
            beta = xraylib.Refractive_Index_Im(matEven,energy,denEven)
            f.write( ("%26.17e  "*2+"\n") % tuple([delta,beta]) )

        for i in range(np):
            energy = 30e0*math.pow(10,elfactor*(istart+i-1)) *1e-3 # in keV!!
            delta = 1e0-xraylib.Refractive_Index_Re(matOdd,energy,denOdd)
            beta = xraylib.Refractive_Index_Im(matOdd,energy,denOdd)
            f.write( ("%26.17e "*2+"\n") % tuple([delta,beta]) )


        #! srio@esrf.eu 2012-06-07 Nevot-Croce ML roughness model implemented.
        #! By convention, starting from the version that includes ML roughness
        #! we set NPAR negative, in order to assure compatibility with old
        #! versions. If NPAR<0, roughness data are read, if NPAR>0 no roughness.
        f.write("%i \n" % -npair)


        for i in range(npair):
            f.write( ("%26.17e "*4+"\n") % tuple([thick[i],gamma1[i],mlroughness1[i],mlroughness2[i]]) )

        f.write("%i \n" % igrade)
        if igrade == 1:
            f.write("%s \n" % fgrade)
        elif igrade == 2:  # igrade=2, coefficients
            f.write("%f  %f  %f  %f\n"%(a0,a1,a2,a3))
        ###--------------------------------------------------------------------------------------
        f.close()
        print("File written to disk: %s" % fileout)
        return None



if __name__ == "__main__":

    # a = MLayer.pre_mlayer()

    filename = "/Users/srio/Oasys/mlayer.dat"

    a = MLayer()
    a.read_preprocessor_file(filename)
    for k in a.pre_mlayer_dict:
        print(k)
