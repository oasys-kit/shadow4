
import numpy
import scipy.constants as codata

tocm = codata.h*codata.c/codata.e*1e2 # 12398.419739640718e-8

class PreRefl(object):

    def __init__(self):

        self.prerefl_dict = None

    def read_preprocessor_file(self,filename):

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


    def get_refraction_index(self,energy1,verbose=False):

        wnum = 2*numpy.pi * energy1 / tocm

        QSTEP = self.prerefl_dict["QSTEP"]
        QMIN = self.prerefl_dict["QMIN"]

        index1 = (wnum - QMIN) / QSTEP

        if index1 > self.prerefl_dict["NREFL"]:
            raise Exception("Error: Photon energy above upper limit.")

        WNUM0  =   QSTEP * int(index1) + QMIN
        DEL_X  =   wnum - WNUM0
        DEL_X  =   DEL_X / QSTEP

        index1 = int(index1)

        myALFA  =   self.prerefl_dict["ZF1"][index1] + (self.prerefl_dict["ZF1"][index1+1]-self.prerefl_dict["ZF1"][index1]) * DEL_X
        myGAMMA  =  self.prerefl_dict["ZF2"][index1] + (self.prerefl_dict["ZF2"][index1+1]-self.prerefl_dict["ZF2"][index1]) * DEL_X

        refraction_index = (1.0 - myALFA / 2) + (myGAMMA / 2)*1j

        if verbose:
            print("   prerefl_test: calculates refraction index for a given energy")
            print("                 using a file created by the prerefl preprocessor.")
            print("   \n")

            print("------------------------------------------------------------------------" )
            print("Inputs: " )
            print("   wavenumber (2 pi/lambda) [cm^-1]:  ",wnum )
            # wnum = 2*numpy.pi * energy1 / tocm
            print("   energy [eV]:                       ",wnum * tocm / (2*numpy.pi))
            print("   wavelength [A]:                    ",(1.0/wnum) * 2 * numpy.pi*1e8 )
            print("Outputs: " )
            print("   refraction index = (1-delta) + i*beta : " )
            print("   delta:                          ",1.0-refraction_index.real )#1.0-rr_ind )
            print("   beta:                           ",refraction_index.imag) #rr_attenuation / (2*wnum) )
            print("   real(n):                        ",refraction_index.real )
            print("   attenuation coef [cm^-1]:       ",2*refraction_index.imag*wnum )
            print("------------------------------------------------------------------------" )


        return refraction_index

if __name__ == "__main__":
    a = PreRefl()

    prerefl_file = "Be5_30.dat"
    a.read_preprocessor_file(prerefl_file)

    refraction_index = a.get_refraction_index(10000.0,verbose=True)

