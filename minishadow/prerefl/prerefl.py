"""

python version of the mirror reflectivity code and refractive index calculations in shadow

TODO: vectorization

"""

import numpy
import scipy.constants as codata

tocm = codata.h * codata.c / codata.e * 1e2 # 12398.419739640718e-8

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

    def preprocessor_info(self,verbose=False):

        print("\n========================================")
        for k in self.prerefl_dict.keys():
            try:
                print(k,self.prerefl_dict[k][0])
            except:
                print(k,self.prerefl_dict[k])

        print("QMIN: %f  EMIN: %f "%(self.prerefl_dict["QMIN"], self.prerefl_dict["QMIN"] * tocm / (2*numpy.pi)))
        print("QMAX: %f  EMAX: %f "%(self.prerefl_dict["QMAX"], self.prerefl_dict["QMAX"] * tocm / (2*numpy.pi)))
        print("========================================")

    def get_refraction_index(self,energy1,verbose=False):

        wnum = 2*numpy.pi * energy1 / tocm

        QSTEP = self.prerefl_dict["QSTEP"]
        QMIN = self.prerefl_dict["QMIN"]

        index1 = (wnum - QMIN) / QSTEP

        if index1 > self.prerefl_dict["NREFL"]:
            raise Exception("Error: Photon energy above upper limit.")

        WNUM0 = QSTEP * int(index1) + QMIN
        DEL_X = wnum - WNUM0
        DEL_X = DEL_X / QSTEP

        index1 = int(index1)



        ALFA = self.prerefl_dict["ZF1"][index1] + (self.prerefl_dict["ZF1"][index1+1]-self.prerefl_dict["ZF1"][index1]) * DEL_X
        GAMMA = self.prerefl_dict["ZF2"][index1] + (self.prerefl_dict["ZF2"][index1+1]-self.prerefl_dict["ZF2"][index1]) * DEL_X

        refraction_index = (1.0 - ALFA / 2) + (GAMMA / 2)*1j

        if verbose:
            print("   prerefl_test: calculates refraction index for a given energy")
            print("                 using a file created by the prerefl preprocessor.")
            print("   \n")

            print("------------------------------------------------------------------------" )
            print("Inputs: " )
            print("   energy [eV]:                       ",wnum * tocm / (2*numpy.pi))
            print("   wavelength [A]:                    ",(1.0/wnum) * 2 * numpy.pi*1e8 )
            print("   wavenumber (2 pi/lambda) [cm^-1]:  ",wnum )
            # wnum = 2*numpy.pi * energy1 / tocm
            print("Outputs: " )
            print("   refraction index = (1-delta) + i*beta : " )
            print("   delta:                          ",1.0-refraction_index.real )#1.0-rr_ind )
            print("   beta:                           ",refraction_index.imag) #rr_attenuation / (2*wnum) )
            print("   real(n):                        ",refraction_index.real )
            print("   attenuation coef [cm^-1]:       ",2*refraction_index.imag*wnum )
            print("------------------------------------------------------------------------" )


        return refraction_index

    def reflectivity_fresnel(self,photon_energy_ev=10000.0,grazing_angle_mrad=3.0,roughness_rms_A=0.0):
        """
        Calculates the reflectivity of an interface using Fresnel formulas.

        Code adapted from XOP and SHADOW

        :param grazing_angle_mrad: scalar with grazing angle in mrad
        :param roughness_rms_A: scalar with roughness rms in Angstroms
        :param photon_energy_ev: scalar or array with photon energies in eV
        :return: (rs,rp,runp) the s-polarized, p-pol and unpolarized reflectivities
        """
        # ;
        # ; calculation of reflectivity (piece of code adapted from shadow/abrefc)
        # ;
        #
        refraction_index = self.get_refraction_index(photon_energy_ev)

        theta1 = grazing_angle_mrad * 1e-3     # in rad
        rough1 = roughness_rms_A

        # ; epsi = 1 - alpha - i gamma
        # alpha = 2.0D0*k*f1
        # gamma = 2.0D0*k*f2

        alpha = 2 * (1.0 - refraction_index.real)
        gamma = 2 * refraction_index.imag

        rho = (numpy.sin(theta1))**2 - alpha
        rho += numpy.sqrt((numpy.sin(theta1)**2 - alpha)**2 + gamma**2)
        rho = numpy.sqrt(rho / 2)

        rs1 = 4 * (rho**2) * (numpy.sin(theta1) - rho)**2 + gamma**2
        rs2 = 4 * (rho**2) * (numpy.sin(theta1) + rho)**2 + gamma**2
        rs = rs1 / rs2

        ratio1 = 4 * rho**2 * (rho * numpy.sin(theta1) - numpy.cos(theta1)**2)**2 + gamma**2 * numpy.sin(theta1)**2
        ratio2 = 4 * rho**2 * (rho * numpy.sin(theta1) + numpy.cos(theta1)**2)**2 + gamma**2 * numpy.sin(theta1)**2
        ratio = ratio1 / ratio2

        rp = rs * ratio
        runp = 0.5 * (rs + rp)

        wavelength_m = codata.h * codata.c / codata.e / photon_energy_ev

        debyewaller = numpy.exp( -(4.0 * numpy.pi * numpy.sin(theta1) * rough1 / (wavelength_m * 1e10))**2 )

        return rs*debyewaller, rp*debyewaller, runp*debyewaller


if __name__ == "__main__":

    from srxraylib.plot.gol import plot

    #
    # refractive index
    #


    a = PreRefl()

    prerefl_file = "Be5_30.dat"
    a.read_preprocessor_file(prerefl_file)

    refraction_index = a.get_refraction_index(10000.0,verbose=True)

    #
    # mirror reflectivity
    #

    a = PreRefl()

    prerefl_file = "Rh1_50.dat"
    a.read_preprocessor_file(prerefl_file)

    Energy = numpy.linspace(5000.0,40000.0,100)

    RS0 = numpy.zeros_like(Energy)
    RS5 = numpy.zeros_like(Energy)

    a.preprocessor_info()

    for ii,ee in enumerate(Energy):

        rs, rp, rav = a.reflectivity_fresnel(grazing_angle_mrad=3.0,photon_energy_ev=ee,roughness_rms_A=0.0)

        RS0[ii] = rs

        rs, rp, rav = a.reflectivity_fresnel(grazing_angle_mrad=3.0,photon_energy_ev=ee,roughness_rms_A=5.0)

        RS5[ii] = rs


    plot(Energy,RS0,Energy,RS5,ylog=True,legend=["no roughness","5A RMS roughness"])

