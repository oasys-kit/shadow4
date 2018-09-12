__author__ = 'srio'

#
import numpy
import scipy.constants as codata
# import math
#
tocm = 12398.419739640718e-8


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
            mylist = lines[index_pointer].strip().split("    ")
            print("///",mylist)
            delta_s[i] = float(mylist[0])
            beta_s[i] = float(mylist[1])

        for i in range(np):
            index_pointer += 1
            mylist = lines[index_pointer].strip().split("    ")
            delta_e[i] = float(mylist[0])
            beta_e[i] = float(mylist[1])

        for i in range(np):
            index_pointer += 1
            mylist = lines[index_pointer].strip().split("    ")
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
            mylist = lines[index_pointer].strip().split("   ")
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

    # !
    # ! PRE_MLAYER_SCAN
    # !
    #
    # !
    # ! this is a simple routine that computes the multilayer reflectivity
    # ! with a multilayer defined in a file created by pre_mlayer.
    # !
    # ! It can be used for testing pre_mlayer, or for simple calculations of
    # ! ML reflectivity.
    # !
    #  !
    # SUBROUTINE pre_mlayer_scan ()
    def scan(self):


        # ! initializa some variables
        F_REFL = 2
        pin = [0.0,0.0,0.0]
        # pin(2) = 0.0
        # pin(3) = 0.0
        fileOut = "mlayer_scan.dat"

        if self.pre_mlayer_dict is None:
            raise Exception("load preprocessor file before!")
        #
        # !
        # ! input section
        # !
        # print *,"   pre_mlayer_scan: calculates reflectivity of a multilayer"
        # print *,"                 using a file created by the pre_mlayer preprocessor."
        # print *,"   "
        # FILE_REFL  = RSTRING("File Name (from pre_mlayer): ")
        # !file_refl = "morawe.dat"
        #
        energyN = 51 # irint(' Number of energy points (1 for angle-scan): ')
        thetaN = 1    # irint(' Number of anglular points (1 for energy-scan): ')
        # !energyN = 1
        # !thetaN = 5000
        #
        # if (energyN .GT. 1) then
        #      !energy1 = 12400.0
        #      !energy2 = 24800.0
        if energyN == 1:
            energy1 = 17800.0
        else:
            energy1 = 5000.0  # rnumber("Photon energy from [eV]: ")
            energy2 = 20000.0 # rnumber("              to [eV]: ")
        # else
        #      !energy1 = 17800.0
        #      energy1 = rnumber("Photon energy [eV]: ")
        # endif
        #
        # if (thetaN .GT. 1) then
        #      !theta1 = 0.05
        #      !theta2 = 3.0
        if thetaN == 1:
            theta1 = 0.75
        else:
            theta1 = 0.03  #  rnumber("Incident grazing angle from [deg]: ")
            theta2 = 3.0  #  rnumber("                       to [deg]: ")
        # else
        #      !theta1 = 0.75
        #      theta1 = rnumber("Incident grazing angle [deg]: ")
        # endif
        #
        #
        # ! calculations
        #
        # ! initialization (read file)
        # k_what = 0
        # sin_ref = 0.0
        # cos_pole = 1.0
        # wnum = twopi*energy1/tocm
        # CALL REFLEC (PIN,WNUM,SIN_REF,COS_POLE,R_P,R_S,PHASEP,PHASES,ABSOR,K_WHAT)
        #
        # !calculate
        k_what = 1
        if (energyN > 1):
            energyS = (energy2-energy1)/float(energyN-1)
        else:
            energyS = 0.0

        if (thetaN > 1):
            thetaS = (theta2-theta1)/float(thetaN-1)
        else:
            thetaS = 0.0

        f = open(fileOut,'w')

        f.write("#F %s\n"%(fileOut))
        f.write("\n")
        f.write("#S 1 pre_mlater_test results\n")
        f.write("#N 4\n" )
        f.write("#L energy[eV]  grazingAngle [deg]  R_S  R_P\n")

        # do i=1,thetaN
        #   do j=1,energyN
        R_S = 0.0
        R_P = 0.0
        for i in range(1,1+thetaN):
            for j in range(1,1+energyN):
                theta = theta1 + float(i-1) * thetaS
                energy = energy1 + float(j-1) * energyS
                sin_ref = numpy.sin ( theta * numpy.pi/180)
                wnum = 2 * numpy.pi * energy / tocm

                theta = theta1 + float(i-1) * thetaS
                energy = energy1+ float(j-1)*energyS
                sin_ref = numpy.sin ( theta * numpy.pi/180)
                wnum = 2 * numpy.pi * energy / tocm

                COS_POLE = 1.0
                R_P,R_S,PHASEP,PHASES,ABSOR = self.reflec(wnum,sin_ref,COS_POLE,k_what)

            #     CALL REFLEC (PIN,WNUM,SIN_REF,COS_POLE,R_P,R_S,PHASEP,PHASES,ABSOR,K_WHAT)
            #
                if ( (thetaN  ==  1) and (energyN == 1) ):
                    print("------------------------------------------------------------------------")
                    print("Inputs: ")
                    print("   for E=",energy1,"eV: ")
                    print("   energy [eV]:                       ",energy)
                    print("   grazing angle [deg]:               ",theta)
                    print("   wavelength [A]:                    ",(1e0/wnum)*2*numpy.pi*1e8)
                    print("   wavenumber (2 pi/lambda) [cm^-1]:  ",wnum)
                    print("Outputs: ")
                    print("   R_S:                          ",R_S)
                    print("   R_P:                          ",R_P)
                    print("------------------------------------------------------------------------")
                else:
                    f.write("%f  %f  %f  %f\n"%(energy,theta,R_S,R_P))

        f.close()
        print("File %s written to disk"%fileOut)


    def reflec(self,WNUM,SIN_REF,COS_POLE,K_WHAT):

        # ! C+++
        # ! C	SUBROUTINE	REFLEC
        # ! C
        # ! C	PURPOSE		To compute the local reflectivity of a mirror or
        # ! C                     multilayer. Also compute filter transmittivity.
        # ! C
        # ! C	FLAGS		k_what:  .eq. 0 --> initialization call. Reads
        # ! C					in data file.
        # ! C			         .gt. 0 --> performs computations.
        # ! C			(see below)
        # ! C
        # ! C	ARGUMENTS	[ I ] PIN	: (x,y,z) of the intercept
        # ! C			[ I ] wnum 	: wavenumber (cm-1)
        # ! C			[ I ] sin_ref	: sine of angle from surface
        # ! C			[ I ] cos_pole  : cosine of angle of normal from pole
        # ! C			[ O ] R_P 	: p-pol reflection coefficient
        # ! C			[ O ] R_S 	: s-pol    "  "
        # ! C			[ I ] ABSOR	: film thickness
        # ! C			[ O ] ABSOR	: absorption coefficient
        # ! C
        # ! C---
        # SUBROUTINE REFLEC (PIN,WNUM,SIN_REF,COS_POLE,R_P,R_S,PHASEP,PHASES,ABSOR,K_WHAT)
        #
        # implicit none
        #
        # real(kind=skr),dimension(3),intent(in)   :: pin
        # real(kind=skr),             intent(in)   :: wnum,sin_ref,cos_pole
        # integer(kind=ski),          intent(in)   :: k_what
        # real(kind=skr),             intent(inout):: absor
        # real(kind=skr),             intent(out)  :: r_s,r_p,phases,phasep
        #
        # integer(kind=ski), parameter  :: dimMLenergy=300
        #
        # real(kind=skr),dimension(:),allocatable   :: zf1,zf2
        # real(kind=skr),dimension(dimMLenergy)    :: ener, &
        #             delta_s,beta_s,delta_e,beta_e,delta_o,beta_o
        #
        # logical,dimension(5)    ::  ele,elo,els
        # character(len=sklen)    ::  file_grade
        #
        # real(kind=skr),dimension(2,101,2,101) :: tspl,gspl
        # real(kind=skr),dimension(101)         :: tx,ty,gx,gy
        # real(kind=skr),dimension(6)           :: pds
        #
        # real(kind=skr)   :: lateral_grade_constant=1D0 !initialization avoids save
        # real(kind=skr)   :: lateral_grade_slope=0D0
        # real(kind=skr)   :: lateral_grade_quadratic=0D0
        # real(kind=skr)   :: lateral_grade_cubic=0D0
        # integer(kind=ski):: i_grade=0
        #
        # real(kind=skr)   :: ab_coeff, cos_ref, del_x, depth0, elfactor, gfact
        # real(kind=skr)   :: qmin, qmax, qstep, ratio, phot_ener, ratio1, ratio2
        # real(kind=skr)   :: rho, rs1, rs2, tfact, tfilm, wnum0, xin, xlam, yin, gamma1
        # integer(kind=ski):: i,j,nrefl,ierr,ier,index1,iunit
        # integer(kind=ski):: ngx, ngy, ntx, nty, nin, npair
        #
        # ! C
        # ! C SAVE the variables that need to be saved across subsequent invocations
        # ! C of this subroutine.
        # ! C
        # SAVE        QMIN, QMAX, QSTEP, DEPTH0, NREFL, TFILM, &
        #             ZF1, ZF2, &
        #             NIN, ENER,  &
        #             DELTA_S, BETA_S,  &
        #             NPAIR, &
        #             DELTA_E,BETA_E, &
        #             DELTA_O,BETA_O, &
        #             TSPL,TX,TY,PDS, &
        #             GSPL,GX,GY, &
        #             NTX, NTY, NGX, NGY  ! added srio@esrf.eu 20130917
        # ! C
        # ! C Initialization call. The ZF1,ZF2 values do NOT correspond to the F1,F2
        # ! C atomic scattering factors, as they contain a more complex form:
        # ! C		ZFi = Fi*RADIUS*4*PI*ATOMS
        # ! C WNUM is the WAVENUMBER (cm-1) of the ray.
        # ! C ALFA and GAMMA are the complex dielectric function
        # ! C		EPSILON	  = 1 - ALFA + i*GAMMA		[ i=sqrt(-1) ]
        # ! C and are dependent ONLY on the material, while the reflectivities
        # ! C depend of course on the geometry too.
        # ! C ALFA and GAMMA may be generated by PREREF.EXE in [CERRINA.ABS],
        # ! C which is based on B.Henke data (see program header for complete reference).
        # ! C
        # ! C Two flags control the execution.
        # ! C 	F_REFL = 0	ZF1,ZF2 are read in as arrays
        # ! C	       = 1	ALFA and GAMMA are defined in the I/O session
        # ! C			and thus wavelength independent
        # ! C	       = 2      d-spacings and optical constants for multilayer
        # ! C			are read in
        # ! C 	K_WHAT = 0	Initialization
        # ! C		 1	Reflection case
        # ! C		 2	Absorption case
        # ! C
        #
        # !todo: compute phases for the usual Fresnel reflection (mirrors). They
        # !      are set to zero here to avoid NaN in Macs. srio@esrf.eu 20150310
        phases = 0.0
        phasep = 0.0
        #
        # IF (K_WHAT.EQ.0) THEN
        #     ELSE IF (F_REFL.EQ.2) THEN  !multilayer
        #         ! C
        #         ! C  this version allows specification of the individual
        #         ! C  layer thicknesses.
        #         ! C
        #         ! C  input parameters:
        #         ! C  npair = no. of layer pairs (npair=0 means an ordinary mirror,
        #         ! C          elo is the mirror
        #         ! C  xlam = wavelength (angstroms)
        #         ! C  elo = odd layer material
        #         ! C  ele = even layer material
        #         ! C  els = substrate material
        #         ! C  1.0 - delo - i*beto = complex refractive index (odd)
        #         ! C  1.0 - dele - i*bete = complex refractive index (even)
        #         ! C  t_oe   = thickness t(odd)+t(even) in Angstroms of each layer pair
        #         ! C  gratio = gamma ratio t(even)/(t(odd)+t(even))  of each layer pair
        #         ! C  phr = grazing angle in radians
        #         ! C
        #         ! C
        #         iunit = 23
        #         ! WARNING: I got sometimes segmentation fault around this point.
        #         !          Problem not identified....  srio@esrf.eu 2010-08-26
        #         open(unit=iunit,FILE=FILE_REFL,status='OLD',IOSTAT=iErr)
        #         ! srio added test
        #         if (iErr /= 0 ) then
        #             print *,"MIRROR: Error: File not found: "//trim(file_refl)
        #             return
        #             ! stop 'File not found. Aborted.'
        #         end if
        #         READ(iunit,*) NIN
        NIN = self.pre_mlayer_dict["np"]


        #         IF (NIN > dimMLenergy) THEN
        #             print *,'REFLEC: Error: In file: '//trim(file_refl)
        #             print *,'               Maximum number of energy points is',dimMLenergy
        #             print *,'               Using number of energy points',NIN
        #             print *,'MIRROR: Error reading file. Aborted.'
        #             return
        #         END IF
        #         READ(iunit,*) (ENER(I), I = 1, NIN)

        #wnum = 2 * numpy.pi * energy / tocm

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

        NPAIR = numpy.abs(self.pre_mlayer_dict["npair"])

        gamma1 = self.pre_mlayer_dict["gamma1"]
        t_oe = self.pre_mlayer_dict["thick"]

        # gamma1 = ratio t(even)/(t(odd)+t(even))  of each layer pair
        t_e = gamma1 * t_oe
        t_o = (1.0 - gamma1) * t_oe

        mlroughness1 = self.pre_mlayer_dict["mlroughness1"]
        mlroughness2 = self.pre_mlayer_dict["mlroughness2"]

        #         DO 13 I=1,NIN
        #             READ(iunit,*) DELTA_S(I),BETA_S(I)
        # 13      CONTINUE
        #         DO 23 I=1,NIN
        #             READ(iunit,*)  DELTA_E(I),BETA_E(I)
        # 23      CONTINUE
        #         DO 33 I=1,NIN
        #             READ(iunit,*) DELTA_O(I),BETA_O(I)
        # 33      CONTINUE
        #         READ(iunit,*) NPAIR
        #         if(npair .lt. 0) then ! if npair<0 roughness data is available
        #             do i = 1, abs(npair)
        #                 read(iunit,*) t_oe(i),gratio(i),mlroughness1(i),mlroughness2(i)
        #             end do
        #         else
        #             do i = 1, npair
        #                 read(iunit,*) t_oe(i),gratio(i)
        #                 mlroughness1(i)=0.0
        #                 mlroughness2(i)=0.0
        #             end do
        #         endif
        #         npair = abs(npair)

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




        # ! C
        # ! C This is the normal calculation part;
        # ! C
        # ! C If F_REFL is 1, ALFA and GAMMA are defined during the input session
        # ! C and are not modified anymore (single line or closely spaced lines case)
        # ! C
        # IF (F_REFL.EQ.0) THEN      !Both absorp and normal
        #     !reflectivity
        #     index1 =   (WNUM - QMIN)/QSTEP + 1

        # wnum = 2 * numpy.pi * energy / tocm
        PHOT_ENER  = WNUM * tocm / (2*numpy.pi)   # eV


        # index1 =   (WNUM - QMIN) / QSTEP #+ 1
        # index1 = (PHOT_ENER - ENER[0]) / (ENER[1]-ENER[0])
        #
        # print(">>>>",index1,ENER.size)
        #
        # print("photon energy: ",PHOT_ENER,ENER[0],ENER[-1],"  q:%g %g %g %g"%(WNUM,QMIN,wnum[-1],QSTEP))

        # if index1 < 0:
        #     raise Exception("Photon energy below lower limit.")
        # if index1 > ENER.size:
        #     raise Exception("Photon energy above upper limit.")

        #     ! see http://ftp.esrf.fr/pub/scisoft/shadow/user_contributions/compilation_fix2008-04-09.txt
        #     IF (index1.LT.1) THEN
        #        index1=1
        #        ! C     ('REFLEC','Photon energy below lower limit.',0)
        #        print *,"REFLEC: Warning: Photon energy below lower limit. Rerun prerefl."
        #     END IF

        #     IF (index1.GT.NREFL) THEN
        #        index1=NREFL-1
        #        ! C     ('REFLEC','Photon energy above upper limit.',0)
        #        print *,"REFLEC: Warning: Photon energy above upper limit. Rerun prerefl."
        #     END IF
        #     IF (index1.EQ.NREFL)  index1  = index1 - 1

        #     WNUM0  =   QSTEP*(index1-1) + QMIN
        #     DEL_X  =   WNUM - WNUM0
        #     DEL_X  =   DEL_X/QSTEP
        #     ALFA  =   ZF1(index1) + (ZF1(index1+1)-ZF1(index1))*DEL_X
        #     gamma1  =   ZF2(index1) + (ZF2(index1+1)-ZF2(index1))*DEL_X



        #     ! D  WRITE (37,1020) WNUM,WNUM0,INDEX,DEL_X,ALFA,GAMMA
        #     ! D1020  FORMAT (1X,2(E15.8,1X),I4,3(E15.8,1X))
        # END IF
        #
        # IF (K_WHAT.EQ.1) THEN
        #     IF (F_REFL.NE.2) THEN
        #         ! C
        #         ! C Computes the optical coefficients.
        #         ! C
        #         COS_REF =  SQRT(1.0D0 - SIN_REF**2)
        #         RHO  =   SIN_REF**2 - ALFA
        #         RHO  =   RHO + SQRT ((SIN_REF**2 - ALFA)**2 + gamma1**2)
        #         RHO  =   SQRT(RHO/2)
        #         ! C
        #         ! C Computes now the reflectivities
        #         ! C
        #         RS1  =   4*(RHO**2)*(ABS(SIN_REF)-RHO)**2 + gamma1**2
        #         RS2  =   4*(RHO**2)*(ABS(SIN_REF)+RHO)**2 + gamma1**2
        #         R_S  =   RS1/RS2
        #         ! C
        #         ! C Computes now the polarization ratio
        #         ! C
        #         RATIO1  =   4*RHO**2*(RHO*ABS(SIN_REF)-COS_REF**2)**2 + &
        #         gamma1**2*SIN_REF**2
        #         RATIO2  =   4*RHO**2*(RHO*ABS(SIN_REF)+COS_REF**2)**2 + &
        #         gamma1**2*SIN_REF**2
        #         RATIO  =   RATIO1/RATIO2
        #         ! C
        #         ! C The reflectivity for p light will be
        #         ! C
        #         R_P  =   R_S*RATIO
        #         R_S  =   SQRT(R_S)
        #         R_P  =   SQRT(R_P)
        #     ELSE


        #         ! C
        #         ! C Multilayers reflectivity.
        #         ! C First interpolate for all the refractive indices.
        #         ! C
        #         XLAM  =   TWOPI/WNUM*1.0D8    ! Angstrom
        #         PHOT_ENER  = WNUM/TWOPI*TOCM  ! eV
        #         ELFACTOR  = LOG10(1.0D04/30.0D0)/300.0D0
        #         index1  = LOG10(PHOT_ENER/ENER(1))/ELFACTOR + 1

        XLAM  =   2*numpy.pi / WNUM * 1.0e8    # Angstrom

        ELFACTOR  = numpy.log10(1.0e4/30.0e0)/300.0e0
        index1  = numpy.log10(PHOT_ENER/ENER[0])/ELFACTOR

        print(PHOT_ENER,index1)

        #         ! C    INDEX  = 96.0*LOG10(PHOT_ENER/ENER(1)) + 1
        #         ! see http://ftp.esrf.fr/pub/scisoft/shadow/user_contributions/compilation_fix2008-04-09.txt
        #         IF (index1.LT.1) index1=1
        #         ! C    ('REFLEC','Photon energy too small.',2)
        #         IF (index1.GT.NIN) index1=NIN-1
        #         ! C    ('REFLEC','Photon energy too large.',2)
        #         DELS  = DELTA_S(index1) + (DELTA_S(index1+1) - DELTA_S(index1)) &
        #              *(PHOT_ENER - ENER(index1))/(ENER(index1+1) - ENER(index1))
        #         BETS  = BETA_S(index1) + (BETA_S(index1+1) - BETA_S(index1)) &
        #              *(PHOT_ENER - ENER(index1))/(ENER(index1+1) - ENER(index1))
        #         DELE  = DELTA_E(index1) + (DELTA_E(index1+1) - DELTA_E(index1)) &
        #              *(PHOT_ENER - ENER(index1))/(ENER(index1+1) - ENER(index1))
        #         BETE  = BETA_E(index1) + (BETA_E(index1+1) - BETA_E(index1)) &
        #              *(PHOT_ENER - ENER(index1))/(ENER(index1+1) - ENER(index1))
        #         DELO  = DELTA_O(index1) + (DELTA_O(index1+1) - DELTA_O(index1)) &
        #              *(PHOT_ENER - ENER(index1))/(ENER(index1+1) - ENER(index1))
        #         BETO  = BETA_O(index1) + (BETA_O(index1+1) - BETA_O(index1)) &
        #              *(PHOT_ENER - ENER(index1))/(ENER(index1+1) - ENER(index1))
        #

        index1 = int(index1)

        DELS  = DELTA_S[index1] + (DELTA_S[index1+1] - DELTA_S[index1]) *(PHOT_ENER - ENER[index1])/(ENER[index1+1] - ENER[index1])
        BETS  =  BETA_S[index1] + ( BETA_S[index1+1] -  BETA_S[index1]) *(PHOT_ENER - ENER[index1])/(ENER[index1+1] - ENER[index1])
        DELE  = DELTA_E[index1] + (DELTA_E[index1+1] - DELTA_E[index1]) *(PHOT_ENER - ENER[index1])/(ENER[index1+1] - ENER[index1])
        BETE  =  BETA_E[index1] + ( BETA_E[index1+1] -  BETA_E[index1]) *(PHOT_ENER - ENER[index1])/(ENER[index1+1] - ENER[index1])
        DELO  = DELTA_O[index1] + (DELTA_O[index1+1] - DELTA_O[index1]) *(PHOT_ENER - ENER[index1])/(ENER[index1+1] - ENER[index1])
        BETO  =  BETA_O[index1] + ( BETA_O[index1+1] -  BETA_O[index1]) *(PHOT_ENER - ENER[index1])/(ENER[index1+1] - ENER[index1])



        #         ! C
        #         ! C CALL FRESNEL (NPAIR,SIN_REF,COS_POLE,XLAM,R_S,R_P,PHASES,PHASEP)
        #         ! C
        #         ! C
        #         ! C If graded, compute the factor for t and gamma at the intercept PIN.
        #         ! C
        #         TFACT       = 1.0D0
        #         GFACT       = 1.0D0
        TFACT = 1.0
        GFACT = 1.0

        #TODO
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
        #
        #
        print(">>>>>>>>>>>>>>> calling fresnel with:",TFACT,GFACT,NPAIR,SIN_REF,COS_POLE,XLAM,)
        print("              ",DELO,DELE,DELS,BETO,BETE,BETS,)

        R_S,R_P,PHASES,PHASEP = self.fresnel(TFACT,GFACT,NPAIR,SIN_REF,COS_POLE,XLAM,
                                             DELO,DELE,DELS,BETO,BETE,BETS,t_o,t_e,mlroughness1,mlroughness2)

        #         CALL FRESNEL  (TFACT,GFACT,NPAIR,SIN_REF,COS_POLE,XLAM, &
        #                          R_S,R_P,PHASES,PHASEP)
        #     END IF
        # ELSE IF(K_WHAT.EQ.2) THEN
        #     ! C
        #     ! C This is the transmission case. SIN_REF is now the incidence angle
        #     ! C onto the filter.
        #     ! C
        #     ! C Computes now the penetration depth. SIN_REF is now the cosine of
        #     ! C the incidence angle of the ray on the screen.
        #     ! C
        #     ! C         DEPTH  =   DEPTH0*GAMMA*SIN_REF/WNUM
        #     AB_COEFF  = WNUM*gamma1/ABS(SIN_REF)
        #     ! C
        #     ! C Computes the film absorption. The thickness is passed at the call with
        #     ! C K_WHAT = 0
        #     ! C
        #     ! C ABSOR is the attenuation of the A vector.
        #     ! C
        #     ABSOR  =   EXP(-TFILM*AB_COEFF/2.0D0)
        # END IF
        # RETURN
        # End Subroutine reflec
        return R_S,R_P,0,phases,phasep

    def fresnel(self,TFACT,GFACT,NPAIR,SIN_REF,COS_POLE,XLAM,
                delo,dele,dels,beto,bete,bets,t_o,t_e,mlroughness1,mlroughness2):
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
        # !c outputs:
        # !c         ans = S polarization  reflectivity
        # !c         anp = P polarization  reflectivity
        # !c         phaseS = change of phase S
        # !c         phaseP = change di phase P
        # !c  other variables:
        # !c        delo,dele,dels = parameter delta odd, even, substrate respectively
        # !c        belo,bele,bels = parametro beta odd, even, substrate respectively
        # !c        1.0 - delo - i*beto = complex refractive index (odd)
        # !c        1.0 - dele - i*bete = complex refractive index (even)
        # !c        t_o = thickness of odd layers (a)
        # !c        t_e = thickness of even layers (a)
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
        #
        #
        # subroutine FRESNEL (tfact,gfact,n,sin_ref,cos_pole,xlam,ans, anp,phaseS,phaseP)
        #
        # implicit none
        #
        # real(kind=skr),  intent(in)  :: tfact,gfact,sin_ref,cos_pole,xlam
        # real(kind=skr),  intent(out) :: ans, anp,phaseS,phaseP
        #
        # real(kind=skr)     :: xmfv,sin_ref2, cos_ref2, pp, qq, refv
        # integer(kind=ski)  :: i,j,n
        #
        # complex(kind=skx)  ::  ci,fo,fe,fv,ffe,ffv,ffvp,ffo,ffep,ffop,re2
        # complex(kind=skx)  ::  ro2,ao,ae,r,rp,fs,ffs,ffsp,rs2
        # real(kind=skr)     ::  gamma,thick,t_e,t_o
        #
        # ! nevot-croce roughness
        # real(kind=skr)     ::  sigma_o2,sigma_e2,sigma_s2,sigma_v2
        # complex(kind=skx)  ::  arg_o,arg_e,arg_s,arg_v
        # complex(kind=skx)  ::  fnevot_o,fnevot_e,fnevot_s,fnevot_v
        # real(kind=skr)     ::  prefact
        #
        # !--------------------------------------
        # ! another way...
        # !complex(kind=skx)  ::  ro,re,rs,sin_tra,cos_tra,qo2,qe2
        # !complex(kind=skx)  ::  sin_s,cos_s,sin_o,cos_o,sin_e,cos_e
        # !real(kind=skr)     ::  cos_ref
        # !--------------------------------------
        #
        #
        # ! C
        #
        # ! "i" opmplex
        # ci=(0.0D0,1.0D0)
        ci = 0.+1.0j
        #
        # ! (refraction index "odd,even,substrate")**2
        # ro2=(1.0D0-delo-ci*beto)**2
        # re2=(1.0D0-dele-ci*bete)**2
        # rs2=(1.0D0-dels-ci*bets)**2

        ro2 = (1.0 - delo - ci * beto)**2
        re2 = (1.0 - dele - ci * bete)**2
        rs2 = (1.0 - dels - ci * bets)**2

        #
        # ! angles
        SIN_REF2 = SIN_REF**2
        COS_REF2 = 1.0 - SIN_REF2

        #
        # ! f(o,e) = sin theta_inc - sin theta_ critical
        # fo = ro2 - COS_REF2
        # fe = re2 - COS_REF2
        # refv=SIN_REF2
        # xmfv=0.0D0


        fo = ro2 - COS_REF2
        fe = re2 - COS_REF2
        refv = SIN_REF2
        xmfv = 0.0

        #
        # fv = Dcmplx(refv,xmfv)
        # fs = rs2 - COS_REF2
        #
        # !fo=cDsqrt(fo)
        # !fe=cDsqrt(fe)
        # !fv=cDsqrt(fv)
        # !fs=cDsqrt(fs)
        # fo=mysqrt(fo)
        # fe=mysqrt(fe)
        # fv=mysqrt(fv)
        # fs=mysqrt(fs)

        fv = refv + ci * xmfv
        fs = rs2 - COS_REF2


        fo = numpy.sqrt(fo) # complex!!
        fe = numpy.sqrt(fe) # complex!!
        fv = numpy.sqrt(fv) # complex!!
        fs = numpy.sqrt(fs) # complex!!


        #
        # ! Fresnel formula "S" (in function of incidence angle and critical angle)
        # ffe=(fe-fo)/(fe+fo)
        # ffo=-ffe
        # ffv=(fv-fo)/(fv+fo)
        # ffs=(fe-fs)/(fe+fs)
        # ! Fresnel formula "P" (in function of incidence angle and critical angle)
        # ffep=(fe/re2-fo/ro2)/(fe/re2+fo/ro2)
        # ffop=-ffep
        # ffvp=(fv-fo/ro2)/(fv+fo/ro2)
        # ffsp=(fe/re2-fs/rs2)/(fe/re2+fs/rs2)
        #

        ffe = (fe-fo)/(fe+fo)
        ffo = -ffe
        ffv = (fv-fo)/(fv+fo)
        ffs = (fe-fs)/(fe+fs)
        # ! Fresnel formula "P" (in function of incidence angle and critical angle)
        ffep = (fe/re2-fo/ro2)/(fe/re2+fo/ro2)
        ffop = -ffep
        ffvp = (fv-fo/ro2)/(fv+fo/ro2)
        ffsp = (fe/re2-fs/rs2)/(fe/re2+fs/rs2)




        # !-----------------------------------
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
        # !-----------------------------------
        #



        # ! reflectivity initialization
        # r=(0.0D0,0.0D0)
        # rp=(0.0D0,0.0D0)
        # prefact=(8.*(PI**2.))/(xlam**2)
        #
        #
        # !c Nevot-Croce roughness
        # !c DO NOT include refraction index in the roughness formula
        # sigma_s2=0.0d0 ! sigma_s**2.0 !roughn. substrate
        # sigma_v2=0.0d0 ! sigma_v**2.0!roughn. vacuum

        r = 0.0 + 0.0j
        rp = 0.0 + 0.0j
        prefact = (8.*(numpy.pi**2.)) / (XLAM**2)


        # !c Nevot-Croce roughness
        # !c DO NOT include refraction index in the roughness formula

        sigma_s2 = 0.0 # ! sigma_s**2.0 !roughn. substrate
        sigma_v2 = 0.0 # ! sigma_v**2.0!roughn. vacuum


        #
        # ! loop over the bilayers
        # ! remember thet "even" is the bottom sublayer
        # do 1 j=1,n   ! n is the number of bilayers
        #          ! C
        #          ! C compute the thickness for the odd and even material :
        #          ! C
        #          THICK = T_OE(J) * TFACT
        #          GAMMA = GRATIO(J) * GFACT
        #          t_e = GAMMA * THICK
        #          t_o = (1.0D0-GAMMA) * THICK
        #          ! C
        #          ao=-ci*(pi*fo*t_o*cos_pole/xlam)
        #          ae=-ci*(pi*fe*t_e*cos_pole/xlam)
        #          ao=cDexp(ao)
        #          ae=cDexp(ae)
        #          if(j.eq.1)go to 6
        #          ! even (botton) sublayer
        #          sigma_e2=mlroughness1(j)**2.0 !roughn. even layer
        #          !arg_e=FO*FE*sigma_e2/(CDSqrt(ro2)*CDSqrt(re2))
        #          arg_e=FO*FE*sigma_e2/(mysqrt(ro2)*mysqrt(re2))
        #          fnevot_e=cdexp(-prefact*arg_e)
        #          r=(ae**4)*(r+ffe*fnevot_e)/(r*ffe*fnevot_e+1.0)
        #          rp=(ae**4)*(rp+ffep*fnevot_e)/(rp*ffep*fnevot_e+1.0)
        #          !r=(ae**4)*(r+ffe)/(r*ffe+1.0D0)
        #          !rp=(ae**4)*(rp+ffep)/(rp*ffep+1.0D0)
        #          go to 7
        # 6        continue
        #          ! layer on top of substrate
        #          !arg_s=FE*FS*sigma_s2/(CDSqrt(re2)*CDSqrt(rs2))
        #          arg_s=FE*FS*sigma_s2/(mysqrt(re2)*mysqrt(rs2))
        #          fnevot_s=cdexp(-prefact*arg_s)
        #          r=(ae**4.0)*(r+ffs*fnevot_s)/(r*ffs*fnevot_s+1.0)
        #          rp=(ae**4.0)*(rp+ffsp*fnevot_s)/(rp*ffsp*fnevot_s+1.0)
        #          !r=(ae**4)*(r+ffs)/(r*ffs+1.0D0)
        #          !rp=(ae**4)*(rp+ffsp)/(rp*ffsp+1.0D0)
        # 7        continue
        #          ! odd layer (top sublayer)
        #          sigma_o2=mlroughness2(j)**2.0 !roughn. odd layer
        #          !arg_o=FO*FE*sigma_o2/(CDSqrt(ro2)*CDSqrt(re2))
        #          arg_o=FO*FE*sigma_o2/(mysqrt(ro2)*mysqrt(re2))
        #          fnevot_o=cdexp(-prefact*arg_o)
        #          r=(ao**4.0)*(r+ffo*fnevot_o)/(r*ffo*fnevot_o+1.0)
        #          rp=(ao**4.0)*(rp+ffop*fnevot_o)/(rp*ffop*fnevot_o+1.0)
        #          !r=(ao**4)*(r+ffo)/(r*ffo+1.0D0)
        #          !rp=(ao**4)*(rp+ffop)/(rp*ffop+1.0D0)
        # 1 continue
        #


        #
        # ! loop over the bilayers
        # ! remember thet "even" is the bottom sublayer

        for j in range(NPAIR): #   =1,n   ! n is the number of bilayers
            # ! C
            # ! C compute the thickness for the odd and even material :
            # ! C
            # THICK = T_OE(J) * TFACT
            # GAMMA = GRATIO(J) * GFACT
            # t_e = GAMMA * THICK
            # t_o = (1.0D0-GAMMA) * THICK
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
        # !arg_v=fo*fv*sigma_v2/CDSqrt(ro2)
        # arg_v=fo*fv*sigma_v2/mysqrt(ro2)
        # fnevot_v=cdexp(-prefact*arg_v)
        # r=(r+ffv*fnevot_v)/(r*ffv*fnevot_v+1.0)
        # !r=(r+ffv)/(r*ffv+1.0D0)
        # !
        # !added srio@esrf.eu 2012-06-07
        # !rp=(rp+ffvp)/(rp*ffvp+1.0)
        # rp=(rp+ffvp*fnevot_v)/(rp*ffvp*fnevot_v+1.0)



        arg_v = fo * fv * sigma_v2 / numpy.sqrt(ro2)
        fnevot_v = numpy.exp(-prefact * arg_v)
        r = (r + ffv * fnevot_v) / (r * ffv * fnevot_v + 1.0)
        rp = (rp + ffvp * fnevot_v) / (rp * ffvp * fnevot_v + 1.0)


        #
        # !
        # ! calculate phases
        # !
        # pp = Dimag(r)
        # qq = Dreal(r)
        pp = r.imag
        qq = r.real
        PHASES = numpy.arctan2(pp,qq)
        ans = numpy.abs(r)  # TODO: check correct function
        # anp=cDabs(rp)

        #
        # PP = DIMAG(RP)
        # QQ = DREAL (RP)
        # CALL ATAN_2(PP,QQ,PHASEP)       ! P phase change in unit
        # s of radians
        # ans=cDabs(r)
        pp = rp.imag
        qq = rp.real
        PHASEP = numpy.arctan2(pp,qq)
        anp = numpy.abs(rp)  # TODO: check correct function
        # anp=cDabs(rp)


        # ! C      ans=ans**2
        # !
        # ! end
        # !
        # return
        #
        # End Subroutine fresnel

        print("                       ",ans,anp,PHASES,PHASEP,pp,qq)
        return ans,anp,PHASES,PHASEP

if __name__ == "__main__":

    # a = MLayer.pre_mlayer()

    filename = "/Users/srio/Oasys/mlayer.dat"

    a = MLayer()
    a.read_preprocessor_file(filename)

    for k in a.pre_mlayer_dict:
        print(k) #,a.pre_mlayer_dict[k])

    a.scan()
