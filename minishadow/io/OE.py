__author__ = 'srio'

import json
import numpy
from gfile import GFile


class OE(object):

    def __init__(self):

        self.FMIRR = 5
        self.F_TORUS = 0
        self.FCYL = 0
        self.F_EXT = 0
        self.FSTAT = 0
        self.F_SCREEN = 0
        self.F_PLATE = 0
        self.FSLIT = 0
        self.FWRITE = 0
        self.F_RIPPLE = 0
        self.F_MOVE = 0
        self.F_THICK = 0
        self.F_BRAGG_A = 0
        self.F_G_S = 0
        self.F_R_RAN = 0
        self.F_GRATING = 0
        self.F_MOSAIC = 0
        self.F_JOHANSSON = 0
        self.F_SIDE = 0
        self.F_CENTRAL = 0
        self.F_CONVEX = 0
        self.F_REFLEC = 0
        self.F_RUL_ABS = 0
        self.F_RULING = 0
        self.F_PW = 0
        self.F_PW_C = 0
        self.F_VIRTUAL = 0
        self.FSHAPE = 1
        self.FHIT_C = 0
        self.F_MONO = 0
        self.F_REFRAC = 0
        self.F_DEFAULT = 1
        self.F_REFL = 0
        self.F_HUNT = 1
        self.F_CRYSTAL = 0
        self.F_PHOT_CENT = 0
        self.F_ROUGHNESS = 0
        self.F_ANGLE = 0
        self.NPOINT = 5000
        self.NCOL = 18
        self.N_SCREEN = 0
        self.ISTAR1 = 0
        self.CIL_ANG = 0.00000000000000
        self.ELL_THE = 0.00000000000000
        self.N_PLATES = 0
        self.IG_SEED = 0
        self.MOSAIC_SEED = 1626261131
        self.ALPHA = 0.00000000000000
        self.SSOUR = 0.00000000000000
        self.THETA = 0.00000000000000
        self.SIMAG = 0.00000000000000
        self.RDSOUR = 0.00000000000000
        self.RTHETA = 0.00000000000000
        self.OFF_SOUX = 0.00000000000000
        self.OFF_SOUY = 0.00000000000000
        self.OFF_SOUZ = 0.00000000000000
        self.ALPHA_S = 0.00000000000000
        self.RLEN1 = 0.00000000000000
        self.RLEN2 = 0.00000000000000
        self.RMIRR = 0.00000000000000
        self.AXMAJ = 0.00000000000000
        self.AXMIN = 0.00000000000000
        self.CONE_A = 0.00000000000000
        self.R_MAJ = 0.00000000000000
        self.R_MIN = 0.00000000000000
        self.RWIDX1 = 0.00000000000000
        self.RWIDX2 = 0.00000000000000
        self.PARAM = 0.00000000000000
        self.HUNT_H = 0.00000000000000
        self.HUNT_L = 0.00000000000000
        self.BLAZE = 0.00000000000000
        self.RULING = 12000.0000000000
        self.ORDER = -1.00000000000000
        self.PHOT_CENT = 5.00000000000000
        self.X_ROT = 0.00000000000000
        self.D_SPACING = 0.00000000000000
        self.A_BRAGG = 0.00000000000000
        self.SPREAD_MOS = 0.00000000000000
        self.THICKNESS = 0.00000000000000
        self.R_JOHANSSON = 0.00000000000000
        self.Y_ROT = 0.00000000000000
        self.Z_ROT = 0.00000000000000
        self.OFFX = 0.00000000000000
        self.OFFY = 0.00000000000000
        self.OFFZ = 0.00000000000000
        self.SLLEN = 0.00000000000000
        self.SLWID = 0.00000000000000
        self.SLTILT = 0.00000000000000
        self.COD_LEN = 0.00000000000000
        self.COD_WID = 0.00000000000000
        self.X_SOUR = 0.00000000000000
        self.Y_SOUR = 0.00000000000000
        self.Z_SOUR = 0.00000000000000
        self.X_SOUR_ROT = 0.00000000000000
        self.Y_SOUR_ROT = 0.00000000000000
        self.Z_SOUR_ROT = 0.00000000000000
        self.R_LAMBDA = 5.00000000000000
        self.THETA_I = 0.00000000000000
        self.ALPHA_I = 0.00000000000000
        self.T_INCIDENCE = 88.0000000000000
        self.T_SOURCE = 10.0000000000000
        self.T_IMAGE = 20.0000000000000
        self.T_REFLECTION = 88.0000000000000
        self.FILE_SOURCE = "begin.dat"
        self.FILE_RIP = ""
        self.FILE_REFL = ""
        self.FILE_MIR = ""
        self.FILE_ROUGH = ""
        self.FZP = 0
        self.HOLO_R1 = 300.000000000000
        self.HOLO_R2 = 300.000000000000
        self.HOLO_DEL = -20.0000000000000
        self.HOLO_GAM = -20.0000000000000
        self.HOLO_W = 4879.86000000000
        self.HOLO_RT1 = 0.00000000000000
        self.HOLO_RT2 = 0.00000000000000
        self.AZIM_FAN = 0.00000000000000
        self.DIST_FAN = 0.00000000000000
        self.COMA_FAC = 0.00000000000000
        self.ALFA = 0.00000000000000
        self.GAMMA = 0.00000000000000
        self.R_IND_OBJ = 1.00000000000000
        self.R_IND_IMA = 1.00000000000000
        self.R_ATTENUATION_OBJ = 0.00000000000000
        self.R_ATTENUATION_IMA = 0.00000000000000
        self.F_R_IND = 0
        self.FILE_R_IND_OBJ = ""
        self.FILE_R_IND_IMA = ""
        self.RUL_A1 = 0.00000000000000
        self.RUL_A2 = 0.00000000000000
        self.RUL_A3 = 0.00000000000000
        self.RUL_A4 = 0.00000000000000
        self.F_POLSEL = 4
        self.F_FACET = 0
        self.F_FAC_ORIENT = 0
        self.F_FAC_LATT = 0
        self.RFAC_LENX = 10.0000000000000
        self.RFAC_LENY = 10.0000000000000
        self.RFAC_PHAX = 0.00000000000000
        self.RFAC_PHAY = 0.00000000000000
        self.RFAC_DELX1 = 0.00000000000000
        self.RFAC_DELX2 = 0.00000000000000
        self.RFAC_DELY1 = 0.00000000000000
        self.RFAC_DELY2 = 0.00000000000000
        self.FILE_FAC = ""
        self.F_SEGMENT = 0
        self.ISEG_XNUM = 1
        self.ISEG_YNUM = 1
        self.FILE_SEGMENT = ""
        self.FILE_SEGP = ""
        self.SEG_LENX = 0.00000000000000
        self.SEG_LENY = 0.00000000000000
        self.F_KOMA = 0
        self.FILE_KOMA = ""
        self.F_EXIT_SHAPE = 0
        self.F_INC_MNOR_ANG = 0
        self.ZKO_LENGTH = 0.00000000000000
        self.RKOMA_CX = 0.00000000000000
        self.RKOMA_CY = 0.00000000000000
        self.F_KOMA_CA = 0
        self.FILE_KOMA_CA = ""
        self.F_KOMA_BOUNCE = 0
        self.X_RIP_AMP = 0.00000000000000
        self.X_RIP_WAV = 0.00000000000000
        self.X_PHASE = 0.00000000000000
        self.Y_RIP_AMP = 0.00000000000000
        self.Y_RIP_WAV = 0.00000000000000
        self.Y_PHASE = 0.00000000000000
        self.N_RIP = 0
        self.ROUGH_X = 0.00000000000000
        self.ROUGH_Y = 0.00000000000000
        self.OE_NUMBER = 0
        self.IDUMMY = 0
        self.DUMMY = 0.00000000000000

        self.CX_SLIT = numpy.zeros(11)
        self.CX_SLIT[1] = 0.00000000000000
        self.CX_SLIT[2] = 0.00000000000000
        self.CX_SLIT[3] = 0.00000000000000
        self.CX_SLIT[4] = 0.00000000000000
        self.CX_SLIT[5] = 0.00000000000000
        self.CX_SLIT[6] = 0.00000000000000
        self.CX_SLIT[7] = 0.00000000000000
        self.CX_SLIT[8] = 0.00000000000000
        self.CX_SLIT[9] = 0.00000000000000
        self.CX_SLIT[10] = 0.00000000000000

        self.CZ_SLIT = numpy.zeros(11)
        self.CZ_SLIT[1] = 0.00000000000000
        self.CZ_SLIT[2] = 0.00000000000000
        self.CZ_SLIT[3] = 0.00000000000000
        self.CZ_SLIT[4] = 0.00000000000000
        self.CZ_SLIT[5] = 0.00000000000000
        self.CZ_SLIT[6] = 0.00000000000000
        self.CZ_SLIT[7] = 0.00000000000000
        self.CZ_SLIT[8] = 0.00000000000000
        self.CZ_SLIT[9] = 0.00000000000000
        self.CZ_SLIT[10] = 0.00000000000000

        self.D_PLATE = numpy.zeros(11)
        self.D_PLATE[1] = 0.00000000000000
        self.D_PLATE[2] = 0.00000000000000
        self.D_PLATE[3] = 0.00000000000000
        self.D_PLATE[4] = 0.00000000000000
        self.D_PLATE[5] = 0.00000000000000
        self.D_PLATE[6] = 0.00000000000000
        self.D_PLATE[7] = 0.00000000000000
        self.D_PLATE[8] = 0.00000000000000
        self.D_PLATE[9] = 0.00000000000000
        self.D_PLATE[10] = 0.00000000000000

        self.FILE_ABS = numpy.chararray(11)
        self.FILE_ABS[1] = ""
        self.FILE_ABS[2] = ""
        self.FILE_ABS[3] = ""
        self.FILE_ABS[4] = ""
        self.FILE_ABS[5] = ""
        self.FILE_ABS[6] = ""
        self.FILE_ABS[7] = ""
        self.FILE_ABS[8] = ""
        self.FILE_ABS[9] = ""
        self.FILE_ABS[10] = ""

        self.FILE_SCR_EXT = numpy.chararray(11)
        self.FILE_SCR_EXT[1] = ""
        self.FILE_SCR_EXT[2] = ""
        self.FILE_SCR_EXT[3] = ""
        self.FILE_SCR_EXT[4] = ""
        self.FILE_SCR_EXT[5] = ""
        self.FILE_SCR_EXT[6] = ""
        self.FILE_SCR_EXT[7] = ""
        self.FILE_SCR_EXT[8] = ""
        self.FILE_SCR_EXT[9] = ""
        self.FILE_SCR_EXT[10] = ""

        self.I_ABS = numpy.zeros(11)
        self.I_ABS[1] = 0
        self.I_ABS[2] = 0
        self.I_ABS[3] = 0
        self.I_ABS[4] = 0
        self.I_ABS[5] = 0
        self.I_ABS[6] = 0
        self.I_ABS[7] = 0
        self.I_ABS[8] = 0
        self.I_ABS[9] = 0
        self.I_ABS[10] = 0

        self.I_SCREEN = numpy.zeros(11)
        self.I_SCREEN[1] = 0
        self.I_SCREEN[2] = 0
        self.I_SCREEN[3] = 0
        self.I_SCREEN[4] = 0
        self.I_SCREEN[5] = 0
        self.I_SCREEN[6] = 0
        self.I_SCREEN[7] = 0
        self.I_SCREEN[8] = 0
        self.I_SCREEN[9] = 0
        self.I_SCREEN[10] = 0

        self.I_SLIT = numpy.zeros(11)
        self.I_SLIT[1] = 0
        self.I_SLIT[2] = 0
        self.I_SLIT[3] = 0
        self.I_SLIT[4] = 0
        self.I_SLIT[5] = 0
        self.I_SLIT[6] = 0
        self.I_SLIT[7] = 0
        self.I_SLIT[8] = 0
        self.I_SLIT[9] = 0
        self.I_SLIT[10] = 0

        self.I_STOP = numpy.zeros(11)
        self.I_STOP[1] = 0
        self.I_STOP[2] = 0
        self.I_STOP[3] = 0
        self.I_STOP[4] = 0
        self.I_STOP[5] = 0
        self.I_STOP[6] = 0
        self.I_STOP[7] = 0
        self.I_STOP[8] = 0
        self.I_STOP[9] = 0
        self.I_STOP[10] = 0

        self.K_SLIT = numpy.zeros(11)
        self.K_SLIT[1] = 0
        self.K_SLIT[2] = 0
        self.K_SLIT[3] = 0
        self.K_SLIT[4] = 0
        self.K_SLIT[5] = 0
        self.K_SLIT[6] = 0
        self.K_SLIT[7] = 0
        self.K_SLIT[8] = 0
        self.K_SLIT[9] = 0
        self.K_SLIT[10] = 0

        self.RX_SLIT = numpy.zeros(11)
        self.RX_SLIT[1] = 0.00000000000000
        self.RX_SLIT[2] = 0.00000000000000
        self.RX_SLIT[3] = 0.00000000000000
        self.RX_SLIT[4] = 0.00000000000000
        self.RX_SLIT[5] = 0.00000000000000
        self.RX_SLIT[6] = 0.00000000000000
        self.RX_SLIT[7] = 0.00000000000000
        self.RX_SLIT[8] = 0.00000000000000
        self.RX_SLIT[9] = 0.00000000000000
        self.RX_SLIT[10] = 0.00000000000000

        self.RZ_SLIT = numpy.zeros(11)
        self.RZ_SLIT[1] = 0.00000000000000
        self.RZ_SLIT[2] = 0.00000000000000
        self.RZ_SLIT[3] = 0.00000000000000
        self.RZ_SLIT[4] = 0.00000000000000
        self.RZ_SLIT[5] = 0.00000000000000
        self.RZ_SLIT[6] = 0.00000000000000
        self.RZ_SLIT[7] = 0.00000000000000
        self.RZ_SLIT[8] = 0.00000000000000
        self.RZ_SLIT[9] = 0.00000000000000
        self.RZ_SLIT[10] = 0.00000000000000

        self.SCR_NUMBER = numpy.zeros(11)
        self.SCR_NUMBER[1] = 0
        self.SCR_NUMBER[2] = 0
        self.SCR_NUMBER[3] = 0
        self.SCR_NUMBER[4] = 0
        self.SCR_NUMBER[5] = 0
        self.SCR_NUMBER[6] = 0
        self.SCR_NUMBER[7] = 0
        self.SCR_NUMBER[8] = 0
        self.SCR_NUMBER[9] = 0
        self.SCR_NUMBER[10] = 0

        self.SL_DIS = numpy.zeros(11)
        self.SL_DIS[1] = 0.00000000000000
        self.SL_DIS[2] = 0.00000000000000
        self.SL_DIS[3] = 0.00000000000000
        self.SL_DIS[4] = 0.00000000000000
        self.SL_DIS[5] = 0.00000000000000
        self.SL_DIS[6] = 0.00000000000000
        self.SL_DIS[7] = 0.00000000000000
        self.SL_DIS[8] = 0.00000000000000
        self.SL_DIS[9] = 0.00000000000000
        self.SL_DIS[10] = 0.00000000000000

        self.THICK = numpy.zeros(11)
        self.THICK[1] = 0.00000000000000
        self.THICK[2] = 0.00000000000000
        self.THICK[3] = 0.00000000000000
        self.THICK[4] = 0.00000000000000
        self.THICK[5] = 0.00000000000000
        self.THICK[6] = 0.00000000000000
        self.THICK[7] = 0.00000000000000
        self.THICK[8] = 0.00000000000000
        self.THICK[9] = 0.00000000000000
        self.THICK[10] = 0.00000000000000

        self.CCC = numpy.zeros(11)
        self.CCC[1] = 0.00000000000000
        self.CCC[2] = 0.00000000000000
        self.CCC[3] = 0.00000000000000
        self.CCC[4] = 0.00000000000000
        self.CCC[5] = 0.00000000000000
        self.CCC[6] = 0.00000000000000
        self.CCC[7] = 0.00000000000000
        self.CCC[8] = 0.00000000000000
        self.CCC[9] = 0.00000000000000
        self.CCC[10] = 0.00000000000000

    def load(self,filename="start.01"):
        self.load_start01(filename=filename)

    def load_start01(self,filename="start.01"):
        a = GFile()
        a.load_gfile((filename))
        dict1 = a.get_as_dictionary()

        for key in dict1.keys():
            if hasattr(self,key):
                # print("assigning: ",key,dict1[key])
                setattr(self, key, dict1[key])
            else: # this is for arrays...
                try:
                    command = "self.%s = %s"%(key,repr(dict1[key]))
                    # print("command: ",command)
                    exec(command)
                except:
                    print("   error executing: ",command)
        # self.THICK[6] = 110.0

if __name__ == "__main__":

    oe1 = OE()

    oe1.load(filename="start.01")

    print(dir(oe1))

    print(">>>>>oe1.RX_SLIT",oe1.RX_SLIT)
    print(">>>>>Y_ROT",oe1.Y_ROT )
    print(">>>>>ISTAR1",oe1.ISTAR1 )
    print(">>>>>THICK",oe1.THICK )

    assert(oe1.THICK[6] == 110.0)



