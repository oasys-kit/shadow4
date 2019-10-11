from shadow4.compatibility.gfile import GFile

class Source(object):

    def __init__(self):
        self.FDISTR = 2
        self.FGRID = 0
        self.FSOUR = 3
        self.FSOURCE_DEPTH = 1
        self.F_COHER = 0
        self.F_COLOR = 1
        self.F_PHOT = 1
        self.F_POL = 3
        self.F_POLAR = 1
        self.F_OPD = 1
        self.F_WIGGLER = 0
        self.F_BOUND_SOUR = 0
        self.F_SR_TYPE = 0
        self.ISTAR1 = 6775431
        self.NPOINT = 5000
        self.NCOL = 18
        self.N_CIRCLE = 0
        self.N_COLOR = 2
        self.N_CONE = 0
        self.IDO_VX = 1
        self.IDO_VZ = 1
        self.IDO_X_S = 1
        self.IDO_Y_S = 1
        self.IDO_Z_S = 1
        self.IDO_XL = 0
        self.IDO_XN = 0
        self.IDO_ZL = 0
        self.IDO_ZN = 0
        self.SIGXL1 = 0.00000000000000
        self.SIGXL2 = 0.00000000000000
        self.SIGXL3 = 0.00000000000000
        self.SIGXL4 = 0.00000000000000
        self.SIGXL5 = 0.00000000000000
        self.SIGXL6 = 0.00000000000000
        self.SIGXL7 = 0.00000000000000
        self.SIGXL8 = 0.00000000000000
        self.SIGXL9 = 0.00000000000000
        self.SIGXL10 = 0.00000000000000
        self.SIGZL1 = 0.00000000000000
        self.SIGZL2 = 0.00000000000000
        self.SIGZL3 = 0.00000000000000
        self.SIGZL4 = 0.00000000000000
        self.SIGZL5 = 0.00000000000000
        self.SIGZL6 = 0.00000000000000
        self.SIGZL7 = 0.00000000000000
        self.SIGZL8 = 0.00000000000000
        self.SIGZL9 = 0.00000000000000
        self.SIGZL10 = 0.00000000000000
        self.CONV_FACT = 0.00000000000000
        self.CONE_MAX = 0.00000000000000
        self.CONE_MIN = 0.00000000000000
        self.EPSI_DX = 0.00000000000000
        self.EPSI_DZ = 0.00000000000000
        self.EPSI_X = 0.00000000000000
        self.EPSI_Z = 0.00000000000000
        self.HDIV1 = 0.500000000000000E-06
        self.HDIV2 = 0.500000000000000E-06
        self.PH1 = 10.0000000000000
        self.PH2 = 1010.00000000000
        self.PH3 = 0.00000000000000
        self.PH4 = 0.00000000000000
        self.PH5 = 0.00000000000000
        self.PH6 = 0.00000000000000
        self.PH7 = 0.00000000000000
        self.PH8 = 0.00000000000000
        self.PH9 = 0.00000000000000
        self.PH10 = 0.00000000000000
        self.RL1 = 0.00000000000000
        self.RL2 = 0.00000000000000
        self.RL3 = 0.00000000000000
        self.RL4 = 0.00000000000000
        self.RL5 = 0.00000000000000
        self.RL6 = 0.00000000000000
        self.RL7 = 0.00000000000000
        self.RL8 = 0.00000000000000
        self.RL9 = 0.00000000000000
        self.RL10 = 0.00000000000000
        self.BENER = 0.00000000000000
        self.POL_ANGLE = 0.00000000000000
        self.POL_DEG = 1.00000000000000
        self.R_ALADDIN = 0.00000000000000
        self.R_MAGNET = 0.00000000000000
        self.SIGDIX = 0.100000000000000E-02
        self.SIGDIZ = 0.100000000000000E-03
        self.SIGMAX = 0.100000000000000E-02
        self.SIGMAY = 0.100000000000000E-02
        self.SIGMAZ = 0.100000000000000E-02
        self.VDIV1 = 0.500000000000000E-05
        self.VDIV2 = 0.500000000000000E-05
        self.WXSOU = 0.100000000000000
        self.WYSOU = 0.200000000000000
        self.WZSOU = 0.200000000000000
        self.PLASMA_ANGLE = 0.00000000000000
        self.FILE_TRAJ = ""
        self.FILE_SOURCE = ""
        self.FILE_BOUND = ""
        self.OE_NUMBER = 0
        self.NTOTALPOINT = 10000000
        self.IDUMMY = 0
        self.DUMMY = 0.00000000000000
        self.F_NEW = 0


    def load(self,filename="start.00"):
        self.load_start00(filename=filename)

    def load_start00(self,filename="start.00"):
        a = GFile()
        a.load_gfile((filename))
        dict1 = a.get_as_dictionary()

        for key in dict1.keys():
            # print("assigning: ",key,dict1[key])
            if hasattr(self,key):
                setattr(self, key, dict1[key])


if __name__ == "__main__":

    oe0 = Source()

    oe0.load(filename="start.00")

    print(dir(oe0))

    assert(oe0.FDISTR == 3)




