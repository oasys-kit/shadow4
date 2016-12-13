
import numpy
import scipy.constants as codata
import json
import os

SHADOW3_BINARY = "/users/srio/OASYS_VE/shadow3/shadow3"


angstroms_to_eV = codata.h*codata.c/codata.e*1e10

def read_json(myinput="xshundul.json"):
    #
    # read inputs from a file created by ShadowVUI ----------------------------
    #
    if isinstance(myinput,str):
        inFileTxt = myinput # "xshundul.json"
        with open(inFileTxt, mode='r') as f1:
            h = json.load(f1)

    elif isinstance(myinput,dict):
        h = myinput
    else:
        raise Exception("Unknown input")

    return h


def info_dictionary(h):

    # list all non-empty keywords
    txt = "-----------------------------------------------------\n"
    for i,j in h.items():
        if (j != None):
            txt += "%s = %s\n" %(i,j)
    txt += "-----------------------------------------------------\n"
    txt += "k: %f \n"%(h['K'])
    return txt

def epath_shadow3(jsn):

    input = "epath\n2\n%f \n%f \n%f \n%d \n1.\nxshundul.par\nxshundul.traj\n1\nxshundul.plt\nexit\n"% \
    (jsn["LAMBDAU"],jsn["K"],jsn["E_ENERGY"],jsn["NG_E"])

    filename = "shadow3_epath.inp"
    f = open(filename,'w')
    f.write(input)
    f.close()

    os.system(SHADOW3_BINARY+" < "+filename)


def epath_compare():
    a = numpy.loadtxt("xshundul.plt").T
    # TTT = ETOFZ(I)/C  # X/c  Y=0
    # BBB = 1.0D0 - EBETAZ(I)
    # WRITE(33,1010) EXOFZ(I), EBETAX(I), EZ(I), BBB, TTT
    # X, BetaX, Z, 1-betaZ, X/c
    from srxraylib.plot.gol import plot
    plot(a[2],a[0],xtitle="Z(along)",ytitle="X(Horizontal)",title="trajectory at the edges??")

def undul_set_shadow3(jsn):

    input = "undul_set\n0\n0\n%d \n%d \n%d \nxshundul.traj\n%d\n%f\n%f\n%f\n%f\n0\n1000\nexit\n"% \
    (jsn["NG_E"],jsn["NG_T"],jsn["NG_P"],
     jsn["NPERIODS"],jsn["EMIN"],jsn["EMAX"],jsn["INTENSITY"],jsn["MAXANGLE"])
    filename = "shadow3_undul_set.inp"
    f = open(filename,'w')
    f.write(input)
    f.close()

    os.system(SHADOW3_BINARY+" < "+filename)

def undul_phot_shadow3():
    input = "undul_phot\nexit\n"
    filename = "shadow3_undul_phot.inp"
    f = open(filename,'w')
    f.write(input)
    f.close()

    os.system(SHADOW3_BINARY+" < "+filename)


def undul_cdf_shadow3():

    input = "undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n"
    filename = "shadow3_undul_cdf.inp"
    f = open(filename,'w')
    f.write(input)
    f.close()

    os.system(SHADOW3_BINARY+" < "+filename)

def input_source_shadow3(jsn):

    input = """input_source
1
0
15000
36255
0
2
xshundul.sha
0.039999999
0.0010000000
4.0000000e-07
0
3.9999999e-09
0
3
1
1
source
systemfile
exit"""

    input = "input source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nsource\nsystemfile\nexit\n"% \
    (jsn["NRAYS"],jsn["SEED"],jsn["SX"],jsn["SZ"],jsn["EX"],0,jsn["EZ"],0,3,1,1)
    filename = "shadow3_input_source.inp"
    f = open(filename,'w')
    f.write(input)
    f.close()

    os.system(SHADOW3_BINARY+" < "+filename)


if __name__ == "__main__":
    jsn = read_json("xshundul.json")
    os.system("rm -f xshundul.plt xshundul.par xshundul.traj")
    print(info_dictionary(jsn))

    epath_shadow3(jsn)
    # epath_compare()

    undul_set_shadow3(jsn)
    undul_phot_shadow3()
    undul_cdf_shadow3()
    input_source_shadow3(jsn)