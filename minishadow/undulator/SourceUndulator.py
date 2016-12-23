
import numpy
import scipy.constants as codata
import json
import os

from numpy.testing import assert_equal, assert_almost_equal
from srxraylib.plot.gol import plot,plot_image,plot_show
import scipy.integrate

import Shadow

import platform
if platform.system() == "Linux":
    SHADOW3_BINARY = "/users/srio/OASYS_VE/shadow3/shadow3"
else:
    SHADOW3_BINARY = "/Users/srio/Oasys/OASYS_VE/shadow3/shadow3"

angstroms_to_eV = codata.h*codata.c/codata.e*1e10

#
# read inputs
#
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


#
# run full sequence
#
def undulator(code='shadow3',force_shadow3_undul_cdf=False,force_srw_undul_phot=False):
    jsn = read_json("xshundul.json")
    os.system("rm -f xshundul.plt xshundul.par xshundul.traj xshundul.info xshundul.sha")
    print(info_dictionary(jsn))

    if code == 'shadow3':

        # epath
        commands = "epath\n2\n%f \n%f \n%f \n%d \n1.\nxshundul.par\nxshundul.traj\n1\nxshundul.plt\nexit\n"% \
        (jsn["LAMBDAU"],jsn["K"],jsn["E_ENERGY"],jsn["NG_E"])
        shadow3_commands(commands=commands,input_file="shadow3_epath.inp")

        # undul_set
        commands = "undul_set\n0\n0\n%d \n%d \n%d \nxshundul.traj\n%d\n%f\n%f\n%f\n%f\n0\n1000\nexit\n"% \
            (jsn["NG_E"],jsn["NG_T"],jsn["NG_P"],
             jsn["NPERIODS"],jsn["EMIN"],jsn["EMAX"],jsn["INTENSITY"],jsn["MAXANGLE"])
        shadow3_commands(commands=commands,input_file="shadow3_undul_set.inp")

        # undul_phot
        shadow3_commands(commands="undul_phot\nexit\n",input_file="shadow3_undul_phot.inp")
        shadow3_commands(commands="undul_phot_dump\nexit\n",input_file="shadow3_undul_phot_dump.inp")

        # undul_cdf
        shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
                        input_file="shadow3_undul_cdf.inp")


        # input source
        commands = "input source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
        (jsn["NRAYS"],jsn["SEED"],jsn["SX"],jsn["SZ"],jsn["EX"],0,jsn["EZ"],0,3,1,1)


        shadow3_commands(commands=commands,input_file="shadow3_input_source.inp")

        # run source
        commands = "source\nsystemfile\nexit\n"
        shadow3_commands(commands=commands,input_file="shadow3_source.inp")

        os.system("cp begin.dat begin_shadow3.dat")
        os.system("cp uphot.dat uphot_shadow3.dat")
        os.system("cp xshundul.sha xshundul_shadow3.sha")

        # make script
        oe0 = Shadow.Source()
        oe0.load("start.00")
        Shadow.ShadowTools.make_python_script_from_list([oe0],script_file="script_undulator.py")

    elif code == 'minishadow':


        # epath
        commands = "epath\n2\n%f \n%f \n%f \n%d \n1.\nxshundul.par\nxshundul.traj\n1\nxshundul.plt\nexit\n"% \
        (jsn["LAMBDAU"],jsn["K"],jsn["E_ENERGY"],jsn["NG_E"])
        shadow3_commands(commands=commands,input_file="shadow3_epath.inp")


        # undul_phot
        if force_srw_undul_phot:
            from undul_phot_srw import undul_phot_srw
            undul_phot_srw(jsn)
        else:
            from undul_phot import undul_phot
            undul_phot(jsn)
        shadow3_commands(commands="undul_phot_dump\nexit\n",input_file="shadow3_undul_phot_dump.inp")

        # undul_cdf
        undul_cdf(file_in="uphot.dat",method='trapz',file_out="xshundul.sha",do_plot=False)


        # input source
        use_shadow3_binary = False
        if use_shadow3_binary:
            commands = "input source\n1\n0\n%d \n%d \n0 \n2 \nxshundul.sha\n%g\n%g\n%g\n%d\n%g\n%d\n%d\n%d\n%d\nexit\n"% \
            (jsn["NRAYS"],jsn["SEED"],jsn["SX"],jsn["SZ"],jsn["EX"],0,jsn["EZ"],0,3,1,1)
            shadow3_commands(commands=commands,input_file="shadow3_input_source.inp")
            # run source
            commands = "source\nsystemfile\nexit\n"
            shadow3_commands(commands=commands,input_file="shadow3_source.inp")
        else:
            iwrite = 1
            # initialize shadow3 source (oe0) and beam
            oe0 = Shadow.Source()
            beam = Shadow.Beam()
            oe0.EPSI_X = jsn["EX"]
            oe0.EPSI_Z = jsn["EZ"]
            oe0.SIGDIX = 0.0
            oe0.SIGDIZ = 0.0
            oe0.SIGMAX = jsn["SX"]
            oe0.SIGMAY = 0.0
            oe0.SIGMAZ = jsn["SZ"]
            oe0.FILE_TRAJ = b'xshundul.sha'
            oe0.ISTAR1 = jsn["SEED"]
            oe0.NPOINT = jsn["NRAYS"]
            oe0.F_WIGGLER = 2
            if iwrite:
                oe0.write("start.00")
            beam.genSource(oe0)
            if iwrite:
                oe0.write("end.00")
                beam.write("begin.dat")

        os.system("cp begin.dat begin_minishadow.dat")
        os.system("cp uphot.dat uphot_minishadow.dat")
        os.system("cp xshundul.sha xshundul_minishadow.sha")

def shadow3_commands(commands="exit\n",input_file="shadow3_tmp.inp"):
    f = open(input_file,'w')
    f.write(commands)
    f.close()
    os.system(SHADOW3_BINARY+" < "+input_file)


#
# load files
#
def load_uphot_dot_dat(file_in="uphot.dat",do_plot=False):
    #
    # read uphot.dat file (like in SHADOW undul_phot_dump)
    #
    f = open(file_in,'r')
    firstline = f.readline()
    f.close()

    NG_E,NG_T,NG_P = numpy.fromstring(firstline,dtype=int,sep=" ")
    print("\nload_uphot_dot_dat: file is %s"%(file_in))
    print("load_uphot_dot_dat: NG_E,NG_T,NG_P  %d  %d  %d"%(NG_E,NG_T,NG_P ))

    tmp = numpy.loadtxt(file_in,skiprows=1)

    if tmp.size != 3*(NG_E*NG_T*NG_P)+NG_E+NG_E*NG_T:
        raise Exception("load_uphot_dot_dat: File not understood")

    E = numpy.zeros(NG_E)
    T = numpy.zeros((NG_E,NG_T))
    P = numpy.zeros((NG_E,NG_T,NG_P))


    # for e in E:
    #     f.write("%g \n"%(e))
    #
    itmp = 0
    for ie in range(NG_E):
        E[ie] = tmp[itmp]
        itmp += 1

    for ie in range(NG_E):
        for it in range(NG_T):
            T[ie,it] = tmp[itmp]
            itmp += 1

    for ie in range(NG_E):
        for it in range(NG_T):
            for ip in range(NG_P):
                P[ie,it,ip] = tmp[itmp]
                itmp += 1

    RN0 = numpy.zeros((NG_E,NG_T,NG_P))
    POL_DEG = numpy.zeros((NG_E,NG_T,NG_P))

    for e in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                RN0[e,t,p] = tmp[itmp]
                itmp += 1

    for e in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                POL_DEG[e,t,p] = tmp[itmp]
                itmp += 1

    TT = T.flatten()[0:NG_T].copy()
    PP = P.flatten()[0:NG_P].copy()
    print("load_uphot_dot_dat: Step in E: %f, Interval in E: %f"%( (E[1]-E[0]), (E[-1]-E[0]) ))
    print("load_uphot_dot_dat: Step in P: %f, Interval in P: %f"%( (PP[1]-PP[0]), (PP[-1]-PP[0]) ))
    print("load_uphot_dot_dat: Step in T: %f, Interval in T: %f"%( (TT[1]-TT[0]), (TT[-1]-TT[0]) ))

    print("load_uphot_dot_dat: RN0 max: %f min: %f"%(RN0.max(),RN0.min()) )
    print("load_uphot_dot_dat: POL_DEG max: %f min: %f"%(POL_DEG.max(),POL_DEG.min()) )

    if do_plot:
        # with abscissas
        # plot_image(RN0[0,:,:],T.flatten()[0:NG_T].copy()*1e6,P.flatten()[0:NG_P].copy(),title="RN0[0]",xtitle="Theta [urad]",ytitle="Phi [rad]",aspect='auto',show=0)
        # plot_image(POL_DEG[0,:,:],T.flatten()[0:NG_T].copy()*1e6,P.flatten()[0:NG_P].copy(),title="POL_DEG[0]",xtitle="Theta [urad]",ytitle="Phi [rad]",aspect='auto',show=1)

        # with indices
        plot_image(RN0[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title=file_in+" RN0[0]",xtitle="Theta [index]",ytitle="Phi [index]",aspect=None,show=0)
        plot_image(POL_DEG[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title=file_in+" POL_DEG[0]",xtitle="Theta [index]",ytitle="Phi [index]",aspect=None,show=0)


    return RN0, POL_DEG, E, TT, PP


def load_xshundun_dot_sha(file_in="xshundul.sha",do_plot=False):
    #
    # read uphot.dat file (like in SHADOW undul_phot_dump)
    #
    f = open(file_in,'r')
    firstline = f.readline()
    f.close()

    NG_E,NG_T,NG_P, IANGLE = numpy.fromstring(firstline,dtype=int,sep=" ")
    print("NG_E,NG_T,NG_P, IANGLE  %d  %d  %d %d \n"%(NG_E,NG_T,NG_P,IANGLE ))

    tmp = numpy.loadtxt(file_in,skiprows=1)

    if tmp.size != 2*(NG_E + NG_E*NG_T + NG_E*NG_T*NG_P) + NG_E*NG_T*NG_P:
        raise Exception("File not understood")

    E = numpy.zeros(NG_E)
    T = numpy.zeros((NG_E,NG_T))
    P = numpy.zeros((NG_E,NG_T,NG_P))


    # for e in E:
    #     f.write("%g \n"%(e))
    #
    itmp = 0
    for ie,e in enumerate(E):
        E[ie] = tmp[itmp]
        itmp += 1

    for ie in range(NG_E):
        for it in range(NG_T):
            T[ie,it] = tmp[itmp]
            itmp += 1

    for ie in range(NG_E):
        for it in range(NG_T):
            for ip in range(NG_P):
                P[ie,it,ip] = tmp[itmp]
                itmp += 1

    TWO = numpy.zeros((NG_E))
    ONE = numpy.zeros((NG_E,NG_T))
    ZERO = numpy.zeros((NG_E,NG_T,NG_P))
    POL_DEGREE = numpy.zeros_like(ZERO)


    for e in range(NG_E):
        TWO[e] = tmp[itmp]
        itmp += 1

    for e in range(NG_E):
        for t in range(NG_T):
            ONE[e,t] = tmp[itmp]
            itmp += 1

    for e in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                ZERO[e,t,p] = tmp[itmp]
                itmp += 1

    for e in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                POL_DEGREE[e,t,p] = tmp[itmp]
                itmp += 1

    if do_plot:
        plot(E,TWO,title="TWO %s"%file_in,xtitle="E",ytitle="TWO",show=0)
        plot_image(ONE,numpy.arange(NG_E),numpy.arange(NG_T),title="ONE %s "%file_in,xtitle="index Energy",ytitle="index Theta",show=0)
        plot_image(ZERO[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="ZERO[0] %s"%file_in,xtitle="index Theta",ytitle="index Phi",show=0)
        # plot_image(POL_DEGREE[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="POL_DEGREE[0]",xtitle="index Theta",ytitle="index Phi",show=0)


    return TWO,ONE,ZERO,E,T,P,POL_DEGREE

#
# undul_cdf
#
def undul_cdf(file_in="uphot.dat",file_out=None,method='trapz',do_plot=False):
    #
    # read uphot.dat file (like in SHADOW undul_phot_dump)
    #
    RN0,POL_DEG,E,T,P = load_uphot_dot_dat(file_in=file_in,do_plot=do_plot)


    NG_E,NG_T,NG_P = RN0.shape
    print("undul_cdf: NG_E,NG_T,NG_P, %d  %d %d \n"%(NG_E,NG_T,NG_P))


    # corrdinates are polar: multiply by sin(theta) to allow dS= r^2 sin(Theta) dTheta dPhi
    YRN0 = numpy.zeros_like(RN0)
    for e in numpy.arange(NG_E):
        for t in numpy.arange(NG_T):
            for p in numpy.arange(NG_P):
                YRN0[e,t,p] = RN0[e,t,p] * numpy.sin(T[t])


    if method == "sum":
        RN1 = YRN0.sum(axis=2) * (P[1] - P[0])             # RN1(e,t)
        RN2 = RN1.sum(axis=1)  * (T[1] - T[0])             # RN2(e)
        ZERO  = numpy.cumsum(RN0,axis=2)   * (P[1] - P[0]) # CDF(e,t,p)
        ONE   = numpy.cumsum(RN1,axis=1)   * (T[1] - T[0]) # CDF(e,t)
        TWO   = numpy.cumsum(RN2)          * (E[1] - E[0]) # CDF(e)
    else:
        RN1 = numpy.trapz(YRN0,axis=2) * (P[1]-P[0])                            # RN1(e,t)
        RN2 = numpy.trapz(RN1,axis=1)  * (T[1]-T[0])                            # RN2(e)
        ZERO  = scipy.integrate.cumtrapz(RN0,initial=0,axis=2)  * (P[1] - P[0]) # CDF(e,t,p)
        ONE   = scipy.integrate.cumtrapz(RN1,initial=0,axis=1)  * (T[1] - T[0]) # CDF(e,t)
        TWO   = scipy.integrate.cumtrapz(RN2,initial=0)         * (E[1] - E[0]) # CDF(e)



    print("undul_cdf: Shadow ZERO,ONE,TWO: ",ZERO.shape,ONE.shape,TWO.shape)
    print("undul_cdf: Total Power emitted in the specified angles is: %g Watts."%( (RN2*E).sum()*(E[1]-E[0])*codata.e) )


    if do_plot:
        plot(E,TWO,title="NEW TWO",xtitle="E",ytitle="TWO",show=0)
        plot_image(ONE,numpy.arange(NG_E),numpy.arange(NG_T),title="NEW ONE",xtitle="index Energy",ytitle="index Theta",show=0)
        plot_image(ZERO[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="NEW ZERO[0]",xtitle="index Theta",ytitle="index Phi",show=0)
        #plot_image(POL_DEGREE[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="POL_DEGREE[0]",xtitle="index Theta",ytitle="index Phi",show=0)



    #
    # create xshundul.sha file (like in SHADOW undul_cdf)
    #
    if file_out is not None:
        f = open(file_out,'w')
        f.write("%d  %d  %d 1 \n"%(NG_E,NG_T,NG_P))

        for e in E:
            f.write("%g \n"%(e))

        for e in E:
            for t in T:
                f.write("%g \n"%t)

        for e in E:
            for t in T:
                for p in P:
                    f.write("%g \n"%p)

        for e in numpy.arange(NG_E):
            f.write("%g \n"%(TWO[e]))

        for e in numpy.arange(NG_E):
            for t in numpy.arange(NG_T):
                f.write("%g \n"%(ONE[e,t]))

        for e in numpy.arange(NG_E):
            for t in numpy.arange(NG_T):
                for p in numpy.arange(NG_P):
                    f.write("%g \n"%(ZERO[e,t,p]))

        for e in numpy.arange(NG_E):
            for t in numpy.arange(NG_T):
                for p in numpy.arange(NG_P):
                    f.write("%g \n"%(POL_DEG[e,t,p]))


        f.close()
        print("undul_cdf: File written to disk: %s"%file_out)




    return TWO,ONE,ZERO,E,T,P



#
# comparisons/tests
#
def epath_compare():
    a = numpy.loadtxt("xshundul.plt").T
    # TTT = ETOFZ(I)/C  # X/c  Y=0
    # BBB = 1.0D0 - EBETAZ(I)
    # WRITE(33,1010) EXOFZ(I), EBETAX(I), EZ(I), BBB, TTT
    # X, BetaX, Z, 1-betaZ, X/c
    plot(a[2],a[0],xtitle="Z(along)",ytitle="X(Horizontal)",title="trajectory at the edges??")


def compare_shadow3_files(file1,file2,do_assert=True):

    Shadow.ShadowTools.plotxy(file1,4,6,nbins=101,nolost=1,title=file1)
    Shadow.ShadowTools.plotxy(file2,4,6,nbins=101,nolost=1,title=file2)

    if do_assert:
        begin1 = Shadow.Beam()
        begin1.load(file1)
        begin2     = Shadow.Beam()
        begin2.load(file2)
        assert_almost_equal(begin1.rays[:,0:6],begin2.rays[:,0:6],3)

def test_undul_cdf(do_plot=True):
    #
    # undul_phot.dat must exist
    #
    shadow3_commands(commands="undul_cdf\n0\n1\nxshundul.sha\nxshundul.info\nexit\n",
                    input_file="shadow3_undul_cdf.inp")

    undul_cdf(file_in="uphot.dat",method='sum',  file_out="xshundul2.sha",do_plot=False)
    undul_cdf(file_in="uphot.dat",method='trapz',file_out="xshundul3.sha",do_plot=False)

    TWO1,ONE1,ZERO1,E1,T1,P1,POL_DEGREE1 = load_xshundun_dot_sha(file_in="xshundul.sha", do_plot=do_plot)
    TWO2,ONE2,ZERO2,E2,T2,P2,POL_DEGREE2 = load_xshundun_dot_sha(file_in="xshundul2.sha",do_plot=do_plot)
    TWO3,ONE3,ZERO3,E3,T3,P3,POL_DEGREE3 = load_xshundun_dot_sha(file_in="xshundul3.sha",do_plot=do_plot)


    tmp = numpy.where(ZERO1 > 0.1*ZERO1.max())
    print("test_undul_cdf: ZERO:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ZERO2[tmp]-ZERO1[tmp])/ZERO1[tmp]) )))
    print("test_undul_cdf: ZERO: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ZERO3[tmp]-ZERO1[tmp])/ZERO1[tmp]) )))

    tmp = numpy.where(ONE1 > 0.1*ONE1.max())
    print(r"test_undul_cdf: ONE:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ONE2[tmp]-ONE1[tmp])/ONE1[tmp]) )))
    print(r"test_undul_cdf: ONE: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((ONE3[tmp]-ONE1[tmp])/ONE1[tmp]) )))

    tmp = numpy.where(TWO1 > 0.1*TWO1.max())
    print("test_undul_cdf: TWO:   sum/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((TWO2[tmp]-TWO1[tmp])/TWO1[tmp]) )))
    print("test_undul_cdf: TWO: trapz/shadow3 %4.2f %%: "%(numpy.average( 100*numpy.abs((TWO3[tmp]-TWO1[tmp])/TWO1[tmp]) )))



if __name__ == "__main__":

    # test_undul_cdf(do_plot=False)


    os.system("rm -f begin*.dat uphot*.dat")
    undulator(code='shadow3')

    # TODO: check polarization
    undulator(code='minishadow',force_shadow3_undul_cdf=False,force_srw_undul_phot=False)




    load_uphot_dot_dat("uphot_shadow3.dat",do_plot=True)
    load_uphot_dot_dat("uphot_minishadow.dat",do_plot=True)
    # load_uphot_dot_dat("uphot.dat",do_plot=True)

    compare_shadow3_files("begin_shadow3.dat","begin_minishadow.dat",do_assert=True)

    # uphot.dat

    # RN0,POL_DEG,E,T,P = load_uphot_dot_dat(do_plot=True)
    # print("RN0,POL_DEG",RN0.shape,POL_DEG.shape)
    # print("E",E)
    # print("T",T)
    # print("P",P)
    #


    # xshundul.sha
    # TWO,ONE,ZERO,E,T,P,POL_DEGREE = load_xshundun_dot_sha(do_plot=True)
    # print("Shadow TWO",TWO.shape)
    # print("Shadow ONE",ONE.shape)
    # print("Shadow ZERO",ZERO.shape)
    #
    # print("Total TWO",  TWO.sum())
    # print("Total ONE",  ONE.sum())
    # print("Total ZERO",ZERO.sum())
    # NG_E,NG_T,NG_P = ZERO.shape

    # new undul_cdf
    # TWO,ONE,ZERO,E,T,P,POL_DEGREE = undul_cdf(do_plot=True)



    plot_show()