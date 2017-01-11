#
# load/write files for the undul_phot and undul_cdf shadow3/undulator preprocessors
#


import numpy
from srxraylib.plot.gol import plot,plot_image,plot_show

def load_file_undul_phot(file_in="uphot.dat",do_plot=False,show=False,verbose=False):
    """
    read uphot.dat file (like in SHADOW undul_phot_dump)

    :param file_in: name of the file to be read
    :param do_plot: flag to plot
    :param verbose: flag to print some info
    :return: a dictionary {'radiation':RN0, 'polarization':POL_DEG, 'photon_energy':E, 'theta':TT, 'phi':PP}
    """

    f = open(file_in,'r')
    firstline = f.readline()
    f.close()

    NG_E,NG_T,NG_P = numpy.fromstring(firstline,dtype=int,sep=" ")
    if verbose:
        print("\nload_uphot_dot_dat: file is %s"%(file_in))
        print("load_uphot_dot_dat: NG_E,NG_T,NG_P  %d  %d  %d"%(NG_E,NG_T,NG_P ))

    tmp = numpy.loadtxt(file_in,skiprows=1)

    if tmp.size != 3*(NG_E*NG_T*NG_P)+NG_E+NG_E*NG_T:
        raise Exception("load_uphot_dot_dat: File not understood")

    E = numpy.zeros(NG_E)
    T = numpy.zeros((NG_E,NG_T))
    P = numpy.zeros((NG_E,NG_T,NG_P))


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
    if verbose:
        if NG_E > 1: print("load_uphot_dot_dat: Step in E: %f, Interval in E: %f"%( (E[1]-E[0]), (E[-1]-E[0]) ))
        print("load_uphot_dot_dat: Step in P: %f, Interval in P: %f"%( (PP[1]-PP[0]), (PP[-1]-PP[0]) ))
        print("load_uphot_dot_dat: Step in T: %f, Interval in T: %f"%( (TT[1]-TT[0]), (TT[-1]-TT[0]) ))

        print("load_uphot_dot_dat: RN0 max: %f min: %f"%(RN0.max(),RN0.min()) )
        print("load_uphot_dot_dat: POL_DEG max: %f min: %f"%(POL_DEG.max(),POL_DEG.min()) )

    if do_plot:
        plot_image(RN0[0,:,:],TT*1e6,PP*180/numpy.pi,title=file_in+" RN0[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
        plot_image(POL_DEG[0,:,:],TT*1e6,PP*180/numpy.pi,title=file_in+" POL_DEG[0]",xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)
        if show: plot_show()
    return {'radiation':RN0, 'polarization':POL_DEG, 'photon_energy':E, 'theta':TT, 'phi':PP}



def write_file_undul_phot(undul_phot_dict,file_out="uphot.dat"):

    Z2      = undul_phot_dict['radiation']
    POL_DEG = undul_phot_dict['polarization']
    e       = undul_phot_dict['photon_energy']
    theta   = undul_phot_dict['theta']
    phi     = undul_phot_dict['phi']

    NG_E = e.size
    NG_T = theta.size
    NG_P = phi.size

    f = open(file_out,'w')
    f.write("%d  %d  %d \n"%(NG_E,NG_T,NG_P))
    for ie in range(NG_E):
        f.write("%20.10f \n"%(e[ie]))

    for ie in range(NG_E):
        for t in range(NG_T):
            f.write("%20.10f \n"%(theta[t]))

    for ie in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                f.write("%20.10f \n"%(phi[p]))


    for ie in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                f.write("%20.10f \n"%Z2[ie,t,p])

    for ie in range(NG_E):
        for t in range(NG_T):
            for p in range(NG_P):
                f.write("%20.10f \n"%(POL_DEG[ie,t,p]))

    f.close()
    print("File written to disk: %s"%file_out)




def load_fule_undul_cdf(file_in="xshundul.sha",do_plot=False,show=False,verbose=True):
    #
    # read uphot.dat file (like in SHADOW undul_phot_dump)
    #
    f = open(file_in,'r')
    firstline = f.readline()
    f.close()

    NG_E,NG_T,NG_P, IANGLE = numpy.fromstring(firstline,dtype=int,sep=" ")
    if verbose: print("NG_E,NG_T,NG_P, IANGLE  %d  %d  %d %d \n"%(NG_E,NG_T,NG_P,IANGLE ))

    tmp = numpy.loadtxt(file_in,skiprows=1)

    if tmp.size != 2*(NG_E + NG_E*NG_T + NG_E*NG_T*NG_P) + NG_E*NG_T*NG_P:
        raise Exception("File not understood")

    E = numpy.zeros(NG_E)
    T = numpy.zeros((NG_E,NG_T))
    P = numpy.zeros((NG_E,NG_T,NG_P))


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
        plot(E,TWO,title="TWO %s"%file_in,xtitle="E",ytitle="TWO",show=False)
        plot_image(ONE,numpy.arange(NG_E),numpy.arange(NG_T),title="ONE %s "%file_in,xtitle="index Energy",ytitle="index Theta",aspect='auto',show=False)
        plot_image(ZERO[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="ZERO[0] %s"%file_in,xtitle="index Theta",ytitle="index Phi",aspect='auto',show=False)
        # plot_image(POL_DEGREE[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),title="POL_DEGREE[0]",xtitle="index Theta",ytitle="index Phi",show=0)
        if show: plot_show()

    return {'cdf_EnergyThetaPhi':TWO,'cdf_EnergyTheta':ONE,'cdf_Energy':ZERO,'energy':E,'theta':T,'phi':P,'polarization':POL_DEGREE}

def write_file_undul_sha(dict,file_out="xshundul.sha"):
    #
    # create xshundul.sha file (like in SHADOW undul_cdf)
    #
    TWO =     dict['cdf_EnergyThetaPhi']
    ONE =     dict['cdf_EnergyTheta']
    ZERO =    dict['cdf_Energy']
    E =       dict['energy']
    T =       dict['theta']
    P =       dict['phi']
    POL_DEG = dict['polarization']

    NG_E = E.size
    NG_T = T.size
    NG_P = P.size

    if file_out is not None:
        f = open(file_out,'w')
        f.write("%d  %d  %d 1 \n"%(NG_E,NG_T,NG_P))

        for e in E:
            f.write("%20.10f \n"%(e))

        for e in E:
            for t in T:
                f.write("%20.10f \n"%t)

        for e in E:
            for t in T:
                for p in P:
                    f.write("%20.10f \n"%p)

        for e in numpy.arange(NG_E):
            f.write("%20.10f \n"%(TWO[e]))

        for e in numpy.arange(NG_E):
            for t in numpy.arange(NG_T):
                f.write("%20.10f \n"%(ONE[e,t]))

        for e in numpy.arange(NG_E):
            for t in numpy.arange(NG_T):
                for p in numpy.arange(NG_P):
                    f.write("%20.10f \n"%(ZERO[e,t,p]))

        for e in numpy.arange(NG_E):
            for t in numpy.arange(NG_T):
                for p in numpy.arange(NG_P):
                    f.write("%20.10f \n"%(POL_DEG[e,t,p]))


        f.close()
        print("File written to disk: %s"%file_out)



# def compare_shadow3_files(file1,file2,do_assert=True):
#
#     Shadow.ShadowTools.plotxy(file1,4,6,nbins=101,nolost=1,title=file1)
#     Shadow.ShadowTools.plotxy(file2,4,6,nbins=101,nolost=1,title=file2)
#
#     if do_assert:
#         begin1 = Shadow.Beam()
#         begin1.load(file1)
#         begin2 = Shadow.Beam()
#         begin2.load(file2)
#         assert_almost_equal(begin1.rays[:,0:6],begin2.rays[:,0:6],3)
