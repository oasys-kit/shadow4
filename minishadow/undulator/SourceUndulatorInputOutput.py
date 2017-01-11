#
# load/write files for the undul_phot and undul_cdf shadow3/undulator preprocessors
#


import numpy


def load_file_undul_phot(file_in="uphot.dat"):
    """
    read uphot.dat file (like in SHADOW undul_phot_dump)

    :param file_in: name of the file to be read
    :return: a dictionary {'radiation':RN0, 'polarization':POL_DEG, 'photon_energy':E, 'theta':TT, 'phi':PP}
    """

    f = open(file_in,'r')
    firstline = f.readline()
    f.close()

    NG_E,NG_T,NG_P = numpy.fromstring(firstline,dtype=int,sep=" ")


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

    return {'radiation':RN0, 'polarization':POL_DEG, 'photon_energy':E, 'theta':TT, 'phi':PP}



def write_file_undul_phot(undul_phot_dict,file_out="uphot.dat"):
    """
    write uphot.dat file from a dictionary with the output of undul_phot
    :param undul_phot_dict: a dictionary {'radiation':RN0, 'polarization':POL_DEG, 'photon_energy':E, 'theta':TT, 'phi':PP}
    :param file_out: name of the output file (Default: uphot.dat)
    :return:
    """

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



def load_file_undul_cdf(file_in="xshundul.sha"):
    """
    Loads a file containing the output of undul_cdf into a dictionary

    :param file_in: Name of the input file (Default: xshundul.sha)
    :return: a dictionary {'cdf_EnergyThetaPhi':TWO,'cdf_EnergyTheta':ONE,'cdf_Energy':ZERO,
            'energy':E,'theta':T,'phi':P,'polarization':POL_DEGREE}
    """

    f = open(file_in,'r')
    firstline = f.readline()
    f.close()

    NG_E,NG_T,NG_P, IANGLE = numpy.fromstring(firstline,dtype=int,sep=" ")


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

    return {'cdf_EnergyThetaPhi':TWO,'cdf_EnergyTheta':ONE,'cdf_Energy':ZERO,'energy':E,'theta':T,'phi':P,'polarization':POL_DEGREE}

def write_file_undul_sha(dict,file_out="xshundul.sha"):
    """
    Create a file (xshundul.sha) with output of undul_cdf

    :param dict: a dictionary as output from undul_cdf
    :param file_out: output file name
    :return:
    """
    #
    #
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

# TODO do these plot directly with matplotlib to avoid dependencies
def plot_undul_cdf(undul_cdf_input,do_show=True):
    #
    # plots output of undul_cdf
    #
    try:
        from srxraylib.plot.gol import plot,plot_image,plot_show
    except:
        print("srxraylib not available: No plot")

    if isinstance(undul_cdf_input,str):
        undul_cdf_dict = load_file_undul_cdf(undul_cdf_input)
    else:
        undul_cdf_dict = undul_cdf_input

    TWO = undul_cdf_dict['cdf_EnergyThetaPhi']
    ONE = undul_cdf_dict['cdf_EnergyTheta']
    ZERO = undul_cdf_dict['cdf_Energy']
    # E = undul_cdf_dict['energy']
    # T = undul_cdf_dict['theta']
    # P = undul_cdf_dict['phi']


    NG_E,NG_T,NG_P = ZERO.shape

    plot(numpy.arange(NG_E),TWO,title="cdf(energy) TWO",xtitle="index Energy",ytitle="cdf(E) TWO",show=0)
    plot_image(ONE,numpy.arange(NG_E),numpy.arange(NG_T),aspect='auto',
               title="cdf(energy,theta) ONE",xtitle="index Energy",ytitle="index Theta",show=0)
    plot_image(ZERO[0,:,:],numpy.arange(NG_T),numpy.arange(NG_P),aspect='auto',
               title="cdf (theta,phi) ZERO[0]",xtitle="index Theta",ytitle="index Phi",show=0)

    if do_show: plot_show()


def plot_undul_phot(undul_phot_input,do_plot_intensity=True,do_plot_polarization=True,do_show=True,title=""):
    #
    # plots the output of undul_phot
    #
    try:
        from srxraylib.plot.gol import plot,plot_image,plot_show
    except:
        print("srxraylib not available: No plot")

    if isinstance(undul_phot_input,str):
        undul_phot_dict = load_file_undul_phot(undul_phot_input)
        title += undul_phot_input
    else:
        undul_phot_dict = undul_phot_input


    if do_plot_intensity: plot_image(undul_phot_dict['radiation'][0,:,:],undul_phot_dict['theta']*1e6,undul_phot_dict['phi']*180/numpy.pi,
               title="INTENS RN0[0] "+title,xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)

    if do_plot_polarization: plot_image(undul_phot_dict['polarization'][0,:,:],undul_phot_dict['theta']*1e6,undul_phot_dict['phi']*180/numpy.pi,
               title="POL_DEG RN0[0] "+title,xtitle="Theta [urad]",ytitle="Phi [deg]",aspect='auto',show=False)


    if do_show: plot_show()

