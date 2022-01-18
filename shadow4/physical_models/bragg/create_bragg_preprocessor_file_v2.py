
import numpy
import xraylib
import scipy.constants as codata

# needed by bragg_calc
from dabax.common_tools import f0_xop

# needed by bragg_calc
from dabax.common_tools import bragg_metrictensor, atomic_symbols

def create_bragg_preprocessor_file_v2(interactive=True,
        DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1, TEMPERATURE_FACTOR=1.0,
        E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
        SHADOW_FILE="bragg.dat",
        material_constants_library=xraylib):

    """
     SHADOW preprocessor for crystals - python+xraylib version

     -"""

    # codata_e2_mc2 = 2.81794032e-15 = Classical electron radius in S.I.
    codata_e2_mc2 = codata.hbar * codata.alpha / codata.m_e / codata.c

    if interactive:
        print("bragg: SHADOW preprocessor for crystals - python+xraylib version")
        fileout = input("Name of output file : ")

        print(" bragg (python) only works now for ZincBlende Cubic structures. ")
        print(" Valid descriptor are: ")
        print("     Si (alternatively Si_NIST, Si2) ")
        print("     Ge")
        print("     Diamond")
        print("     GaAs, GaSb, GaP")
        print("     InAs, InP, InSb")
        print("     SiC")

        descriptor = input("Name of crystal descriptor : ")

        print("Miller indices of crystal plane of reeflection.")
        miller = input("H K L: ")
        miller = miller.split()
        hh = int(miller[0])
        kk = int(miller[1])
        ll = int(miller[2])

        temper = input("Temperature (Debye-Waller) factor (set 1 for default): ")
        temper = float(temper)

        emin = input("minimum photon energy (eV): ")
        emin = float(emin)
        emax = input("maximum photon energy (eV): ")
        emax = float(emax)
        estep = input("energy step (eV): ")
        estep = float(estep)

    else:
        fileout    = SHADOW_FILE
        descriptor = DESCRIPTOR
        hh         = int(H_MILLER_INDEX)
        kk         = int(K_MILLER_INDEX)
        ll         = int(L_MILLER_INDEX)
        temper     = float(TEMPERATURE_FACTOR)
        emin       = float(E_MIN)
        emax       = float(E_MAX)
        estep      = float(E_STEP)


    #
    # end input section, start calculations
    #

    out_dict = bragg_calc(descriptor=descriptor,
                          hh=hh,kk=kk,ll=ll,temper=temper,
                          emin=emin,emax=emax,estep=estep,fileout=fileout,
                          material_constants_library=material_constants_library)

    # dump_bragg_preprocessor_file_v1(out_dict, fileout=fileout)

    return out_dict


def bragg_calc(descriptor="Si",hh=1,kk=1,ll=1,temper=1.0,emin=5000.0,emax=15000.0,estep=100.0,fileout=None,
               material_constants_library=xraylib):
    """
    Preprocessor for Structure Factor (FH) calculations. It calculates the basic ingredients of FH.

    :param descriptor: crystal name (as in xraylib)
    :param hh: miller index H
    :param kk: miller index K
    :param ll: miller index L
    :param temper: temperature factor (scalar <=1.0 )
    :param emin:  photon energy minimum
    :param emax: photon energy maximum
    :param estep: photon energy step
    :param fileout: name for the output file (default=None, no output file)
    :return: a dictionary with all ingredients of the structure factor.
    """

    output_dictionary = {}

    codata_e2_mc2 = codata.e**2 / codata.m_e / codata.c**2 / (4*numpy.pi*codata.epsilon_0)  # in m

    # f = open(fileout,'w')

    txt = ""
    txt += "# Bragg version, Data file type\n"
    txt += "2.5 1\n"

    cryst = material_constants_library.Crystal_GetCrystal(descriptor)
    volume = cryst['volume']

    #test crystal data - not needed
    itest = 0
    if itest:

        print ("  Unit cell dimensions are %f %f %f" % (cryst['a'],cryst['b'],cryst['c']))
        print ("  Unit cell angles are %f %f %f" % (cryst['alpha'],cryst['beta'],cryst['gamma']))
        print ("  Unit cell volume is %f A^3" % volume )
        print ("  Atoms at:")
        print ("     Z  fraction    X        Y        Z")
        for i in range(cryst['n_atom']):
            atom =  cryst['atom'][i]
            print ("    %3i %f %f %f %f" % (atom['Zatom'], atom['fraction'], atom['x'], atom['y'], atom['z']) )
        print ("  ")

    volume = volume*1e-8*1e-8*1e-8 # in cm^3
    dspacing = material_constants_library.Crystal_dSpacing(cryst, hh, kk, ll)
    rn = (1e0/volume)*(codata_e2_mc2*1e2)
    dspacing *= 1e-8 # in cm

    txt += "# RN = (e^2/(m c^2))/V) [cm^-2], d spacing [cm]\n"
    txt += "%e %e \n" % (rn , dspacing)

    output_dictionary["rn"] = rn
    output_dictionary["dspacing"] = dspacing

    atom = cryst['atom']
    list_Zatom = [ atom[i]['Zatom'] for i in range(len(atom))]
    number_of_atoms = len(list_Zatom)
    list_fraction = [ atom[i]['fraction'] for i in range(len(atom))]
    try:
        list_charge = [atom[i]['charge'] for i in range(len(atom))]
    except:
        list_charge = [0.0] * number_of_atoms
    list_x = [ atom[i]['x'] for i in range(len(atom))]
    list_y = [ atom[i]['y'] for i in range(len(atom))]
    list_z = [ atom[i]['z'] for i in range(len(atom))]

    # creates an is that contains Z, occupation and charge, that will
    # define the different sites.
    IDs = []
    number_of_atoms = len(list_Zatom)
    for i in range(number_of_atoms):
        IDs.append("Z:%2d-F:%g-C:%g" % (list_Zatom[i],list_fraction[i], list_charge[i]))

    # calculate indices of uniqte Id's sorted by Z
    unique_indexes1 = numpy.unique(IDs, return_index=True) [1]
    unique_Zatom1 = [list_Zatom[i] for i in unique_indexes1]
    # sort by Z
    ii = numpy.argsort(unique_Zatom1)
    unique_indexes = unique_indexes1[ii]

    unique_Zatom = [list_Zatom[i] for i in unique_indexes]
    unique_charge = [list_charge[i] for i in unique_indexes]
    unique_scattering_electrons = []
    for i, Zi in enumerate(unique_Zatom):
        unique_scattering_electrons.append(Zi - unique_charge[i])

    nbatom = (len(unique_Zatom))

    txt += "# Number of different element-sites in unit cell NBATOM:\n%d \n" % nbatom
    output_dictionary["nbatom"] = nbatom

    txt += "# for each element-site, the number of scattering electrons (Z_i + charge_i)\n"
    for i in unique_Zatom:
        txt += "%d "%i
    txt += "\n"
    output_dictionary["atnum"] = list(unique_scattering_electrons)

    txt += "# for each element-site, the occupation factor\n"
    unique_fraction = []
    for i in range(len(unique_indexes)):
        unique_fraction.append(list_fraction[unique_indexes[i]])
        txt += "%g "%(unique_fraction[i])
    txt += "\n"
    output_dictionary["fraction"] = unique_fraction


    txt += "# for each element-site, the temperature factor\n" # temperature parameter
    list_temper = []
    for i in range(len(unique_indexes)):
        txt += "%5.3f "%temper
        list_temper.append(temper)
    txt += "\n"
    output_dictionary["temper"] = list_temper

    #
    # Geometrical part of structure factor:  G and G_BAR
    #
    txt += "# for each type of element-site, COOR_NR=G_0\n"
    list_multiplicity = []
    for i in range(len(unique_indexes)):
        id = IDs[unique_indexes[i]]
        txt += "%d "%IDs.count(id)
        list_multiplicity.append(IDs.count(id))
    txt += "\n"
    output_dictionary["G_0"] = list_multiplicity

    txt += "# for each type of element-site, G and G_BAR (both complex)\n"
    list_g = []
    list_g_bar = []
    for i in range(len(unique_indexes)):
        id = IDs[unique_indexes[i]]
        ga = 0.0 + 0j
        for i,zz in enumerate(IDs):
            if zz == id:
                ga += numpy.exp(2j*numpy.pi*(hh*list_x[i]+kk*list_y[i]+ll*list_z[i]))
        txt += "(%g,%g) \n"%(ga.real,ga.imag)
        txt += "(%g,%g) \n"%(ga.real,-ga.imag)
        list_g.append(ga)
        list_g_bar.append(ga.conjugate())
    output_dictionary["G"] = list_g
    output_dictionary["G_BAR"] = list_g_bar

    #
    # F0 part
    #
    txt += "# for each type of element-site, the number of f0 coefficients followed by them\n"
    list_f0 = []
    for i in range(len(unique_indexes)):
        zeta = list_Zatom[unique_indexes[i]]
        tmp = f0_xop(zeta)
        txt += ("11 "+"%g "*11+"\n")%(tuple(tmp))
        list_f0.append(tmp.tolist())
    output_dictionary["f0coeff"] = list_f0


    npoint  = int( (emax - emin)/estep + 1 )
    txt += "# The number of energy points NPOINT: \n"
    txt +=  ("%i \n") % npoint
    output_dictionary["npoint"] = npoint
    txt += "# for each energy point, energy, F1(1),F2(1),...,F1(nbatom),F2(nbatom)\n"
    list_energy = []
    out_f1 = numpy.zeros( (len(unique_indexes),npoint), dtype=float)
    out_f2 = numpy.zeros( (len(unique_indexes),npoint), dtype=float)
    out_fcompton = numpy.zeros( (len(unique_indexes),npoint), dtype=complex)
    for i in range(npoint):
        energy = (emin+estep*i)
        txt += ("%20.11e \n") % (energy)
        list_energy.append(energy)

        for j in range(len(unique_indexes)):
            zeta = list_Zatom[unique_indexes[j]]
            f1a =  material_constants_library.Fi(int(zeta),energy*1e-3)
            f2a = -material_constants_library.Fii(int(zeta),energy*1e-3) # TODO: check the sign!!
            txt +=  (" %20.11e %20.11e 1.000 \n")%(f1a, f2a)
            out_f1[j,i] = f1a
            out_f2[j,i] = f2a
            out_fcompton[j,i] = 1.0

    output_dictionary["energy"] = list_energy
    output_dictionary["f1"] = out_f1
    output_dictionary["f2"] = out_f2
    output_dictionary["fcompton"] = out_fcompton

    if fileout != None:
        with open(fileout,"w") as f:
            f.write(txt)
            print("File written to disk: %s" % fileout)

    return output_dictionary

def crystal_fh(input_dictionary,phot_in,theta=None,forceratio=0):
    """

    :param input_dictionary: as resulting from bragg_calc()
    :param phot_in: photon energy in eV
    :param theta: incident angle (half of scattering angle) in rad
    :return: a dictionary with structure factor
    """

    # outfil    = input_dictionary["outfil"]
    # fract     = input_dictionary["fract"]
    rn        = input_dictionary["rn"]
    dspacing  = numpy.array(input_dictionary["dspacing"])
    nbatom    = numpy.array(input_dictionary["nbatom"])
    atnum     = numpy.array(input_dictionary["atnum"])
    temper    = numpy.array(input_dictionary["temper"])
    G_0       = numpy.array(input_dictionary["G_0"])
    G         = numpy.array(input_dictionary["G"])
    G_BAR     = numpy.array(input_dictionary["G_BAR"])
    f0coeff   = numpy.array(input_dictionary["f0coeff"])
    npoint    = numpy.array(input_dictionary["npoint"])
    energy    = numpy.array(input_dictionary["energy"])
    fp        = numpy.array(input_dictionary["f1"])
    fpp       = numpy.array(input_dictionary["f2"])
    fraction = numpy.array(input_dictionary["fraction"])



    phot_in = numpy.array(phot_in,dtype=float).reshape(-1)

    toangstroms = codata.h * codata.c / codata.e * 1e10


    itheta = numpy.zeros_like(phot_in)
    for i,phot in enumerate(phot_in):

        if theta is None:
            itheta[i] = numpy.arcsin(toangstroms*1e-8/phot/2/dspacing)
        else:
            itheta[i] = theta

        # print("energy= %g eV, theta = %15.13g deg"%(phot,itheta[i]*180/numpy.pi))
        if phot < energy[0] or phot > energy[-1]:
            raise Exception("Photon energy %g eV outside of valid limits [%g,%g]"%(phot,energy[0],energy[-1]))

        if forceratio == 0:
            ratio = numpy.sin(itheta[i]) / (toangstroms / phot)
        else:
            ratio = 1 / (2 * dspacing * 1e8)
        # print("Ratio: ",ratio)

        F0 = numpy.zeros(nbatom)
        F000 = numpy.zeros(nbatom)
        for j in range(nbatom):
            #icentral = int(f0coeff.shape[1]/2)
            #F0[j] = f0coeff[j,icentral]
            icentral = int(len(f0coeff[j])/2)
            F0[j] = f0coeff[j][icentral]
            F000[j] = F0[j]
            for i in range(icentral):
                #F0[j] += f0coeff[j,i] * numpy.exp(-1.0*f0coeff[j,i+icentral+1]*ratio**2)
                F0[j] += f0coeff[j][i] * numpy.exp(-1.0*f0coeff[j][i+icentral+1]*ratio**2)
                #srio F000[j] += f0coeff[j][i]  #actual number of electrons carried by each atom, X.J. Yu, slsyxj@nus.edu.sg
            F000[j] = atnum[j] # srio
        # ;C
        # ;C Interpolate for the atomic scattering factor.
        # ;C
        for j,ienergy in enumerate(energy):
            if ienergy > phot:
                break
        nener = j - 1


        F1 = numpy.zeros(nbatom,dtype=float)
        F2 = numpy.zeros(nbatom,dtype=float)
        F = numpy.zeros(nbatom,dtype=complex)

        for j in range(nbatom):
            F1[j] = fp[j,nener] + (fp[j,nener+1] - fp[j,nener]) * \
            (phot - energy[nener]) / (energy[nener+1] - energy[nener])
            F2[j] = fpp[j,nener] + (fpp[j,nener+1] - fpp[j,nener]) * \
            (phot - energy[nener]) / (energy[nener+1] - energy[nener])

        r_lam0 = toangstroms * 1e-8 / phot
        for j in range(nbatom):
            F[j] = F0[j] + F1[j] + 1j * F2[j]
            # print("F",F)


        F_0 = 0.0 + 0.0j
        FH = 0.0 + 0.0j
        FH_BAR = 0.0 + 0.0j
        FHr = 0.0 + 0.0j
        FHi = 0.0 + 0.0j
        FH_BARr = 0.0 + 0.0j
        FH_BARi = 0.0 + 0.0j


        TEMPER_AVE = 1.0
        for j in range(nbatom):
            FH  += fraction[j] * (G[j] *   F[j] * 1.0) * temper[j]
            FHr += fraction[j] * (G[j] * (F0[j] + F1[j])* 1.0) * temper[j]
            FHi += fraction[j] * (G[j] *  F2[j] * 1.0) * temper[j]
            FN = F000[j] + F1[j] + 1j * F2[j]
            F_0 += fraction[j] * (G_0[j] *  FN  * 1.0)
            # TEMPER_AVE *= (temper[j])**(G_0[j]/(G_0.sum()))

            FH_BAR  += fraction[j] * ((G_BAR[j] * F[j] * 1.0)) * temper[j]
            FH_BARr += fraction[j] * ((G_BAR[j] * (F0[j]  + F1[j]) *1.0)) * temper[j]
            FH_BARi += fraction[j] * ((G_BAR[j] *  F2[j] * 1.0)) * temper[j]
            # print("TEMPER_AVE: ",TEMPER_AVE)


        # ;C
        # ;C multiply by the average temperature factor
        # ;C


        # FH      *= TEMPER_AVE
        # FHr     *= TEMPER_AVE
        # FHi     *= TEMPER_AVE
        # FH_BAR  *= TEMPER_AVE
        # FH_BARr *= TEMPER_AVE
        # FH_BARi *= TEMPER_AVE

        STRUCT = numpy.sqrt(FH * FH_BAR)

        # ;C
        # ;C   PSI_CONJ = F*( note: PSI_HBAR is PSI at -H position and is
        # ;C   proportional to fh_bar but PSI_CONJ is complex conjugate os PSI_H)
        # ;C


        psi_over_f = rn * r_lam0**2 / numpy.pi
        psi_h      = rn * r_lam0**2 / numpy.pi * FH
        psi_hr     = rn * r_lam0**2 / numpy.pi * FHr
        psi_hi     = rn * r_lam0**2 / numpy.pi * FHi
        psi_hbar   = rn * r_lam0**2 / numpy.pi * FH_BAR
        psi_hbarr  = rn * r_lam0**2 / numpy.pi * FH_BARr
        psi_hbari  = rn * r_lam0**2 / numpy.pi * FH_BARi
        psi_0      = rn * r_lam0**2 / numpy.pi * F_0
        psi_conj   = rn * r_lam0**2 / numpy.pi * FH.conjugate()

        # ;
        # ; Darwin width
        # ;
        # print(rn,r_lam0,STRUCT,itheta)
        ssvar = rn * (r_lam0**2) * STRUCT / numpy.pi / numpy.sin(2.0*itheta)
        spvar = ssvar * numpy.abs((numpy.cos(2.0*itheta)))
        ssr = ssvar.real
        spr = spvar.real

        # ;C
        # ;C computes refractive index.
        # ;C ([3.171] of Zachariasen's book)
        # ;C
        REFRAC = (1.0+0j) - r_lam0**2 * rn * F_0 / 2/ numpy.pi
        DELTA_REF = 1.0 - REFRAC.real
        ABSORP = 4.0 * numpy.pi * (-REFRAC.imag) / r_lam0


        txt = ""
        txt += '\n******************************************************'
        txt += '\n       at energy    = '+repr(phot)+' eV'
        txt += '\n                    = '+repr(r_lam0*1e8)+' Angstroms'
        txt += '\n       and at angle = '+repr(itheta*180.0/numpy.pi)+' degrees'
        txt += '\n                    = '+repr(itheta)+' rads'
        txt += '\n******************************************************'

        for j in range(nbatom):
            txt += '\n  '
            txt += '\nFor atom '+repr(j+1)+':'
            txt += '\n       fo + fp+ i fpp = '
            txt += '\n        '+repr(F0[j])+' + '+ repr(F1[j].real)+' + i'+ repr(F2[j])+" ="
            txt += '\n        '+repr(F0[j] + F1[j] + 1j * F2[j])
            txt += '\n       Z = '+repr(atnum[j])
            txt += '\n       Temperature factor = '+repr(temper[j])
        txt += '\n  '
        txt += '\n Structure factor F(0,0,0) = '+repr(F_0)
        txt += '\n Structure factor FH = '      +repr(FH)
        txt += '\n Structure factor FH_BAR = '  +repr(FH_BAR)
        txt += '\n Structure factor F(h,k,l) = '+repr(STRUCT)
        txt += '\n  '
        txt += '\n Psi_0  = '   +repr(psi_0)
        txt += '\n Psi_H  = '   +repr(psi_h)
        txt += '\n Psi_HBar  = '+repr(psi_hbar)
        txt += '\n  '
        txt += '\n Psi_H(real) Real and Imaginary parts = '   + repr(psi_hr)
        txt += '\n Psi_H(real) Modulus  = '                   + repr(numpy.abs(psi_hr))
        txt += '\n Psi_H(imag) Real and Imaginary parts = '   + repr(psi_hi)
        txt += '\n Psi_H(imag) Modulus  = '                   + repr(abs(psi_hi))
        txt += '\n Psi_HBar(real) Real and Imaginary parts = '+ repr(psi_hbarr)
        txt += '\n Psi_HBar(real) Modulus  = '                + repr(abs(psi_hbarr))
        txt += '\n Psi_HBar(imag) Real and Imaginary parts = '+ repr(psi_hbari)
        txt += '\n Psi_HBar(imag) Modulus  = '                + repr(abs(psi_hbari))
        txt += '\n  '
        txt += '\n Psi/F factor = '                           + repr(psi_over_f)
        txt += '\n  '
        txt += '\n Average Temperature factor = '             + repr(TEMPER_AVE)
        txt += '\n Refraction index = 1 - delta - i*beta'
        txt += '\n            delta = '                       + repr(DELTA_REF)
        txt += '\n             beta = '                       + repr(1.0e0*REFRAC.imag)
        txt += '\n Absorption coeff = '                       + repr(ABSORP)+' cm^-1'
        txt += '\n  '
        txt += '\n e^2/(mc^2)/V = '                           + repr(rn)+' cm^-2'
        txt += '\n d-spacing = '                              + repr(dspacing*1.0e8)+' Angstroms'
        txt += '\n SIN(theta)/Lambda = '                      + repr(ratio)
        txt += '\n  '
        txt += '\n Darwin width for symmetric s-pol [microrad] = ' + repr(2.0e6*ssr)
        txt += '\n Darwin width for symmetric p-pol [microrad] = ' + repr(2.0e6*spr)

    return {"PHOT":phot, "WAVELENGTH":r_lam0*1e-2 ,"THETA":itheta, "F_0":F_0, "FH":FH, "FH_BAR":FH_BAR,
	        "STRUCT":STRUCT, "psi_0":psi_0, "psi_h":psi_h, "psi_hbar":psi_hbar,
        	"DELTA_REF":DELTA_REF, "REFRAC":REFRAC, "ABSORP":ABSORP, "RATIO":ratio,
        	"ssr":ssr, "spr":spr, "psi_over_f":psi_over_f, "info":txt}


#
#
#
#
#
#
#
#
#
#

def bragg_calc2(descriptor="YB66", hh=1, kk=1, ll=1, temper=1.0, emin=5000.0, emax=15000.0, estep=100.0, ANISO_SEL=0,
                fileout=None,
                sourceCryst=2, # 0=xraylib, 1=dabax, 2=auto
                sourceF0=2,    # 0=xraylib, 1=dabax, 2=auto
                do_not_prototype=0, # 0=use site groups (recommended), 1=use all individual sites
                verbose=True,
                dabax_lib=None,
                ):
    """
    Preprocessor for Structure Factor (FH) calculations. It calculates the basic ingredients of FH.

    :param descriptor: crystal name (as in xraylib)
    :param hh: miller index H
    :param kk: miller index K
    :param ll: miller index L
    :param temper: temperature factor (scalar <=1.0 )
    :param emin:  photon energy minimum
    :param emax: photon energy maximum
    :param estep: photon energy step
    :param ANISO_SEL: source of temperature factor:
                0: use scalar value defined in temper
                1: use isotropic value calculated from keyword UNIANISO_COFF in dabax Crystal.dat file
                2: use anisotropic value calculated from keyword UNIANISO_COFF in dabax Crystal.dat file
    :param fileout: name for the output file (default=None, no output file)
    :return: a dictionary with all ingredients of the structure factor.
    """

    output_dictionary = {}

    codata_e2_mc2 = codata.e ** 2 / codata.m_e / codata.c ** 2 / (4 * numpy.pi * codata.epsilon_0)  # in m

    # f = open(fileout,'w')

    txt = ""
    txt += "# Bragg version, Data file type\n"
    txt += "2.6 1\n"

    if sourceCryst == 0:
        cryst = xraylib.Crystal_GetCrystal(descriptor)
        if cryst is None:
            raise Exception("Crystal descriptor %s not found in xraylib" % descriptor)
    elif sourceCryst == 1:
        cryst = dabax_lib.Crystal_GetCrystal(entry_name=descriptor)
    elif sourceCryst == 2:
        try:
            cryst = xraylib.Crystal_GetCrystal(descriptor)
            if cryst is None:
                raise Exception("Crystal descriptor %s not found in xraylib" % descriptor)
            sourceCryst = 0
        except:
            try:
                cryst = dabax_lib.Crystal_GetCrystal(entry_name=descriptor)
                sourceCryst = 1
            except:
                raise Exception("Crystal descriptor %s not found in xraylib nor in dabax" % descriptor)

    volume = cryst['volume']

    # test crystal data - not needed
    itest = 0
    if itest:

        print("  Unit cell dimensions are %f %f %f" % (cryst['a'], cryst['b'], cryst['c']))
        print("  Unit cell angles are %f %f %f" % (cryst['alpha'], cryst['beta'], cryst['gamma']))
        print("  Unit cell volume is %f A^3" % volume)
        print("  Atoms at:")
        print("     Z  fraction    X        Y        Z")
        for i in range(cryst['n_atom']):
            atom = cryst['atom'][i]
            print("    %3i %f %f %f %f" % (atom['Zatom'], atom['fraction'], atom['x'], atom['y'], atom['z']))
        print("  ")

    volume = volume * 1e-8 * 1e-8 * 1e-8  # in cm^3
    rn = (1e0 / volume) * (codata_e2_mc2 * 1e2)

    dspacing = bragg_metrictensor(cryst['a'], cryst['b'], cryst['c'], cryst['alpha'], cryst['beta'], cryst['gamma'], HKL=[hh, kk, ll])
    dspacing *= 1e-8  # in cm

    txt += "# RN = (e^2/(m c^2))/V) [cm^-2], d spacing [cm]\n"
    txt += "%e %e \n" % (rn, dspacing)

    output_dictionary["rn"] = rn
    output_dictionary["dspacing"] = dspacing

    atom = cryst['atom']
    number_of_atoms = len(atom)
    list_Zatom = [atom[i]['Zatom'] for i in range(len(atom))]

    list_fraction = [atom[i]['fraction'] for i in range(number_of_atoms)]
    try:
        list_charge = [atom[i]['charge'] for i in range(number_of_atoms)]
    except:
        list_charge = [0.0] * number_of_atoms
    list_x = [atom[i]['x'] for i in range(number_of_atoms)]
    list_y = [atom[i]['y'] for i in range(number_of_atoms)]
    list_z = [atom[i]['z'] for i in range(number_of_atoms)]

    # calculate array of temperature factor for all atoms
    #
    # Consider anisotropic temperature factor
    # X.J. Yu, slsyxj@nus.edu.sg
    # A dummy dictionary Aniso with start =0 if no aniso temperature factor input
    # start
    if 'Aniso' in cryst.keys() and cryst['Aniso'][0]['start'] > 0:  # most crystals have no Anisotropic input
        TFac = TemperFactor(1.0 / (2.0 * dspacing * 1e8), cryst['Aniso'], Miller={'h': hh, 'k': kk, 'l': ll}, \
                            cell={'a': cryst['a'], 'b': cryst['b'], 'c': cryst['c']}, n=len(atom))
        B_TFac = 1
    else:
        B_TFac = 0
    #
    #
    #
    list_temper = []
    list_temper_label = []
    if ANISO_SEL == 0:
        for i in range(number_of_atoms):
            list_temper.append(temper)
            list_temper_label.append(-1)
    elif ANISO_SEL == 1:
        if B_TFac:
            for i in range(number_of_atoms):
                list_temper.append(TFac[0, i])
                list_temper_label.append(TFac[2, i])
        else:
            raise Exception("No crystal data to calculate isotropic temperature factor for crystal %s" % descriptor)
    elif ANISO_SEL == 2:
        if B_TFac:
            for i in range(number_of_atoms):
                list_temper.append(TFac[1, i])
                list_temper_label.append(TFac[2, i])
        else:
            raise Exception("No crystal data to calculate anisotropic temperature factor for crystal %s" % descriptor)

    list_AtomicName = []
    for i in range(number_of_atoms):
        s = atomic_symbols()[atom[i]['Zatom']]
        if sourceCryst == 1: # charge is not available in xraylib
            if atom[i]['charge'] != 0.0:  # if charge is 0, s is symbol only, not B0, etc
                s = s + f'%+.6g' % atom[i]['charge']
        list_AtomicName.append(s)

    # identify the prototypical atoms
    labels_prototypical = []
    for i in range(number_of_atoms):
        labels_prototypical.append("Z=%d C=%g F=%g T=%g" % (list_Zatom[i], list_charge[i], list_fraction[i], list_temper_label[i]))

    if do_not_prototype:
        indices_prototypical = numpy.arange(number_of_atoms)  # different with diff_pat for complex crystal
    else:
        indices_prototypical = numpy.unique(labels_prototypical, return_index=True)[1]

    number_of_prototypical_atoms = len(indices_prototypical)

    # for i in range(number_of_prototypical_atoms):
    #     print("   >>> ", i, indices_prototypical[i], labels_prototypical[indices_prototypical[i]])
    #
    # for i in indices_prototypical:
    #     print("   >>>>> ", i, labels_prototypical[i])
    #
    # print(">>>>  list_labels", len(labels_prototypical), len(indices_prototypical), labels_prototypical)

    #
    # get f0 coefficients
    #

    f0coeffs = []
    if sourceF0 == 0:
        for i in indices_prototypical:
            f0coeffs.append(f0_xop(atom[i]['Zatom']))
    elif sourceF0 == 1:
        for i in indices_prototypical:
                f0coeffs.append(dabax_lib.f0_with_fractional_charge(atom[i]['Zatom'], atom[i]['charge']) )
    elif sourceF0 == 2:
        total_charge_flag = numpy.abs(numpy.array(list_charge)).sum() # note the abs(): to be used as flag...

        if total_charge_flag != 0: # Use dabax
            for i in indices_prototypical:
                f0coeffs.append(dabax_lib.f0_with_fractional_charge(atom[i]['Zatom'], atom[i]['charge']))
        else: # use xraylib
            if 'AtomicName' not in atom[0].keys():
                for i in indices_prototypical:  #normal case come in here
                    f0coeffs.append(f0_xop(atom[i]['Zatom']))
            else:   #for case with like 'Y3+' entries in f0_xop
                import re
                for i in indices_prototypical:
                    x = atom[i]['AtomicName']
                    tmp_x = re.search('(^[a-zA-Z]*)',x)
                    if tmp_x.group(0) == x:
                        f0coeffs.append(f0_xop(atom[i]['Zatom']))  #neutral atom
                    else:
                        f0coeffs.append(f0_xop(0,AtomicName=x))    #charged atom

    txt += "# Number of different element-sites in unit cell NBATOM:\n%d \n" % number_of_prototypical_atoms
    output_dictionary["nbatom"] = number_of_prototypical_atoms

    txt += "# for each element-site, the number of scattering electrons (Z_i + charge_i)\n"
    atnum_list = []
    for i in indices_prototypical:
        txt += "%f " % (list_Zatom[i] + list_charge[i])
        atnum_list.append(list_Zatom[i] + list_charge[i])
    txt += "\n"
    output_dictionary["atnum"] = atnum_list


    txt += "# for each element-site, the occupation factor\n"
    unique_fraction = [list_fraction[i] for i in indices_prototypical]
    for z in unique_fraction:
        txt += "%g " % (z)
    txt += "\n"
    output_dictionary["fraction"] = unique_fraction

    txt += "# for each element-site, the temperature factor\n"  # temperature parameter
    unique_temper = []
    for i in indices_prototypical:
        txt += "%g " % list_temper[i]
        unique_temper.append(list_temper[i])
    txt += "\n"
    output_dictionary["temper"] = unique_temper

    #
    # Geometrical part of structure factor:  G and G_BAR
    #
    txt += "# for each type of element-site, COOR_NR=G_0\n"
    list_multiplicity = []
    for i in indices_prototypical:
        # zz = list_AtomicName[i]
        # fraction = list_fraction[i]
        # temper = list_temper[i]
        # count = 0
        # for j in range(len(list_Zatom)):
        #     if (list_AtomicName[j] == zz) and (list_fraction[j] == fraction) and (list_temper[j] == temper): count += 1

        if do_not_prototype:
            txt += "%d " % 1
            list_multiplicity.append(1)
        else:
            count = 0
            for j in range(number_of_atoms):
                if labels_prototypical[j] == labels_prototypical[i]: count += 1
            txt += "%d " % count
            list_multiplicity.append(count)
    txt += "\n"
    output_dictionary["G_0"] = list_multiplicity


    txt += "# for each type of element-site, G and G_BAR (both complex)\n"
    list_g = []
    list_g_bar = []
    for i in indices_prototypical:

        if do_not_prototype:
            # # ga_item = numpy.exp(2j * numpy.pi * (hh * list_x[i] + kk * list_y[i] + ll * list_z[i]))
            # ga += ga_item
            ga = numpy.exp(2j * numpy.pi * (hh * list_x[i] + kk * list_y[i] + ll * list_z[i]))
        else:
            ga = 0.0 + 0j
            for j in range(number_of_atoms):
                if labels_prototypical[j] == labels_prototypical[i]:
                # if list_AtomicName[j] == zz and list_fraction[j] == ff and list_temper[j] == tt:
                    ga_item = numpy.exp(2j * numpy.pi * (hh * list_x[j] + kk * list_y[j] + ll * list_z[j]))
                    ga += ga_item

        txt += "(%g,%g) \n" % (ga.real, ga.imag)
        txt += "(%g,%g) \n" % (ga.real, -ga.imag)
        list_g.append(ga)
        list_g_bar.append(ga.conjugate())
    output_dictionary["G"] = list_g
    output_dictionary["G_BAR"] = list_g_bar

    #
    # F0 part
    #
    txt += "# for each type of element-site, the number of f0 coefficients followed by them\n"
    for f0coeffs_item in f0coeffs:
        txt += "%d " % len(f0coeffs_item)
        for cc in f0coeffs_item:
            txt += "%g " % cc
        txt += "\n"
    output_dictionary["f0coeff"] = f0coeffs


    # X.J. Yu, use ceil to round up, otherwise we may get actual max energy less than emax
    npoint = int(numpy.ceil(((emax - emin) / estep + 1)))
    txt += "# The number of energy points NPOINT: \n"
    txt += ("%i \n") % npoint
    output_dictionary["npoint"] = npoint

    txt += "# for each energy point, energy, F1(1),F2(1),...,F1(nbatom),F2(nbatom)\n"
    list_energy = []
    out_f1 = numpy.zeros((len(indices_prototypical), npoint), dtype=float)
    out_f2 = numpy.zeros((len(indices_prototypical), npoint), dtype=float)
    out_fcompton = numpy.zeros((len(indices_prototypical), npoint), dtype=complex)

    for i in range(npoint):
        energy = (emin + estep * i)
        txt += ("%20.11e \n") % (energy)
        list_energy.append(energy)

        for j,jj in enumerate(indices_prototypical):
            f1a = xraylib.Fi(list_Zatom[jj], energy * 1e-3)
            f2a = -xraylib.Fii(list_Zatom[jj], energy * 1e-3)  # TODO: check the sign!!
            txt += (" %20.11e %20.11e 1.000 \n") % (f1a, f2a)
            out_f1[j, i] = f1a
            out_f2[j, i] = f2a
            out_fcompton[j, i] = 1.0

    output_dictionary["energy"] = list_energy
    output_dictionary["f1"] = out_f1
    output_dictionary["f2"] = out_f2
    output_dictionary["fcompton"] = out_fcompton

    if fileout != None:
        with open(fileout, "w") as f:
            f.write(txt)
        if verbose: print("File written to disk: %s" % fileout)

    return output_dictionary

def TemperFactor(sinTheta_lambda,anisos,Miller={'h':1,'k':1,'l':1},cell={'a':23.44,'b':23.44,'c':23.44},n=1936):
    '''
    #+
    # Singapore Synchrotron Light Source (SSLS)
    # :Author: X.J. Yu, slsyxj@nus.edu.sg
    # :Name:  TemperFactor
    # :Purpose: Calculation isotropic & anisotropic temerature factors
    # :Input:
    #     Miller: Miller indices
    #     cell:  dictionary of lattice [a,b,c] in units of Angstrom
    #     sinTheta_lambda: Sin(theta)/lambda, lambda in units of Angstrom
    #     n: number of atomic sites
    #     anisos: array of dictionary containing anisotropic coefficients
    #     Out: output results in a 2-elements list: [[sotropic],[anisotropic]]
    #-
    '''
    #0: isotropic, 1: anisotropic temerature factors
    # results = numpy.zeros([2,n])
    results = numpy.zeros([3,n]) # srio adds "start"

    for i,aniso in enumerate(anisos):
        s = aniso['start']-1
        e = aniso['end']
        if aniso['beta11'] >= 1:
            print(">>>>>>>>>>>>>>>>> WHY AM I HERE? ")
            #if beta11>=1, then beta22 is Beq, the other fields are unused
            #if Beq specified, anisotropic temperature factor same as isotropic
            Beq = aniso['beta22']
            results[1,s:e] = numpy.exp(-sinTheta_lambda*sinTheta_lambda*Beq)
        else:
            Beq = 4.0/3.0*( aniso['beta11']*cell['a']*cell['a']+aniso['beta22']*cell['b']*cell['b']+ \
                aniso['beta33']*cell['c']*cell['c'] ) # this is true only for cubic, tetragonal and orthorhombic Giacovazzo pag 188
            results[1,s:e] = numpy.exp(-(aniso['beta11']*Miller['h']*Miller['h'] + \
                  aniso['beta22']*Miller['k']*Miller['k'] + aniso['beta33']*Miller['l']*Miller['l'] + \
                  2.0*Miller['h']*Miller['k']*aniso['beta12'] + 2.0*Miller['h']*Miller['l']*aniso['beta13'] + 2.0*Miller['k']*Miller['l']*aniso['beta23']))
        results[0,s:e] = numpy.exp(-sinTheta_lambda*sinTheta_lambda*Beq)

        results[2, s:e] = s

    return results


if __name__ == "__main__":

    # method = 1 # 0 = xraylib, 1=Dabax
    #
    #
    # if method == 0:
    #     tmp = create_bragg_preprocessor_file_v2(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1,
    #                                                            K_MILLER_INDEX=1, L_MILLER_INDEX=1,
    #                                                            TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0,
    #                                                            E_STEP=100.0, SHADOW_FILE="bragg0_v2.dat",
    #                                                            material_constants_library=xraylib)
    # else:
    #     from dabax.dabax_xraylib import DabaxXraylib
    #     import socket
    #     if socket.getfqdn().find("esrf") >= 0:
    #         dx = DabaxXraylib(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
    #     else:
    #         dx = DabaxXraylib()
    #
    #     print(dx.info())
    #     tmp = create_bragg_preprocessor_file_v2(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1,
    #                                             K_MILLER_INDEX=1, L_MILLER_INDEX=1,
    #                                             TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0,
    #                                             E_STEP=100.0, SHADOW_FILE="bragg1_v2.dat",
    #                                             material_constants_library=dx)
    #
    # for key in tmp.keys():
    #     print("---------------", key)
    #     # print(out[key] - tmp[key])
    #
    # tmp_fh = crystal_fh(tmp, 8000.0)
    # print(tmp_fh['info'])

    from dabax.dabax_xraylib import DabaxXraylib
    import socket
    if socket.getfqdn().find("esrf") >= 0:
        dx = DabaxXraylib(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
    else:
        dx = DabaxXraylib()

    tmp = bragg_calc2(descriptor="YB66", hh=1, kk=1, ll=1, temper=1.0,
                emin=5000.0, emax=15000.0, estep=100.0, ANISO_SEL=0,
                fileout=None,
                sourceCryst=2, # 0=xraylib, 1=dabax, 2=auto
                sourceF0=2,    # 0=xraylib, 1=dabax, 2=auto
                do_not_prototype=0, # 0=use site groups (recommended), 1=use all individual sites
                verbose=True,
                dabax_lib=dx,)
