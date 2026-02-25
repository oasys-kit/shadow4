import numpy
from dabax.dabax_xraylib import DabaxXraylib
from dabax.common_tools import f0_xop, f0_xop_with_fractional_charge
from dabax.common_tools import bragg_metrictensor, lorentz, atomic_symbols

import scipy.constants as codata

def create_bragg_preprocessor_file_v2(
        interactive=False, # back compatibility, not used
        DESCRIPTOR="Si",
        H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
        TEMPERATURE_FACTOR=1.0,
        E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0,
        SHADOW_FILE="bragg.dat",
        material_constants_library=None):

    """
    Preprocessor for Crystal Structure Factor calculations.
    It calculates the SHADOW bragg preprocessor file version 2.

    Note that the created file is read using crystalpy:  DiffractionSetupShadowPreprocessorV2()
    SHADOW4 can read version 1 file, but it cannot create it.


    Parameters
    ----------
    descriptor: str, optional
        crystal name (as in dabax, xraylib)
    H_MILLER_INDEX: int, optional
        miller index H
    K_MILLER_INDEX: int, optional
        miller index K
    L_MILLER_INDEX: int, optional
        miller index L
    TEMPERATURE_FACTOR: float, optional
        temperature factor (scalar <=1.0 )
    E_MIN: float, optional
        photon energy minimum in eV
    E_MAX: float, optional
        photon energy maximum in eV
    E_STEP: float, optional
        photon energy step in eV
    SHADOW_FILE: None or str, optional
        name for the output file (default=None, no output file)
    material_constants_library: xraylib or instance of DabaxXraylib, optional
        The pointer to the material library to be used to retrieve scattering data.

    Returns
    -------
    dict
        a dictionary with all ingredients of the structure factor.
    """

    if material_constants_library is None:
        try:    material_constants_library = xraylib
        except: material_constants_library = DabaxXraylib()
    # codata_e2_mc2 = 2.81794032e-15 = Classical electron radius in S.I.
    # codata_e2_mc2 = codata.hbar * codata.alpha / codata.m_e / codata.c

    return bragg_calc2(
                descriptor= DESCRIPTOR,
                hh        = int(H_MILLER_INDEX),
                kk        = int(K_MILLER_INDEX),
                ll        = int(L_MILLER_INDEX),
                temper=     float(TEMPERATURE_FACTOR),
                emin      = float(E_MIN),
                emax      = float(E_MAX),
                estep     = float(E_STEP),
                fileout   = SHADOW_FILE,
                material_constants_library=material_constants_library,
                )


def bragg_calc2(descriptor="YB66", hh=1, kk=1, ll=1, temper=1.0,
                emin=5000.0, emax=15000.0, estep=100.0, ANISO_SEL=0,
                fileout=None,
                do_not_prototype=0, # 0=use site groups (recommended), 1=use all individual sites
                material_constants_library=None,
                verbose=False
                ):
    """
    Preprocessor for Structure Factor (FH) calculations. It calculates the basic ingredients of FH.


    Parameters
    ----------
    descriptor: str, optional
        crystal name (as in xraylib)
    hh: int, optional
        miller index H
    kk: int, optional
        miller index K
    ll: int, optional
        miller index L
    temper: float, optional
        temperature factor (scalar <=1.0 )
    emin: float, optional
        photon energy minimum in eV
    emax: float, optional
        photon energy maximum in eV
    estep: float, optional
        photon energy step in eV
    ANISO_SEL: int, optional
        0: Do not use anisotropy.
        1: Use anisotropy.
    fileout: None or str, optional
        name for the output file (default=None, no output file)
    do_not_prototype: int, optional
        for computing the structure factor, 0=sum over site groups (recommended), 1=sum over each individual sites
    verbose: int, optional
        Set to 1 for verbose output.
    material_constants_library: xraylib or instance of DabaxXraylib, optional
        The pointer to the material library to be used to retrieve scattering data.

    Returns
    -------
    dict
        a dictionary with all ingredients of the structure factor.
    """

    if material_constants_library is None:
        material_constants_library = DabaxXraylib()

    output_dictionary = {}

    codata_e2_mc2 = codata.e ** 2 / codata.m_e / codata.c ** 2 / (4 * numpy.pi * codata.epsilon_0)  # in m

    # f = open(fileout,'w')

    version = "2.6"
    output_dictionary["version"] = version


    txt = ""
    txt += "# Bragg version, Data file type\n"
    txt += "%s\n" % version

    cryst = material_constants_library.Crystal_GetCrystal(descriptor)

    if cryst is None:
        raise Exception("Crystal descriptor %s not found in material constants library" % descriptor)

    volume = cryst['volume']

    # test crystal data - not needed
    icheck= 0
    if icheck:
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
        TFac = _temper_factor(1.0 / (2.0 * dspacing * 1e8), cryst['Aniso'], Miller={'h': hh, 'k': kk, 'l': ll}, \
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
        # if sourceCryst == 1: # charge is not available in xraylib
        try: # charge is not available in xraylib
            if atom[i]['charge'] != 0.0:  # if charge is 0, s is symbol only, not B0, etc
                s = s + f'%+.6g' % atom[i]['charge']
        except:
            pass
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

    f0coeffs = []
    for i in indices_prototypical:
        try:
            charge = atom[i]['charge']
        except:
            charge = 0.0
        f0coeffs.append(f0_xop_with_fractional_charge(atom[i]['Zatom'], charge))


    txt += "# Number of different element-sites in unit cell NBATOM:\n%d \n" % number_of_prototypical_atoms
    output_dictionary["nbatom"] = number_of_prototypical_atoms

    txt += "# for each element-site, the number of scattering electrons (Z_i + charge_i)\n"
    atnum_list = []
    for i in indices_prototypical:
        txt += "%f " % (list_Zatom[i] - list_charge[i])
        atnum_list.append(list_Zatom[i] - list_charge[i])
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
    out_fcompton = numpy.zeros((len(indices_prototypical), npoint), dtype=float) # todo: is complex?

    if isinstance(material_constants_library, DabaxXraylib ):
        # vectorize with DABAX
        energies = numpy.zeros(npoint)
        for i in range(npoint):
            energies[i] = (emin + estep * i)

        DABAX_F_RESULTS = []
        for j, jj in enumerate(indices_prototypical):
            DABAX_F_RESULTS.append(numpy.array( material_constants_library.FiAndFii(list_Zatom[jj], energies * 1e-3)))

        for i in range(npoint):
            energy = (emin + estep * i)
            txt += ("%20.11e \n") % (energy)
            list_energy.append(energy)

            for j, jj in enumerate(indices_prototypical):
                f1a =  (DABAX_F_RESULTS[j])[0, i]  # material_constants_library.Fi(list_Zatom[jj], energy * 1e-3)
                f2a = -(DABAX_F_RESULTS[j])[1, i]  # -material_constants_library.Fii(list_Zatom[jj], energy * 1e-3)
                txt += (" %20.11e %20.11e 1.000 \n") % (f1a, f2a)
                out_f1[j, i] = f1a
                out_f2[j, i] = f2a
                out_fcompton[j, i] = 1.0
    else:
        # make a simple loop with xraylib (fast)
        for i in range(npoint):
            energy = (emin + estep * i)
            txt += ("%20.11e \n") % (energy)
            list_energy.append(energy)

            for j,jj in enumerate(indices_prototypical):
                f1a = material_constants_library.Fi(list_Zatom[jj], energy * 1e-3)
                f2a = -material_constants_library.Fii(list_Zatom[jj], energy * 1e-3)
                txt += (" %20.11e %20.11e 1.000 \n") % (f1a, f2a)
                out_f1[j, i] = f1a
                out_f2[j, i] = f2a
                out_fcompton[j, i] = 1.0


    output_dictionary["energy"] = list_energy
    output_dictionary["f1"] = out_f1
    output_dictionary["f2"] = out_f2
    output_dictionary["fcompton"] = out_fcompton

    if fileout != None:
        _bragg_preprocessor_file_v2_write(output_dictionary, fileout)

    return output_dictionary

def _temper_factor(sinTheta_lambda, anisos, Miller={'h':1,'k':1,'l':1}, cell={'a':23.44,'b':23.44,'c':23.44}, n=1936):
    """
    Calculation isotropic & anisotropic temperature factors.

    Parameters
    ----------
    sinTheta_lambda: float
        Sin(theta)/lambda, lambda in units of Angstrom.
    anisos: numpy array
        array of dictionary containing anisotropic coefficients.
    Miller: dict
        The miller indices, example: {'h':1,'k':1,'l':1}.
    cell: dict
        The cell a,b,c parameters, example: {'a':23.44,'b':23.44,'c':23.44}
    n: int, optional
        number of atomic sites.

    Returns
    -------
    list
        output results in a 2-elements list: [[isotropic],[anisotropic]].

    """

    #+
    # Singapore Synchrotron Light Source (SSLS)
    # :Author: X.J. Yu, slsyxj@nus.edu.sg
    # :Name:  _temper_factor
    # :Purpose: Calculation isotropic & anisotropic temerature factors
    # :Input:
    #     Miller: Miller indices
    #     cell:  dictionary of lattice [a,b,c] in units of Angstrom
    #     sinTheta_lambda: Sin(theta)/lambda, lambda in units of Angstrom
    #     n: number of atomic sites
    #     anisos: array of dictionary containing anisotropic coefficients
    #     Out: output results in a 2-elements list: [[isotropic],[anisotropic]]
    #-

    #0: isotropic, 1: anisotropic temerature factors
    # results = numpy.zeros([2,n])
    results = numpy.zeros([3,n]) # srio adds "start"

    for i, aniso in enumerate(anisos):
        s = aniso['start'] - 1
        e = aniso['end']
        if aniso['beta11'] >= 1:
            #if beta11>=1, then beta22 is Beq, the other fields are unused
            #if Beq specified, anisotropic temperature factor same as isotropic
            Beq = aniso['beta22']
            results[1, s:e] = numpy.exp(-sinTheta_lambda * sinTheta_lambda*Beq)
        else:
            Beq = 4.0 / 3.0 * ( aniso['beta11'] * cell['a'] * cell['a'] + aniso['beta22'] * cell['b'] * cell['b']+ \
                aniso['beta33'] * cell['c'] * cell['c'] ) # this is true only for cubic, tetragonal and orthorhombic Giacovazzo pag 188
            results[1,s:e] = numpy.exp(-(aniso['beta11'] * Miller['h'] * Miller['h'] + \
                  aniso['beta22'] * Miller['k'] * Miller['k'] + aniso['beta33'] * Miller['l'] * Miller['l'] + \
                  2.0 * Miller['h'] * Miller['k'] * aniso['beta12'] + 2.0 * Miller['h'] * Miller['l'] * aniso['beta13'] + \
                                         2.0 * Miller['k'] * Miller['l'] * aniso['beta23']))
        results[0,s:e] = numpy.exp(-sinTheta_lambda * sinTheta_lambda * Beq)

        results[2, s:e] = s

    return results

def _bragg_preprocessor_file_v2_write(output_dictionary, fileout=None):

    txt = ""
    txt += "# Bragg version, Data file type\n"
    txt += "%s 1\n" % output_dictionary["version"]
    txt += "# RN = (e^2/(m c^2))/V) [cm^-2], d spacing [cm]\n"
    txt += "%e %e \n" % (output_dictionary["rn"] , output_dictionary["dspacing"])

    txt += "# Number of different element-sites in unit cell NBATOM:\n%d \n" % output_dictionary["nbatom"]


    txt += "# for each element-site, the number of scattering electrons (Z_i + charge_i)\n"
    for i in output_dictionary['atnum']:
        txt += "%g "%i
    txt += "\n"

    txt += "# for each element-site, the occupation factor\n"
    for i in output_dictionary["fraction"]:
        txt += "%g "%i
    txt += "\n"

    txt += "# for each element-site, the temperature factor\n" # temperature parameter
    for i in output_dictionary["temper"]:
        txt += "%g "%i
    txt += "\n"

    #
    # Geometrical part of structure factor:  G and G_BAR
    #
    txt += "# for each type of element-site, COOR_NR=G_0\n"


    for i in output_dictionary["G_0"]:
        txt += "%g "%i
    txt += "\n"

    #
    txt += "# for each type of element-site, G and G_BAR (both complex)\n"

    for i,ga in enumerate(output_dictionary["G"]):
        ga_bar = output_dictionary["G_BAR"][i]
        txt += "(%g,%g) \n"%(ga.real,ga.imag)
        txt += "(%g,%g) \n"%(ga_bar.real,ga_bar.imag)

    #
    # F0 part
    #
    txt += "# for each type of element-site, the number of f0 coefficients followed by them\n"
    for i in range(len(output_dictionary['f0coeff'])):
        coeff = output_dictionary['f0coeff'][i]
        nn = len(coeff)
        txt += ("%d "%(nn)+"%g "*nn+"\n")%(tuple(coeff))


    txt += "# The number of energy points NPOINT: \n"
    txt +=  ("%i \n") % output_dictionary["npoint"]
    txt += "# for each energy point, energy, F1(1),F2(1),...,F1(nbatom),F2(nbatom)\n"

    for i in range(output_dictionary["npoint"]):
        txt += ("%20.11e \n") % (output_dictionary["energy"][i])

        for j in range(output_dictionary['nbatom']):
            f1a = output_dictionary['f1'][j,i]
            f2a = output_dictionary['f2'][j,i]
            fcompton = output_dictionary['fcompton'][j,i]
            txt +=  (" %20.11e %20.11e %20.11e \n")%(f1a, f2a, fcompton)


    if fileout != None:
        with open(fileout,"w") as f:
            f.write(txt)
            print("File written to disk: %s" % fileout)

    return txt

if __name__ == "__main__":

    for method in [0]:
        if method == 0:
            SHADOW_FILE = "bragg_v2_dabax.dat"
            tmp = create_bragg_preprocessor_file_v2(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                  TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0, SHADOW_FILE=SHADOW_FILE,
                  material_constants_library=DabaxXraylib())

        else:
            import xraylib
            SHADOW_FILE = "bragg_v2_xraylib.dat"
            tmp = create_bragg_preprocessor_file_v2(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1,
                                                                   K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                                                                   TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0,
                                                                   E_STEP=100.0, SHADOW_FILE=SHADOW_FILE,
                                                                   material_constants_library=xraylib)

