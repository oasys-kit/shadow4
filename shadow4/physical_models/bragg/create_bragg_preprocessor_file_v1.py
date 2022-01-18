
import numpy
import xraylib
import scipy.constants as codata



def dump_bragg_preprocessor_file_v1(out_dict, fileout=""):
    if fileout != "":
        f = open(fileout, 'wt')

        f.write("%i " % out_dict["i_latt"])  # flag ZincBlende
        # f.write("%e " % ((1e0 / volume_in_cm3) * (codata_e2_mc2 * 1e2))) # 1/V*electronRadius
        f.write("%e " % (out_dict["one_over_volume_times_electron_radius_in_cm"]))
        f.write("%e " % out_dict["dspacing_in_cm"])
        f.write("\n")
        f.write("%i " % out_dict["zeta_a"])
        f.write("%i " % out_dict["zeta_b"])
        f.write("%e " % out_dict["temper"])  # temperature parameter
        f.write("\n")
        f.write("(%20.11e,%20.11e ) \n" % (out_dict["ga.real"], out_dict["ga.imag"]))
        f.write("(%20.11e,%20.11e ) \n" % (out_dict["ga_bar.real"], out_dict["ga_bar.imag"]))
        f.write("(%20.11e,%20.11e ) \n" % (out_dict["gb.real"], out_dict["gb.imag"]))
        f.write("(%20.11e,%20.11e ) \n" % (out_dict["gb_bar.real"], out_dict["gb_bar.imag"]))
        f.write("%e %e %e  \n" % (out_dict["fit_a"][0], out_dict["fit_a"][1], out_dict["fit_a"][2]  ))
        f.write("%e %e %e  \n" % (out_dict["fit_b"][0], out_dict["fit_b"][1], out_dict["fit_b"][2]  ))
        f.write(("%i \n") % out_dict["npoint"])
        for i in range(out_dict["npoint"]):
            f.write(("%20.11e %20.11e %20.11e \n %20.11e %20.11e \n") % ( \
                out_dict["Energy"][i],
                out_dict["F1a"][i],
                out_dict["F2a"][i],
                out_dict["F1b"][i],
                out_dict["F2b"][i]
            ))
        f.close()
        print("File written to disk: %s" % fileout)


def create_bragg_preprocessor_file_v1(interactive=True,
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
    cryst = material_constants_library.Crystal_GetCrystal(descriptor)
    volume = cryst['volume']
    volume_in_cm3 = volume * 1e-8 * 1e-8 * 1e-8  # in cm^3
    one_over_volume_times_electron_radius_in_cm = (1e0 / volume_in_cm3) * (codata_e2_mc2 * 1e2)

    # test crystal data - not needed

    if (cryst == None):
        raise Exception("Undefined crystal")
    print("  Unit cell dimensions are %f %f %f" % (cryst['a'], cryst['b'], cryst['c']))
    print("  Unit cell angles are %f %f %f" % (cryst['alpha'], cryst['beta'], cryst['gamma']))
    print("  Unit cell volume is %f A^3" % volume)
    print("  Atoms at:")
    print("     Z  fraction    X        Y        Z")
    for i in range(cryst['n_atom']):
        atom = cryst['atom'][i]
        print("    %3i %f %f %f %f" % (atom['Zatom'], atom['fraction'], atom['x'], atom['y'], atom['z']))
    print("  ")

    # dspacing
    dspacing = material_constants_library.Crystal_dSpacing(cryst, hh, kk, ll)
    dspacing_in_cm = dspacing * 1e-8
    atom = cryst['atom'] # Z's

    ga = (1e0 + 0j) + numpy.exp(1j * numpy.pi * (hh + kk)) \
         + numpy.exp(1j * numpy.pi * (hh + ll)) \
         + numpy.exp(1j * numpy.pi * (kk + ll))
    gb = ga * numpy.exp(1j * numpy.pi * 0.5 * (hh + kk + ll))
    ga_bar = ga.conjugate()
    gb_bar = gb.conjugate()
    zetas = numpy.array([atom[0]["Zatom"], atom[7]["Zatom"]])
    zeta_a = zetas[0]
    zeta_b = zetas[1]
    npoint = int((emax - emin) / estep + 1)
    for i,zeta in enumerate(zetas):
        xx01 = 1e0 / 2e0 / dspacing
        xx00 = xx01 - 0.1
        xx02 = xx01 + 0.1
        yy00 = material_constants_library.FF_Rayl(int(zeta), xx00)
        yy01 = material_constants_library.FF_Rayl(int(zeta), xx01)
        yy02 = material_constants_library.FF_Rayl(int(zeta), xx02)
        xx = numpy.array([xx00, xx01, xx02])
        yy = numpy.array([yy00, yy01, yy02])
        fit = numpy.polyfit(xx, yy, 2)
        if i == 0:
            fit_a = fit[::-1] # (tuple(fit[::-1].tolist()))
        elif i == 1:
            fit_b = fit[::-1] # (tuple(fit[::-1].tolist()))
        else:
            raise Exception("Unknown crystal structure")

    Energy = numpy.zeros(npoint)
    F1a = numpy.zeros(npoint)
    F1b = numpy.zeros(npoint)
    F2a = numpy.zeros(npoint)
    F2b = numpy.zeros(npoint)


    for i in range(npoint):
        energy = (emin + estep * i)
        Energy[i] = energy
        F1a[i] =     material_constants_library.Fi(int(zetas[0]), energy * 1e-3)
        F2a[i] = abs(material_constants_library.Fii(int(zetas[0]), energy * 1e-3))
        F1b[i] =     material_constants_library.Fi(int(zetas[1]), energy * 1e-3)
        F2b[i] = abs(material_constants_library.Fii(int(zetas[1]), energy * 1e-3))

    out_dict = {}
    # inputs
    out_dict["fileout"]    = fileout
    out_dict["descriptor"] = descriptor
    out_dict["hh"]         = hh
    out_dict["kk"]         = kk
    out_dict["ll"]         = ll
    out_dict["temper"]     = temper
    out_dict["emin"]       = emin
    out_dict["emax"]       = emax
    out_dict["estep"]      = estep
    # outputs
    out_dict["i_latt"] = 0
    out_dict["one_over_volume_times_electron_radius_in_cm"] = one_over_volume_times_electron_radius_in_cm
    out_dict["dspacing_in_cm"] = dspacing_in_cm
    out_dict["zeta_a"] = zeta_a
    out_dict["zeta_b"] = zeta_b
    out_dict["temper"] = temper
    out_dict["ga.real"] = ga.real
    out_dict["ga.imag"] = ga.imag
    out_dict["ga_bar.real"] = ga_bar.real
    out_dict["ga_bar.imag"] = ga_bar.imag

    out_dict["gb.real"] = gb.real
    out_dict["gb.imag"] = gb.imag
    out_dict["gb_bar.real"] = gb_bar.real
    out_dict["gb_bar.imag"] = gb_bar.imag

    out_dict["fit_a"] = fit_a
    out_dict["fit_b"] = fit_b
    out_dict["npoint"] = npoint
    out_dict["Energy"] = Energy
    out_dict["F1a"] = F1a
    out_dict["F2a"] = F2a
    out_dict["F1b"] = F1b
    out_dict["F2b"] = F2b

    dump_bragg_preprocessor_file_v1(out_dict, fileout=fileout)

    return out_dict


if __name__ == "__main__":

    method = 1 # 0 = xraylib, 1=Dabax


    if method == 0:
        tmp = create_bragg_preprocessor_file_v1(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1,
                                                               K_MILLER_INDEX=1, L_MILLER_INDEX=1,
                                                               TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0,
                                                               E_STEP=100.0, SHADOW_FILE="bragg0.dat",
                                                               material_constants_library=xraylib)
    else:
        from dabax.dabax_xraylib import DabaxXraylib
        import socket
        if socket.getfqdn().find("esrf") >= 0:
            dx = DabaxXraylib(dabax_repository="http://ftp.esrf.fr/pub/scisoft/DabaxFiles/")
        else:
            dx = DabaxXraylib()

        print(dx.info())
        tmp = create_bragg_preprocessor_file_v1(interactive=False, DESCRIPTOR="Si", H_MILLER_INDEX=1, K_MILLER_INDEX=1, L_MILLER_INDEX=1,
              TEMPERATURE_FACTOR=1.0, E_MIN=5000.0, E_MAX=15000.0, E_STEP=100.0, SHADOW_FILE="bragg1.dat",
              material_constants_library=dx)

    for key in tmp.keys():
        print("---------------", key)
        # print(out[key] - tmp[key])

