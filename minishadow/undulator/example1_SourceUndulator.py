#
# examples for SourceUndulator
#

import os
from SourceUndulator import SourceUndulator


#
# switch on/off plots
#
DO_PLOT = True

#
# Tests
#

if __name__ == "__main__":

    import Shadow
    import numpy

    # some inputs
    E_ENERGY=6.04
    INTENSITY=0.2
    SX=0.04
    SZ=0.001
    SXP=10e-6
    SZP=4e-6
    FLAG_EMITTANCE=0
    LAMBDAU=0.032
    NPERIODS=50
    K=0.25
    EMIN= 10498.0000
    EMAX= 10499.0000
    NG_E=11
    MAXANGLE=0.1
    NG_T=51
    NG_P=11
    N_J=20
    SEED=36255
    NRAYS=15000


    u = SourceUndulator()

    u.set_from_keywords(
        E_ENERGY = E_ENERGY,
        N_J = N_J,
        INTENSITY = INTENSITY,
        SX = SX,
        SZ = SZ,
        SXP = SXP,
        SZP = SZP,
        FLAG_EMITTANCE = FLAG_EMITTANCE,
        LAMBDAU = LAMBDAU,
        NPERIODS = NPERIODS,
        K = K,
        EMIN = EMIN,
        EMAX = EMAX,
        NG_E = NG_E,
        MAXANGLE = MAXANGLE,
        NG_T = NG_T,
        NG_P = NG_P,
        SEED = SEED,
        NRAYS = NRAYS,
        )

    # print(u.info(debug=True))

    #
    # reset harmonic
    #
    harmonic_number=3 # zero for no changes

    if harmonic_number != 0:
        u.set_energy_monochromatic_at_resonance(harmonic_number=harmonic_number)
    #
    # u.EMAX = u.EMIN  + 5
    # u.NG_E = 1
    # #


    print(u.info(debug=True))

    os.system("rm begin.dat start.00 uphot.dat")
    beam = u.calculate_beam(code_undul_phot='internal',dump_undul_phot_file=True,dump_start_files=True)
    beam.write("begin.dat")


    tkt = Shadow.ShadowTools.plotxy("begin.dat",4,6,nbins=150)
    print("FWHM: h histo: %f, v histo: %f urad"%(1e6*tkt['fwhm_h'],1e6*tkt['fwhm_h']))

    if harmonic_number!= 0:
        print("FWHM: calc1: %f, calc2: %f urad"%(1e6*2.35*0.69*u.get_resonance_central_cone(harmonic_number),
              1e6*2.35*0.69*numpy.sqrt(u.get_resonance_wavelength(harmonic_number)/(u.NPERIODS*u.LAMBDAU))))


    from SourceUndulatorInputOutput import plot_undul_cdf,plot_undul_phot
    # plot_undul_cdf("xshundul.sha")
    plot_undul_phot("uphot.dat")