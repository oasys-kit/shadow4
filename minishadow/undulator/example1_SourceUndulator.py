#
# examples for SourceUndulator
#


import Shadow
import numpy
import os
from SourceUndulatorInputOutput import plot_undul_cdf,plot_undul_phot
from SourceUndulator import SourceUndulator


def main1():
    u = SourceUndulator()

    u.set_from_keywords(E_ENERGY=6.04,INTENSITY=0.2,FLAG_EMITTANCE=1,SX=0.04,SZ=0.001,SXP=10e-6,SZP=4e-6,
                        LAMBDAU=0.032,NPERIODS=50,K=0.25,
                        EMIN= 10500.0000,EMAX= 10550.0000,
                        NG_E=101,NG_T=51,NG_P=11,N_J=20,
                        MAXANGLE=0.1,
                        SEED=36255,NRAYS=15000)


    #
    # overwrite energy for selected harmonic
    #
    harmonic_number = 1 # zero for do nothing

    if harmonic_number != 0:
        u.set_energy_monochromatic_at_resonance(harmonic_number=harmonic_number)


    print(u.info(debug=True))

    os.system("rm begin.dat start.00 uphot.dat")

    method = 1 # 0=direct, 1=get intermediate radiation
    if method == 0:
        # direct calculation
        beam = u.calculate_shadow3_beam(code_undul_phot='internal',dump_undul_phot_file=True,dump_start_files=True)
        plot_undul_phot("uphot.dat",do_plot_intensity=True,do_plot_polarization=False)
    else:
        # via intermediate radiation (without writing uphot.dat)
        dict_radiation = u.calculate_radiation(code_undul_phot='srw') # internal, pysru, srw
        plot_undul_phot(dict_radiation,do_plot_intensity=True,do_plot_polarization=False)
        beam = u.calculate_shadow3_beam(use_existing_undul_phot_output=dict_radiation,dump_undul_phot_file=False,dump_start_files=True)

    beam.write("begin.dat")

    #
    # tkt = Shadow.ShadowTools.plotxy("begin.dat",4,6,nbins=150)
    # print("FWHM: h histo: %f, v histo: %f urad"%(1e6*tkt['fwhm_h'],1e6*tkt['fwhm_h']))

    if harmonic_number!= 0:
        print("FWHM: calculated: %f (or %f) urad"%(1e6*2.35*0.69*u.get_resonance_central_cone(harmonic_number),
              1e6*2.35*0.69*numpy.sqrt(u.get_resonance_wavelength(harmonic_number)/(u.NPERIODS*u.LAMBDAU))))


    # plot_undul_cdf("xshundul.sha")




if __name__ == "__main__":

    main1()


    tkt = Shadow.ShadowTools.plotxy("begin.dat",4,6,nbins=150)
    print("FWHM: h histo: %f, v histo: %f urad"%(1e6*tkt['fwhm_h'],1e6*tkt['fwhm_h']))