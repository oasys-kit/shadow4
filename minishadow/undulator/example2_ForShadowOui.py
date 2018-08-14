#
# examples for SourceUndulator to be used in ShadowOui
#

import Shadow
import numpy
import os
# from SourceUndulatorInputOutput import plot_undul_cdf,plot_undul_phot
# from SourceUndulatorInputOutput import write_file_undul_phot_h5, write_file_undul_cdf_h5
# from SourceUndulator import SourceUndulator

from syned.storage_ring.electron_beam import ElectronBeam
from syned.storage_ring.magnetic_structures.undulator import Undulator
# from syned.storage_ring.light_source import LightSource
from SampleUndulator import SampleUndulator
if __name__ == "__main__":

    #
    # syned
    #
    su = Undulator.initialize_as_vertical_undulator(K=0.25,period_length=0.032,periods_number=50)

    ebeam = ElectronBeam(energy_in_GeV=6.04,
                 energy_spread = 0.0,
                 current = 0.2,
                 number_of_bunches = 400,
                 moment_xx=(400e-6)**2,
                 moment_xxp=0.0,
                 moment_xpxp=(10e-6)**2,
                 moment_yy=(10e-6)**2,
                 moment_yyp=0.0,
                 moment_ypyp=(4e-6)**2 )

    sampleundulator = SampleUndulator(name="test",syned_electron_beam=ebeam,syned_undulator=su,
                    # FLAG_EMITTANCE=0,EMIN=10000.0,EMAX=10650.0,NG_E=13,
                    FLAG_EMITTANCE=0,EMIN=10490.0,EMAX=10510.0,NG_E=13,
                    MAXANGLE=0.015,NG__T=51,NG_P=11,NG_J=20,
                    SEED=36255,NRAYS=15000)




    # sampleundulator.set_energy_monochromatic_at_resonance(1)

    print(sampleundulator.info())


    #
    # Calculate Radiation
    #

    # # via intermediate radiation (without writing uphot.dat)
    #
    # dict_radiation = sampleundulator.calculate_radiation(code_undul_phot='internal') # internal, pysru, srw
    #
    # # plot_undul_phot(dict_radiation)
    # # write_file_undul_phot_h5(dict_radiation,file_out="tmp.h5",mode='w',entry_name="radiation")
    # # write_file_undul_phot_h5(dict_radiation,mode='a',entry_name="radiation2")
    #
    # for k in dict_radiation.keys():
    #     try:
    #         print(k,dict_radiation[k].shape)
    #     except:
    #         pass
    #
    # dict_cdf = sampleundulator.calculate_cdf()
    # #
    #
    #
    # for k in dict_radiation.keys():
    #     print(k)
    #     try:
    #         print("    >>>",k,dict_cdf[k].shape)
    #     except:
    #         pass
    #
    # # write_file_undul_cdf_h5(dict_cdf,file_out="tmp.h5",mode='a',entry_name="cdf")
    # # plot_undul_cdf(dict_cdf)
    #
    # sampleundulator.write_file_h5("tmp.h5")
    #
    # # plot_undul_phot(dict_radiation,do_plot_intensity=True,do_plot_polarization=False)
    beam = sampleundulator.calculate_shadow3_beam()
    #
    # beam.write("begin.dat")



