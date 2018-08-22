#
# examples for SourceUndulator to be used in ShadowOui
#

import os
from syned.storage_ring.electron_beam import ElectronBeam
from syned.storage_ring.magnetic_structures.undulator import Undulator
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
                    FLAG_EMITTANCE=0,FLAG_SIZE=0,
                    EMIN=10490.0,EMAX=10510.0,NG_E=13,
                    MAXANGLE=0.015,NG_T=51,NG_P=11,NG_J=20,
                    SEED=36255,NRAYS=15000,
                    code_undul_phot="internal")




    # sampleundulator.set_energy_monochromatic_at_resonance(0.98)

    print(sampleundulator.info())


    beam = sampleundulator.calculate_shadow3_beam()

    # print(sampleundulator.info())


    os.system("rm -f begin.dat start.00 end.00")
    beam.write("begin.dat")
    print("File written to disk: begin.dat")




