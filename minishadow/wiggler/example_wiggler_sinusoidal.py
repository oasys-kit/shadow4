
import numpy


from srxraylib.plot.gol import plot, plot_scatter, plot_image

from syned.storage_ring.magnetic_structures.wiggler import Wiggler
from syned.storage_ring.electron_beam import ElectronBeam

from minishadow.wiggler.source_wiggler import SourceWiggler

if __name__ == "__main__":


    e_min = 5000.0 # 70490.0 #
    e_max = 100000.0 # 70510.0 #
    e_min = 70490.0 #
    e_max = 70510.0 #
    NRAYS = 5000
    use_emittances=True



    wigFile = "xshwig.sha"
    inData = ""

    nPer = 5 # 50
    nTrajPoints = 501
    ener_gev = 6.04
    per = 0.040
    kValue = 7.85
    trajFile = "tmp.traj"
    shift_x_flag = 0
    shift_x_value = 0.0
    shift_betax_flag = 0
    shift_betax_value = 0.0


    sw = SourceWiggler()

    #
    # syned
    #
    syned_wiggler = Wiggler(K_vertical=kValue,K_horizontal=0.0,period_length=per,number_of_periods=nPer)

    syned_electron_beam = ElectronBeam(energy_in_GeV=6.04,
                 energy_spread = 0.0,
                 current = 0.2,
                 number_of_bunches = 400,
                 moment_xx=(400e-6)**2,
                 moment_xxp=0.0,
                 moment_xpxp=(10e-6)**2,
                 moment_yy=(10e-6)**2,
                 moment_yyp=0.0,
                 moment_ypyp=(4e-6)**2 )

    sourcewiggler = SourceWiggler(name="test",syned_electron_beam=syned_electron_beam,
                    syned_wiggler=syned_wiggler,
                    flag_emittance=use_emittances,
                    emin=e_min,emax=e_max,ng_e=10, ng_j=nTrajPoints)



    print(sourcewiggler.info())


    rays = sourcewiggler.calculate_rays(NRAYS=NRAYS)

    plot_scatter(rays[:,1],rays[:,0],title="trajectory",show=False)
    plot_scatter(rays[:,0],rays[:,2],title="real space",show=False)
    plot_scatter(rays[:,3],rays[:,5],title="divergence space")
