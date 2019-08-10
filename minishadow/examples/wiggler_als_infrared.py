
import numpy


from srxraylib.plot.gol import plot,plot_scatter


from syned.storage_ring.magnetic_structures.wiggler import Wiggler
from minishadow.wiggler.magnetic_structure_1D_field import MagneticStructure1DField

from syned.storage_ring.electron_beam import ElectronBeam

from minishadow.wiggler.source_wiggler import SourceWiggler

from scipy.interpolate import interp1d


if __name__ == "__main__":

    use_emittances=True
    e_min = 0.4
    e_max = 0.4
    NRAYS = 5000




    wigFile = "xshwig.sha"
    inData = ""

    nTrajPoints = 501
    ener_gev = 1.90
    per = 0.5
    kValue = 4
    trajFile = "tmp.traj"
    shift_x_flag = 0
    shift_x_value = 0.0
    shift_betax_flag = 0
    shift_betax_value = 0.0


    sw = SourceWiggler()

    #
    # syned
    #


    syned_electron_beam = ElectronBeam(energy_in_GeV=1.9,current=0.4,
                                       moment_xx=(39e-6)**2,
                                       moment_xpxp=(2000e-12 / 51e-6)**2,
                                       moment_yy=(31e-6)**2,
                                       moment_ypyp=(30e-12 / 31e-6)**2,
                                       )

    # conventional wiggler
    # syned_wiggler = Wiggler(K_vertical=kValue,K_horizontal=0.0,period_length=per,number_of_periods=nPer)

    # B from file
    # filename = "/home/manuel/Oasys/BM_smooth.b"
    filename = "/home/manuel/Oasys/BM_multi.b"
    syned_wiggler = MagneticStructure1DField.initialize_from_file(filename)
    syned_wiggler.add_spatial_shift(-0.478)
    # syned_wiggler.flip_B()



    if e_min == e_max:
        ng_e = 1
    else:
        ng_e = 10

    sourcewiggler = SourceWiggler(name="test",
                    syned_electron_beam=syned_electron_beam,
                    syned_wiggler=syned_wiggler,
                    flag_emittance=use_emittances,
                    emin=e_min,
                    emax=e_max,
                    ng_e=ng_e,
                    ng_j=nTrajPoints)

    # sourcewiggler.set_electron_initial_conditions_by_label(position_label="maximum",velocity_label="half_excursion")
    sourcewiggler.set_electron_initial_conditions_by_label(position_label="maximum",
                                                           velocity_label="user_value", velocity_value=0.094,)


    print(sourcewiggler.info())


    rays = sourcewiggler.calculate_rays(NRAYS=NRAYS)

    plot_scatter(rays[:,0]*1e6,rays[:,2]*1e6,xtitle="X um",ytitle="Z um")
    plot_scatter(rays[:,1],rays[:,0]*1e6,xtitle="Y m",ytitle="X um")
    plot_scatter(rays[:,1],rays[:,2]*1e6,xtitle="Y m",ytitle="Z um")
    plot_scatter(rays[:,3]*1e6,rays[:,5]*1e6,xtitle="X' urad",ytitle="Z' urad")


