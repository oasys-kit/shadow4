
import numpy

from srxraylib.plot.gol import plot,plot_scatter, set_qt

from shadow4.syned.magnetic_structure_1D_field import MagneticStructure1DField

from syned.storage_ring.electron_beam import ElectronBeam

from shadow4.sources.wiggler.source_wiggler import SourceWiggler

from shadow4.beam.beam import Beam

from scipy.ndimage import gaussian_filter1d


def create_file_with_magnetic_field(filename):

    L = 1605.0 #mm

    y = numpy.linspace(0,L, 2000)

    B = y * 0.0


    for i in range(y.size):
        if y[i] > 75 and y[i] < 575: B[i] = -0.876
        if y[i] > 650 and y[i] < 975: B[i] = 0.16
        if y[i] > 1030 and y[i] < 1530: B[i] = -0.8497

    # plot(y, B)



    B2 = gaussian_filter1d(B, 2.5)

    y -= y[y.size//2]
    y *= 1e-3

    plot(y, B, y, B2, legend=["original","smoothed"],xtitle="y / m",ytitle="B / T")

    f = open(filename, "w")
    for i in range(y.size):
        f.write("%f  %f\n" % (y[i], B2[i]))
    f.close()
    print("File written to disk: %s"%filename)




if __name__ == "__main__":

    set_qt()

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
    trajFile = ""
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

    filename = "BM_multi.b"
    create_file_with_magnetic_field(filename)
    syned_wiggler = MagneticStructure1DField.initialize_from_file(filename)
    # syned_wiggler.add_spatial_shift(-0.478)
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


    sourcewiggler.set_electron_initial_conditions_by_label(position_label="value_at_zero",
                                                           velocity_label="value_at_zero")


    print(sourcewiggler.info())


    rays = sourcewiggler.calculate_rays(NRAYS=NRAYS)

    plot_scatter(rays[:,0]*1e6,rays[:,2]*1e6,xtitle="X um",ytitle="Z um")
    plot_scatter(rays[:,1],rays[:,0]*1e6,xtitle="Y m",ytitle="X um")
    plot_scatter(rays[:,1],rays[:,2]*1e6,xtitle="Y m",ytitle="Z um")
    plot_scatter(rays[:,3]*1e6,rays[:,5]*1e6,xtitle="X' urad",ytitle="Z' urad")

    Beam.initialize_from_array(rays).write("begin.h5")


