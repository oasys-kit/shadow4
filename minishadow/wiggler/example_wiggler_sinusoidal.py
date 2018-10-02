
import numpy
import srxraylib.sources.srfunc as srfunc

from srxraylib.plot.gol import plot

from syned.storage_ring.magnetic_structures.wiggler import Wiggler
from syned.storage_ring.electron_beam import ElectronBeam

from source_wiggler import SourceWiggler

def run_python_preprocessors():

    wigFile = "xshwig.sha"
    inData = ""

    nPer = 50
    nTrajPoints = 501
    ener_gev = 6.04
    per = 0.040
    kValue = 7.85
    trajFile = "tmp.traj"
    shift_x_flag = 0
    shift_x_value = 0.0
    shift_betax_flag = 0
    shift_betax_value = 0.0

    e_min = 5000.0
    e_max = 100000.0

    # ("Calculate electron trajectory"
    (traj, pars) = srfunc.wiggler_trajectory(b_from=0,
                                             inData=inData,
                                             nPer=nPer,
                                             nTrajPoints=nTrajPoints,
                                             ener_gev=ener_gev,
                                             per=per,
                                             kValue=kValue,
                                             trajFile="tmp.traj",
                                             shift_x_flag=shift_x_flag,
                                             shift_x_value=shift_x_value,
                                             shift_betax_flag=shift_betax_flag,
                                             shift_betax_value=shift_betax_value)


    # traj[0,ii] = yx[i]
    # traj[1,ii] = yy[i]+j * per - start_len
    # traj[2,ii] = 0.0
    # traj[3,ii] = betax[i]
    # traj[4,ii] = betay[i]
    # traj[5,ii] = 0.0
    # traj[6,ii] = curv[i]
    # traj[7,ii] = bz[i]

    # plot(traj[1,:],traj[0,:])
    # print(pars)

    #
    # calculate cdf and write file for Shadow/Source
    #

    srfunc.wiggler_cdf(traj,
                       enerMin=e_min,
                       enerMax=e_max,
                       enerPoints=1001,
                       outFile=wigFile,
                       elliptical=False)

    print("CDF written to file %s \n"%(str(wigFile)))

def run_shadow3_source(ener_gev=6.04,use_emittances=True):
    import Shadow

    oe0 = Shadow.Source()
    beam = Shadow.Beam()


    oe0.BENER = ener_gev
    oe0.CONV_FACT = 100.0
    oe0.FDISTR = 0
    oe0.FILE_TRAJ = b'xshwig.sha'
    oe0.FSOUR = 0
    oe0.FSOURCE_DEPTH = 0
    oe0.F_COLOR = 0
    oe0.F_PHOT = 0
    oe0.F_WIGGLER = 1
    oe0.HDIV1 = 1.0
    oe0.HDIV2 = 1.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 5676561
    oe0.NCOL = 0
    oe0.N_COLOR = 0
    oe0.PH1 = 5000.0
    oe0.PH2 = 100000.0
    oe0.POL_DEG = 0.0

    if use_emittances:
        oe0.SIGMAX = 0.0078
        oe0.SIGMAY = 0.0
        oe0.SIGMAZ = 0.0036
        oe0.EPSI_X = 3.8e-07
        oe0.EPSI_Z = 3.8e-09
    else:
        oe0.SIGMAX = 0.0
        oe0.SIGMAY = 0.0
        oe0.SIGMAZ = 0.0
        oe0.EPSI_X = 0.0
        oe0.EPSI_Z = 0.0

    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0

    beam.genSource(oe0)

    Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")



if __name__ == "__main__":

    #
    # current way
    #

    # run_python_preprocessors()
    # run_shadow3_source(ener_gev=6.04,use_emittances=False)

    #
    # new way
    #

    wigFile = "xshwig.sha"
    inData = ""

    nPer = 50
    nTrajPoints = 501
    ener_gev = 6.04
    per = 0.040
    kValue = 7.85
    trajFile = "tmp.traj"
    shift_x_flag = 0
    shift_x_value = 0.0
    shift_betax_flag = 0
    shift_betax_value = 0.0

    e_min = 5000.0
    e_max = 100000.0

    sw = SourceWiggler()

    #
    # syned
    #
    syned_wiggler = Wiggler()

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
                    flag_emittance=1,
                    emin=10490.0,emax=10510.0,ng_e=3)


    sourcewiggler.set_energy_monochromatic(5000.0)

    print(sourcewiggler.info())

