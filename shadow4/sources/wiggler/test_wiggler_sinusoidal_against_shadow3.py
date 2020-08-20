
import numpy


from srxraylib.plot.gol import plot, plot_scatter, plot_image

from syned.storage_ring.magnetic_structures.wiggler import Wiggler
from syned.storage_ring.electron_beam import ElectronBeam

from shadow4.sources.wiggler.source_wiggler import SourceWiggler

def run_python_preprocessors(e_min=1000.0,e_max=10000.0 ):

    import srxraylib.sources.srfunc as srfunc

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

     # 5000.0
     # 100000.0

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

def run_shadow3_source(ener_gev=6.04,use_emittances=True,EMIN=10000.0,EMAX=11000.0,NRAYS=500):
    import Shadow

    oe0 = Shadow.Source()
    beam = Shadow.Beam()

    oe0.BENER = ener_gev
    oe0.CONV_FACT = 100.0
    oe0.FDISTR = 0
    oe0.FILE_TRAJ = b'xshwig.sha' # b'xshwig.sha'
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
    oe0.PH1 = EMIN
    oe0.PH2 = EMAX
    oe0.POL_DEG = 0.0
    oe0.NPOINT = NRAYS

    if use_emittances:
        moment_xx=(400e-6)**2
        moment_xxp=0.0
        moment_xpxp=(10e-6)**2
        moment_yy=(10e-6)**2
        moment_yyp=0.0
        moment_ypyp=(4e-6)**2
        oe0.SIGMAX = 1e2 * numpy.sqrt(moment_xx)
        oe0.SIGMAY = 0.0
        oe0.SIGMAZ = 1e2 * numpy.sqrt(moment_yy)
        oe0.EPSI_X = 1e2 * numpy.sqrt(moment_xx) * numpy.sqrt(moment_xpxp)
        oe0.EPSI_Z = 1e2 * numpy.sqrt(moment_yy) * numpy.sqrt(moment_ypyp)
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


    oe0.write("start.00")
    beam.genSource(oe0)
    beam.write("begin.dat")


    return beam

def compare_rays_with_shadow3_beam(raysnew,beam,do_plot=True,do_assert=False):

    import Shadow
    # Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")

    # a = Shadow.Beam()
    # a.load("begin.dat")
    # rays = a.rays

    rays = beam.rays

    if do_plot:
        plot_scatter(rays[:,1],rays[:,0],title="Trajectory shadow3",show=False)
        plot_scatter(raysnew[:,1],raysnew[:,0],title="Trajectory new")


        plot_scatter(rays[:,3],rays[:,5],title="Divergences shadow3",show=False)
        plot_scatter(raysnew[:,3],raysnew[:,5],title="Divergences new")

        plot_scatter(rays[:,0],rays[:,2],title="Real Space shadow3",show=False)
        plot_scatter(raysnew[:,0],raysnew[:,2],title="Real Space new")

        b = Shadow.Beam()
        b.rays = raysnew
        Shadow.ShadowTools.histo1(beam_shadow3,11,ref=23,nolost=1)
        Shadow.ShadowTools.histo1(b,11,ref=23,nolost=1)

    if do_assert:
        print("Comparing...")
        for i in range(6):
            if i <= 2:
                fact = 100
            else:
                fact = 1.0
            m0 = (raysnew[:,i]*fact).mean()
            m1 = beam.rays[:,i].mean()
            print("\ncol %d, mean: "%(i+1),m0,m1,numpy.abs(m0-m1)/numpy.abs(m1))
            std0 = (raysnew[:,i]*fact).std()
            std1 = beam.rays[:,i].std()
            print("col %d, std: "%(i+1),std0,std1,numpy.abs(std0-std1)/numpy.abs(std1))

            assert((numpy.abs(m0-m1)/numpy.abs(m1)) < 15.0)
            assert((numpy.abs(std0-std1)/numpy.abs(std1)) < 0.05)


if __name__ == "__main__":
    import Shadow
    import os

    e_min = 5000.0 # 70490.0 #
    e_max = 100000.0 # 70510.0 #
    e_min = 70490.0 #
    e_max = 70510.0 #
    NRAYS = 5000
    use_emittances=True


    #
    # current way
    #
    # os.system("rm begin.dat start.00 xshwig.sha")
    # os.system("pwd")
    run_python_preprocessors(e_min=e_min,e_max=e_max)
    beam_shadow3 = run_shadow3_source(ener_gev=6.04,EMIN=e_min,EMAX=e_max,NRAYS=NRAYS,use_emittances=use_emittances)
    # Shadow.ShadowTools.histo1("begin.dat",11)

    #
    # new way
    #

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


    # sourcewiggler.set_energy_monochromatic(5000.0)

    print(sourcewiggler.info())

    # sourcewiggler.calculate_radiation()

    rays = sourcewiggler.calculate_rays(NRAYS=NRAYS)

    compare_rays_with_shadow3_beam(rays,beam_shadow3,do_plot=False,do_assert=True)



    #
    # check new sync_f function
    #

    from numpy.testing import assert_almost_equal
    from shadow4.sources.wiggler.source_wiggler import sync_f_sigma_and_pi
    from srxraylib.sources.srfunc import sync_f
    rAngle = numpy.array([1,2,3,4])
    rEnergy = 1.3
    s,p = sync_f_sigma_and_pi(rAngle,rEnergy)
    s0 = sync_f(rAngle, rEnergy, polarization=1)
    p0 = sync_f(rAngle, rEnergy, polarization=2)
    assert_almost_equal(s,s0[:,0])
    assert_almost_equal(p, p0[:, 0])
    print(">>>>s: ",s,s0[:,0])
    print(">>>>p: ",p,p0[:,0])
