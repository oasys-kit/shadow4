"""
Tools to compare beams from shadow3 and Shadow4
"""
import numpy
from srxraylib.plot.gol import plot_scatter
import Shadow
from shadow4.beam.beam import Beam

def compare_shadow3_and_shadow4_beams(beam3,beam4,do_plot=True,do_assert=False,assert_value=1e-2):

    return compare_shadow3_and_shadow4_rays(beam3.rays,beam4.rays,do_plot=do_plot,do_assert=do_assert,
                                            assert_value=assert_value)


def compare_shadow3_and_shadow4_rays(rays3,rays4,do_plot=True,do_assert=False,assert_value=1e-2):

    raysnew = rays4
    rays = rays3

    if do_plot:
        # plot_scatter(rays[:,1],rays[:,0],title="Trajectory shadow3",show=False)
        # plot_scatter(raysnew[:,1],raysnew[:,0],title="Trajectory shadow4")


        plot_scatter(rays[:,3],rays[:,5],title="Divergences shadow3",show=False)
        plot_scatter(raysnew[:,3],raysnew[:,5],title="Divergences shadow4")

        plot_scatter(rays[:,0],rays[:,2],title="Real Space shadow3",show=False)
        plot_scatter(raysnew[:,0],raysnew[:,2],title="Real Space shadow4")

        #
        b3 = Shadow.Beam()
        b3.rays = rays

        b4 = Shadow.Beam()
        b4.rays = raysnew
        Shadow.ShadowTools.histo1(b3,11,ref=23,nolost=1)
        Shadow.ShadowTools.histo1(b4,11,ref=23,nolost=1)



    print("Comparing...")
    for i in range(6):
        if i <= 2:
            fact = 1e-3 # shadow3 units to m
        else:
            fact = 1.0
        m0 = (raysnew[:,i]).mean()
        m1 = (rays[:,i]*fact).mean()
        print("\ncol %d, mean sh3, sh4, |sh4-sh3|: %10g  %10g  %10g"%(i+1,m1,m0,numpy.abs(m0-m1)))
        std0 = raysnew[:,i].std()
        std1 = (rays[:,i]*fact).std()
        print("col %d, stdv sh3, sh4, |sh4-sh3|: %10g  %10g  %10g"%(i+1,std1,std0,numpy.abs(std0-std1)))

        if do_assert:
            assert(numpy.abs(m0-m1) < assert_value)
            assert(numpy.abs(std0-std1) < assert_value)