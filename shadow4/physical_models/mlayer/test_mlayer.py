
from srxraylib.plot.gol import plot
from shadow4.physical_models.mlayer.mlayer import MLayer
import numpy
from numpy.testing import assert_almost_equal

import os

DO_PLOTS = True
COMPARE_WITH_SHADOW3  = True
SHADOW3_PATH = "/home/manuel/OASYS1.2/shadow3/shadow3"

def test_xoppy_defauls():

    a = MLayer.pre_mlayer(
        interactive=False,
        FILE="pre_mlayer.dat",
        E_MIN=5000.0, E_MAX=10000.0,
        O_DENSITY=2.33, O_MATERIAL="Si",  # odd: closer to vacuum
        E_DENSITY=19.3, E_MATERIAL="W",  # even: closer to substrate
        S_DENSITY=2.33, S_MATERIAL="Si",  # substrate
        GRADE_DEPTH=0,
        N_PAIRS=50,
        THICKNESS=50.0,
        GAMMA=0.5,  # gamma ratio  =  t(even) / (t(odd) + t(even))")
        ROUGHNESS_EVEN=0.0,
        ROUGHNESS_ODD=0.0,
        FILE_DEPTH="myfile_depth.dat",
        GRADE_SURFACE=0,
        FILE_SHADOW="mlayer1.sha",
        FILE_THICKNESS="mythick.dat",
        FILE_GAMMA="mygamma.dat",
        AA0=1.0, AA1=0.0, AA2=0.0, AA3=0.0)


    #
    # theta scan
    #
    rs, rp, e, t = a.scan(fileOut=None,  # "pre_mlayer_scan.dat",
                          energyN=1, energy1=8050.0, energy2=8050.0,
                          thetaN=600, theta1=0.0, theta2=6.0)

    print(rs.shape, rp.shape, e.shape, t.shape)

    if DO_PLOTS:
        plot(t, rs[0], xtitle="angle [deg]", ytitle="Reflectivity", title="Default xoppy",  ylog=False)
        # a.plot_optical_constants()


    if COMPARE_WITH_SHADOW3:

        f = open("shadow3.inp",'w')
        f.write("pre_mlayer_scan\npre_mlayer.dat\n1\n600\n8050\n0.0\n6.0\nexit\n")
        f.close()

        os.system(SHADOW3_PATH+" < shadow3.inp")

        s3 = numpy.loadtxt("pre_mlayer_scan.dat",skiprows=5)
        print(s3.shape)
        if DO_PLOTS:
            plot(s3[:,1],s3[:,2],t, rs[0], xtitle="angle [deg]", ytitle="Reflectivity", title="Default xoppy",
                 legend = ["shadow3 (pre_mlayer_scan)","shadow4 (MLayer)"])

        print(s3[:,1].shape,s3[:,2].shape,t.shape, rs[0].shape)
        assert_almost_equal(s3[:,2], rs[0])

def test_aw1():


    a = MLayer.pre_mlayer(
        interactive=False,
        FILE="pre_mlayer.dat",
        E_MIN=100.0, #100.0,
        E_MAX=5000.0, #1000.0,
        O_DENSITY=7.19, O_MATERIAL="Cr",  # odd: closer to vacuum
        E_DENSITY=3.00, E_MATERIAL="Sc",  # even: closer to substrate
        S_DENSITY=2.33, S_MATERIAL="Si",  # substrate
        GRADE_DEPTH=0,
        N_PAIRS=50,
        THICKNESS=22.0,
        GAMMA=10.0/22.0,  #  gamma ratio  =  t(even) / (t(odd) + t(even))")
        ROUGHNESS_EVEN=0.0,
        ROUGHNESS_ODD=0.0,
        FILE_DEPTH="myfile_depth.dat",
        GRADE_SURFACE=0,
        FILE_SHADOW="mlayer1.sha",
        FILE_THICKNESS="mythick.dat",
        FILE_GAMMA="mygamma.dat",
        AA0=1.0,AA1=0.0,AA2=0.0,AA3=0.0)

    #
    # energy scan
    #
    rs, rp, e, t = a.scan(h5file="", #"pre_mlayer_scan.dat",
            energyN=300,
            energy1=300.0, #300.0,
            energy2=500.0, #500.0,
            thetaN=1,
            theta1=45.0,
            theta2=45.0,
            )

    print(rs.shape,rp.shape,e.shape,t.shape)

    if DO_PLOTS:
        plot(e,rs[:,0],xtitle="Photon energy [eV]",ytitle="Reflectivity",title="Example AW",)

    if COMPARE_WITH_SHADOW3:

        f = open("shadow3.inp",'w')
        f.write("pre_mlayer_scan\npre_mlayer.dat\n300\n1\n300.0\n500.0\n45.0\nexit\n")
        f.close()

        os.system(SHADOW3_PATH+" < shadow3.inp")

        s3 = numpy.loadtxt("pre_mlayer_scan.dat",skiprows=5)
        print(s3.shape)
        if DO_PLOTS:
            plot(s3[:,0],s3[:,2],e, rs[:,0], xtitle="energy [eV]", ytitle="Reflectivity", title="Example AW",
                 legend = ["shadow3 (pre_mlayer_scan)","shadow4 (MLayer)"])

        print(s3[:,0].shape,s3[:,2].shape,t.shape, rs[:,0].shape)
        assert_almost_equal(s3[:,2],rs[:,0],6)



def test_substrate():


    a = MLayer.pre_mlayer(
        interactive=False,
        FILE="pre_mlayer.dat",
        E_MIN=100.0, #100.0,
        E_MAX=10000.0, #1000.0,
        O_DENSITY=7.19, O_MATERIAL="Cr",  # odd: closer to vacuum
        E_DENSITY=3.00, E_MATERIAL="Sc",  # even: closer to substrate
        S_DENSITY=2.33, S_MATERIAL="Si",  # substrate
        GRADE_DEPTH=0,
        N_PAIRS=0,
        THICKNESS=22.0,
        GAMMA=10.0/22.0,  #  gamma ratio  =  t(even) / (t(odd) + t(even))")
        ROUGHNESS_EVEN=0.0,
        ROUGHNESS_ODD=0.0,
        FILE_DEPTH="myfile_depth.dat",
        GRADE_SURFACE=0,
        FILE_SHADOW="mlayer1.sha",
        FILE_THICKNESS="mythick.dat",
        FILE_GAMMA="mygamma.dat",
        AA0=1.0,AA1=0.0,AA2=0.0,AA3=0.0)

    # a.read_preprocessor_file("pre_mlayer.dat")
    #
    # energy scan
    #
    rs, rp, e, t = a.scan(h5file="", #"pre_mlayer_scan.dat",
            energyN=1000,
            energy1=300.0, #300.0,
            energy2=6000.0, #500.0,
            thetaN=1,
            theta1=10e-3*180.0/numpy.pi,
            theta2=10e-3*180.0/numpy.pi,
            )

    print(rs.shape,rp.shape,e.shape,t.shape)

    if DO_PLOTS:
        plot(e,rs[:,0],xtitle="Photon energy [eV]",ytitle="Reflectivity",title="Only substrate")


    assert( rs[:,0][0] > 0.9)
    assert (rs[:, 0][-1] < 0.1)


def test_no_preprocessor():
    out = MLayer.initialize_from_bilayer_stack(
            material_S="Si", density_S=None, roughness_S=0.0,  #  2.33
            material_E="W" , density_E=None, roughness_E=0.0,  #  19.3
            material_O="Si", density_O=None, roughness_O=0.0,  #  2.33
            bilayer_pairs=50,
            bilayer_thickness=50.0,
            bilayer_gamma=0.5,
    )

    for key in out.pre_mlayer_dict.keys():
        print(key,out.pre_mlayer_dict[key])
    #
    # theta scan
    #

    #
    rs, rp, e, t = out.scan(h5file="",  # "pre_mlayer_scan.dat",
                          energyN=1, energy1=8050.0, energy2=8050.0,
                          thetaN=600, theta1=0.0, theta2=6.0)
    #
    print(rs.shape, rp.shape, e.shape, t.shape)
    #
    if DO_PLOTS:
        plot(t, rs[0], xtitle="angle [deg]", ytitle="Reflectivity", title="Default xoppy - no preprocessor (direct calling xraylib)", ylog=False)
        # a.plot_optical_constants()

    if COMPARE_WITH_SHADOW3:
        a = MLayer.pre_mlayer(
            interactive=False,
            FILE="pre_mlayer.dat",
            E_MIN=5000.0, E_MAX=10000.0,
            O_DENSITY=2.33, O_MATERIAL="Si",  # odd: closer to vacuum
            E_DENSITY=19.3, E_MATERIAL="W",  # even: closer to substrate
            S_DENSITY=2.33, S_MATERIAL="Si",  # substrate
            GRADE_DEPTH=0,
            N_PAIRS=50,
            THICKNESS=50.0,
            GAMMA=0.5,  # gamma ratio  =  t(even) / (t(odd) + t(even))")
            ROUGHNESS_EVEN=0.0,
            ROUGHNESS_ODD=0.0,
            FILE_DEPTH="myfile_depth.dat",
            GRADE_SURFACE=0,
            FILE_SHADOW="mlayer1.sha",
            FILE_THICKNESS="mythick.dat",
            FILE_GAMMA="mygamma.dat",
            AA0=1.0, AA1=0.0, AA2=0.0, AA3=0.0)

        f = open("shadow3.inp",'w')
        f.write("pre_mlayer_scan\npre_mlayer.dat\n1\n600\n8050\n0.0\n6.0\nexit\n")
        f.close()

        os.system(SHADOW3_PATH+" < shadow3.inp")

        s3 = numpy.loadtxt("pre_mlayer_scan.dat",skiprows=5)
        print(s3.shape)
        if DO_PLOTS:
            plot(s3[:,1],s3[:,2],t, rs[0], xtitle="angle [deg]", ytitle="Reflectivity", title="Default xoppy",
                 legend = ["shadow3 (pre_mlayer_scan)","shadow4 (MLayer - no preprocessor)",])

        print(s3[:,1].shape,s3[:,2].shape,t.shape, rs[0].shape)
        assert_almost_equal(s3[:,2], rs[0],3)

def test_xoppy():
    # from orangecontrib.xoppy.util.mlayer import MLayer

    out = MLayer.initialize_from_bilayer_stack(
        material_S="Si", density_S=None, roughness_S=0.0,
        material_E="W", density_E="10", roughness_E=0.0,
        material_O="Si", density_O=None, roughness_O=0.0,
        bilayer_pairs=50,
        bilayer_thickness=50.0,
        bilayer_gamma=0.5,
    )

    for key in out.pre_mlayer_dict.keys():
        print(key, out.pre_mlayer_dict[key])
    #
    rs, rp, e, t = out.scan(h5file="",
                            energyN=1, energy1=8050.0, energy2=15000.0,
                            thetaN=600, theta1=0.0, theta2=6.0)

    #
    # plot (example)
    #
    myscan = 0
    from srxraylib.plot.gol import plot, plot_image

    if myscan == 0:  # angle scan
        plot(t, rs[0], xtitle="angle [deg]", ytitle="Reflectivity-s", title="")
    elif myscan == 1:  # energy scan
        plot(e, rs[:, 0], xtitle="Photon energy [eV]", ytitle="Reflectivity-s", title="")
    elif myscan == 2:  # double scan
        plot_image(rs, e, t, xtitle="Photon energy [eV]", ytitle="Grazing angle [deg]", title="Reflectivity-s",
                   aspect="auto")


if __name__ == "__main__":

    #
    # test_xoppy_defauls()
    #
    test_aw1()
    #
    test_substrate()

    test_no_preprocessor()

    test_xoppy()

