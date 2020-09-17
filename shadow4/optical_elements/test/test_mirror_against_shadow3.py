#
# this test corresponds to the branch4 (using prerefl) in  test_shadow4_mirrors.ows in shadow4/oasys_workspaces
#


import numpy

from syned.beamline.optical_elements.mirrors.mirror import Mirror as SyMirror

from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.syned.shape import Plane, Sphere, Ellipsoid, Paraboloid, Hyperboloid # TODO from syned.beamline.shape

from shadow4.beam.beam import Beam


from Shadow.ShadowTools import plotxy
from shadow4.compatibility.beam3 import Beam3

from shadow4.syned.shape import MultiplePatch

from shadow4.syned.shape import Rectangle, Ellipse, TwoEllipses # TODO from syned.beamline.shape

from shadow4.optical_elements.s4_mirror import S4Mirror, S4MirrorElement

from numpy.testing import assert_almost_equal

def run_shadow3():
    #
    # Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
    #
    import Shadow
    import numpy

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 0

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    source = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 3
    oe0.F_COHER = 1
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 5676561
    oe0.NPOINT = 100000
    oe0.PH1 = 5.0
    oe0.POL_DEG = 0.5
    oe0.SIGDIX = 1e-06
    oe0.SIGDIZ = 1e-06
    oe0.SIGMAX = 0.0
    oe0.SIGMAZ = 0.0
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0

    oe1.DUMMY = 100.0
    oe1.FHIT_C = 1
    oe1.FILE_REFL = b'/Users/srio/Oasys/SiC.dat'
    oe1.FWRITE = 0
    oe1.F_REFLEC = 0 #1
    oe1.RLEN1 = 5e-05
    oe1.RLEN2 = 5e-05
    oe1.RWIDX1 = 2e-05
    oe1.RWIDX2 = 2e-05
    oe1.T_IMAGE = 6.0
    oe1.T_INCIDENCE = 88.8
    oe1.T_REFLECTION = 88.8

    # Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    source.genSource(oe0)
    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    #
    # run optical element 1
    #
    print("    Running optical element: %d" % (1))
    if iwrite:
        oe1.write("start.01")

    beam.traceOE(oe1, 1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    return source, beam, oe1

if __name__ == "__main__":

    from srxraylib.plot.gol import set_qt
    set_qt()
    do_plot = False
    do_assert = True

    OASYS_HOME = "/Users/srio/Oasys/"

    source3, beam3, oe1 = run_shadow3()


    beam0 = Beam.initialize_from_array(source3.rays)

    #
    # syned definitopns
    #

    # surface shape
    surface_shape = Plane()

    # boundaries
    rlen1 = 5e-05
    rlen2 = 5e-05
    rwidx1 = 2e-05
    rwidx2 = 2e-05
    boundary_shape = Rectangle(x_left=-rwidx2,x_right=rwidx1,y_bottom=-rlen2,y_top=rlen1)

    symirror1 = SyMirror(
                name="M1",
                surface_shape=surface_shape,
                boundary_shape=boundary_shape,
                coating=None, #"%s/SiC.dat" % OASYS_HOME,
                coating_thickness=None)

    coordinates_syned = ElementCoordinates(p = 10.0,
                                           q = 6.0,
                                           angle_radial = 88.8 * numpy.pi / 180,)


    #
    # shadow definitions
    #
    mirror1 = S4Mirror(name="M1",
                surface_shape=surface_shape,
                boundary_shape=boundary_shape)
                # element_coordinates_syned=coordinates_syned)
    print(mirror1.info())

    mirror1e = S4MirrorElement(optical_element=mirror1, coordinates=coordinates_syned)
    #
    # run
    #
    beam1, mirr1 = mirror1e.trace_beam(beam0, flag_lost_value=-11000.0)
    print(mirror1e.info())

    #
    # check
    #

    if do_plot:
        plotxy(beam3, 1, 3, title="Image 1 shadow3", nbins=101, nolost=1)
        beam1s3 = Beam3.initialize_from_shadow4_beam(beam1)
        plotxy(beam1s3, 1, 3, title="Image 1 shadow4", nbins=101, nolost=1)


    mirr3 = Beam3(N=beam0.rays.shape[0])
    mirr3.load("mirr.01")

    print("\ncol#   m-shadow4  m-shadow3  source")
    for i in range(18):
        print("col%d   %20.10f  %20.10f  %20.10f  " % (i+1, mirr1.rays[10,i], mirr3.rays[10,i], source3.rays[10,i]))
        if do_assert:
            assert_almost_equal (mirr1.rays[:,i], mirr3.rays[:,i], 1)

    print("\ncol#   shadow4  shadow3  source")
    for i in range(18):
        print("col%d   %20.10f  %20.10f  %20.10f  " % (i+1, beam1.rays[10,i], beam3.rays[10,i], source3.rays[10,i]))
        if do_assert:
            assert_almost_equal (beam1.rays[:,i], beam3.rays[:,i], 1)