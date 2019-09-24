import numpy
import platform

from numpy.testing import assert_equal, assert_almost_equal

#
# shadow3
#
import Shadow

#
# minishadow
#
from shadow4.beam.beam import Beam
from shadow4.optical_surfaces.conic import Conic



def run_shadow3_from_start_files(iwrite=0):
    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe0_before_run = Shadow.Source()


    # TODO: this is a turn-around for the Linux bug...
    if platform.system() == "Linux":
        str = open('start.00', 'r').read()
        lines = str.split("\n")
        for line in lines:
            command = "oe0."+line
            try:
                exec(command)
            except:
                print("run_shadow3_from_start_files: Failed to exec: %s"%command)

            command = "oe0_before_run."+line
            try:
                exec(command)
            except:
                print("run_shadow3_from_start_files: Failed to exec: %s"%command)
    else:
        oe0.load("start.00")
        oe0_before_run.load("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")


    beam_source = beam.duplicate()
    # return beam,beam.rays.T.copy()

    oe1 = Shadow.OE()
    oe1_before_run = Shadow.OE()


    if platform.system() == "Linux":
        str = open('start.01', 'r').read()
        lines = str.split("\n")
        for line in lines:
            command1 = "oe1."+line.replace("(1)","[0]")\
            .replace("(2)","[1]")\
            .replace("(3)","[2]")\
            .replace("(4)","[3]")\
            .replace("(5)","[4]")\
            .replace("(6)","[5]")\
            .replace("(7)","[6]")\
            .replace("(8)","[7]")\
            .replace("(9)","[8]")\
            .replace("(10)","[9]")
            command2 = "oe1_before_run."+line.replace("(1)","[0]")\
            .replace("(2)","[1]")\
            .replace("(3)","[2]")\
            .replace("(4)","[3]")\
            .replace("(5)","[4]")\
            .replace("(6)","[5]")\
            .replace("(7)","[6]")\
            .replace("(8)","[7]")\
            .replace("(9)","[8]")\
            .replace("(10)","[9]")

            try:
                exec(command1)
            except:
                print("run_shadow3_from_start_files: Failed to exec: %s"%command1)

            try:
                exec(command2)
            except:
                print("run_shadow3_from_start_files: Failed to exec: %s"%command2)


    else:
        oe1.load("start.01")
        oe1_before_run.load("start.01")

    beam.traceOE(oe1,1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    return beam_source,beam,oe0_before_run,oe1_before_run



def compare_results(do_assert=True):
    Shadow.ShadowTools.plotxy("minimirr.01",2,1,nbins=101,nolost=1,title="Mirror (Python)")
    Shadow.ShadowTools.plotxy("mirr.01",2,1,nbins=101,nolost=1,title="Mirror (SHADOW)")

    Shadow.ShadowTools.plotxy("ministar.01",1,3,nbins=101,nolost=1,title="Image (Python)")
    Shadow.ShadowTools.plotxy("star.01",1,3,nbins=101,nolost=1,title="Image (SHADOW)")


    if do_assert:

        minimirr = Shadow.Beam()
        minimirr.load("minimirr.01")
        mirr     = Shadow.Beam()
        mirr.load("mirr.01")
        assert_almost_equal(minimirr.rays[:,0:6],mirr.rays[:,0:6])


        ministar = Shadow.Beam()
        ministar.load("ministar.01")
        star     = Shadow.Beam()
        star.load("star.01")
        assert_almost_equal(ministar.rays[:,0:6],star.rays[:,0:6])


def minishadow_run_conic_mirror():

    # ;
    # ; ray tracing of a single conic mirror using minishadow
    # ; results are compared with shadow3
    # ;
    #

    # ;
    # ; Runs shadow3
    #
    shadow3_beam_source,shadow3_beam,oe0,oe1 = run_shadow3_from_start_files(iwrite=1)


    # copy source to new Beam object
    newbeam = Beam.initialize_from_array(shadow3_beam_source.rays.copy())

    # ;
    # ; INPUTS
    # ;
    #
    fmirr         = oe1.FMIRR    # 1
    p             = oe1.T_SOURCE # 1000.0       # source-mirror
    q             = oe1.T_IMAGE  # 300.0        # mirror-image
    alpha         = oe1.ALPHA    # 0.0      # mirror orientation angle
    theta_grazing = (90.0-oe1.T_INCIDENCE) * numpy.pi / 180  # 5e-3     # grazing angle, rad
    fcyl          = oe1.FCYL
    f_convex      = oe1.F_CONVEX

    print("fmirr = %s, p=%f, q=%f, alpha=%f, theta_grazing=%f rad, fcyl=%d"%\
              (fmirr,p,q,alpha,theta_grazing,fcyl))

    # ccc = SurfaceConic()
    # ccc.set_sphere_from_focal_distances(p,q,theta_grazing,itype=fmirr,cylindrical=fcyl)


    if fmirr == 1: # sphere
        ccc = Conic.initialize_as_sphere_from_focal_distances(p,q,theta_grazing,cylindrical=fcyl,switch_convexity=f_convex)
    elif fmirr == 2: # Ellipsoid
        ccc = Conic.initialize_as_ellipsoid_from_focal_distances(p,q,theta_grazing,cylindrical=fcyl,switch_convexity=f_convex)
    elif fmirr == 3: # Toroidal
        raise Exception("fmirr=3 (toroidal) is NOT A CONIC surface")
    elif fmirr == 4: # Paraboloid
        ccc = Conic.initialize_as_paraboloid_from_focal_distances(p,q,theta_grazing,cylindrical=fcyl,switch_convexity=f_convex)
    elif fmirr == 5:
        ccc = Conic.initialize_as_plane()
    elif fmirr == 6: # Codling slit
        raise Exception("fmirr=6 (Codling slit) is NOT A CONIC surface")
    elif fmirr == 7: # Hyperboloid
        ccc = Conic.initialize_as_hyperboloid_from_focal_distances(p,q,theta_grazing,cylindrical=fcyl,switch_convexity=f_convex)
    elif fmirr == 8: # Cone
        raise Exception("fmirr=8 (Cone) is NOT A CONIC surface")
    else:
        raise Exception("fmirr invalid")

    #
    # put beam in mirror reference system
    #
    # TODO: calculate rotation matrices? Invert them for putting back to the lab system?

    newbeam.rotate(alpha,axis=2)
    newbeam.rotate(theta_grazing,axis=1)
    newbeam.translation([0.0,-p*numpy.cos(theta_grazing),p*numpy.sin(theta_grazing)])

    #
    # reflect beam in the mirror surface and dump mirr.01
    #
    newbeam = ccc.apply_specular_reflection_on_beam(newbeam)
    newbeam.dump_shadow3_file('minimirr.01')

    #
    # put beam in lab frame and compute image
    #
    newbeam.rotate(theta_grazing,axis=1)
    # TODO what about alpha?
    newbeam.retrace(q,resetY=True)
    newbeam.dump_shadow3_file('ministar.01')

if __name__ == "__main__":
    minishadow_run_conic_mirror()
    compare_results()
