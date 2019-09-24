#
# uses example 21 of tutorials
#
import numpy
import platform
from numpy.testing import assert_equal, assert_almost_equal

import Shadow
# from Shadow.ShadowPreprocessorsXraylib import bragg

#
# minishadow
#
from shadow4.beam.beam import Beam
from shadow4.optical_surfaces.mesh import Mesh

def write_bragg_preprocessor_file():
    bragg(interactive=False, DESCRIPTOR="Si",H_MILLER_INDEX=1,K_MILLER_INDEX=1,L_MILLER_INDEX=1,
          TEMPERATURE_FACTOR=1.0,E_MIN=5000.0,E_MAX=55000.0,E_STEP=100.0,SHADOW_FILE="si5_55.111")


def write_shadow_surface(s,xx,yy,outFile='presurface.dat'):
    """
      write_shadowSurface: writes a mesh in the SHADOW/presurface format
      SYNTAX:
           out = write_shadowSurface(z,x,y,outFile=outFile)
      INPUTS:
           z - 2D array of heights
           x - 1D array of spatial coordinates along mirror width.
           y - 1D array of spatial coordinates along mirror length.

      OUTPUTS:
           out - 1=Success, 0=Failure
           outFile - output file in SHADOW format. If undefined, the
                     file is names "presurface.dat"

    """
    out = 1

    try:
       fs = open(outFile, 'w')
    except IOError:
       out = 0
       print ("Error: can\'t open file: "+outFile)
       return
    else:
        # dimensions
        fs.write( repr(xx.size)+" "+repr(yy.size)+" \n" )
        # y array
        for i in range(yy.size):
            fs.write(' ' + repr(yy[i]) )
        fs.write("\n")
        # for each x element, the x value and the corresponding z(y)
        # profile
        for i in range(xx.size):
            tmps = ""
            for j in range(yy.size):
                tmps = tmps + "  " + repr(s[j,i])
            fs.write(' ' + repr(xx[i]) + " " + tmps )
            fs.write("\n")
        fs.close()
        print ("write_shadow_surface: File for SHADOW "+outFile+" written to disk.")


def create_gaussian_bump(do_plot=True):
    # calculate a Gaussian bump
    # create an array of 2 cm length
    npoints = 51
    length = 2.0

    x = numpy.linspace(-0.5*length,0.5*length,npoints)
    y = numpy.linspace(-0.5*length,0.5*length,npoints)

    # create a surface with pixel value its distance to the center

    x1 = numpy.outer(x,numpy.ones(y.size))
    y1 = numpy.outer(numpy.ones(x.size),x)

    r = numpy.sqrt( (x1)**2 + (y1)**2 )

    # define bump FWHM
    bump_fwhm = 0.5 # cm

    #pizel sizes
    pixel = length / (npoints-1)

    # sigma value corresponding to FWHM
    sigma = (bump_fwhm / pixel) / ( 2*numpy.sqrt(2*numpy.log(2)) )
    print("sigma: ",sigma)

    # evaluate the 2D Gaussian
    z = numpy.exp(-((r/pixel)/2/sigma)**2)
    # give a heigth of 1 microns
    z = z * 1e-4

    #write file for SHADOW
    write_shadow_surface(z,x,y,outFile='bump.dat')

    if do_plot:
        #
        #plot
        #
        from matplotlib import pylab as plt
        plt.figure(1)
        plt4 = plt.imshow(z.T*1e4,extent=[-0.5*length,0.5*length,-0.5*length,0.5*length])
        plt.title('Gaussian Bump')
        plt.xlabel('H [cm]')
        plt.ylabel('V [cm]')
        cbar = plt.colorbar(plt4 , format="%.2f")
        cbar.ax.set_ylabel('Deformation [um]')

        plt.show()


def create_start_files():

    # write (1) or not (0) SHADOW files start.xx end.xx star.xx
    iwrite = 1

    #
    # initialize shadow3 source (oe0) and beam
    #
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    oe0.FDISTR = 3
    oe0.F_COLOR = 3
    oe0.F_PHOT = 0
    oe0.HDIV1 = 0.0
    oe0.HDIV2 = 0.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.NPOINT = 5000
    oe0.PH1 = 9990.0
    oe0.PH2 = 10010.0
    oe0.SIGDIX = 8.84999972e-05
    oe0.SIGDIZ = 7.1999998e-06
    oe0.SIGMAX = 5.70000011e-05
    oe0.SIGMAZ = 1.04000001e-05
    oe0.VDIV1 = 0.0
    oe0.VDIV2 = 0.0

    oe1.DUMMY = 1.0
    oe1.FHIT_C = 1
    oe1.FILE_REFL = b'si5_55.111'
    oe1.FILE_RIP = b'bump.dat'
    oe1.F_CENTRAL = 1
    oe1.F_CRYSTAL = 1
    oe1.F_G_S = 2
    oe1.F_RIPPLE = 1
    oe1.PHOT_CENT = 10000.0
    oe1.RLEN1 = 1.0
    oe1.RLEN2 = 1.0
    oe1.RWIDX1 = 1.0
    oe1.RWIDX2 = 1.0
    oe1.R_LAMBDA = 5000.0
    oe1.T_IMAGE = 1000.0
    oe1.T_INCIDENCE = 78.595143
    oe1.T_REFLECTION = 78.595143
    oe1.T_SOURCE = 3000.0


    oe0.write("start.00")
    oe1.write("start.01")
    print("Files written to disk: start.00 start.01")
    return oe0,oe1


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

    # oe1 = Shadow.OE()
    # oe1_before_run = Shadow.OE()


    # if platform.system() == "Linux":
    #     str = open('start.01', 'r').read()
    #     lines = str.split("\n")
    #     for line in lines:
    #         command1 = "oe1."+line.replace("(1)","[0]")\
    #         .replace("(2)","[1]")\
    #         .replace("(3)","[2]")\
    #         .replace("(4)","[3]")\
    #         .replace("(5)","[4]")\
    #         .replace("(6)","[5]")\
    #         .replace("(7)","[6]")\
    #         .replace("(8)","[7]")\
    #         .replace("(9)","[8]")\
    #         .replace("(10)","[9]")
    #         command2 = "oe1_before_run."+line.replace("(1)","[0]")\
    #         .replace("(2)","[1]")\
    #         .replace("(3)","[2]")\
    #         .replace("(4)","[3]")\
    #         .replace("(5)","[4]")\
    #         .replace("(6)","[5]")\
    #         .replace("(7)","[6]")\
    #         .replace("(8)","[7]")\
    #         .replace("(9)","[8]")\
    #         .replace("(10)","[9]")
    #
    #         try:
    #             exec(command1)
    #         except:
    #             print("run_shadow3_from_start_files: Failed to exec: %s"%command1)
    #
    #         try:
    #             exec(command2)
    #         except:
    #             print("run_shadow3_from_start_files: Failed to exec: %s"%command2)
    #
    #
    # else:
    #     oe1.load("start.01")
    #     oe1_before_run.load("start.01")


    #TODO this gives error in Mac
    # oe1.load("start.01")
    # oe1_before_run.load("start.01")
    oe0,oe1 = create_start_files()
    oe0,oe1_before_run = create_start_files()
    #
    beam.traceOE(oe1,1)
    #
    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    return beam_source,beam,oe0_before_run,oe1_before_run



def compare_results(do_plot=True,do_assert=True):

    if do_plot:
        Shadow.ShadowTools.plotxy("minimirr.01",2,1,nbins=101,nolost=1,title="Mirror (Python)",ref=0)
        Shadow.ShadowTools.plotxy("mirr.01",2,1,nbins=101,nolost=1,title="Mirror (SHADOW)",ref=0)

        Shadow.ShadowTools.plotxy("ministar.01",1,3,nbins=101,nolost=1,title="Image (Python)",ref=0)
        Shadow.ShadowTools.plotxy("star.01",1,3,nbins=101,nolost=1,title="Image (SHADOW)",ref=0)


    if do_assert:
        print("Comparing files mirr.01 and minimirr.01")
        minimirr = Shadow.Beam()
        minimirr.load("minimirr.01")
        mirr     = Shadow.Beam()
        mirr.load("mirr.01")
        assert_almost_equal(minimirr.rays[:,0:6],mirr.rays[:,0:6],2)

        print("Comparing files star.01 and ministar.01")
        ministar = Shadow.Beam()
        ministar.load("ministar.01")
        star     = Shadow.Beam()
        star.load("star.01")
        assert_almost_equal(ministar.rays[:,0:6],star.rays[:,0:6],2)



def minishadow_run_mesh_mirror():

    # ;
    # ; ray tracing of a surface defined with a mesh using minishadow
    # ; results are compared with shadow3
    # ;

    # ;
    # ; Runs shadow3
    # ;
    shadow3_beam_source,shadow3_beam,oe0,oe1 = run_shadow3_from_start_files(iwrite=1)


    # copy source to new Beam object
    newbeam = Beam.initialize_from_array(shadow3_beam_source.rays.copy())


    # ;
    # ; INPUTS
    # ;

    p             = oe1.T_SOURCE # 1000.0       # source-mirror
    q             = oe1.T_IMAGE  # 300.0        # mirror-image
    alpha         = oe1.ALPHA    # 0.0      # mirror orientation angle
    theta_grazing = (90.0-oe1.T_INCIDENCE) * numpy.pi / 180  # 5e-3     # grazing angle, rad

    print("p=%f, q=%f, alpha=%f, theta_grazing=%f rad"%(p,q,alpha,theta_grazing))

    mm = Mesh()
    mm.load_file("bump.dat")

    newbeam.rotate(alpha,axis=2)
    newbeam.rotate(theta_grazing,axis=1)
    newbeam.translation([0.0,-p*numpy.cos(theta_grazing),p*numpy.sin(theta_grazing)])


    # #
    # # reflect beam in the mirror surface and dump mirr.01
    # #

    newbeam,t,x1,v1,x2,v2 = mm.apply_specular_reflection_on_beam(newbeam)

    newbeam.dump_shadow3_file('minimirr.01')

    # #
    # # put beam in lab frame and compute image
    # #
    newbeam.rotate(theta_grazing,axis=1)
    # TODO what about alpha?
    newbeam.retrace(q,resetY=True)
    newbeam.dump_shadow3_file('ministar.01')

if __name__ == "__main__":

    create_start_files()
    create_gaussian_bump(do_plot=False)
    # write_bragg_preprocessor_file()
    minishadow_run_mesh_mirror()
    compare_results(do_plot=True,do_assert=True)

