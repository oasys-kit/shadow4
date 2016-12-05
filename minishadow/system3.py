import numpy

# from conicset import conicset
# from conicintercept import conicintercept
import platform


import Shadow
from Beam import Beam
from SurfaceConic import SurfaceConic

def compare_results():
    Shadow.ShadowTools.plotxy("minimirr.01",2,1,nbins=101,nolost=1,title="Mirror (Python)")
    Shadow.ShadowTools.plotxy("mirr.01",2,1,nbins=101,nolost=1,title="Mirror (SHADOW)")

    Shadow.ShadowTools.plotxy("ministar.01",1,3,nbins=101,nolost=1,title="Image (Python)")
    Shadow.ShadowTools.plotxy("star.01",1,3,nbins=101,nolost=1,title="Image (SHADOW)")

def dump_shadow_file(a0,file):

    ncol,npoint = a0.shape
    if ncol != 18:
        raise Exception("dump_shadow_file: Not a beam (must have [18,:] ).")
    # print("ncol = ",ncol,"; npoint = ",npoint)

    beam = Shadow.Beam(N=npoint)

    beam.rays = a0.T.copy()

    beam.write(file)
    print("File %s written to disk. "%file)

    return beam



def calculate_source_from_start_file(iwrite=0):
    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()


    # TODO: this is a turn-around for the Linux bug...
    if platform.system() == "Linux":
        str = open('start.00', 'r').read()
        lines = str.split("\n")
        for line in lines:
            command = "oe0."+line
            try:
                exec(command)
            except:
                print("calculate_source_from_start_file: Failed to exec: %s"%command)
    else:
        oe0.load("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")


    return beam,beam.rays.T.copy()


def calculate_system_from_start_file(beamIn,iwrite=0):
    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = beamIn.duplicate()
    oe1 = Shadow.OE()
    oe1bis = Shadow.OE()


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
            command2 = "oe1bis."+line.replace("(1)","[0]")\
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
                print("calculate_system_from_start_file: Failed to exec: %s"%command1)

            try:
                exec(command2)
            except:
                print("calculate_system_from_start_file: Failed to exec: %s"%command2)


    else:
        oe1.load("start.01")
        oe1bis.load("start.01")

    beam.traceOE(oe1,1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    return oe1bis



# def vector_reflection(v1,normal):
#     tmp = v1 * normal
#     tmp2 = tmp[0,:] + tmp[1,:] + tmp[2,:]
#     tmp3 = normal.copy()
#
#     for jj in (0,1,2):
#         tmp3[jj,:] = tmp3[jj,:] * tmp2
#
#     v2 = v1 - 2 * tmp3
#     v2mod = numpy.sqrt(v2[0,:]**2 + v2[1,:]**2 + v2[2,:]**2)
#     v2 /= v2mod
#
#     return v2
#
# def ccc_normal(ccc,x2):
#     # ;
#     # ; Calculates the normal at each intercept [see shadow's normal.F]
#     # ;
#     normal = numpy.zeros_like(x2)
#
#     normal[0,:] = 2 * ccc[1-1] * x2[0,:] + ccc[4-1] * x2[1,:] + ccc[6-1] * x2[2,:] + ccc[7-1]
#     normal[1,:] = 2 * ccc[2-1] * x2[1,:] + ccc[4-1] * x2[0,:] + ccc[5-1] * x2[2,:] + ccc[8-1]
#     normal[2,:] = 2 * ccc[3-1] * x2[2,:] + ccc[5-1] * x2[1,:] + ccc[6-1] * x2[0,:] + ccc[9-1]
#
#     normalmod =  numpy.sqrt( normal[0,:]**2 + normal[1,:]**2 + normal[2,:]**2 )
#     normal[0,:] /= normalmod
#     normal[1,:] /= normalmod
#     normal[2,:] /= normalmod
#
#     return normal
#
#
# def ccc_reflection_beam(ccc,newbeam):
#     # ;
#     # ; TRACING...
#     # ;
#
#     x1 =   newbeam.get_columns([1,2,3]) # numpy.array(a3.getshcol([1,2,3]))
#     v1 =   newbeam.get_columns([4,5,6]) # numpy.array(a3.getshcol([4,5,6]))
#     flag = newbeam.get_column(10)        # numpy.array(a3.getshonecol(10))
#
#
#     t,iflag = conicintercept(ccc,x1,v1)
#     x2 = x1 + v1 * t
#     for i in range(flag.size):
#         if iflag[i] < 0: flag[i] = -100
#
#
#     # ;
#     # ; Calculates the normal at each intercept [see shadow's normal.F]
#     # ;
#
#     normal = ccc_normal(ccc,x2)
#
#     # ;
#     # ; reflection
#     # ;
#
#     v2 = vector_reflection(v1,normal)
#
#     # ;
#     # ; writes the mirr.XX file
#     # ;
#
#     newbeam.set_column(1, x2[0])
#     newbeam.set_column(2, x2[1])
#     newbeam.set_column(3, x2[2])
#     newbeam.set_column(4, v2[0])
#     newbeam.set_column(5, v2[1])
#     newbeam.set_column(6, v2[2])
#     newbeam.set_column(10, flag )
#
#     return newbeam


def minishadow_run(iwrite=1):



    # ;
    # ; 1) reads shadow source 
    #
    shadow_beam,a = calculate_source_from_start_file(iwrite=iwrite)
    newbeam = Beam.initialize_from_array(a)


    #
    # runs the system with shadow3
    oe1 = calculate_system_from_start_file(shadow_beam,iwrite=iwrite)




    # ;
    # ; ray tracing of a single conic mirror calculating the ray intersection
    # ; analytically
    # ;
    #


    # ;
    # ; INPUTS
    # ;
    #
    fmirr   = oe1.FMIRR  # 1
    p       = oe1.T_SOURCE  # 1000.0       # source-mirror
    q       = oe1.T_IMAGE  # 300.0        # mirror-image
    alpha   = oe1.ALPHA  # 0.0      # mirror orientation angle
    theta   = (90.0-oe1.T_INCIDENCE) * numpy.pi / 180  # 5e-3     # grazing angle, rad
    fcyl    = oe1.FCYL

    print("fmirr = %s, p=%f, q=%f, alpha=%f, theta=%f rad, fcyl=%d"%\
              (fmirr,p,q,alpha,theta,fcyl))

    ccc = conicset(p,q,theta,itype=fmirr,cylindrical=fcyl)


    #
    # put beam in mirror reference system
    #
    newbeam.rotate(alpha,axis=2)
    newbeam.rotate(theta,axis=1)
    newbeam.translation([0.0,-p*numpy.cos(theta),p*numpy.sin(theta)])


    #
    # reflect beam in the mirror surface and dump mirr.01
    #
    newbeam = ccc_reflection_beam(ccc,newbeam)
    dump_shadow_file(newbeam.get_rays(),'minimirr.01')

    #
    # put beam in lab frame and compute image
    #
    newbeam.rotate(theta,axis=1)
    # TODO what about alpha?
    newbeam.retrace(q,resetY=True)
    dump_shadow_file(newbeam.get_rays(),'ministar.01')

if __name__ == "__main__":
    minishadow_run()
    compare_results()
