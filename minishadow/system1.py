import numpy

from conicset import conicset
from shrot import shrot
from shtranslation import shtranslation
from conicintercept import conicintercept

import Shadow


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

def retrace(a0,dist,resetY=False):
    a1 = a0.copy()
    try:
        tof = (-a0[1,:] + dist)/a0[4,:]
        a1[0,:] += tof * a0[3,:]
        a1[1,:] += tof * a0[4,:]
        a1[2,:] += tof * a0[5,:]

        if resetY:
            a1[1,:] = 0.0

    except AttributeError:
        print ('retrace: No rays')

    return a1


def calculate_source_from_start_file(iwrite=0):
    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()

    oe0.load("start.00")

    # str = open('start.00', 'r').read()
    # lines = str.split("\n")
    # for line in lines:
    #     command = "oe0."+line
    #     try:
    #         exec(command)
    #     except:
    #         print("calculate_source_from_start_file: Failed to exec: %s"%command)

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

    oe1.load("start.01")
    oe1bis.load("start.01")
    # str = open('start.01', 'r').read()
    # lines = str.split("\n")
    # for line in lines:
    #     command = "oe1."+line.replace("(1)","[0]")\
    #     .replace("(2)","[1]")\
    #     .replace("(3)","[2]")\
    #     .replace("(4)","[3]")\
    #     .replace("(5)","[4]")\
    #     .replace("(6)","[5]")\
    #     .replace("(7)","[6]")\
    #     .replace("(8)","[7]")\
    #     .replace("(9)","[8]")\
    #     .replace("(10)","[9]")
    #
    #     try:
    #         exec(command)
    #     except:
    #         print("calculate_system_from_start_file: Failed to exec: %s"%command)

    beam.traceOE(oe1,1)

    if iwrite:
        oe1.write("end.01")
        beam.write("star.01")

    return oe1bis

def minishadow_run(iwrite=1):

    fileM='minimirr.01'
    fileI='ministar.01'

    # ;
    # ; 1) reads shadow source and puts it in the mirror reference system
    # ;
    # ;
    # ; get source from a SHADOW file and put it in the frame of mirror
    # ; get the starting points x1 and directions v1
    # ;
    beam,a = calculate_source_from_start_file(iwrite=iwrite)
    oe1 = calculate_system_from_start_file(beam,iwrite=iwrite)


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

    print(">>> T_INCIDENCE",oe1.T_INCIDENCE)
    ccc = conicset(p,q,theta,itype=fmirr,cylindrical=fcyl)


    #
    #
    #


    a1 = shrot(a,alpha,axis=2,rad=1)
    a2 = shrot(a1,theta,axis=1,rad=1)
    a3 = shtranslation(a2,[0.0,-p*numpy.cos(theta),p*numpy.sin(theta)])


    # ;
    # ; TRACING...
    # ;
    x1 =   a3[[0,1,2],:] # numpy.array(a3.getshcol([1,2,3]))
    v1 =   a3[[3,4,5],:] # numpy.array(a3.getshcol([4,5,6]))
    flag = a3[9,:]       # numpy.array(a3.getshonecol(10))

    x2 = numpy.zeros_like(x1)

    for kk in range(flag.size):
        t,iflag = conicintercept(ccc,x1[:,kk],v1[:,kk])
        x2[:,kk] = x1[:,kk] + v1[:,kk] * t
        if iflag < 0:
            flag[kk] = iflag

    # ;
    # ; Calculates the normal at each intercept [see shadow's normal.F]
    # ;
    normal = numpy.zeros_like(x1)

    normal[0,:] = 2 * ccc[1-1] * x2[0,:] + ccc[4-1] * x2[1,:] + ccc[6-1] * x2[2,:] + ccc[7-1]
    normal[1,:] = 2 * ccc[2-1] * x2[1,:] + ccc[4-1] * x2[0,:] + ccc[5-1] * x2[2,:] + ccc[8-1]
    normal[2,:] = 2 * ccc[3-1] * x2[2,:] + ccc[5-1] * x2[1,:] + ccc[6-1] * x2[0,:] + ccc[9-1]

    normalmod =  numpy.sqrt( normal[0,:]**2 + normal[1,:]**2 + normal[2,:]**2 )
    normal[0,:] /= normalmod
    normal[1,:] /= normalmod
    normal[2,:] /= normalmod


    # ;
    # ; reflection
    # ;
    tmp = v1 * normal
    tmp2 = tmp[0,:] + tmp[1,:] + tmp[2,:]
    tmp3 = normal.copy()

    for jj in (0,1,2):
        tmp3[jj,:] = tmp3[jj,:] * tmp2

    v2 = v1 - 2 * tmp3
    v2mod = numpy.sqrt(v2[0,:]**2 + v2[1,:]**2 + v2[2,:]**2)
    v2 /= v2mod

    # ;
    # ; writes the mirr.XX file
    # ;
    a4 = a3.copy()
    a4[[0,1,2],:] = x2
    a4[[3,4,5],:] = v2
    a4[9,:] = flag

    #
    # image
    #
    a5 = shrot(a4,theta,axis=1,rad=1)
    a6 = retrace(a5,q,resetY=True)



    dump_shadow_file(a4,fileM)
    dump_shadow_file(a6,fileI)


if __name__ == "__main__":
    minishadow_run()
    compare_results()
