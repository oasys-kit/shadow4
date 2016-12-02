import numpy

from conicset import conicset
from shrot import shrot
from shtranslation import shtranslation
from conicintercept import conicintercept

import Shadow

def dump_shadow_file(a0,file):

    ncol,npoint = a0.shape
    if ncol != 18:
        raise Exception("dump_shadow_file: Not a beam (must have [18,:] ).")
    print("ncol = ",ncol,"npoint = ",npoint)

    beam = Shadow.Beam(N=npoint)

    beam.rays = a0.T

    beam.write("file")
    print("File %s written to disk. "%file)

    return beam

def retrace(a0,dist):
    a1 = a0.copy()
    try:
        tof = (-a0[1,:] + dist)/a0[4,:]
        a1[0,:] += tof * a0[3,:]
        a1[1,:] += tof * a0[4,:]
        a1[2,:] += tof * a0[5,:]
    except AttributeError:
        print ('retrace: No rays')

    return a1


def calculate_source_from_start_file(iwrite=0):
    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()


    str = open('start.00', 'r').read()
    lines = str.split("\n")
    for line in lines:
        command = "oe0."+line
        try:
            exec(command)
        except:
            print("calculate_source_from_start_file: Failed to exec: %s"%command)

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")


    return beam.rays.T.copy()


# ;
# ; ray tracing of a single conic mirror calculating the ray intersection
# ; analytically
# ;
#
# ;
# ; INPUTS
# ;
#
fmirr = 1
# ;
p = 1000.0       # source-mirror
q = 300.0        # mirror-image
alpha = 0.0      # mirror orientation angle
theta = 5e-3     # grazing angle, rad


ccc = conicset(p,q,theta,itype=fmirr,cylindrical=0)
print("ccc: ",ccc)
#
#
iwrite = 1 # for source
fileM='minimirr.01'
fileI='ministar.01'
#
#
# ;
# ; 1) reads shadow source and puts it in the mirror reference system
# ;
# ;
# ; get source from a SHADOW file and put it in the frame of mirror
# ; get the starting points x1 and directions v1
# ;
#
# a = readsh(fileS)

a = calculate_source_from_start_file(iwrite=iwrite)


print(">>>>>>>>>>>>>>>>>>>>",alpha)
a1 = shrot(a,alpha,axis=2,rad=1)


# a2= shrot(a1,theta,axis=1,/rad)
print(">>>>>>>>>>>>>>>>>>>>",theta)
a2= shrot(a1,theta,axis=1,rad=1)

# a3= shtranslation(a2,[0D,-p*cos(theta),p*sin(theta)])

a3= shtranslation(a2,[0.0,-p*numpy.cos(theta),p*numpy.sin(theta)])


# ;
# ; TRACING...
# ;


x1 =   a3[[0,1,2],:] # numpy.array(a3.getshcol([1,2,3]))
v1 =   a3[[3,4,5],:] # numpy.array(a3.getshcol([4,5,6]))
flag = a3[9,:]       # numpy.array(a3.getshonecol(10))

print('\n\n>>>>>>>>>x1',x1.shape)
print('\n\n>>>>>>>>>v1',v1.shape)


x2 = numpy.zeros_like(x1)
v2 = v1.copy()


print(">>>>>>>>>>x2, v2",x2.shape,v2.shape)

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



# normalmod = numpy.zeros_like(normal)
# normalmod[0,:] =  numpy.sqrt( normal[0,:]**2 + normal[1,:]**2 + normal[2,:]**2 )
# normalmod[1,:] =  normalmod[0,:]
# normalmod[2,:] =  normalmod[0,:]


normalmod =  numpy.sqrt( normal[0,:]**2 + normal[1,:]**2 + normal[2,:]**2 )
normal[0,:] /= normalmod
normal[1,:] /= normalmod
normal[2,:] /= normalmod
# print("<><><>normal",normal[0,0:10],normal[1,0:10],normal[2,0:10])
# print("<><><>===================",normal.shape,normalmod.shape)

# for i in range(10):
#     # print("<><> i: ",i,x1[:,i],v1[:,i],x2[:,i],v2[:,i])
#     print("<><> i: ",i,normal[0,i],normal[1,i],normal[2,i])

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
# a4=a3

a4 = a3.copy()

# print(">>>>>>>>>>>>>>>>",a4.shape,x2.shape,v2.shape,flag.shape)
#TODO MIRROR COL 5????


a4[[0,1,2],:] = x2
a4[[3,4,5],:] = v2
a4[9,:] = flag

#
# image
#


a5 = shrot(a4,theta,axis=1,rad=1)
# beam2 = dump_shadow_file(a5,"shrot.01")
#
# beam2.retrace(q)
# beam2.rays[:,1] = 0.0
# beam2.write("ministar2.01")


a6 = retrace(a5,q)
#reset Y
a6[1,:] = 0.0


# for i in range(10):
#     print("<><> i: ",i,v2[2,i],tmp3[:,i])

print("\n\n")

for i in range(10):
    print("<><> i: ",i,a6[0:5,i])


dump_shadow_file(a4,fileM)
dump_shadow_file(a5,"shrot.01")
dump_shadow_file(a6,fileI)
#
# files
#
# if fileM != None: a4.write(fileM)
# if fileI != None: a6.write(fileI)
#
# print("<><><>",a.rays.shape)

