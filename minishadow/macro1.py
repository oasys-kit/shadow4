import numpy

from conicset import conicset
from shrot import shrot
from shtranslation import shtranslation
from conicintercept import conicintercept

def calculate_source(iwrite=0):
    import Shadow

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()

    #
    # Define variables. See meaning of variables in:
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #
    oe0.FSOUR = 0
    oe0.F_PHOT = 0
    oe0.HDIV1 = 5e-06
    oe0.HDIV2 = 5e-06
    oe0.NPOINT = 100
    oe0.PH1 = 1000.0
    oe0.ISTAR1 =  6775431


    #Run SHADOW to create the source

    if iwrite:
        oe0.write("start.00")

    beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
        beam.write("begin.dat")

    return beam

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
#
# ; radius to focus
# ; conic coeficients
# ;ccc = DblArr(10)
#
# ;radius=2D0*p*q/(p+q)/sin(theta)
# ccc=conicset(p,q,theta,type=fmirr,cylindrical=0)
# ;print,'ccc: ',ccc
# ;print,'p,q,theta,radius: ',p,q,theta,radius

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

a = calculate_source(iwrite=iwrite)


# a1= shrot(a,alpha,axis=2,/rad)
print(">>>>>>>>>>>>>>>>>>>>",alpha)
a1 = shrot(a,alpha,axis=2,rad=1)


# a2= shrot(a1,theta,axis=1,/rad)
print(">>>>>>>>>>>>>>>>>>>>",theta)
a2= shrot(a1,theta,axis=1,rad=1)

# a3= shtranslation(a2,[0D,-p*cos(theta),p*sin(theta)])

a3= shtranslation(a2,[0.0,-p*numpy.cos(theta),p*numpy.sin(theta)])

# ;
#
# ;
# ; TRACING...
# ;
# x1= GetShCol(a3,[1,2,3])
# v1= GetShCol(a3,[4,5,6])
# flag= GetShCol(a3,10)

x1 =   numpy.array(a3.getshcol([1,2,3]))
v1 =   numpy.array(a3.getshcol([4,5,6]))
flag = numpy.array(a3.getshonecol(10))

print('\n\n>>>>>>>>>x1',x1.shape)
print('\n\n>>>>>>>>>v1',v1.shape)
# for i in range(100):
#     print("i,xi: ",i,x1[0,i],x1[1,i],x1[2,i])
#
# for i in range(100):
#     print("i,vi: ",i,v1[0,i],v1[1,i],v1[2,i])

# ;
#
# ; calculates x2,v2,flag
# x2=x1
# v2=v1
# FOR kk=0L,N_Elements(flag)-1 DO BEGIN
#   t = conicintercept(ccc,x1[*,kk],v1[*,kk],iflag)
#   ;print,'kk,t,iflag: ',kk,t,iflag
#   x2[*,kk] = x1[*,kk]+v1[*,kk]*t
#   flag[kk]=iflag
# ENDFOR
#


x2 = numpy.zeros_like(x1)
v2 = v1.copy()


print(">>>>>>>>>>x2, v2",x2.shape,v2.shape)

for kk in range(flag.size):
    t,iflag = conicintercept(ccc,x1[:,kk],v1[:,kk])
    # print ('kk,t,iflag: ',kk,t,flag)
    x2[:,kk] = x1[:,kk] + v1[:,kk] * t
    if iflag < 0:
        flag[kk] = iflag




# ;
# ; Calculates the normal at each intercept [see shadow's normal.F]
# ;
# normal = x1*0D0
#
# normal[0,*]=2*CCC(1-1)*X2[0,*] + CCC(4-1)*x2[1,*] + CCC(6-1)*x2[2,*] + CCC(7-1)
# normal[1,*]=2*CCC(2-1)*X2[1,*] + CCC(4-1)*x2[0,*] + CCC(5-1)*x2[2,*] + CCC(8-1)
# normal[2,*]=2*CCC(3-1)*X2[2,*] + CCC(5-1)*x2[1,*] + CCC(6-1)*x2[0,*] + CCC(9-1)
# normalmod = normal*0
# normalmod[0,*] =  sqrt( normal[0,*]^2+normal[1,*]^2+normal[2,*]^2 )
# normalmod[1,*] =  normalmod[0,*]
# normalmod[2,*] =  normalmod[0,*]
# ;help,normal,normalmod
# normal=normal/normalmod
# ;help,normal
#


normal = numpy.zeros_like(x1)

normal[0,:] = 2 * ccc[1-1] * x2[0,:] + ccc[4-1] * x2[1,:] + ccc[6-1] * x2[2,:] + ccc[7-1]
normal[1,:] = 2 * ccc[2-1] * x2[1,:] + ccc[4-1] * x2[0,:] + ccc[5-1] * x2[2,:] + ccc[8-1]
normal[2,:] = 2 * ccc[3-1] * x2[2,:] + ccc[5-1] * x2[1,:] + ccc[6-1] * x2[0,:] + ccc[9-1]

normalmod = numpy.zeros_like(x1)
normalmod[0,:] =  numpy.sqrt( normal[0,:]**2 + normal[1,:]**2 + normal[2,:]**2 )
normalmod[1,:] =  normalmod[0,:]
normalmod[2,:] =  normalmod[0,:]

normal = normal / normalmod


# ;print,normal[*,0:10]
# ;
# ; reflection
# ;
# tmp = v1*normal
# tmp2 = Reform(tmp[0,*]+tmp[1,*]+tmp[2,*])
# tmp3 = normal
# ;help,tmp2,tmp3
# FOR jj=0,2 DO tmp3[jj,*]=tmp3[jj,*]*tmp2
# v2=v1-2D*tmp3
# v2mod=sqrt(v2[0,*]^2+v2[1,*]^2+v2[2,*]^2)
# FOR jj=0,2 DO v2[jj,*]=v2[jj,*]/v2mod
#


tmp = v1 * normal
tmp2 = tmp[0,:] + tmp[1,:] + tmp[2,:]
tmp3 = normal

# ;help,tmp2,tmp3
#
for jj in (0,2):
    tmp3[jj,:] = tmp3[jj,:] * tmp2

v2 = v1 - 2 * tmp3
v2mod = numpy.sqrt(v2[0,:]**2 + v2[1,:]**2 + v2[2,:]**2)
for jj in (0,2):
    v2[jj,:] = v2[jj,:] / v2mod


for i in range(10):
    # print("<><> i: ",i,x1[:,i],v1[:,i],x1[:,i],v1[:,i])
    print("<><> i: ",i,x1[:,i],v1[:,i],x2[:,i],v2[:,i])

# ;
# ; writes the mirr.XX file
# ;
# a4=a3

a4 = a3.duplicate()

print(">>>>>>>>>>>>>>>>",a4.rays.shape,x2.shape,v2.shape,flag.shape)
a4.rays[:,[0,1,2]] = x2.T
a4.rays[:,[3,4,5]] = v2.T
a4.rays[:,9] = flag

#
# image
#

# a5 = shrot(a4,theta,axis=1,/rad)
# putrays,a5,"shrot.01"
# a6 = retrace(a5,dist=q,/resetY)
# putrays,a6,fileI

a5 = shrot(a4,theta,axis=1,rad=1)
a5.write("shrot.01")

a6 = a5.duplicate()
a6.retrace(q)
#reset Y
a6.rays[:,1] = 0.0


#
# files
#
if fileM != None: a4.write(fileM)
if fileI != None: a6.write(fileI)


# ray=a4.ray
# ray([0,1,2],*)=x2
# ray([3,4,5],*)=v2
# ray(9,*)=flag
# a4.ray=ray
# putrays,a4,fileM
# ;
# ; writes the star.XX file
# ;

#
# END

