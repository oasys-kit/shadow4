
import numpy

def shrot(a0,theta1,axis=1,file=None,rad=0):

# FUNCTION shrot,a0,theta1,axis=axis,file=file,rad=rad
#
# ; alpha: alpha [deg]
# ; axis: 1:x, 2:y, 3:z
#
# ;+
# ; NAME:
# ;	SHROT
# ; PURPOSE:
# ;	This function rotates a Shadow file/structure and angle
# ;	a0 around a given axis
# ; CATEGORY:
# ;	SHADOW's utilities.
# ; CALLING SEQUENCE:
# ;	result = shRot(shadow_in,thetaIn)
# ; INPUTS:
# ;	shadow_in: The IDL structure with SHADOW data
# ;	thetaIn:  The rotation angle in degreed (default=0)
# ; OUTPUTS:
# ;	 Result: The IDL/Shadow structure of the rotated beam
# ; KEYWORD PARAMETERS (INPUT):
# ;	AXIS = The axis number (Shadow's column) for the rotation
# ;              (i.e, 1:x (default), 2:y, 3:z)
# ;	RAD = set this kewyword when  shadow_in is in radiants
# ; KEYWORD PARAMETERS (OUTPUT):
# ;	FILE = set this keyword to a string with the file name for
# ;       writing the rotated beam.
# ;
# ; MODIFICATION HISTORY:
# ;	M. Sanchez del Rio. ESRF. Grenoble May 2007
# ;
# ;-
# on_error,2
#
# IF N_ELEMENTS(theta1) EQ 0 THEN theta=0.0 ELSE theta=theta1 ; t_incidence[deg]
# IF N_ELEMENTS(axis) EQ 0 THEN axis = 1
#
# ;
# ;
# ;
# IF Not(Keyword_Set(rad)) THEN theta=theta*!dpi/180  ; rads

    if not rad:
        theta1 = theta1 * numpy.pi / 180
#
# a0=readsh(a0)
# a=a0
# ray0=a0.ray
# ray1=a.ray

    a = a0.duplicate()
    ray0 = a0.rays
    ray1 = a.rays

#
#
# torot = IndGen(3)+1
# torot = torot[Where(torot NE axis)]
# print,'torot: ',torot
    if axis == 1:
        torot = [2,3]
    elif axis == 2:
        torot = [1,3]
    elif axis == 3:
        torot = [1,2]

# costh = cos(theta)
# sinth = sin(theta)
#
# tstart = [1,4,7,16]

    costh = numpy.cos(theta1)
    sinth = numpy.sin(theta1)

    tstart = numpy.array([1,4,7,16])

#
# FOR i=0,N_Elements(tstart)-1 DO BEGIN
#   tstarti=tstart-1
#   cols = IndGen(3) + tstart[i]
#   colsi=cols-1
#   newaxis= axis + tstart[i] -1
#   newaxisi = newaxis-1
#   newtorot = torot + tstart[i] -1
#   newtoroti = newtorot -1
# ;print,'>>> '
# ;print,'>>> colsi             : ',colsi
# ;print,'>>> rotatingi(fix,1,2): ',newaxisi,newtoroti
#
# ray1[newtoroti[0],*] =  $
#    ray0[newtoroti[0],*]*costh + ray0[newtoroti[1],*]*sinth
# ray1[newtoroti[1],*] = $
#   -ray0[newtoroti[0],*]*sinth + ray0[newtoroti[1],*]*costh
# ray1[newaxisi,*]     =  ray0[newaxisi,*]
#
# ENDFOR


    for i in range(len(tstart)):

        newaxis = axis + tstart[i] - 1
        newaxisi = newaxis - 1
        newtorot = torot + tstart[i] - 1
        newtoroti = newtorot -1

        print(ray0.shape,ray1.shape,"rotating axes %d, %d around axis %d"%(1+newtoroti[0],1+newtoroti[1],1+newaxisi))
        ray1[:,newtoroti[0]] =  ray0[:,newtoroti[0]] * costh + ray0[:,newtoroti[1]] * sinth
        ray1[:,newtoroti[1]] = -ray0[:,newtoroti[0]] * sinth + ray0[:,newtoroti[1]] * costh
        ray1[:,newaxisi]     =  ray0[:,newaxisi]


    # a.rays = ray1

    if file is not None:
        a.write(file)


    return a