
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


    if not rad:
        theta1 = theta1 * numpy.pi / 180

    a1 = a0.copy()

    if axis == 1:
        torot = [2,3]
    elif axis == 2:
        torot = [1,3]
    elif axis == 3:
        torot = [1,2]


    costh = numpy.cos(theta1)
    sinth = numpy.sin(theta1)

    tstart = numpy.array([1,4,7,16])

    for i in range(len(tstart)):

        newaxis = axis + tstart[i] - 1
        newaxisi = newaxis - 1
        newtorot = torot + tstart[i] - 1
        newtoroti = newtorot -1

        print(a0.shape,a1.shape,"rotating axes %d, %d around axis %d"%(1+newtoroti[0],1+newtoroti[1],1+newaxisi))
        a1[newtoroti[0],:] =  a0[newtoroti[0],:] * costh + a0[newtoroti[1],:] * sinth
        a1[newtoroti[1],:] = -a0[newtoroti[0],:] * sinth + a0[newtoroti[1],:] * costh
        a1[newaxisi]     =  a0[newaxisi,:]


    return a1