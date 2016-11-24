FUNCTION shrot,a0,theta1,axis=axis,file=file,rad=rad

; alpha: alpha [deg]
; axis: 1:x, 2:y, 3:z

;+
; NAME:
;	SHROT
; PURPOSE:
;	This function rotates a Shadow file/structure and angle
;	a0 around a given axis
; CATEGORY:
;	SHADOW's utilities.
; CALLING SEQUENCE:
;	result = shRot(shadow_in,thetaIn)
; INPUTS:
;	shadow_in: The IDL structure with SHADOW data
;	thetaIn:  The rotation angle in degreed (default=0)
; OUTPUTS:
;	 Result: The IDL/Shadow structure of the rotated beam
; KEYWORD PARAMETERS (INPUT):
;	AXIS = The axis number (Shadow's column) for the rotation
;              (i.e, 1:x (default), 2:y, 3:z)
;	RAD = set this kewyword when  shadow_in is in radiants
; KEYWORD PARAMETERS (OUTPUT):
;	FILE = set this keyword to a string with the file name for
;       writing the rotated beam.
;
; MODIFICATION HISTORY:
;	M. Sanchez del Rio. ESRF. Grenoble May 2007
;
;-
on_error,2

IF N_ELEMENTS(theta1) EQ 0 THEN theta=0.0 ELSE theta=theta1 ; t_incidence[deg]
IF N_ELEMENTS(axis) EQ 0 THEN axis = 1

;
;
;
IF Not(Keyword_Set(rad)) THEN theta=theta*!dpi/180  ; rads

a0=readsh(a0)
a=a0
ray0=a0.ray
ray1=a.ray


torot = IndGen(3)+1
torot = torot[Where(torot NE axis)]
print,'torot: ',torot
costh = cos(theta)
sinth = sin(theta)

tstart = [1,4,7,16]

FOR i=0,N_Elements(tstart)-1 DO BEGIN
  tstarti=tstart-1
  cols = IndGen(3) + tstart[i] 
  colsi=cols-1
  newaxis= axis + tstart[i] -1
  newaxisi = newaxis-1
  newtorot = torot + tstart[i] -1 
  newtoroti = newtorot -1
;print,'>>> '
;print,'>>> colsi             : ',colsi
;print,'>>> rotatingi(fix,1,2): ',newaxisi,newtoroti

ray1[newtoroti[0],*] =  $
   ray0[newtoroti[0],*]*costh + ray0[newtoroti[1],*]*sinth
ray1[newtoroti[1],*] = $
  -ray0[newtoroti[0],*]*sinth + ray0[newtoroti[1],*]*costh
ray1[newaxisi,*]     =  ray0[newaxisi,*]

ENDFOR

a.ray=ray1

IF Keyword_Set(file) THEN putrays,a,file
RETURN,a

END



