import numpy

# IF N_Elements(iType) EQ 0 THEN iType=5
# IF N_Elements(p) EQ 0 THEN p=3000D0
# IF N_Elements(q) EQ 0 THEN q=1000D0
# IF N_Elements(theta1) EQ 0 THEN theta=5D-3 ELSE theta=theta1
# IF N_Elements(cyl) EQ 0 THEN cyl=0
# IF N_Elements(cylAngle) EQ 0 THEN cylAngle=0.0
# IF N_Elements(convex) EQ 0 THEN convex=0
# IF N_Elements(anglenordeg) EQ 0 THEN anglenordeg=0


def conicset(p, q, theta1, itype=1, cylindrical=0, cylangle=0.0, convex=0, anglenordeg=0):



# FUNCTION conicset,p,q,theta1,type=itype,  $
#   cylindrical=cylindrical, cylangle=cylangle, $ ; in degrees
#   convex=convex,txt=txt,anglenordeg=anglenordeg
#
# ;+
# ;
# ;       NAME:
# ;               CONICSET
# ;       PURPOSE:
# ;               This function Calculates the 10 conic
# ;               coefficients from the focusing conditions
# ;       CATEGORY:
# ;               SHADOW tools
# ;       CALLING SEQUENCE:
# ;               c = conicset(p,q,theta)
# ;
# ; 	INPUTS:
# ;		p: distance source-surface [cm]
# ;		q: distance surface-image [cm]
# ;		theta: grazing andgle in rad (for angle with respect to
# ;                      the normal, set the ANGLENORDEG keyword)
# ;
# ; 	KEYWORD PARAMETERS
# ;	        ANGLENORDEG: Set this keyword to indicate that the
# ;		  angle is in Degrees with respect to the surface normal
# ;		TYPE:
# ;        	'0'
# ;        	'Spherical-1'
# ;        	'Ellipsoidal-2'
# ;        	'Toroidal-3'
# ;        	'Paraboloid-4'
# ;        	'Plane-5'
# ;        	'6'
# ;        	'Hyperboloid-7'
# ;        	'Cone-8??'
# ;        	'Polynomial-9??
# ;		CYLINDRICAL: set to 1 for cylindrical symmetry
# ;		CYLANGLE: if cylindrical is set, the angle between
# ;		  the cylinder axis and the X axis in deg.
# ;		CONVEX: set to 1 for invering convexity
# ;		TXT: set to a named variable to receive some
# ;		     info on the calculations.
# ;
# ;	ALGORITHM:
# ;		 Translated from SHADOW/MSETUP
# ;
# ;		 Equation of the conic:
# ;
# ;	         c[0]*X^2 + c[1]*Y^2 + c[2]*Z^2 +
# ;                c[3]*X*Y + c[4]*Y*Z + c[5]*X*Z  +
# ;                c[6]*X + c[7]*Y + c[8]*Z + c[9] = 0
#
# ;
# ;	AUTHOR:
# ;		M. Sanchez del Rio srio@esrf.eu, 2006
# ;
# ;	MODIFICATION HISTORY:
# ;		2007-11-08 added doc.
# ;
# ;-
# ;C+++
# ;C
# ;C	SUBROUTINE	MSETUP
# ;C
# ;C	PURPOSE		To compute the parameters specifying a given
# ;C			mirror.
# ;C
# ;C	OUTPUT		a) to common block
# ;C			b) MSETxx, xx being the OE number
# ;C
# ;C---
#
# IF N_Elements(iType) EQ 0 THEN iType=5
# IF N_Elements(p) EQ 0 THEN p=3000D0
# IF N_Elements(q) EQ 0 THEN q=1000D0
# IF N_Elements(theta1) EQ 0 THEN theta=5D-3 ELSE theta=theta1
# IF N_Elements(cyl) EQ 0 THEN cyl=0
# IF N_Elements(cylAngle) EQ 0 THEN cylAngle=0.0
# IF N_Elements(convex) EQ 0 THEN convex=0
# IF N_Elements(anglenordeg) EQ 0 THEN anglenordeg=0
#
#
# SSOUR=Double(p)
# SIMAG=Double(q)
# THETA=Double(THETA)

    ssour = float(p)
    simag = float(q)
    theta = float(theta1)

#
#
    oeType=['0',
            'Spherical-1',
            'Ellipsoidal-2',
            'Toroidal-3',
            'Paraboloid-4',
            'Plane-5',
            '6',
            'Hyperboloid-7',
            'Cone-8??',
            'Polynomial-9??']
#
#
    if anglenordeg == 0:
        angleTxt='Grazing angle [rad]: '
    else:
        angleTxt='Angle to normal [deg]: '


    txt = ""
#
# 	  txt = ['************ conicset *************',$
#              'Computation of conic equation of the form: ',$
#              '    c[0]*X^2 + c[1]*Y^2 + c[2]*Z^2 + ',$
#              '    c[3]*X*Y + c[4]*Y*Z + c[5]*X*Z  + ',$
#              '    c[6]*X + c[7]*Y + c[8]*Z + c[9] = 0  ',$
#              '','Inputs: ',$
#                 String('p[cm]: ',p,Format='(A40,G20.15)'),$
#                 String('q[cm]: ',q,Format='(A40,G20.15)'),$
#                 String(angleNorDeg,theta, $
#                        Format='(A40,G20.15)'),$
#                 String('Conic type: ', $
#                        oeType[iType], $
#                        Format='(2A40)'),$
#                 String('Invert convexity flag: ', $
#                        convex, $
#                        Format='(A40,I5)'),$
#                 String('Cylindrical symmetry flag: ', $
#                        cyl, $
#                        Format='(A40,I5)'),$
#                 String('Cylindrial axes to X [deg]: ', $
#                        cylAngle, $
#                        Format='(A40,G20.15)'), $
#              '','','','Outputs: ']
#
# ;COSTHE=sin(theta) ; For shadow, theta is measured with respect to the normal
# ;SINTHE=cos(theta)



# ;
# ; theta in rad with respect to the normal
# ;
#
# IF angleNorDeg EQ 0 THEN BEGIN
#    theta = (!dpi/2)-theta ; *1D-3
# ENDIF ELSE BEGIN
#    theta = theta*!dpi/180.0
# ENDELSE
# ; print,'Angle with respect to the surface normal [rad]:',theta
#


    if anglenordeg == 0:
       theta = (numpy.pi/2)-theta
    else:
       theta = theta * numpy.pi / 180.0

    print('Angle with respect to the surface normal [rad]:',theta)

    COSTHE = numpy.cos(theta)
    SINTHE = numpy.sin(theta)


#
#
# ccc = DblArr(10)
    ccc = numpy.zeros(10)
#
    if itype == 0:
        pass
    elif itype == 1:  # Spherical
        # 		;rmirr=2D0*p*q/(p+q)/sin(theta)
        rmirr = ssour * simag * 2 / COSTHE / (ssour + simag)

        ccc[1-1] =  1.0	        # X^2  # = 0 in cylinder case
        ccc[2-1] =  1.0	        # Y^2
        ccc[3-1] =  1.0	        # Z^2
        ccc[4-1] =   .0	        # X*Y   # = 0 in cylinder case
        ccc[5-1] =   .0	        # Y*Z
        ccc[6-1] =   .0	        # X*Z   # = 0 in cylinder case
        ccc[7-1] =   .0	        # X     # = 0 in cylinder case
        ccc[8-1] =   .0	        # Y
        ccc[9-1] = -2 * rmirr	# Z
        ccc[10-1]  =   .0       # G
        txt += "%s Spherical radius: %f \n"%(" "*40,rmirr)
    else:
        raise Exception("itype Not yet implemented")

# 	   END
# 	2: BEGIN ; Ellipsoidal
# 		AXMAJ 	=  ( SSOUR + SIMAG )/2
# 		AXMIN 	=  SQRT(SIMAG*SSOUR)*COSTHE
# 		AFOCI 	=  SQRT( AXMAJ^2-AXMIN^2 )
# 		ECCENT 	=  AFOCI/AXMAJ
# 		;C
# 		;C The center is computed on the basis of the object and image positions
# 		;C
# 		YCEN  = (SSOUR-SIMAG)*0.5D0/ECCENT
# 		ZCEN  = -SQRT(1-YCEN^2/AXMAJ^2)*AXMIN
# 		;C
# 		;C Computes now the normal in the mirror center.
# 		;C
# 		rncen=DblArr(3)
# 		RNCEN(1-1)  =   .0D0
# 		RNCEN(2-1)  = - 2*YCEN/AXMAJ^2
# 		RNCEN(3-1)  = - 2*ZCEN/AXMIN^2
# 		;CALL NORM(RNCEN,RNCEN)
# 		rncen = rncen/Sqrt(Total(rncen^2))
# 		;C
# 		;C Computes the tangent versor in the mirror center.
# 		;C
# 		rtcen=DblArr(3)
# 		RTCEN(1-1)  =  .0D0
# 		RTCEN(2-1)  =   RNCEN(3-1)
# 		RTCEN(3-1)  = - RNCEN(2-1)
# 		;C
# 		txt = [txt,  $
#                        String('Rev Ellipsoid a: ', $
#                        AXMAJ, Format='(A40,G20.15)'), $
#                        String('Rev Ellipsoid b: ', $
#                        AXMIN, Format='(A40,G20.15)'), $
#                        String('Rev Ellipsoid c=sqrt(a^2-b^2): ', $
#                        AFOCI, Format='(A40,G20.15)'), $
#                        String('Rev Ellipsoid focal discance c^2: ', $
#                        AFOCI^2, Format='(A40,G20.15)'), $
#                        String('Rev Ellipsoid excentricity: ', $
#                        ECCENT, Format='(A40,G20.15)'),$
#                        'Mirror center at: '+vect2string([0,YCEN,ZCEN]), $
#                        'Mirror normal: '+vect2string(RNCEN), $
#                        'Mirror tangent: '+vect2string(RTCEN) ]
# 		;C Computes now the quadric coefficient with the mirror center
# 		;C located at (0,0,0) and normal along (0,0,1)
# 		;C
# 		A 	=  1/AXMIN^2
# 		B 	=  1/AXMAJ^2
# 		C 	=  A
# 		CCC(1-1) 	=  A
# 		CCC(2-1) 	=  B*RTCEN(2-1)^2 + C*RTCEN(3-1)^2
# 		CCC(3-1) 	=  B*RNCEN(2-1)^2 + C*RNCEN(3-1)^2
# 		CCC(4-1) 	=  .0D0
# 		CCC(5-1) 	=  2*(B*RNCEN(2-1)*RTCEN(2-1)+C*RNCEN(3-1)*RTCEN(3-1))
# 		CCC(6-1) 	=  .0D0
# 		CCC(7-1) 	=  .0D0
# 		CCC(8-1) 	=  .0D0
# 		CCC(9-1) 	=  2*(B*YCEN*RNCEN(2-1)+C*ZCEN*RNCEN(3-1))
# 		CCC(10-1) =  .0D0
# 	   END
# 	3: BEGIN ; Toroidal ?????????????????
# 		R_MAJ	=   SSOUR*SIMAG*2/COSTHE/(SSOUR + SIMAG)
# 		R_MIN	=   SSOUR*SIMAG*2*COSTHE/(SSOUR + SIMAG)
# 		txt = [txt, 'Toroid is not a conic!', $
#                        String('Major radius: ', $
#                        R_MAJ, Format='(A40,G20.15)'), $
#                        String('Minor radius: ', $
#                        R_MIN, Format='(A40,G20.15)')]
# 		;C
# 		;C NOTE : The major radius is the in reality the radius of the torus
# 		;C max. circle. The true major radius is then
# 		;C
# 		R_MAJ	=   R_MAJ - R_MIN
# 	   END
# 	4: BEGIN ; Paraboloid
# 		;C
# 		;C Computes the parabola
# 		;C
# 		IF SSOUR LT SIMAG THEN BEGIN
# 		  PARAM	=   2*SSOUR*COSTHE^2
#      		  YCEN	= - SSOUR*SINTHE^2
#      		  ZCEN	= - 2*SSOUR*SINTHE*COSTHE
# 		  fact = -1.0D0
# 		ENDIF ELSE BEGIN
# 		  PARAM	=   2*SIMAG*COSTHE^2
#      		  YCEN	= - SIMAG*SINTHE^2
#      		  ZCEN	= - 2*SIMAG*SINTHE*COSTHE
# 		  fact = 1.0D0
# 		ENDELSE
# 		txt = [txt, $
#                        String('Parabolois p: ', $
#                        PARAM, Format='(A40,G20.15)')]
# 		CCC(1-1)	=   1.0D0
# 		CCC(2-1)	=   COSTHE^2
# 		CCC(3-1)	=   SINTHE^2
# 		CCC(4-1)  =    .0D0
# 		CCC(5-1)	=   2*fact*COSTHE*SINTHE
# 		CCC(6-1)	=    .0D0
# 		CCC(7-1)	=    .0D0
# 		CCC(8-1)	=    .0D0
# 		CCC(9-1)	=   2*ZCEN*SINTHE - 2*PARAM*COSTHE
# 		CCC(10-1) =    .0D0
# 	   END
# 	5: BEGIN  ; Plane
# 		;C
# 		;C The sign of CCC(9) is < 0 to keep consistency with the other cases
# 		;C normal.
# 		;C
# 		CCC(9-1) = - 1.0D0
# 		txt = [txt, 'Plane is a degenerated conic']
# 	   END
#         6:
# 	7: BEGIN  ; Hyperbolical
# 		AXMAJ 	=  ( SSOUR - SIMAG )/2
# 		;C
# 		;C If AXMAJ > 0, then we are on the left branch of the hyp. Else we
# 		;C are onto the right one. We have to discriminate between the two cases
# 		;C In particular, if AXMAJ.LT.0 then the hiperb. will be convex.
# 		;C
# 		AFOCI = 0.5D0*SQRT( SSOUR^2 + SIMAG^2 + 2*SSOUR*SIMAG*COS(2*THETA) )
# ;; why this works better?
# ;;		AFOCI = 0.5D0*SQRT( SSOUR^2 + SIMAG^2 - 2*SSOUR*SIMAG*COS(2*THETA) )
# 		AXMIN  =  SQRT( AFOCI^2 - AXMAJ^2 )
#
# 		ECCENT 	=  AFOCI/ABS( AXMAJ )
#
#                 BRANCH=-1.0D0 ; branch=+1,-1
# 		;C
# 		;C Computes the center coordinates in the hiperbola RF.
# 		;C
# 		;IF AXMAJ GT 0.0D0 THEN BEGIN
# 		;  YCEN	=   ( AXMAJ - SSOUR )/ECCENT			; < 0
# 		;ENDIF ELSE BEGIN
# 		;  YCEN	=   ( SSOUR - AXMAJ )/ECCENT			; > 0
# 		;ENDELSE
#
#                 YCEN =   ABS( SSOUR - AXMAJ )/ECCENT*BRANCH
#
# 		ZCEN_ARG = ABS( YCEN^2/AXMAJ^2 - 1.0D0)
# 		IF ZCEN_ARG GT 1.0D-14  THEN BEGIN
# 		  ZCEN	= - AXMIN * SQRT(ZCEN_ARG)			; < 0
# 		ENDIF ELSE BEGIN
# 		  ZCEN  = 0.0D0
# 		ENDELSE
# ;
# ; THIS GIVES BETTER LOOKING HYPERBOLA BUT WORSE TRACING. WHY?
# ;YCEN=ABS(YCEN)
# ;ZCEN=ABS(ZCEN)
# 		;C
# 		;C Computes now the normal in the same RF. The signs are forced to
# 		;C suit our RF.
# 		;C
# 		rncen = DblArr(3)
# 		RNCEN (1-1) =   0.0D0
# 		RNCEN (2-1) = - ABS( YCEN )/AXMAJ^2			; < 0
# 		RNCEN (3-1) = - ZCEN/AXMIN^2				; > 0
# ;Csrio		RNCEN (2-1) =   YCEN/AXMAJ^2			; < 0
# ;Csrio		RNCEN (3-1) = - ZCEN/AXMIN^2				; > 0
#
# 		;CALL 	NORM	(RNCEN,RNCEN)
# 		rncen = rncen/Sqrt(Total(rncen^2))
# 		;C
# 		;C Computes the tangent in the same RF
# 		;C
# 		rtcen = DblArr(3)
# 		RTCEN (1-1) =   0.0D0
# 		RTCEN (2-1) = - RNCEN(3-1)   ; > 0
# 		RTCEN (3-1) =   RNCEN(2-1)   ; > 0
# ;;		RTCEN (2-1) = ABS(RNCEN(3-1))
# ;;		RTCEN (3-1) = ABS(RNCEN(2-1))
#
#
# 		txt = [txt,  $
#                        String('Rev Hyperboloid a: ', $
#                        AXMAJ, Format='(A40,G20.15)'), $
#                        String('Rev Hyperboloid b: ', $
#                        AXMIN, Format='(A40,G20.15)'), $
#                        String('Rev Hyperboloid c: ', $
#                        AFOCI, Format='(A40,G20.15)'), $
#                        String('Rev Hyperboloid focal discance c^2: ', $
#                        AFOCI^2, Format='(A40,G20.15)'), $
#                        String('Rev Hyperboloid excentricity: ', $
#                        ECCENT, Format='(A40,G20.15)'), $
#                        'Mirror BRANCH: '+String(branch), $
#                        'Mirror center at: '+vect2string([0,YCEN,ZCEN]), $
#                        'Mirror normal: '+vect2string(RNCEN), $
#                        'Mirror tangent: '+vect2string(RTCEN) ]
# 		;C
# 		;C Coefficients of the canonical form
# 		;C
# 		A	= - 1/AXMIN^2
# 		B	=   1/AXMAJ^2
# 		C	=   A
# 		;C
# 		;C Rotate now in the mirror RF. The equations are the same as for the
# 		;C ellipse case.
# 		;C
# 		CCC(1-1) 	=  A
# 		CCC(2-1) 	=  B*RTCEN(2-1)^2 + C*RTCEN(3-1)^2
# 		CCC(3-1) 	=  B*RNCEN(2-1)^2 + C*RNCEN(3-1)^2
# 		CCC(4-1) 	=  .0D0
# 		CCC(5-1) 	= 2*(B*RNCEN(2-1)*RTCEN(2-1)+C*RNCEN(3-1)*RTCEN(3-1))
# 		CCC(6-1) 	=  .0D0
# 		CCC(7-1) 	=  .0D0
# 		CCC(8-1) 	=  .0D0
# 		CCC(9-1) 	= 2*(B*YCEN*RNCEN(2-1)+C*ZCEN*RNCEN(3-1))
# 		CCC(10-1) =  .0D0
# 	   END
# 	8: BEGIN ; Ice-cream cone
# 		CCC(1-1)	=  1.0D0
# 		CCC(2-1)	=  1.0D0
# 		;CCC(3-1)	=  -(TAN (!dpi/180*CONE_A))^2
# 		CCC(3-1)	=  -(TAN (THETA))^2
# 		txt = [txt, 'Cone is an Ice-cream cone']
# 	   END
#         9:
# 	else: BEGIN
#            print,'Optical surface type not found: '+oeType[iType]
# 	   END
# ENDCASE
#
# ;
# ; cylindrical case
# ;
# ;C
# ;C Set to zero the coeff. involving X for the cylinder case, after
# ;C projecting the surface on the desired plane.
# ;C
# 	IF Keyword_Set(cylindrical) THEN BEGIN
# 	  IF N_Elements(cylangle) EQ 0 THEN CIL_ANG=0.0D0 ELSE CIL_ANG=cylAngle

    if cylindrical:
        CIL_ANG = numpy.pi / 180 * CIL_ANG
        COS_CIL = numpy.cos(CIL_ANG)
        SIN_CIL = numpy.sin(CIL_ANG)
#
        A_1	 =   ccc[1-1]
        A_2	 =   ccc[2-1]
        A_3	 =   ccc[3-1]
        A_4	 =   ccc[4-1]
        A_5	 =   ccc[5-1]
        A_6	 =   ccc[6-1]
        A_7	 =   ccc[7-1]
        A_8	 =   ccc[8-1]
        A_9	 =   ccc[9-1]
        A_10 =   ccc[10-1]

#      	  CCC(1-1) =  A_1*SIN_CIL^4 + A_2*COS_CIL^2*SIN_CIL^2 - $   ;X^2
#      		    A_4*COS_CIL*SIN_CIL^3
#      	  CCC(2-1) =  A_2*COS_CIL^4 + A_1*COS_CIL^2*SIN_CIL^2 - $ ;Y^2
#                     A_4*COS_CIL^3*SIN_CIL
#      	  CCC(3-1) =  A_3						 ;Z^2
#      	  CCC(4-1) =  - 2*A_1*COS_CIL*   SIN_CIL^3 - $
# 		      2*A_2*COS_CIL^3*SIN_CIL + $
# 		      2*A_4*COS_CIL^2*SIN_CIL^2		 ;X Y
#      	  CCC(5-1) =  A_5*COS_CIL^2 - A_6*COS_CIL*SIN_CIL	 ;Y Z
#      	  CCC(6-1) =  A_6*SIN_CIL^2 - A_5*COS_CIL*SIN_CIL	 ;X Z
#      	  CCC(7-1) =  A_7*SIN_CIL^2 - A_8*COS_CIL*SIN_CIL	 ;X
#      	  CCC(8-1) =  A_8*COS_CIL^2 - A_7*COS_CIL*SIN_CIL	 ;Y
#      	  CCC(9-1) =  A_9						 ;Z
#      	  CCC(10-1)=  A_10

        ccc[1-1] =  A_1 * SIN_CIL**4 + A_2 * COS_CIL**2 * SIN_CIL**2 - A_4 * COS_CIL * SIN_CIL**3
        ccc[2-1] =  A_2 * COS_CIL**4 + A_1 * COS_CIL**2 * SIN_CIL**2 - A_4 * COS_CIL**3 * SIN_CIL
        ccc[3-1] =  A_3						 # Z^2
        ccc[4-1] =  - 2*A_1 * COS_CIL * SIN_CIL**3 - 2 * A_2 * COS_CIL**3 * SIN_CIL + 2 * A_4 * COS_CIL**2 *SIN_CIL**2 # X Y
        ccc[5-1] =  A_5 * COS_CIL**2 - A_6 * COS_CIL * SIN_CIL	 # Y Z
        ccc[6-1] =  A_6 * SIN_CIL**2 - A_5 * COS_CIL * SIN_CIL	 # X Z
        ccc[7-1] =  A_7 * SIN_CIL**2 - A_8 * COS_CIL * SIN_CIL	 # X
        ccc[8-1] =  A_8 * COS_CIL**2 - A_7 * COS_CIL * SIN_CIL	 # Y
        ccc[9-1] =  A_9						 # Z
        ccc[10-1]=  A_10



#      	ENDIF
# ;C
# ;C Set the correct mirror convexity, i.e., Z->-Z
# ;C
#      	IF CONVEX EQ 1 THEN BEGIN

#      	ENDIF

        if convex:
            ccc[5-1]  = - ccc[5-1]
            ccc[6-1]  = - ccc[6-1]
            ccc[9-1]  = - ccc[9-1]

#
#       txt=[txt,'','']
#       FOR i=0,9 DO txt=[txt,String(' c['+StrCompress(i,/Rem)+']=',$
#          ccc[i],Format='(A6,G20.10)')]
#
    print(txt)
    return ccc