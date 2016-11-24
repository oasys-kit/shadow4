FUNCTION conicintercept,ccc,xIn1,vIn1,iflag,keep=keep


;+
;
;       NAME:
;               CONICINTERCEPT
;
;       PURPOSE:
;               This function Calculates the intersection of a 
;               conic (defined by its 10 coefficients in ccc) 
;               with a straight line, defined by a point xIn and
;               an unitary direction vector vIn
;
;       CATEGORY:
;               SHADOW tools
;
;       CALLING SEQUENCE:
;               t = conicIntercept(ccc,xIn,vIn,iFlag)
;
; 	INPUTS:
;		ccc: the array with the 10 coefficients defining the 
;                    conic.
;		xIn: a vector DblArr(3) or stack of vectors DblArr(3,nvectors)
;		vIn: a vector DblArr(3) or stack of vectors DblArr(3,nvectors)
;
;       OUTPUTS
;		t the "travelled" distance between xIn and the surface
;	
; 	OUTPUT KEYWORD PARAMETERS
;		IFLAG: A flag (negative if no intersection)
;
; 	KEYWORD PARAMETERS
;               keep: 0 [default] keep the max t from both solutions
;                     1 keep the MIN t from both solutions
;                     2 keep the first solution
;                     3 keep the second solution
;	ALGORITHM: 
;		 Adapted from SHADOW/INTERCEPT
;	
;		 Equation of the conic: 
;
;	         c[0]*X^2 + c[1]*Y^2 + c[2]*Z^2 +
;                c[3]*X*Y + c[4]*Y*Z + c[5]*X*Z  +
;                c[6]*X + c[7]*Y + c[8]*Z + c[9] = 0 
;
;       NOTE that the vectors, that are usually DblArr(3) can be 
;            stacks of vectors DblArr(3,nvectors). In such a case, 
;            the routine returns t
;
;
;	AUTHOR: 
;		M. Sanchez del Rio srio@esrf.eu, Sept. 29, 2009
;	
;	MODIFICATION HISTORY:
;
;-


;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
;C
;C	subroutine	intercept	( xin, vin, tpar, iflag)
;C
;C	purpose		computes the intercepts onto the mirror surface
;C
;C	arguments	xin	ray starting position     mirror RF
;C			vin	ray direction		  mirror RF
;C			tpar	distance from start of
;C				intercept
;C			iflag   input		1	ordinary case
;C					       -1	ripple case
;C			iflag	output		0	success
;C					       -1       complex sol.
;C
;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

on_error,2

; 
IF N_Elements(keep) EQ 0 THEN keep=0

IF ( (N_Elements(xIn) EQ 3) AND ( (Size(xIn))[0] EQ 1)) $
  THEN xIn=Reform(xIn,3,1) ELSE xIn=xIn1

IF ( (N_Elements(vIn) EQ 3) AND ( (Size(vIn))[0] EQ 1)) $
  THEN vIn=Reform(vIn,3,1) ELSE vIn=vIn1

;C
;C conic mirrors
;C

	AA 	=   CCC[1-1]*VIN[1-1,*]^2   $
      		  + CCC[2-1]*VIN[2-1,*]^2   $
      		  + CCC[3-1]*VIN[3-1,*]^2  $
      		  + CCC[4-1]*VIN[1-1,*]*VIN[2-1,*]   $
      		  + CCC[5-1]*VIN[2-1,*]*VIN[3-1,*]  $
                  + CCC[6-1]*VIN[1-1,*]*VIN[3-1,*]
	BB 	=   CCC[1-1]*XIN[1-1,*]*VIN[1-1,*]*2  $
      		  + CCC[2-1]*XIN[2-1,*]*VIN[2-1,*]*2   $
       		  + CCC[3-1]*XIN[3-1,*]*VIN[3-1,*]*2   $
      		  + CCC[4-1]*(XIN[2-1,*]*VIN[1-1,*]    $
      			  + XIN[1-1,*]*VIN[2-1,*])    $
        	  + CCC[5-1]*(XIN[3-1,*]*VIN[2-1,*]   $
      			  + XIN[2-1,*]*VIN[3-1,*])  $
        	  + CCC[6-1]*(XIN[1-1,*]*VIN[3-1,*]   $
      		  	  + XIN[3-1,*]*VIN[1-1,*])   $
        	  + CCC[7-1]*VIN[1-1,*]   $
      		  + CCC[8-1]*VIN[2-1,*]   $
      		  + CCC[9-1]*VIN[3-1,*]   ;$
;      		  + CCC[10-1]
	CC 	=   CCC[1-1]*XIN[1-1,*]^2   $
      		  + CCC[2-1]*XIN[2-1,*]^2   $
      		  + CCC[3-1]*XIN[3-1,*]^2    $
        	  + CCC[4-1]*XIN[2-1,*]*XIN[1-1,*]   $
                  + CCC[5-1]*XIN[2-1,*]*XIN[3-1,*]   $
      		  + CCC[6-1]*XIN[1-1,*]*XIN[3-1,*]  $
          	  + CCC[7-1]*XIN[1-1,*]  $
      		  + CCC[8-1]*XIN[2-1,*]  $
      		  + CCC[9-1]*XIN[3-1,*]  $
      		  + CCC[10-1]
;C
;C Solve now the second deg. equation **
;C
;help,aa,bb,cc,xin,vin
;print,'>> AA,BB,CC: ',AA[0:10],BB[0:10],CC[0:10]

; initialize output arrays
         DENOM = AA*0.0D0
         DETER = AA*0.0D0
         TPAR =  AA*0.0D0
         TPAR1 = AA*0.0D0
         TPAR2 = AA*0.0D0
         IFLAG = Long(AA*0)+1L

; normal case
         iTest1 = where(ABS(AA) GT 1.0D-15)
         IF itest1[0] NE -1 THEN BEGIN
;print,'>>conicIntercept: conic ',N_Elements(iTest1)
	   DENOM[iTest1] = 0.5D0/AA[iTest1]
	 
           DETER[iTest1] = BB[iTest1]^2 - CC[iTest1]*AA[iTest1]*4

           TMP = DETER[iTest1] 
           iBad = Where(TMP LT 0)
           IF iBad[0] NE -1 THEN BEGIN
      	     Print,'CONICINTERCEPT: Warning. Discriminant LT 0',N_ELEMENTS(iBad)
	     iFlag[itest1[iBad]] = -1
           ENDIF
           
           iGood = Where(TMP GE 0)
           IF iGood[0] NE -1 THEN BEGIN
             iTmp = iTest1[iGood] 
	     TPAR1[iTmp] = -(BB[iTmp] + SQRT(DETER[iTmp]))*DENOM[iTmp]
	     TPAR2[iTmp] = -(BB[iTmp] - SQRT(DETER[iTmp]))*DENOM[iTmp]
           
;help,keep
             CASE keep OF
     	       0: TPAR	=   TPAR1>TPAR2
     	       1: TPAR	=   TPAR1<TPAR2
     	       2: TPAR	=   TPAR1
     	       3: TPAR	=   TPAR2
               else: TPAR = TPAR1
             ENDCASE
           ENDIF

         ENDIF

; linear equation
         iTest1 = where(AA EQ 0.0D0 AND BB NE 0.0D0)
         IF itest1[0] NE -1 THEN BEGIN
print,'>>conicIntercept: Straight line ',N_Elements(iTest1)
help,cc,bb,tpar
	   tmp = -CC[itest1]/BB[iTest1]
	   TPAR[iTest1] = tmp
	   ;TPAR[iTest1] = -CC[itest1]/BB[iTest1]
         ENDIF

; bad coefficients
         iTest1 = where( (abs(AA)+abs(BB)+abs(CC)) EQ 0.0D0)
         IF itest1[0] NE -1 THEN BEGIN
print,'>>conicIntercept: bad coefficients ',N_Elements(iTest1)
	   IFLAG[iTest1] = -1
         ENDIF


        TPAR = Reform(TPAR)

        IF N_Elements(TPAR) EQ 1 THEN TPAR=TPAR[0]
        IF N_Elements(IFLAG) EQ 1 THEN IFLAG=IFLAG[0]

     	RETURN,TPAR

     	END
