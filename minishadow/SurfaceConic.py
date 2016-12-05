
import numpy

from numpy.testing import assert_equal, assert_almost_equal

# OE surface in form of conic equation:
#      ccc[0]*X^2 + ccc[1]*Y^2 + ccc[2]*Z^2 +
#      ccc[3]*X*Y + ccc[4]*Y*Z + ccc[5]*X*Z  +
#      ccc[6]*X   + ccc[7]*Y   + ccc[8]*Z + ccc[9] = 0


class SurfaceConic(object):

    def __init__(self, ccc=numpy.zeros(10)):

        if ccc is not None:
            self.ccc = ccc.copy()
        else:
            self.ccc = numpy.zeros(10)

    @classmethod
    def initialize_from_coefficients(cls, ccc):
        if numpy.array(ccc).size != 10:
            raise Exception("Invalid coefficients (dimension must be 10)")
        return SurfaceConic(ccc=ccc)

    @classmethod
    def initialize_as_plane(cls):
        return SurfaceConic(numpy.array([0,0,0,0,0,0,0,0,-1.,0]))


    def duplicate(self):
        return SurfaceConic.initialize_from_coefficients(self.ccc.copy())

    #
    # getters
    #

    def get_coefficients(self):
        return self.ccc.copy()


    #
    # setters
    #

    def set_coefficients(self,ccc):
        if numpy.array(ccc).size != 10:
            raise Exception("Invalid coefficients (dimension must be 10)")
        self.ccc = ccc


    def vector_reflection(self,v1,normal):
        tmp = v1 * normal
        tmp2 = tmp[0,:] + tmp[1,:] + tmp[2,:]
        tmp3 = normal.copy()

        for jj in (0,1,2):
            tmp3[jj,:] = tmp3[jj,:] * tmp2

        v2 = v1 - 2 * tmp3
        v2mod = numpy.sqrt(v2[0,:]**2 + v2[1,:]**2 + v2[2,:]**2)
        v2 /= v2mod

        return v2

    def get_normal(self,x2):
        # ;
        # ; Calculates the normal at each intercept [see shadow's normal.F]
        # ;
        normal = numpy.zeros_like(x2)
        ccc = self.ccc

        normal[0,:] = 2 * ccc[1-1] * x2[0,:] + ccc[4-1] * x2[1,:] + ccc[6-1] * x2[2,:] + ccc[7-1]
        normal[1,:] = 2 * ccc[2-1] * x2[1,:] + ccc[4-1] * x2[0,:] + ccc[5-1] * x2[2,:] + ccc[8-1]
        normal[2,:] = 2 * ccc[3-1] * x2[2,:] + ccc[5-1] * x2[1,:] + ccc[6-1] * x2[0,:] + ccc[9-1]

        normalmod =  numpy.sqrt( normal[0,:]**2 + normal[1,:]**2 + normal[2,:]**2 )
        normal[0,:] /= normalmod
        normal[1,:] /= normalmod
        normal[2,:] /= normalmod

        return normal


    def calculate_reflected_beam(self,newbeam):
        # ;
        # ; TRACING...
        # ;

        x1 =   newbeam.get_columns([1,2,3]) # numpy.array(a3.getshcol([1,2,3]))
        v1 =   newbeam.get_columns([4,5,6]) # numpy.array(a3.getshcol([4,5,6]))
        flag = newbeam.get_column(10)        # numpy.array(a3.getshonecol(10))


        t,iflag = self.calculate_intercept(x1,v1)
        x2 = x1 + v1 * t
        for i in range(flag.size):
            if iflag[i] < 0: flag[i] = -100


        # ;
        # ; Calculates the normal at each intercept [see shadow's normal.F]
        # ;

        normal = ccc_normal(ccc,x2)

        # ;
        # ; reflection
        # ;

        v2 = vector_reflection(v1,normal)

        # ;
        # ; writes the mirr.XX file
        # ;

        newbeam.set_column(1, x2[0])
        newbeam.set_column(2, x2[1])
        newbeam.set_column(3, x2[2])
        newbeam.set_column(4, v2[0])
        newbeam.set_column(5, v2[1])
        newbeam.set_column(6, v2[2])
        newbeam.set_column(10, flag )

        return newbeam

    def calculate_intercept(self,XIN,VIN,keep=0):
    # FUNCTION conicintercept,ccc,xIn1,vIn1,iflag,keep=keep
    #
    #
    # ;+
    # ;
    # ;       NAME:
    # ;               CONICINTERCEPT
    # ;
    # ;       PURPOSE:
    # ;               This function Calculates the intersection of a
    # ;               conic (defined by its 10 coefficients in ccc)
    # ;               with a straight line, defined by a point xIn and
    # ;               an unitary direction vector vIn
    # ;
    # ;       CATEGORY:
    # ;               SHADOW tools
    # ;
    # ;       CALLING SEQUENCE:
    # ;               t = conicIntercept(ccc,xIn,vIn,iFlag)
    # ;
    # ; 	INPUTS:
    # ;		ccc: the array with the 10 coefficients defining the
    # ;                    conic.
    # ;		xIn: a vector DblArr(3) or stack of vectors DblArr(3,nvectors)
    # ;		vIn: a vector DblArr(3) or stack of vectors DblArr(3,nvectors)
    # ;
    # ;       OUTPUTS
    # ;		t the "travelled" distance between xIn and the surface
    # ;
    # ; 	OUTPUT KEYWORD PARAMETERS
    # ;		IFLAG: A flag (negative if no intersection)
    # ;
    # ; 	KEYWORD PARAMETERS
    # ;               keep: 0 [default] keep the max t from both solutions
    # ;                     1 keep the MIN t from both solutions
    # ;                     2 keep the first solution
    # ;                     3 keep the second solution
    # ;	ALGORITHM:
    # ;		 Adapted from SHADOW/INTERCEPT
    # ;
    # ;		 Equation of the conic:
    # ;
    # ;	         c[0]*X^2 + c[1]*Y^2 + c[2]*Z^2 +
    # ;                c[3]*X*Y + c[4]*Y*Z + c[5]*X*Z  +
    # ;                c[6]*X + c[7]*Y + c[8]*Z + c[9] = 0
    # ;
    # ;       NOTE that the vectors, that are usually DblArr(3) can be
    # ;            stacks of vectors DblArr(3,nvectors). In such a case,
    # ;            the routine returns t
    # ;
    # ;
    # ;	AUTHOR:
    # ;		M. Sanchez del Rio srio@esrf.eu, Sept. 29, 2009
    # ;
    # ;	MODIFICATION HISTORY:
    # ;
    # ;-
    #
    #
    # ;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    # ;C
    # ;C	subroutine	intercept	( xin, vin, tpar, iflag)
    # ;C
    # ;C	purpose		computes the intercepts onto the mirror surface
    # ;C
    # ;C	arguments	xin	ray starting position     mirror RF
    # ;C			vin	ray direction		  mirror RF
    # ;C			tpar	distance from start of
    # ;C				intercept
    # ;C			iflag   input		1	ordinary case
    # ;C					       -1	ripple case
    # ;C			iflag	output		0	success
    # ;C					       -1       complex sol.
    # ;C
    # ;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    #

        CCC = self.ccc

        if XIN.shape==(3,):
            XIN.shape = (3,1)
        if VIN.shape==(3,):
            VIN.shape = (3,1)

        AA 	=       CCC[1-1]*VIN[1-1,:]**2  \
                        + CCC[2-1]*VIN[2-1,:]**2  \
                        + CCC[3-1]*VIN[3-1,:]**2  \
                        + CCC[4-1]*VIN[1-1,:]*VIN[2-1,:]  \
                        + CCC[5-1]*VIN[2-1,:]*VIN[3-1,:]  \
                        + CCC[6-1]*VIN[1-1,:]*VIN[3-1,:]


        BB 	=       CCC[1-1] *  XIN[1-1,:] * VIN[1-1,:]*2    \
                        + CCC[2-1] *  XIN[2-1,:] * VIN[2-1,:]*2    \
                        + CCC[3-1] *  XIN[3-1,:] * VIN[3-1,:]*2    \
                        + CCC[4-1] * (XIN[2-1,:] * VIN[1-1,:]    \
                        + XIN[1-1,:] * VIN[2-1,:])    \
                        + CCC[5-1]*(XIN[3-1,:]*VIN[2-1,:]    \
                        + XIN[2-1,:]*VIN[3-1,:])    \
                        + CCC[6-1]*(XIN[1-1,:]*VIN[3-1,:]    \
                        + XIN[3-1,:]*VIN[1-1,:])    \
                        + CCC[7-1] * VIN[1-1,:]    \
                        + CCC[8-1] * VIN[2-1,:]    \
                        + CCC[9-1] * VIN[3-1,:]

        CC 	=             CCC[1-1] * XIN[1-1,:]**2    \
                        + CCC[2-1] * XIN[2-1,:]**2    \
                        + CCC[3-1] * XIN[3-1,:]**2    \
                        + CCC[4-1] * XIN[2-1,:] * XIN[1-1,:]    \
                        + CCC[5-1] * XIN[2-1,:] * XIN[3-1,:]    \
                        + CCC[6-1] * XIN[1-1,:]*XIN[3-1,:]    \
                        + CCC[7-1] * XIN[1-1,:]    \
                        + CCC[8-1] * XIN[2-1,:]    \
                        + CCC[9-1] * XIN[3-1,:]    \
                        + CCC[10-1]


    # ;C
    # ;C Solve now the second deg. equation **
    # ;C


        DENOM = AA*0.0
        DETER = AA*0.0
        TPAR1 = AA*0.0
        TPAR2 = AA*0.0
        IFLAG = numpy.ones(AA.size) # int(AA*0)+1

        itest1 = numpy.argwhere( numpy.abs(AA) > 1e-15)

        if len(itest1) > 0:

            DENOM[itest1] = 0.5 / AA[itest1]
            DETER[itest1] = BB[itest1]**2 - CC[itest1] * AA[itest1] * 4

            TMP = DETER[itest1]

            ibad = numpy.argwhere(TMP < 0)
            if len(ibad) == 0:
                IFLAG[itest1[ibad]] = -1

            igood = numpy.argwhere(TMP >= 0)
            if len(igood) > 0:
                itmp = itest1[igood]
                TPAR1[itmp] = -(BB[itmp] + numpy.sqrt(DETER[itmp])) * DENOM[itmp]
                TPAR2[itmp] = -(BB[itmp] - numpy.sqrt(DETER[itmp])) * DENOM[itmp]

                if keep == 0:
                    TPAR = numpy.maximum(TPAR1,TPAR2)
                elif keep == 1:
                    TPAR = numpy.minimum(TPAR1,TPAR2)
                elif keep == 2:
                    TPAR = TPAR1
                elif keep == 3:
                    TPAR = TPAR2
                else:
                    TPAR = TPAR1


        return TPAR,IFLAG

    def set_cylindrical(self,CIL_ANG):

        COS_CIL = numpy.cos(CIL_ANG)
        SIN_CIL = numpy.sin(CIL_ANG)

        A_1	 =   self.ccc[1-1]
        A_2	 =   self.ccc[2-1]
        A_3	 =   self.ccc[3-1]
        A_4	 =   self.ccc[4-1]
        A_5	 =   self.ccc[5-1]
        A_6	 =   self.ccc[6-1]
        A_7	 =   self.ccc[7-1]
        A_8	 =   self.ccc[8-1]
        A_9	 =   self.ccc[9-1]
        A_10 =   self.ccc[10-1]


        self.ccc[1-1] =  A_1 * SIN_CIL**4 + A_2 * COS_CIL**2 * SIN_CIL**2 - A_4 * COS_CIL * SIN_CIL**3
        self.ccc[2-1] =  A_2 * COS_CIL**4 + A_1 * COS_CIL**2 * SIN_CIL**2 - A_4 * COS_CIL**3 * SIN_CIL
        self.ccc[3-1] =  A_3						 # Z^2
        self.ccc[4-1] =  - 2*A_1 * COS_CIL * SIN_CIL**3 - 2 * A_2 * COS_CIL**3 * SIN_CIL + 2 * A_4 * COS_CIL**2 *SIN_CIL**2 # X Y
        self.ccc[5-1] =  A_5 * COS_CIL**2 - A_6 * COS_CIL * SIN_CIL	 # Y Z
        self.ccc[6-1] =  A_6 * SIN_CIL**2 - A_5 * COS_CIL * SIN_CIL	 # X Z
        self.ccc[7-1] =  A_7 * SIN_CIL**2 - A_8 * COS_CIL * SIN_CIL	 # X
        self.ccc[8-1] =  A_8 * COS_CIL**2 - A_7 * COS_CIL * SIN_CIL	 # Y
        self.ccc[9-1] =  A_9						 # Z
        self.ccc[10-1]=  A_10



    def set_sphere_from_focal_distances(p, q, theta1, itype=1, cylindrical=0, cylangle=0.0, convex=0, anglenordeg=0):

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



    #
    # info
    #
    def info(self):
        """

        :return:
        """
        txt = ""

        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        txt += "OE surface in form of conic equation: \n"
        txt += "  ccc[0]*X^2 + ccc[1]*Y^2 + ccc[2]*Z^2  \n"
        txt += "  ccc[3]*X*Y + ccc[4]*Y*Z + ccc[5]*X*Z  \n"
        txt += "  ccc[6]*X   + ccc[7]*Y   + ccc[8]*Z + ccc[9] = 0 \n"
        txt += " with \n"
        txt += " c[0] = %f \n "%self.ccc[0]
        txt += " c[1] = %f \n "%self.ccc[1]
        txt += " c[2] = %f \n "%self.ccc[2]
        txt += " c[3] = %f \n "%self.ccc[3]
        txt += " c[4] = %f \n "%self.ccc[4]
        txt += " c[5] = %f \n "%self.ccc[5]
        txt += " c[6] = %f \n "%self.ccc[6]
        txt += " c[7] = %f \n "%self.ccc[7]
        txt += " c[8] = %f \n "%self.ccc[8]
        txt += " c[9] = %f \n "%self.ccc[9]
        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'n"

        return txt



def tests():
    #
    # initializers
    #
    a = SurfaceConic()
    print(a.info())

    a = SurfaceConic.initialize_from_coefficients([1,1,1,1,1,1,1,1,1,1])
    print(a.info())

    a = SurfaceConic.initialize_as_plane()
    print(a.info())

if __name__ == "__main__":
    tests()
