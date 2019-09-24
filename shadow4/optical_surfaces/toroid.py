
import numpy

from numpy.testing import assert_equal, assert_almost_equal

class Toroid(object):

    def __init__(self, coeff=None):

        self.r_maj = 1e10 # initialize as plane
        self.r_min = 1e10 # initialize as plane

        if coeff is not None:
            self.coeff = coeff.copy()
        else:
            self.coeff = numpy.zeros(5)

        self.f_torus = 0        # - for fmirr=3; mirror pole location:
                                #  lower/outer (concave/concave) (0),
                                # lower/inner (concave/convex) (1),
                                # upper/inner (convex/concave) (2),
                                # upper/outer (convex/convex) (3).


    @classmethod
    def initialize_from_coefficients(cls, coeff):
        if numpy.array(coeff).size != 5:
            raise Exception("Invalid coefficients (dimension must be 5)")
        return Toroid(coeff=coeff)
    #
    # @classmethod
    # def initialize_as_plane(cls):
    #     return Conic(numpy.array([0,0,0,0,0,0,0,0,-1.,0]))

    #
    # initializers from focal distances
    #

    # @classmethod
    # def initialize_as_sphere_from_focal_distances(cls,p, q, theta1, cylindrical=0, cylangle=0.0, switch_convexity=0):
    #     ccc = Conic()
    #     ccc.set_sphere_from_focal_distances(p,q,theta1)
    #     if cylindrical:
    #         ccc.set_cylindrical(cylangle)
    #     if switch_convexity:
    #         ccc.switch_convexity()
    #     return ccc



    #
    # initializars from surface parameters
    #

    # @classmethod
    # def initialize_as_sphere_from_curvature_radius(cls, radius, cylindrical=0, cylangle=0.0, switch_convexity=0):
    #     ccc = Conic()
    #     ccc.set_sphere_from_curvature_radius(radius)
    #     if cylindrical:
    #         ccc.set_cylindrical(cylangle)
    #     if switch_convexity:
    #         ccc.switch_convexity()
    #     return ccc

    def duplicate(self):
        return Toroid.initialize_from_coefficients(self.coeff.copy())

    #
    # getters
    #

    def get_coefficients(self):
        return self.coeff.copy()


    #
    # setters
    #
    def set_from_focal_distances(self, ssour, simag, theta_grazing):

        theta = (numpy.pi/2) - theta_grazing

        R_TANGENTIAL = ssour * simag * 2 / numpy.cos(theta) / (ssour + simag)
        R_SAGITTAL  = ssour * simag * 2 * numpy.cos(theta) / (ssour + simag)
        # ! C
        # ! C NOTE : The major radius is the in reality the radius of the torus
        # ! C max. circle. The true major radius is then
        # ! C
        print(">>>>> RTAN, RSAG: ",R_TANGENTIAL,R_SAGITTAL)
        self.r_maj = R_TANGENTIAL - R_SAGITTAL
        self.r_min = R_SAGITTAL
        #TODO

    def set_toroid_radii(self,r_maj,r_min):
        self.r_maj = r_maj
        self.r_min = r_min

    def set_tangential_and_sagittal_radii(self,rtan,rsag):
        self.r_min = rsag
        self.r_maj = rtan - rsag

    def get_toroid_radii(self):
        return self.r_maj,self.r_min

    def get_tangential_and_sagittal_radii(self):
        return self.r_maj+self.r_min, self.r_min

    def set_coefficients(self,coeff):
        if numpy.array(coeff).size != 5:
            raise Exception("Invalid coefficients (dimension must be 5)")
        self.coeff = coeff


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
        # ; Calculates the normal at intercept points x2 [see shadow's normal.F]
        # ;

        normal = numpy.zeros_like(x2)

        X_IN = x2[0]
        Y_IN = x2[1]
        Z_IN = x2[2]

        # normal[0,:] = 2 * self.ccc[1-1] * x2[0,:] + self.ccc[4-1] * x2[1,:] + self.ccc[6-1] * x2[2,:] + self.ccc[7-1]
        # normal[1,:] = 2 * self.ccc[2-1] * x2[1,:] + self.ccc[4-1] * x2[0,:] + self.ccc[5-1] * x2[2,:] + self.ccc[8-1]
        # normal[2,:] = 2 * self.ccc[3-1] * x2[2,:] + self.ccc[5-1] * x2[1,:] + self.ccc[6-1] * x2[0,:] + self.ccc[9-1]
        #
        # normalmod =  numpy.sqrt( normal[0,:]**2 + normal[1,:]**2 + normal[2,:]**2 )


        # ! ** Torus case. The z coordinate is offsetted due to the change in
        # ! ** ref. frame for this case.
        # 	  IF (F_TORUS.EQ.0) THEN
        #      	    Z_IN 	= Z_IN - R_MAJ - R_MIN
        # 	  ELSE IF (F_TORUS.EQ.1) THEN
        #      	    Z_IN 	= Z_IN - R_MAJ + R_MIN
        # 	  ELSE IF (F_TORUS.EQ.2) THEN
        #      	    Z_IN 	= Z_IN + R_MAJ - R_MIN
        # 	  ELSE IF (F_TORUS.EQ.3) THEN
        #      	    Z_IN 	= Z_IN + R_MAJ + R_MIN
        # 	  END IF
        #
        if self.f_torus == 0:
            Z_IN = Z_IN - self.r_maj - self.r_min
        elif self.f_torus == 1:
            Z_IN = Z_IN - self.r_maj + self.r_min
        elif self.f_torus == 2:
            Z_IN = Z_IN + self.r_maj - self.r_min
        elif self.f_torus == 3:
            Z_IN = Z_IN + self.r_maj + self.r_min


        #      	  PART	= X_IN**2 + Y_IN**2 + Z_IN**2
        #
        #      	  VOUT(1)  = 4*X_IN*(PART + R_MAJ**2 - R_MIN**2)
        #      	  VOUT(2)  = 4*Y_IN*(PART - (R_MAJ**2 + R_MIN**2))
        #      	  VOUT(3)  = 4*Z_IN*(PART - (R_MAJ**2 + R_MIN**2))

        PART = X_IN**2 + Y_IN**2 + Z_IN**2

        normal[0,:] = 4*X_IN*(PART +  self.r_maj**2 - self.r_min**2)
        normal[1,:] = 4*Y_IN*(PART - (self.r_maj**2 + self.r_min**2))
        normal[2,:] = 4*Z_IN*(PART - (self.r_maj**2 + self.r_min**2))

        n2 = numpy.sqrt(normal[0,:]**2 + normal[1,:]**2 + normal[2,:]**2)

        normal[0,:] /= n2
        normal[1,:] /= n2
        normal[2,:] /= n2

        return normal



    def calculate_intercept(self,XIN,VIN,keep=0):

        P1 = XIN[0,:]
        P2 = XIN[1,:]
        P3 = XIN[2,:]

        V1 = VIN[0,:]
        V2 = VIN[1,:]
        V3 = VIN[2,:]


        # ! C
        # ! C move the ref. frame to the torus one.
        # ! C

        if self.f_torus == 0:
            P3 = P3 - self.r_maj - self.r_min
        elif self.f_torus == 1:
            P3 = P3 - self.r_maj + self.r_min
        elif self.f_torus == 2:
            P3 = P3 + self.r_maj - self.r_min
        elif self.f_torus == 3:
            P3 = P3 + self.r_maj + self.r_min



        #P1[-1],P2[-1],P3[-1],V1[-1],V2[-1],V3[-1]=-8.5306017543434476E-003,-2999.9940033165662   ,    -749991.61550754297      ,     6.4050812398434073E-005 , 0.99999800361779412,-1.9971624670828548E-003


        #     ! ** Evaluates the quartic coefficients **

        A	= self.r_maj**2 - self.r_min**2
        B	= - (self.r_maj**2 + self.r_min**2)

        AA	= P1*V1**3 + P2*V2**3 + P3*V3**3 + \
            V1*V2**2*P1 + V1**2*V2*P2 + \
            V1*V3**2*P1 + V1**2*V3*P3 + \
            V2*V3**2*P2 + V2**2*V3*P3
        AA	= 4*AA

        BB	= 3*P1**2*V1**2 + 3*P2**2*V2**2 +  \
            3*P3**2*V3**2 + \
            V2**2*P1**2 + V1**2*P2**2 + \
            V3**2*P1**2 + V1**2*P3**2 + \
            V3**2*P2**2 + V2**2*P3**2 + \
            A*V1**2 + B*V2**2 + B*V3**2 + \
            4*V1*V2*P1*P2 +  \
            4*V1*V3*P1*P3 +  \
            4*V2*V3*P2*P3
        BB	= 2*BB

        CC	= P1**3*V1 + P2**3*V2 + P3**3*V3 + \
            P2*P1**2*V2 + P1*P2**2*V1 + \
            P3*P1**2*V3 + P1*P3**2*V1 + \
            P3*P2**2*V3 + P2*P3**2*V2 + \
            A*V1*P1 + B*V2*P2 + B*V3*P3
        CC	= 4*CC


        # TODO check DD that is the result of adding something like:
        # DD0 + A**2 = -3.16397160937e+23 + 3.16397160937e+23 = 23018340352.0
        # In fortran I get:
        #  -3.1639716093723415E+023  + 3.1639716093725710E+023 =  22951231488.000000
        DD	= P1**4 + P2**4 + P3**4 + \
            2*P1**2*P2**2 + 2*P1**2*P3**2 + \
            2*P2**2*P3**2 + \
            2*A*P1**2 + 2*B*P2**2 + 2*B*P3**2 + \
            A**2

        # for i in range(AA.size):
        #     print("R,r,AA,BB,CC,DD",self.r_maj,self.r_min,AA[i],BB[i],CC[i],DD[i])
        #     print("              P,V",P1[i],P2[i],P3[i],V1[i],V2[i],V3[i])
        #     print("              A,B,DD",A,B,DD[i],A**2,DD[i]+A**2)



        # print(coeff,coeff.shape)

        AA.shape = -1
        BB.shape = -1
        CC.shape = -1
        DD.shape = -1
        # print(">>>>",AA.shape,BB.shape,CC.shape,DD.shape)

        i_res = numpy.ones_like(AA)
        answer = numpy.ones_like(AA)
        for k in range(AA.size):
            # print("coeff: ",i,1.0,AA[i],BB[i],CC[i],DD[i])
            coeff = numpy.array([1.0,AA[k],BB[k],CC[k],DD[k]])
            # print("coeff: ",i,coeff.shape,coeff)
            h_output = numpy.roots(coeff)
            # print(i,h_output)






            # test1 = h_output.imag
            # test2 = numpy.zeros_like(test1)
            # # print(test1)

            if h_output.imag.prod() != 0:
                print("all the solutions are complex")
                i_res[k] = -1
                answer[k] = 0.0
            else:
                Answers = []

                for i in range(4):
                    if h_output[i].imag == 0:
                        Answers.append(h_output[i].real)

                #! C
                #! C Sort the real intercept in ascending order.
                #! C

                Answers = numpy.sort(numpy.array(Answers))

                # ! C
                # ! C Pick the output according to F_TORUS.
                # ! C

                # TODO check correctness of  indices not shifted
                if self.f_torus == 0:
                    answer[k] = Answers[-1]
                elif self.f_torus == 1:
                    if len(Answers) > 1:
                        answer[k] = Answers[-1]
                    else:
                        i_res[k] = -1
                elif self.f_torus == 2:
                    if len(Answers) > 1:
                        answer[k] = Answers[1]
                    else:
                        i_res[k] = -1
                elif self.f_torus == 3:
                    answer[k] = Answers[0]


        return answer,i_res




    def set_cylindrical(self,CIL_ANG):
        pass



    def switch_convexity(self):
        pass





    #
    # info
    #
    def info(self):
        """

        :return:
        """
        txt = ""

        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        txt += " Toroid major radius (tangential-sagittal) = %f \n"%(self.r_maj)
        txt += " Toroid minor radius (sagittal) = %f \n\n"%(self.r_min)
        txt += " Toroid tangential (Rmaj+Rmin) = %f \n"%(self.r_maj+self.r_min)
        txt += " Toroid sagittal (Rmin) = %f \n\n"%(self.r_min)

        txt += "OE surface in form of toroidal equation: \n"
        # txt += "  ccc[0]*X^2 + ccc[1]*Y^2 + ccc[2]*Z^2  \n"
        # txt += "  ccc[3]*X*Y + ccc[4]*Y*Z + ccc[5]*X*Z  \n"
        # txt += "  ccc[6]*X   + ccc[7]*Y   + ccc[8]*Z + ccc[9] = 0 \n"
        txt += " with \n"
        # txt += " c[0] = %f \n "%self.ccc[0]
        # txt += " c[1] = %f \n "%self.ccc[1]
        # txt += " c[2] = %f \n "%self.ccc[2]
        # txt += " c[3] = %f \n "%self.ccc[3]
        # txt += " c[4] = %f \n "%self.ccc[4]
        # txt += " c[5] = %f \n "%self.ccc[5]
        # txt += " c[6] = %f \n "%self.ccc[6]
        # txt += " c[7] = %f \n "%self.ccc[7]
        # txt += " c[8] = %f \n "%self.ccc[8]
        # txt += " c[9] = %f \n "%self.ccc[9]
        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'n"

        return txt

    def apply_specular_reflection_on_beam(self,newbeam):
        # ;
        # ; TRACING...
        # ;

        x1 =   newbeam.get_columns([1,2,3]) # numpy.array(a3.getshcol([1,2,3]))
        v1 =   newbeam.get_columns([4,5,6]) # numpy.array(a3.getshcol([4,5,6]))
        flag = newbeam.get_column(10)        # numpy.array(a3.getshonecol(10))

        t,iflag = self.calculate_intercept(x1,v1)

        # print(">>>>>",x1,t)
        # for i in range(t.size):
        #     print(">>>>",x1[0:3,i],t[i],iflag[i])

        x2 = x1 + v1 * t
        for i in range(flag.size):
            if iflag[i] < 0: flag[i] = -100


        # ;
        # ; Calculates the normal at each intercept [see shadow's normal.F]
        # ;

        normal = self.get_normal(x2)

        # for i in range(t.size):
        #     print(">>>>",t[i],normal[:,i])

        # ;
        # ; reflection
        # ;

        v2 = self.vector_reflection(v1,normal)

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


if __name__ == "__main__":

    t = Toroid()

    print(t.info())