import numpy
from shadow4.optical_surfaces.s4_optical_surface import S4OpticalSurface
from shadow4.tools.arrayofvectors import vector_refraction, vector_scattering
from shadow4.tools.arrayofvectors import vector_cross, vector_dot, vector_multiply_scalar, vector_sum, vector_diff
from shadow4.tools.arrayofvectors import vector_modulus_square, vector_modulus, vector_norm, vector_rotate_around_axis

class S4Toroid(S4OpticalSurface):
    def __init__(self,
                 r_maj=1e10,# initialize as plane  Nota bene: r_maj is not the tangential radius!!!
                 r_min=1e10,# initialize as plane
                 f_torus=0, # - for fmirr=3; mirror pole location:
                                # lower/outer (concave/concave) (0),
                                # lower/inner (concave/convex) (1),
                                # upper/inner (convex/concave) (2),
                                # upper/outer (convex/convex) (3).
                 ):

        self.r_maj = r_maj
        self.r_min = r_min
        self.f_torus = f_torus

    def set_from_focal_distances(self, ssour, simag, theta_grazing):

        theta = (numpy.pi/2) - theta_grazing
        R_TANGENTIAL = ssour * simag * 2 / numpy.cos(theta) / (ssour + simag)
        R_SAGITTAL  = ssour * simag * 2 * numpy.cos(theta) / (ssour + simag)

        print(">>>>> RTAN, RSAG: ",R_TANGENTIAL,R_SAGITTAL)
        self.r_maj = R_TANGENTIAL - R_SAGITTAL
        self.r_min = R_SAGITTAL

    def set_toroid_radii(self, r_maj, r_min):
        self.r_maj = r_maj
        self.r_min = r_min

    def set_tangential_and_sagittal_radii(self, rtan, rsag):
        self.r_min = rsag
        self.r_maj = rtan - rsag

    def get_toroid_radii(self):
        return self.r_maj, self.r_min

    def get_tangential_and_sagittal_radii(self):
        return self.r_maj + self.r_min, self.r_min

    def set_f_torus(self, f_torus=0):
        self.f_torus = f_torus

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

        if self.f_torus == 0:
            Z_IN = Z_IN - self.r_maj - self.r_min
        elif self.f_torus == 1:
            Z_IN = Z_IN - self.r_maj + self.r_min
        elif self.f_torus == 2:
            Z_IN = Z_IN + self.r_maj - self.r_min
        elif self.f_torus == 3:
            Z_IN = Z_IN + self.r_maj + self.r_min

        PART = X_IN**2 + Y_IN**2 + Z_IN**2

        normal[0,:] = 4*X_IN*(PART +  self.r_maj**2 - self.r_min**2)
        normal[1,:] = 4*Y_IN*(PART - (self.r_maj**2 + self.r_min**2))
        normal[2,:] = 4*Z_IN*(PART - (self.r_maj**2 + self.r_min**2))

        n2 = numpy.sqrt(normal[0,:]**2 + normal[1,:]**2 + normal[2,:]**2)

        normal[0,:] /= n2
        normal[1,:] /= n2
        normal[2,:] /= n2

        return normal

    def surface_height(self, X, Y, solution_index=3, method=0):

        if method == 0: # fast
            Rt, Rs = self.get_tangential_and_sagittal_radii()
            if solution_index == 0:
                return Rt + Rs + numpy.sqrt(Rt ** 2 - Y ** 2) + numpy.sqrt(Rs ** 2 - X ** 2)
            elif solution_index == 1:
                return Rt + Rs + numpy.sqrt(Rt ** 2 - Y ** 2) - numpy.sqrt(Rs ** 2 - X ** 2)
            elif solution_index == 2:
                return Rt + Rs - numpy.sqrt(Rt ** 2 - Y ** 2) + numpy.sqrt(Rs ** 2 - X ** 2)
            elif solution_index == 3:
                return Rt + Rs - numpy.sqrt(Rt**2 - Y**2) - numpy.sqrt(Rs**2 - X**2)
        else: # using calculate_intercept
            xx = X.flatten()
            yy = Y.flatten()
            zz = numpy.zeros_like(xx)

            XIN = numpy.vstack((xx, yy, zz))
            VIN = numpy.vstack((zz, zz, numpy.ones_like(xx)))

            HEIGHT, i_res = self.calculate_intercept_and_choose_solution(XIN, VIN, return_all_solutions=True)

            height = HEIGHT[solution_index,:]
            height.shape = X.shape
            return height


    def calculate_intercept(self, XIN, VIN):

        P1 = XIN[0,:]
        P2 = XIN[1,:]
        P3 = XIN[2,:]

        V1 = VIN[0,:]
        V2 = VIN[1,:]
        V3 = VIN[2,:]


        #
        # r_min and r_maj are like in shadow3
        #
        r_min = self.r_min
        r_maj = self.r_maj

        if self.f_torus == 0:
            P3 = P3 - r_maj - r_min
        elif self.f_torus == 1:
            P3 = P3 - r_maj + r_min
        elif self.f_torus == 2:
            P3 = P3 + r_maj - r_min
        elif self.f_torus == 3:
            P3 = P3 + r_maj + r_min


        #     ! ** Evaluates the quartic coefficients **

        A	=   r_maj**2 - r_min**2
        B	= -(r_maj**2 + r_min**2)

        AA	= P1 * V1**3 + P2 * V2**3 + P3 * V3**3 + \
            V1 * V2**2 * P1 + V1**2 * V2 * P2 + \
            V1 * V3**2 * P1 + V1**2 * V3 * P3 + \
            V2 * V3**2 * P2 + V2**2 * V3 * P3
        AA	= 4*AA

        BB	= 3 * P1**2 * V1**2 + 3 * P2**2 * V2**2 +  \
            3 * P3**2 * V3**2 + \
            V2**2 * P1**2 + V1**2 * P2**2 + \
            V3**2 * P1**2 + V1**2 * P3**2 + \
            V3**2 * P2**2 + V2**2 * P3**2 + \
            A * V1**2 + B * V2**2 + B * V3**2 + \
            4 * V1 * V2 * P1 * P2 +  \
            4 * V1 * V3 * P1 * P3 +  \
            4 * V2 * V3 * P2 * P3
        BB	= 2 * BB

        CC	= P1**3 * V1 + P2**3 * V2 + P3**3 * V3 + \
            P2 * P1**2 * V2 + P1 * P2**2 * V1 + \
            P3 * P1**2 * V3 + P1 * P3**2 * V1 + \
            P3 * P2**2 * V3 + P2 * P3**2 * V2 + \
            A * V1 * P1 + B * V2 * P2 + B * V3 * P3
        CC	= 4 * CC


        # TODO check DD that is the result of adding something like:
        # DD0 + A**2 = -3.16397160937e+23 + 3.16397160937e+23 = 23018340352.0
        # In fortran I get:
        #  -3.1639716093723415E+023  + 3.1639716093725710E+023 =  22951231488.000000
        DD	= P1**4 + P2**4 + P3**4 + \
            2 * P1**2 * P2**2 + 2 * P1**2 * P3**2 + \
            2 * P2**2 * P3**2 + \
            2 * A * P1**2 + 2 * B * P2**2 + 2 * B * P3**2 + \
            A**2


        AA.shape = -1
        BB.shape = -1
        CC.shape = -1
        DD.shape = -1


        t0 = numpy.zeros_like(AA, dtype=complex)
        t1 = numpy.zeros_like(AA, dtype=complex)
        t2 = numpy.zeros_like(AA, dtype=complex)
        t3 = numpy.zeros_like(AA, dtype=complex)
        for k in range(AA.size):
            h_output2 = numpy.polynomial.polynomial.polyroots([DD[k], CC[k], BB[k], AA[k], 1.0])
            t0[k] = h_output2[0]
            t1[k] = h_output2[1]
            t2[k] = h_output2[2]
            t3[k] = h_output2[3]
            # print(">>>> solutions: ", k, t0[k], t1[k], t2[k], t3[k])

        return t0, t1, t2, t3

    def calculate_intercept_and_choose_solution(self, x1, v1, reference_distance=0.0):

        t0, t1, t2, t3 = self.calculate_intercept(x1, v1)
        out = self.choose_solution(t0, t1, t2, t3)
        return out


    def choose_solution(self, t0, t1, t2, t3):
        i_res  = numpy.ones( t0.size )
        answer = numpy.ones( t0.size )

        for k in range(t0.size):

            h_output = numpy.array([t0[k], t1[k], t2[k], t3[k]])

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
                if self.f_torus == 0:
                    answer[k] = Answers[-1]
                elif self.f_torus == 1:
                    if len(Answers) > 1:
                        answer[k] = Answers[-2]
                    else:
                        i_res[k] = -1
                elif self.f_torus == 2:
                    if len(Answers) > 1:
                        answer[k] = Answers[1]
                    else:
                        i_res[k] = -1
                elif self.f_torus == 3:
                    answer[k] = Answers[0]
        return answer, i_res

    def set_cylindrical(self, CIL_ANG):
        raise Exception("Cannot set_cylindrical() in a Toroid")

    def switch_convexity(self):
        raise Exception("Cannot switch_convexity() in a Toroid")


    #
    # info
    #
    def info(self):
        txt = ""

        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        txt += " Toroid major radius (tangential-sagittal) = %f \n"%(self.r_maj)
        txt += " Toroid minor radius (sagittal) = %f \n\n"%(self.r_min)
        txt += " Toroid tangential (Rmaj+Rmin) = %f \n"%(self.r_maj+self.r_min)
        txt += " Toroid sagittal (Rmin) = %f \n\n"%(self.r_min)
        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'n"

        return txt

    # todo: move the apply_* methods to the parent class
    def apply_specular_reflection_on_beam(self,newbeam):
        # ;
        # ; TRACING...
        # ;

        x1 =   newbeam.get_columns([1,2,3]) # numpy.array(a3.getshcol([1,2,3]))
        v1 =   newbeam.get_columns([4,5,6]) # numpy.array(a3.getshcol([4,5,6]))
        flag = newbeam.get_column(10)        # numpy.array(a3.getshonecol(10))
        optical_path = newbeam.get_column(13)

        t, iflag = self.calculate_intercept_and_choose_solution(x1, v1)

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
        newbeam.set_column(13, optical_path + t)

        return newbeam, normal

    #
    # grating routines
    #
    def apply_grating_diffraction_on_beam(self, beam, ruling=[0.0], order=0, f_ruling=0):

        newbeam = beam.duplicate()

        x1 = newbeam.get_columns([1, 2, 3])  # numpy.array(a3.getshcol([1,2,3]))
        v1 = newbeam.get_columns([4, 5, 6])  # numpy.array(a3.getshcol([4,5,6]))
        flag = newbeam.get_column(10)  # numpy.array(a3.getshonecol(10))
        kin = newbeam.get_column(11) * 1e2 # in m^-1
        optical_path = newbeam.get_column(13)
        nrays = flag.size

        # t1, t2, iflag = self.calculate_intercept(x1, v1)
        reference_distance = -newbeam.get_column(2).mean() + newbeam.get_column(3).mean()
        # t = self.choose_solution(t1, t2, reference_distance=reference_distance)
        t, iflag = self.calculate_intercept_and_choose_solution(x1, v1, reference_distance=reference_distance)

        x2 = x1 + v1 * t
        for i in range(flag.size):
            if iflag[i] < 0: flag[i] = -100

        # ;
        # ; Calculates the normal at each intercept [see shadow's normal.F]
        # ;

        normal = self.get_normal(x2)

        # ;
        # ; reflection
        # ;
        # v2 =  v1.T - 2 * vector_multiply_scalar(normal.T, vector_dot(v1.T, normal.T))
        # V_OUT = v2.copy()
        # v2 = v2.T

        # ;
        # ; grating scattering
        # ;
        if True:
            DIST = x2[1]
            RDENS = 0.0
            for n in range(len(ruling)):
                RDENS += ruling[n] * DIST**n

            PHASE = optical_path + 2 * numpy.pi * order * DIST * RDENS / kin
            G_MOD = 2 * numpy.pi * RDENS * order


            # capilatized vectors are [:,3] as required for vector_* operations
            VNOR = normal.T
            VNOR = vector_multiply_scalar(VNOR, -1.0) # outward normal


            # print(">>>> VNOR: (%20.18g,%20.18g,%20.18f) mod: %20.18f" % (VNOR[-1, 0], VNOR[-1, 1], VNOR[-1, 2],
            #                                          (VNOR[-1, 0]**2 + VNOR[-1, 1]**2 + VNOR[-1, 2]**2)))

            # versors
            X_VRS = numpy.zeros((nrays,3))
            X_VRS[:,0] = 1
            Y_VRS = numpy.zeros((nrays, 3))
            Y_VRS[:,1] = 1

            if f_ruling == 0:
                G_FAC = vector_dot(VNOR, Y_VRS)
                G_FAC = numpy.sqrt(1 - G_FAC**2)
            elif f_ruling == 1:
                G_FAC = 1.0
            elif f_ruling == 5:
                G_FAC = vector_dot(VNOR, Y_VRS)
                G_FAC = numpy.sqrt(1 - G_FAC**2)

            G_MODR = G_MOD * G_FAC


            K_IN = vector_multiply_scalar(v1.T, kin)
            K_IN_NOR = vector_multiply_scalar(VNOR, vector_dot(K_IN, VNOR) )
            K_IN_PAR = vector_diff(K_IN, K_IN_NOR)


            VTAN = vector_cross(VNOR, X_VRS)
            GSCATTER = vector_multiply_scalar(VTAN, G_MODR)


            K_OUT_PAR = vector_sum(K_IN_PAR, GSCATTER)
            K_OUT_NOR = vector_multiply_scalar(VNOR,  numpy.sqrt(kin**2 - vector_modulus_square(K_OUT_PAR)))
            K_OUT = vector_sum(K_OUT_PAR, K_OUT_NOR)
            V_OUT = vector_norm(K_OUT)

        # ;
        # ; writes the mirr.XX file
        # ;

        newbeam.set_column(1, x2[0])
        newbeam.set_column(2, x2[1])
        newbeam.set_column(3, x2[2])
        newbeam.set_column(4, V_OUT.T[0])
        newbeam.set_column(5, V_OUT.T[1])
        newbeam.set_column(6, V_OUT.T[2])
        newbeam.set_column(10, flag)
        newbeam.set_column(13, optical_path + t)

        return newbeam, normal

if __name__ == "__main__":

    t = S4Toroid()
    t.set_from_focal_distances(30.0, 10.0, 0.003)
    print(t.info())

    print("Rmaj, Rmin", t.get_toroid_radii())
    print("Rtan, Rsag", t.get_tangential_and_sagittal_radii())

    from srxraylib.plot.gol import plot_surface

    x = numpy.linspace(-0.01, 0.01, 100)
    y = numpy.linspace(-0.03, 0.03, 200)
    X = numpy.outer(x,numpy.ones_like(y))
    Y = numpy.outer(numpy.ones_like(x),y)

    Z = t.surface_height(X, Y)
    plot_surface(Z, x, y, xtitle="x")

    x2 = numpy.zeros((3, 10))
    x2 = numpy.random.rand(30).reshape((3, 10))
    print("normal: ", t.get_normal(x2))
