"""

Defines the shadow4 Toroid class to deal with a toroidal surface (Quartic equation).

"""
import numpy
from shadow4.optical_surfaces.s4_optical_surface import S4OpticalSurface

from shadow4.tools.arrayofvectors import vector_cross, vector_dot, vector_multiply_scalar, vector_sum, vector_diff
from shadow4.tools.arrayofvectors import vector_modulus_square, vector_modulus, vector_norm, vector_rotate_around_axis
from shadow4.tools.arrayofvectors import vector_reflection
from shadow4.tools.arrayofvectors import vector_refraction, vector_scattering
from shadow4.tools.logger import is_verbose, is_debug
from numpy import sqrt as Sqrt

class S4Toroid(S4OpticalSurface):
    """
    Class to manage toroidal optical surfaces [expressed as a quartic polynomial].

    Parameters
    ----------
    r_maj : float, optional
        Toroid major radius in m. Note that this **is not** the tangential radius mut the radius of
         the toroidal axis, therefore for the usual case of concave surface it is the
         tangential radius minus the sagittal radius.
    r_min : float, optional
        Toroid minor radius (sagittal) in m.
    f_torus : int, optional
        A flag to indicate the mirror pole location within the toroid:
        - (0): lower/outer (concave/concave),
        - (1): lower/inner (concave/convex),
        - (2): upper/inner (convex/concave),
        - (3): upper/outer (convex/convex).

    References
    ----------
    See a graphic in Pag 26 of
    https://github.com/srio/shadow3-docs/blob/master/doc/shadow-trace.pdf

    """

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

    #
    # setters + getters
    #
    def set_from_focal_distances(self, ssour, simag, theta_grazing):
        """
        Sets the toroid radii from factory parameters (p, q, theta).

        Parameters
        ----------
        p : float
            The distance from the source to the mirror pole in m.
        q : float
            The distance from the mirror pole to the image in m.
        theta_grazing : float
            The grazing angle in deg.
        cylindrical : int, optional
            Flag:  0=the surface is curved in both directions. 1=the surface is flat in one direction.
        cylangle : float, optional
            For cylindrical=1, the angle of the cylinder axis with the X axis (CCW).
        switch_convexity : int, optional
            Flag to indicate that the convexity os inverted.

        Returns
        -------
        instance of S4Conic
        """
        theta = (numpy.pi/2) - theta_grazing
        R_TANGENTIAL = ssour * simag * 2 / numpy.cos(theta) / (ssour + simag)
        R_SAGITTAL  = ssour * simag * 2 * numpy.cos(theta) / (ssour + simag)

        self.r_maj = R_TANGENTIAL - R_SAGITTAL
        self.r_min = R_SAGITTAL

        if is_verbose():
            print("R_TANGENTIAL, R_SAGITTAL: ", R_TANGENTIAL, R_SAGITTAL)
            print("Toroid r_maj, r_min: ", self.r_maj, self.r_min)

    def set_toroid_radii(self, r_maj, r_min):
        """
        Sets the toroid radii.

        Parameters
        ----------
        r_maj : float
            The toroi major radius in m (This is **not* the tangential radius).
        r_min : float
            The toroi minur radius in m.
        """
        self.r_maj = r_maj
        self.r_min = r_min

    def set_tangential_and_sagittal_radii(self, rtan, rsag):
        """
        Sets the toroid radii from the tangential and sagittal radius.

        Parameters
        ----------
        rtan : float
            The surface tangential radius in m.
        rsag : float
            The surface sagittal radius in m.
        """
        self.r_min = rsag
        self.r_maj = rtan - rsag

    def set_f_torus(self, f_torus=0):
        """
        Sets the flag for the selected section of the toroid.

        Parameters
        ----------
        f_torus : int, optional
            A flag to indicate the mirror pole location within the toroid:
            - (0): lower/outer (concave/concave),
            - (1): lower/inner (concave/convex),
            - (2): upper/inner (convex/concave),
            - (3): upper/outer (convex/convex).
        """
        self.f_torus = f_torus

    def get_toroid_radii(self):
        """
        Gets the toroid radii r_maj, r_min.
        Note that r_maj is **not** the tangential radius!

        Returns
        -------
        tuple
        (r_maj, r_min) Note that r_maj is not the tangential radius!
        """
        return self.r_maj, self.r_min

    def get_tangential_and_sagittal_radii(self, signed=0):
        """
        Gets the radii of the focusing surface.

        Parameters
        ----------
        signed : int, optional
            Flag (0): returns all positive values, (1): return negative value if radius is convex

        Returns
        -------
        tuple
        (r_rangential, r_sagittal) in m.
        """
        if self.f_torus == 0:
            rtan = self.r_maj + self.r_min
            return rtan, self.r_min
        elif self.f_torus == 1:
            rtan = self.r_maj - self.r_min
            if signed:
                return rtan, -self.r_min
            else:
                return rtan, self.r_min
        elif self.f_torus == 2:
            rtan = self.r_maj - self.r_min
            if signed:
                return -rtan, self.r_min
            else:
                return rtan, self.r_min
        elif self.f_torus == 3:
            rtan = self.r_maj + self.r_min
            if signed:
                return -rtan, -self.r_min
            else:
                return rtan, self.r_min

    def set_cylindrical(self, CIL_ANG):
        raise Exception("Cannot set_cylindrical() in a Toroid")

    def _set_cylindrical(self, CIL_ANG):
        raise Exception("Cannot set_cylindrical() in a Toroid")

    #
    # overloaded methods
    #

    def info(self):
        """
        Creates an info text.

        Returns
        -------
        str
        """
        if self.f_torus == 0:
            rtan = self.r_maj + self.r_min
            rr = "(r_maj+r_min)"
            ff = "(Lower/Outer = Concave/Concave)"
        elif self.f_torus == 1:
            rtan = self.r_maj - self.r_min
            rr = "(r_maj-r_min)"
            ff = "(Lower/Inner = Concave/Convex)"
        elif self.f_torus == 2:
            rtan = self.r_maj - self.r_min
            rr = "(r_maj-r_min)"
            ff = "(Upper/Inner = Convex/Concave)"
        elif self.f_torus == 3:
            rtan = self.r_maj + self.r_min
            rr = "(r_maj+r_min)"
            ff = "(Upper/Outer = Convex/Convex)"

        txt = ""

        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        txt += " Toroid major radius (not optical) r_maj = %f \n"%(self.r_maj)
        txt += " Toroid minor radius (sagittal) r_min = %f \n\n"%(self.r_min)
        txt += " Toroid tangential (optical) %s = %f \n"%(rr, rtan)
        txt += " Toroid sagittal (Rmin) = %f \n\n"%(self.r_min)
        txt += " Toroid f_torus = %d %s \n\n" % (self.f_torus, ff)
        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'n"

        return txt

    def duplicate(self):
        """
        Duplicates an instance of S4Toroid

        Returns
        -------
        instance of S4Toroid.
        """
        return S4Toroid(r_maj=self.r_maj, r_min=self.r_min, f_torus=self.f_torus)

    def get_normal(self, x2):
        """
        Calculates the normal vector (or stack of vectors) at a point on the surface.

        Parameters
        ----------
        x2 : numpy array
            The coordinates vector(s) of shape [3, NRAYS].

        Returns
        -------
        numpy array
            The normal vector(s) of shape [3, NRAYS].

        """

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

        PART = X_IN ** 2 + Y_IN ** 2 + Z_IN ** 2

        normal[0, :] = 4 * X_IN * (PART + self.r_maj ** 2 - self.r_min ** 2)
        normal[1, :] = 4 * Y_IN * (PART - (self.r_maj ** 2 + self.r_min ** 2))
        normal[2, :] = 4 * Z_IN * (PART - (self.r_maj ** 2 + self.r_min ** 2))

        n2 = numpy.sqrt(normal[0, :] ** 2 + normal[1, :] ** 2 + normal[2, :] ** 2)

        normal[0, :] /= n2
        normal[1, :] /= n2
        normal[2, :] /= n2

        return normal

    def calculate_intercept_and_choose_solution(self, x1, v1, reference_distance=0.0, method=0, use_newton_solution=0):
        """

        Calculates the intercept point (or stack of points) for a given ray or stack of rays,
        given a point XIN and director vector VIN.
        Uses the exact solutions of the quartic equation (intersection torus ray).

        Parameters
        ----------
        XIN : numpy array
            The coordinates of a point of origin of the ray: shape [3, NRAYS].
        VIN : numpy array
            The coordinates of a director vector the ray: shape [3, NRAYS].
        reference_distance : float, optional
            Not used in S4Toroid.
        method : int, optional
            Not used in S4Toroid.
        use_newton_solution : int, optional
            Calculate solution using the exact (0) or approximated Newton method (1). In this case, reference_distance
            is used as a first solution guess.

        Returns
        -------
        numpy array
            The selected solution (time or flight path).

        """
        if use_newton_solution == 0:
            t0, t1, t2, t3 = self.calculate_intercept(x1, v1)
            out = self.choose_solution(t0, t1, t2, t3)
            return out
        else:
            return self.calculate_intercept_and_choose_solution_newton(x1, v1,
                                                                       reference_distance=reference_distance,
                                                                       method=method)

    def calculate_intercept_and_choose_solution_newton(self, x1, v1, reference_distance=0.0, method=0):
        """

        Calculates the intercept point (or stack of points) for a given ray or stack of rays,
        given a point XIN and director vector VIN.
        Uses the Newton approximated solution.

        Parameters
        ----------
        XIN : numpy array
            The coordinates of a point of origin of the ray: shape [3, NRAYS].
        VIN : numpy array
            The coordinates of a director vector the ray: shape [3, NRAYS].
        reference_distance : float, optional
            The initial solution guess.
        method : int, optional
            Not used in S4Toroid.

        Returns
        -------
        numpy array
            The selected solution (time or flight path).

        """
        print("Using NEWTON approximated solution.")
        AA, BB, CC, DD = self._calculate_quartic_coefficients(x1, v1)

        if False: # iteration
            t = numpy.zeros_like(AA)
            i_res = numpy.zeros_like(AA)
            epsilon = 1e-10
            for i in range(AA.size):
                p  = lambda x: self._pol4 (x, ABCD=[AA[i], BB[i], CC[i], DD[i]])
                Dp = lambda x: self._dpol4(x, ABCD=[AA[i], BB[i], CC[i], DD[i]])
                t_approx, flag_approx = self._newton(p, Dp, reference_distance, epsilon=epsilon)
                t[i] = t_approx
                i_res[i] = flag_approx

        else: # vectorized
            p = lambda x: self._pol4(x, ABCD=[AA, BB, CC, DD])
            Dp = lambda x: self._dpol4(x, ABCD=[AA, BB, CC, DD])
            epsilon = 1e-10
            iterations = 10
            t, i_res = self._newton_vectorized(p, Dp, reference_distance, iterations=iterations)
            if (numpy.abs(self._pol4(t, ABCD=[AA, BB, CC, DD])) > epsilon).any():
                print("Warning: Lack of precision in Newton method.")

        return t, i_res


    def calculate_intercept(self, XIN, VIN, vectorize=0, do_test_solution=0): #todo vectorized=1,2 fail - search another solution...
        """
        Calculates the intercept point (or stack of points) for a given ray or stack of rays,
        given a point XIN and director vector VIN.

        Parameters
        ----------
        XIN : numpy array
            The coordinates of a point of origin of the ray: shape [3, NRAYS].
        VIN : numpy array
            The coordinates of a director vector the ray: shape [3, NRAYS].
        vectorize : int, optional
            Flag:
                * 2 use S4 new way to compute iterative-solutions,
                * 1 use S4 new way to compute vectorized-solutions,
                * 0 use numpy polynomial.polyroots iteratively to compute solutions.
        do_test_solution : int, optional
            Flag to test the solutions (substitute in the polynomial and check if result is zero): 0=No, 1=Yes.

        Returns
        -------
        tuple
            (t0, t1, t2, t3) The four solutions of the quartic equation (time or flight path).
        """


        AA, BB, CC, DD = self._calculate_quartic_coefficients(XIN, VIN, method=vectorize)

        # get the four solutions of the quartic equation: t^4 + AA t^3 + BB t^2 + CC t + DD = 0
        if vectorize==0: # traditional method (numpy)

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

            for k in range(AA.size):
                z = t0[k]
                result0 = z ** 4 + AA[k] * z ** 3 + BB[k] * z ** 2 + CC[k] * z + DD[k]
                z = t1[k]
                result1 = z ** 4 + AA[k] * z ** 3 + BB[k] * z ** 2 + CC[k] * z + DD[k]
                z = t2[k]
                result2 = z ** 4 + AA[k] * z ** 3 + BB[k] * z ** 2 + CC[k] * z + DD[k]
                z = t3[k]
                result3 = z ** 4 + AA[k] * z ** 3 + BB[k] * z ** 2 + CC[k] * z + DD[k]
                if numpy.abs(result0) > 1e-3: print("Error in solution 0: ray, Delta, A, B, C, D: ", k, numpy.abs(result0), AA[k], BB[k], CC[k], DD[k])
                if numpy.abs(result1) > 1e-3: print("Error in solution 1: ray, Delta, A, B, C, D: ", k, numpy.abs(result1), AA[k], BB[k], CC[k], DD[k])
                if numpy.abs(result2) > 1e-3: print("Error in solution 2: ray, Delta, A, B, C, D: ", k, numpy.abs(result2), AA[k], BB[k], CC[k], DD[k])
                if numpy.abs(result3) > 1e-3: print("Error in solution 3: ray, Delta, A, B, C, D: ", k, numpy.abs(result3), AA[k], BB[k], CC[k], DD[k])

            return t0, t1, t2, t3

        elif vectorize >= 1:  # new method
            if vectorize == 1: # new method, vectorized
                t0, t1, t2, t3 = self._solve_quartic_vectorized(AA, BB, CC, DD)
            elif vectorize == 2: # new method, iterative
                t0 = numpy.zeros_like(AA, dtype=complex)
                t1 = numpy.zeros_like(AA, dtype=complex)
                t2 = numpy.zeros_like(AA, dtype=complex)
                t3 = numpy.zeros_like(AA, dtype=complex)
                for k in range(AA.size):
                    h_output2 = self._solve_quartic(AA[k], BB[k], CC[k], DD[k])
                    t0[k] = h_output2[0]
                    t1[k] = h_output2[1]
                    t2[k] = h_output2[2]
                    t3[k] = h_output2[3]
            else:
                raise Exception("Invalid value of vectorize=" % vectorize)

            if do_test_solution:
                for k in range(AA.size):
                    z = t0[k]
                    result0 = z ** 4 + AA[k] * z ** 3 + BB[k] * z ** 2 + CC[k] * z + DD[k]
                    z = t1[k]
                    result1 = z ** 4 + AA[k] * z ** 3 + BB[k] * z ** 2 + CC[k] * z + DD[k]
                    z = t2[k]
                    result2 = z ** 4 + AA[k] * z ** 3 + BB[k] * z ** 2 + CC[k] * z + DD[k]
                    z = t3[k]
                    result3 = z ** 4 + AA[k] * z ** 3 + BB[k] * z ** 2 + CC[k] * z + DD[k]
                    if numpy.abs(result0) > 1e-3: print("Error in solution 0: ray, Delta, A, B, C, D: ", k, numpy.abs(result0), AA[k], BB[k], CC[k], DD[k])
                    if numpy.abs(result1) > 1e-3: print("Error in solution 1: ray, Delta, A, B, C, D: ", k, numpy.abs(result1), AA[k], BB[k], CC[k], DD[k])
                    if numpy.abs(result2) > 1e-3: print("Error in solution 2: ray, Delta, A, B, C, D: ", k, numpy.abs(result2), AA[k], BB[k], CC[k], DD[k])
                    if numpy.abs(result3) > 1e-3: print("Error in solution 3: ray, Delta, A, B, C, D: ", k, numpy.abs(result3), AA[k], BB[k], CC[k], DD[k])

            return t0, t1, t2, t3

        # elif vectorize == 2:  # vectorised (fqs modified) IT DOES NOT WORK!
        #
        #     # calculate solutions array
        #     #         Input data are coefficients of the Quartic polynomial of the form:
        #     #
        #     #             p[0]*x^4 + p[1]*x^3 + p[2]*x^2 + p[3]*x + p[4] = 0
        #     P = numpy.zeros((AA.size, 5))
        #     P[:, 0] = numpy.ones_like(AA)
        #     P[:, 1] = AA.copy()
        #     P[:, 2] = BB.copy()
        #     P[:, 3] = CC.copy()
        #     P[:, 4] = DD.copy()
        #
        #     from srxraylib.profiles.diaboloid.fqs import quartic_roots
        #     SOLUTION = quartic_roots(P, modified=1, zero_below=1e-6)
        #
        #     SOLUTION_T = SOLUTION.T
        #     return SOLUTION_T[0], SOLUTION_T[1], SOLUTION_T[2], SOLUTION_T[3]

    def choose_solution(self, t0, t1, t2, t3, vectorize=0, zero_below=1e-6):
        """
        Selects the wanted single solution from the total of solutions.

        Parameters
        ----------
        t0 : numpy array
            The array with the first solution.
        t1 : numpy array
            The array with the second solution.
        t2 : numpy array
            The array with the thirs solution.
        t3 : numpy array
            The array with the fourth solution.
        reference_distance : float, optional
            A reference distance. The selected solution is the closer to this reference distance.
        vectorize : int, optional
            0: iterative (loop) method,
            1: use vectorized method.
        zero_below : float, optional
            A level of zero for intermediate step in finding the solutions using vectorized=1.

        Returns
        -------
        numpy array
            The chosen solution.
        """
        i_res  = numpy.ones( t0.size )
        answer = numpy.ones( t0.size )

        if vectorize:
            h_output = numpy.stack((t0, t1, t2, t3))

            # set small imaginary part to zero (to be considered as real)
            mask_precision = numpy.abs(h_output.imag) < zero_below
            h_output[mask_precision] = numpy.real(h_output[mask_precision])

            # create a list of solutions
            ANSWERS = list(h_output[:, k] for k in range(t0.size))

            # remove the imaginary solutions
            ANSWERS = list(ANSWERS[k][ANSWERS[k].imag == 0] for k in range(t0.size))

            # sort solutions
            ANSWERS = list(numpy.sort(ANSWERS[k].real) for k in range(t0.size))

            if self.f_torus == 0:
                answer = [item[-1] for item in ANSWERS]
            elif self.f_torus == 1:
                answer = [item[-2] for item in ANSWERS]
            elif self.f_torus == 2:
                answer = [item[1] for item in ANSWERS]
            elif self.f_torus == 3:
                answer = [item[0] for item in ANSWERS]

            # find and set the bad rays (all solutions are complex)
            mask_all_imag = h_output.imag.prod(axis=0) > zero_below # != 0
            if mask_all_imag.sum() > 0:
                print("all the solutions are complex for %d rays: " % mask.sum())
                i_res[mask_all_imag] = -1
                answer[mask_all_imag] = 0.0
        else:
            for k in range(t0.size):
                h_output = numpy.array([t0[k], t1[k], t2[k], t3[k]])
                mask = numpy.abs(h_output.imag) < zero_below
                h_output[mask] = numpy.real(h_output[mask])

                if h_output.imag.prod() != 0:
                    print("all the solutions are complex for ray index: ", k, h_output)
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

    # todo: method=1 does not work
    #       method=0 does not give enough precision,
    #       (now default sis method=3).
    def surface_height(self, X, Y, solution_index=-1, method=2):
        """
        Calculates a 2D mesh array with the surface heights.

        Parameters
        ----------
        X : numpy 2D array (mesh)
            The x coordinate(s).
        Y : numpy 2D array
            The y coordinate(s).
        method : int, optional
            0 = use a fast direct method to solve the quartic equation,
            1 = use self.calculate_intercept (TODO: most times it does not work, problems in finding the roots).
            2 = Just add the tangential and sagittal circular profiles.
        solution_index : int, optional
            For method=0,1, the index of the solution used for creating the height map:
            * -1 = take the solution according with f_torus,
            * 0 = get the first solution,
            * 1 = get the second solution,
            * 2 = get the third solution,
            * 3 = get the fourth solution.

        Returns
        -------
        2D numpy array
            the height mesh.
        """
        if numpy.abs(self.r_min) < numpy.abs(X).max():
            raise Exception("Cannot calculate toroid height for points with X > r_min. Please check input parameters.")

        if numpy.abs(self.r_maj) < numpy.abs(Y).max():
            raise Exception("Cannot calculate toroid height for points with Y > r_maj. Please check input parameters.")

        if method == 0: # fast
            # memorandum: torus equation (x^2 + y^2 + z^2 + R^2 - r^2)^2 = 4 R^2 (y^2 + z^2)
            # solve the equation for x=y=0 to sort the solutions and get the shift value.
            RR = self.r_maj ** 2 - self.r_min ** 2
            B = 2 * RR - 4 * self.r_maj ** 2
            C = RR ** 2
            t1 = 0.5 * (-B + numpy.sqrt(B ** 2 - 4 * C))
            t2 = 0.5 * (-B - numpy.sqrt(B ** 2 - 4 * C))
            z1 = numpy.sqrt(t1)
            z2 = -numpy.sqrt(t1)
            z3 = numpy.sqrt(t2)
            z4 = -numpy.sqrt(t2)
            ZZ0 = numpy.array([z1, z2, z3, z4])
            if solution_index == -1: # automatic
                indices_sorted = numpy.argsort(ZZ0)
                solution_index = indices_sorted[self.f_torus]

            # now for all values of X, Y
            J = X**2 + Y**2 + self.r_maj**2 - self.r_min**2
            B = 2 * J - 4 * self.r_maj**2
            C = (J**2 - 4 * self.r_maj**2 * Y**2)

            t1 = 0.5 * (-B + numpy.sqrt(B**2 - 4 * C))
            t2 = 0.5 * (-B - numpy.sqrt(B**2 - 4 * C))
            z1 =  numpy.sqrt(t1)
            z2 = -numpy.sqrt(t1)
            z3 =  numpy.sqrt(t2)
            z4 = -numpy.sqrt(t2)
            ZZ = numpy.array([z1, z2, z3, z4])

            Zout = (ZZ[solution_index] - ZZ0[solution_index])
            # from srxraylib.plot.gol import plot, plot_image
            # plot_image(Zout, X[:,0], Y[0,:], aspect='auto')
            return Zout

        elif method == 1: # using calculate_intercept
            xx = X.flatten()
            yy = Y.flatten()
            zz = numpy.zeros_like(xx) + 10000.0

            #
            # r_min and r_maj are like in shadow3
            #

            if self.f_torus == 0:
                zz0 =  - r_maj - r_min
            elif self.f_torus == 1:
                zz0 =  - r_maj + r_min
            elif self.f_torus == 2:
                zz0 =  + r_maj - r_min
            elif self.f_torus == 3:
                zz0 =  + r_maj + r_min


            XIN = numpy.vstack((xx, yy, zz))
            VIN = numpy.vstack((numpy.zeros_like(xx), numpy.zeros_like(xx), numpy.ones_like(xx)))

            t0, t1, t2, t3 = self.calculate_intercept(XIN, VIN)
            nn = t0.size
            ZZ0 = numpy.array([t0[nn//2], t1[nn//2], t2[nn//2], t3[nn//2]])
            if solution_index == -1: # automatic
                indices_sorted = numpy.argsort(ZZ0)
                solution_index = indices_sorted[self.f_torus]

            ZZ = numpy.array([t0, t1, t2, t3])

            height = zz + ZZ[solution_index]
            height.shape = X.shape
            return height

        elif method == 2:

            if self.f_torus == 0:
                Rt =  self.r_maj + self.r_min
                sg = 1.0
            elif self.f_torus == 1:
                Rt =  self.r_maj - self.r_min
                sg = -1.0
            elif self.f_torus == 2:
                Rt =  self.r_maj - self.r_min
                sg = 1.0
            elif self.f_torus == 3:
                Rt =  self.r_maj + self.r_min
                sg = -1.0
            Rs = self.r_min

            x = X[:, 0]
            y = Y[0, :]

            height_tangential = Rt - numpy.sqrt(Rt ** 2 - y ** 2)
            height_sagittal = Rs - numpy.sqrt(Rs ** 2 - x ** 2)

            Z = numpy.zeros((x.size, y.size))

            for i in range(x.size):
                Z[i, :] = height_tangential

            for i in range(y.size):
                Z[:, i] += height_sagittal

            return Z * sg

    #
    # calculations
    #

    def switch_convexity(self):
        raise Exception("Cannot switch_convexity() in a Toroid. Select the adequated f_torus.")

    def _newton_vectorized(self, f, Df, x0, iterations=10):
        # modified from https://patrickwalls.github.io/mathematicalpython/root-finding/newton/
        xn = x0
        for n in range(0, iterations):
            fxn = f(xn)
            Dfxn = Df(xn)
            xn = xn - fxn / Dfxn
        return xn, numpy.ones_like(xn)

    def _newton(self, f, Df, x0, epsilon=1e-10, max_iter=10):
        # modified from https://patrickwalls.github.io/mathematicalpython/root-finding/newton/
        # '''Approximate solution of f(x)=0 by Newton's method.
        #
        # Parameters
        # ----------
        # f : function
        #     Function for which we are searching for a solution f(x)=0.
        # Df : function
        #     Derivative of f(x).
        # x0 : number
        #     Initial guess for a solution f(x)=0.
        # epsilon : number
        #     Stopping criteria is abs(f(x)) < epsilon.
        # max_iter : integer
        #     Maximum number of iterations of Newton's method.
        #
        # Returns
        # -------
        # xn : number
        #     Implement Newton's method: compute the linear approximation
        #     of f(x) at xn and find x intercept by the formula
        #         x = xn - f(xn)/Df(xn)
        #     Continue until abs(f(xn)) < epsilon and return xn.
        #     If Df(xn) == 0, return None. If the number of iterations
        #     exceeds max_iter, then return None.
        #
        # Examples
        # --------
        # >>> f = lambda x: x**2 - x - 1
        # >>> Df = lambda x: 2*x - 1
        # >>> newton(f,Df,1,1e-8,10)
        # Found solution after 5 iterations.
        # 1.618033988749989
        # '''
        xn = x0
        for n in range(0, max_iter):
            fxn = f(xn)
            if numpy.abs(fxn) < epsilon:
                print('Found solution after', n, 'iterations.')
                return xn, 1
            Dfxn = Df(xn)
            if Dfxn == 0:
                print('Zero derivative. No solution found.')
                return x0, -1
            xn = xn - fxn / Dfxn
        print('Exceeded maximum iterations. No solution found.')
        return x0, -1

    def _pol4(self, z0, ABCD=None): # quartic polynomial
        return z0 ** 4 + ABCD[0] * z0 ** 3 + ABCD[1] * z0 ** 2 + ABCD[2] * z0 + ABCD[3]

    def _dpol4(self, z0, ABCD=None): # derivative of the quartic polynomial
        return 4 * z0 ** 3 + 3 * ABCD[0] * z0 ** 2 + 2 * ABCD[1] * z0 + ABCD[2]

    @classmethod
    def _solve_quartic(cls, b, c, d, e):
        # See General Formula for Roots in https://en.wikipedia.org/wiki/Quartic_function
        e += 0j

        D1 = 2 * c**3 - 9 * b * c * d + 27 * b**2 * e + 27 * d**2 - 72 * c * e
        D0 = c**2 - 3 * b * d + 12 * e # Delta0 in https://en.wikipedia.org/wiki/Quartic_function
        D = (D1**2 - 4 * D0**3) / (-27)


        k = (8 * c - 3 * b**2) / 8


        Q = numpy.power(2, -1.0/3) * numpy.power((D1 + numpy.sqrt(D1**2 - 4 * D0**3)), 1.0/3)
        S = 0.5 * numpy.sqrt((1.0/3) * (Q + D0 / Q) - 2 * k / 3)
        if numpy.abs(S) < 1e-3: #
            # print(">>>> changed sign of sqrt", numpy.abs(S), numpy.abs(Q))
            Q = numpy.power(2, -1.0/3) * numpy.power((D1 - numpy.sqrt(D1**2 - 4 * D0**3)), 1.0/3)
            S = 0.5 * numpy.sqrt((1.0/3) * (Q + D0 / Q) - 2 * k / 3)

        if numpy.abs(S) < 1e-6: # If after changing sign is still zero, the depressed quartic is biquadratic
            # print(">>>> S~0: ", S)
            # if numpy.abs(Q) < 1e-6: print("**** Q~0: ", Q)

            # depressed quartic
            p = (8 * c - 3 * b**2) / 8
            q = (b**3 - 4 * b * c + 8 * d) / 8
            r = (- 3 * b**4 + 256 * e - 64 * b * d + 16 * b**2 * c) / 256
            # print(">>>> Depressed quartic y**4 + p y**2 + q y + r = 0: ")
            # print(">>>>   p: ", p)
            # print(">>>>   q: ", q)
            # print(">>>>   r: ", r)

            if True: # numpy.abs(q) < 1e-2: # supposed zero, biquadratic depressed equation
                DD = p**2 - 4 * r
                zz1 = 0.5 * (-p + numpy.sqrt(DD))
                zz2 = 0.5 * (-p - numpy.sqrt(DD))
                z1 =  numpy.sqrt(zz1) - b / 4
                z2 = -numpy.sqrt(zz1) - b / 4
                z3 =  numpy.sqrt(zz2) - b / 4
                z4 = -numpy.sqrt(zz2) - b / 4
            else:
                print(">>>>   S: ", S)
                print(">>>>   q: ", q)
                raise Exception("Cannot treat S=0 and q != 0")

        else:
            m = (b**3 - 4 * b * c + 8 * d) / 8

            z1 = -b / 4 - S + 0.5 * numpy.sqrt(-4 * S**2 - 2 * k + m / S)
            z2 = -b / 4 - S - 0.5 * numpy.sqrt(-4 * S**2 - 2 * k + m / S)
            z3 = -b / 4 + S + 0.5 * numpy.sqrt(-4 * S**2 - 2 * k - m / S)
            z4 = -b / 4 + S - 0.5 * numpy.sqrt(-4 * S**2 - 2 * k - m / S)

        return z2, z3, z4, z1

    @classmethod
    def _solve_quartic_mathematica(cls, b, c, d, e):
        # Solve[ x^4 + b x^3 + c x^2 + d x + e == 0, x]  // FortranForm

        b += 0j
        c += 0j
        d += 0j
        e += 0j

        onethird = 1. / 3
        z1 = -b/4. - Sqrt(b**2/4. - (2*c)/3. + (2**onethird*(c**2 - 3*b*d + 12*e))/
              (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) +
             (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
              (3.*2**onethird))/2. - Sqrt(b**2/2. - (4*c)/3. - (2**onethird*(c**2 - 3*b*d + 12*e))/
              (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) -
             (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
              (3.*2**onethird) - (-b**3 + 4*b*c - 8*d)/
              (4.*Sqrt(b**2/4. - (2*c)/3. + (2**onethird*(c**2 - 3*b*d + 12*e))/
                   (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) +
                  (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
                   (3.*2**onethird))))/2.

        z2 = -b/4. - Sqrt(b**2/4. - (2*c)/3. +
              (2**onethird*(c**2 - 3*b*d + 12*e))/
               (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) +
              (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
               (3.*2**onethird))/2. + Sqrt(b**2/2. - (4*c)/3. - (2**onethird*(c**2 - 3*b*d + 12*e))/
               (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) -
              (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
               (3.*2**onethird) - (-b**3 + 4*b*c - 8*d)/
               (4.*Sqrt(b**2/4. - (2*c)/3. + (2**onethird*(c**2 - 3*b*d + 12*e))/
                    (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) +
                   (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
                    (3.*2**onethird))))/2.

        z3 = -b/4. + Sqrt(b**2/4. - (2*c)/3. +
             (2**onethird*(c**2 - 3*b*d + 12*e))/
              (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) +
             (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
              (3.*2**onethird))/2. - Sqrt(b**2/2. - (4*c)/3. - (2**onethird*(c**2 - 3*b*d + 12*e))/
              (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) -
             (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
              (3.*2**onethird) + (-b**3 + 4*b*c - 8*d)/
              (4.*Sqrt(b**2/4. - (2*c)/3. + (2**onethird*(c**2 - 3*b*d + 12*e))/
                   (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) +
                  (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
                   (3.*2**onethird))))/2.


        z4 = -b/4. + Sqrt(b**2/4. - (2*c)/3. +
             (2**onethird*(c**2 - 3*b*d + 12*e))/
              (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) +
             (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
              (3.*2**onethird))/2. + Sqrt(b**2/2. - (4*c)/3. - (2**onethird*(c**2 - 3*b*d + 12*e))/
              (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) -
             (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
              (3.*2**onethird) + (-b**3 + 4*b*c - 8*d)/
              (4.*Sqrt(b**2/4. - (2*c)/3. + (2**onethird*(c**2 - 3*b*d + 12*e))/
                   (3.*(2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird) +
                  (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e + Sqrt(-4*(c**2 - 3*b*d + 12*e)**3 + (2*c**3 - 9*b*c*d + 27*d**2 + 27*b**2*e - 72*c*e)**2))**onethird/
                   (3.*2**onethird))))/2.


        return z1, z4, z3, z2

    @classmethod
    def _solve_quartic_vectorized(cls, b, c, d, e):

        e = e + numpy.zeros_like(e, dtype=complex)

        D1 = 2 * c**3 - 9 * b * c * d + 27 * b**2 * e + 27 * d**2 - 72 * c * e
        D0 = c**2 - 3 * b * d + 12 * e
        # D = (D1**2 - 4 * D0**3) / (-27)

        k = (8 * c - 3 * b**2) / 8

        Q = numpy.power(2, -1.0/3) * numpy.power((D1 + numpy.sqrt(D1**2 - 4 * D0**3)), 1.0/3)
        S = 0.5 * numpy.sqrt((1.0/3) * (Q + D0 / Q) - 2 * k / 3)
        mask = numpy.abs(S) < 1e-3
        Q[mask] = numpy.power(2, -1.0 / 3) * numpy.power((D1[mask] - numpy.sqrt(D1[mask] ** 2 - 4 * D0[mask] ** 3)), 1.0 / 3)
        S[mask] = 0.5 * numpy.sqrt((1.0 / 3) * (Q[mask] + D0[mask] / Q[mask]) - 2 * k[mask] / 3)


        m = (b**3 - 4 * b * c + 8 * d) / 8

        z1 = -b / 4 - S + 0.5 * numpy.sqrt(-4 * S**2 - 2 * k + m / S)
        z2 = -b / 4 - S - 0.5 * numpy.sqrt(-4 * S**2 - 2 * k + m / S)
        z3 = -b / 4 + S + 0.5 * numpy.sqrt(-4 * S**2 - 2 * k - m / S)
        z4 = -b / 4 + S - 0.5 * numpy.sqrt(-4 * S**2 - 2 * k - m / S)

        # fix the cases when the depressed quartic is biquadratic
        newmask = numpy.abs(S) < 1e-6
        p = (8 * c - 3 * b ** 2) / 8
        q = (b ** 3 - 4 * b * c + 8 * d) / 8
        r = (- 3 * b ** 4 + 256 * e - 64 * b * d + 16 * b ** 2 * c) / 256

        DD = p ** 2 - 4 * r
        zz1 = 0.5 * (-p + numpy.sqrt(DD))
        zz2 = 0.5 * (-p - numpy.sqrt(DD))
        z1[newmask] =  numpy.sqrt(zz1[newmask]) - b[newmask] / 4
        z2[newmask] = -numpy.sqrt(zz1[newmask]) - b[newmask] / 4
        z3[newmask] =  numpy.sqrt(zz2[newmask]) - b[newmask] / 4
        z4[newmask] = -numpy.sqrt(zz2[newmask]) - b[newmask] / 4

        return z2, z3, z4, z1

    def _calculate_quartic_coefficients(self, XIN, VIN, method=1):
        #calculates the coefficients of the quartic polynomial resulting from
        #the intersection of the torus with a ray.
        # For the new equations, see https://arxiv.org/pdf/2301.03191.pdf but pay attention that the
        # equations are full of typos...

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


        if method==0: # shadow3
            #     ! ** Evaluates the quartic coefficients **
            # z^4 + AA z^3 + BB z^2 + CC z + DD = 0
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

            DD	= P1**4 + P2**4 + P3**4 + \
                2 * P1**2 * P2**2 + 2 * P1**2 * P3**2 + \
                2 * P2**2 * P3**2 + \
                2 * A * P1**2 + 2 * B * P2**2 + 2 * B * P3**2 + \
                A**2

            AA.shape = -1
            BB.shape = -1
            CC.shape = -1
            DD.shape = -1

        else: # shadow4 (much simpler, the same result...)
            A = r_maj ** 2 - r_min ** 2

            PdotV = P1 * V1 + P2 * V2 + P3 * V3
            PdotP = P1 ** 2 + P2 ** 2 + P3 ** 2

            AA = 4 * PdotV
            BB = 4 * PdotV ** 2 + 2 * (A + PdotP) - 4 * r_maj ** 2 * (V2 ** 2 + V3 ** 2)
            CC = 4 * PdotV * (A + PdotP) - 8 * r_maj ** 2 * (P2 * V2 + P3 * V3)
            DD = (A + PdotP) ** 2 - 4 * r_maj ** 2 * (P2 ** 2 + P3 ** 2)

        return AA, BB, CC, DD

if __name__ == "__main__":

    # t = S4Toroid(r_maj=5000.0, r_min=100, f_torus=3)
    # # t.set_from_focal_distances(30.0, 10.0, 0.003)
    # print(t.info())
    #
    # print("Rmaj, Rmin", t.get_toroid_radii())
    # print("Rtan, Rsag", t.get_tangential_and_sagittal_radii())
    # print("Rtan, Rsag (signed)", t.get_tangential_and_sagittal_radii(signed=1))
    #
    # from srxraylib.plot.gol import plot_surface
    #
    # x = numpy.linspace(-0.01, 0.01, 100)
    # y = numpy.linspace(-0.03, 0.03, 200)
    # X = numpy.outer(x,numpy.ones_like(y))
    # Y = numpy.outer(numpy.ones_like(x),y)
    #
    # Z = t.surface_height(X, Y, method=0, solution_index=-1)
    # plot_surface(Z, x, y, xtitle="x")
    #
    # x2 = numpy.zeros((3, 10))
    # x2 = numpy.random.rand(30).reshape((3, 10))
    # print("normal: ", t.get_normal(x2))

    BB, CC, DD, EE = 20.246392743580742, 993.7703042220965, 9022.71585643664, -4285.150145292282
    print("using numpy: ", numpy.polynomial.polynomial.polyroots([EE, DD, CC, BB, 1.0]))
    print("using new (itemized): ", S4Toroid._solve_quartic(BB, CC, DD, EE))
    print("using new (vectorized): ", S4Toroid._solve_quartic_vectorized(numpy.array([BB]), numpy.array([CC]), numpy.array([DD]), numpy.array([EE])))
    print("using mathematica: ", S4Toroid._solve_quartic_mathematica(BB, CC, DD, EE))

