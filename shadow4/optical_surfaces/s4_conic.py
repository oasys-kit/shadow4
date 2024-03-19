"""

Defines the shadow4 Conic class to deal with conic surfaces (plane, sphere, ellipsoid, paraboloid and hyperboloid).

The conic surface is expressed as F(x,y,z)=0, with F a quadratic finction of x, y, and z. In other words

      F(x,y,z) =
      ccc[0]*X^2 + ccc[1]*Y^2 + ccc[2]*Z^2 +
      ccc[3]*X*Y + ccc[4]*Y*Z + ccc[5]*X*Z  +
      ccc[6]*X   + ccc[7]*Y   + ccc[8]*Z + ccc[9] = 0



"""
import numpy

from shadow4.optical_surfaces.s4_optical_surface import S4OpticalSurface
from shadow4.tools.arrayofvectors import vector_refraction, vector_scattering
from shadow4.tools.arrayofvectors import vector_cross, vector_dot, vector_multiply_scalar, vector_sum, vector_diff
from shadow4.tools.arrayofvectors import vector_modulus_square, vector_modulus, vector_norm, vector_rotate_around_axis

from shadow4.tools.arrayofvectors import vector_reflection

class S4Conic(S4OpticalSurface):
    """
    Class to manage conic optical surfaces [expressed as a quadratoc polynomial].

    Parameters
    ----------
    ccc : None, list or numpy array.
        Input for the 10 conic coefficients.

    """

    def __init__(self, ccc=None):
        self.set_coefficients(ccc)

    @classmethod
    def initialize_from_coefficients(cls, ccc):
        """
        Create an instance of S4Conic from the coefficients.

        Parameters
        ----------
        ccc : None, list or numpy array.
            Input for the 10 conic coefficients.

        Returns
        -------
        instance of S4Conic
        """
        return S4Conic(ccc=ccc) # errors are taken care in set_coefficients

    #
    # getters / setters
    #

    def get_coefficients(self):
        """
        Returns the conic coefficients.

        Returns
        -------
        numpy array:
            the conic coefficients (copy).
        """
        return self.ccc.copy()

    def set_coefficients(self, ccc):
        """
        Sets the connic coefficients.

        Parameters
        ----------
        ccc : None, list or numpy array.
            Input for the 10 conic coefficients.

        """
        if ccc is not None:
            if isinstance(ccc, list): ccc = numpy.array(ccc)
            self.ccc = ccc.copy()
        else:
            self.ccc = numpy.zeros(10)

    #
    # initializers from surface external parameters
    #

    @classmethod
    def initialize_as_plane(cls):
        """
        Create an instance of S4Conic representing a plane.

        Returns
        -------
        instance of S4Conic
        """
        return S4Conic([0, 0, 0, 0, 0, 0, 0, 0, -1., 0])

    @classmethod
    def initialize_as_sphere_from_external_parameters(cls, radius, cylindrical=0, cylangle=0.0, switch_convexity=0):
        """
        Create an instance of S4Conic representing a sphere (defined from external parameters, i.e. the radius).

        Parameters
        ----------
        radius : float
            The sphere radius in m.
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
        conic = S4Conic()
        conic.set_sphere_from_external_parameters(radius)
        return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)

    @classmethod
    def initialize_as_ellipsoid_from_external_parameters(cls, AXMAJ, AXMIN, ELL_THE,
                                                         cylindrical=0, cylangle=0.0, switch_convexity=0):
        """
        Create an instance of S4Conic representing a ellipsoid from external parameters.

        Parameters
        ----------
        AXMAJ : float
            The major axis of the ellipsoid in m (a).
        AXMIN : float
            The minor axis of the ellipsoid in m (b).
        ELL_THE : float
            The angle from the line joining the center of the ellipsoid with the mirror pole to the major axis (CCW).
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
        conic = S4Conic()
        conic.set_ellipsoid_from_external_parameters(AXMAJ, AXMIN, ELL_THE)
        return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)

    @classmethod
    def initialize_as_hyperboloid_from_external_parameters(cls, AXMAJ, AXMIN, ELL_THE,
                                                         cylindrical=0, cylangle=0.0, switch_convexity=0):
        """
        Create an instance of S4Conic representing a hyperboloid from external parameters.

        Parameters
        ----------
        AXMAJ : float
            The major axis of the ellipsoid in m (a).
        AXMIN : float
            The minor axis of the ellipsoid in m (b).
        ELL_THE : float
            The angle from the line joining the center of the ellipsoid with the mirror pole to the major axis (CCW).
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
        conic = S4Conic()
        conic.set_hyperboloid_from_external_parameters(AXMAJ, AXMIN, ELL_THE)
        return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)

    @classmethod
    def initialize_as_paraboloid_from_external_parameters(cls, parabola_parameter,
                                                         cylindrical=0, cylangle=0.0, switch_convexity=0):
        """
        Create an instance of S4Conic representing a paraboloid from external parameters.

        Parameters
        ----------
        parabola_parameters : float
            The parabola parameter in m.
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
        conic = S4Conic()
        conic.set_paraboloid_from_external_parameters(parabola_parameter)
        return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)


    #
    # define coefficients from surface external parameters
    #
    def set_sphere_from_external_parameters(self, rmirr):
        """
        Sets the conic coefficients for a sphere.

        Parameters
        ----------
        rmirr : float
            The radius of the sphere in m.
        """
        self.ccc[0] =  1.0        # X^2  # = 0 in cylinder case
        self.ccc[1] =  1.0        # Y^2
        self.ccc[2] =  1.0        # Z^2
        self.ccc[3] =  0.0        # X*Y   # = 0 in cylinder case
        self.ccc[4] =  0.0        # Y*Z
        self.ccc[5] =  0.0        # X*Z   # = 0 in cylinder case
        self.ccc[6] =  0.0        # X     # = 0 in cylinder case
        self.ccc[7] =  0.0        # Y
        self.ccc[8] = -2 * rmirr  # Z
        self.ccc[9] = 0.0         # G


    def set_ellipsoid_from_external_parameters(self, AXMAJ, AXMIN, ELL_THE): # todo: remove? or change to new nomenclature?
        """
        Sets the conic coefficients for an ellipsoid given the external parameters (as defined in SHADOW3).

        Parameters
        ----------
        AXMAJ : float
            The major axis of the ellipsoid in m (a).
        AXMIN : float
            The minor axis of the ellipsoid in m (b).
        ELL_THE : float
            The angle from the line joining the center of the ellipsoid with the mirror pole to the major axis (CCW).
        """

        # tan(ELL_THE) = ZCEN / YCEN
        # therefore ZCEN = YCEN tan(ELL_THE)
        # and using the ellipse equation (YCEN / a)**2 + (ZCEN / b)**2 = 1
        # we obtain YCEN = a b / sqrt(b**2 + a**2 tan(ELL_THE)**2)
        #
        YCEN  = AXMAJ * AXMIN
        YCEN  = YCEN / numpy.sqrt(AXMIN**2 + AXMAJ**2 * numpy.tan(ELL_THE)**2)
        ZCEN  = YCEN * numpy.tan(ELL_THE)
        ZCEN  = - numpy.abs(ZCEN)
        if (numpy.cos(ELL_THE) < 0): YCEN = -numpy.abs(YCEN)
        else:                        YCEN =  numpy.abs(YCEN)

        # ;C
        # ;C Computes now the normal in the mirror center.
        # ;C
        RNCEN = numpy.zeros(3)
        RNCEN[1-1] =  0.0
        RNCEN[2-1] = -2 * YCEN / AXMAJ**2
        RNCEN[3-1] = -2 * ZCEN / AXMIN**2
        RNCEN = RNCEN / numpy.sqrt((RNCEN**2).sum())
        # ;C
        # ;C Computes the tangent versor in the mirror center.
        # ;C
        RTCEN = numpy.zeros(3)
        RTCEN[1-1] =  0.0
        RTCEN[2-1] =  RNCEN[3-1]
        RTCEN[3-1] = -RNCEN[2-1]

        # ;C Computes now the quadric coefficient with the mirror center
        # ;C located at (0,0,0) and normal along (0,0,1)

        A = 1 / AXMIN ** 2
        B = 1 / AXMAJ ** 2
        C = A
        self.ccc[0] = A
        self.ccc[1] = B * RTCEN[2 - 1] ** 2 + C * RTCEN[3 - 1] ** 2
        self.ccc[2] = B * RNCEN[2 - 1] ** 2 + C * RNCEN[3 - 1] ** 2
        self.ccc[3] = 0.0
        self.ccc[4] = 2 * (B * RNCEN[2 - 1] * RTCEN[2 - 1] + C * RNCEN[3 - 1] * RTCEN[3 - 1])
        self.ccc[5] = 0.0
        self.ccc[6] = 0.0
        self.ccc[7] = 0.0
        self.ccc[8] = 2 * (B * YCEN * RNCEN[2 - 1] + C * ZCEN * RNCEN[3 - 1])
        self.ccc[9] = 0.0


    def set_hyperboloid_from_external_parameters(self, AXMAJ, AXMIN, ELL_THE): # todo: change to new nomenclature.
        raise NotImplementedError("To be implemented...") # todo

    def set_paraboloid_from_external_parameters(self, parabola_parameter):
        raise NotImplementedError("To be implemented...") # todo

    #
    # initializers from focal distances and angle (factory parameters)
    #
    @classmethod
    def initialize_as_sphere_from_focal_distances(cls, p, q, theta_grazing,
                                                  cylindrical=0, cylangle=0.0, switch_convexity=0, verbose=1):
        """
        Creates an instance of S4Connic representing a sphere from dactory parameters (p, q, theta).

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
        verbose : int, optional
            Flag for verbose output.

        Returns
        -------
        instance of S4Conic
        """
        # todo: implement also sagittal bending?

        theta = (numpy.pi / 2) - theta_grazing
        rmirr = p * q * 2 / numpy.cos(theta) / (p + q)

        if verbose:
            txt = ""
            txt += "p=%f, q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (p, q, theta_grazing, theta)
            txt += "Radius= %f \n" % (rmirr)
            print(txt)

        ccc = [
            1.0,         # X^2  # = 0 in cylinder case
            1.0,         # Y^2
            1.0,         # Z^2
            .0 ,         # X*Y   # = 0 in cylinder case
            .0 ,         # Y*Z
            .0 ,         # X*Z   # = 0 in cylinder case
            .0 ,         # X     # = 0 in cylinder case
            .0 ,         # Y
            -2 * rmirr,  # Z
            .0           # G
            ]

        conic = S4Conic(ccc)
        return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)

    @classmethod
    def initialize_as_ellipsoid_from_focal_distances(cls, p, q, theta1,
                                                     cylindrical=0, cylangle=0.0, switch_convexity=0,
                                                     method=1, verbose=1):
        """
        Creates an instance of S4Conic representing an ellipsoid from factory parameters (p, q, theta).

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
        verbose : int, optional
            Flag for verbose output.
        method : int, optional
            See reference, 0: use Table 5, 1: use table 4.

        Returns
        -------
        Instance of S4Conic

        References
        ----------
        "Conic Surfaces and Transformations for X-Ray Beamline Optics Modeling", M. Sanchez del Rio and K. Goldberg,
        To be published (2024).

        """
        if method == 0: # See Table 5 in Sanchez del Rio and Goldberg
            ccc = [1,
                    numpy.sin(theta1) ** 2,
                    1 - (numpy.sin(theta1) * (p - q) / (p + q)) ** 2,
                    0,
                    -2 * numpy.sin(theta1) * numpy.cos(theta1) * (q - p) / (p + q),
                    0,
                    0,
                    0,
                    -4 * numpy.sin(theta1) * p * q / (p + q),
                    0,
                    ]
        else: # See Table 4 in Sanchez del Rio and Goldberg
            a = (p + q) / 2
            b = numpy.sqrt(p * q) * numpy.sin(theta1)
            c = numpy.sqrt(a ** 2 - b ** 2)
            Yc = (p ** 2 - q ** 2) / (4 * c)
            X = numpy.array([0, Yc, -b * numpy.sqrt(1 - (Yc / a) ** 2)])
            N = numpy.array([0, -2 * Yc / a ** 2, -2 * X[2] / b ** 2])
            n = N / numpy.sqrt(N[0] ** 2 + N[1] ** 2 + N[2] ** 2)

            b_over_a = b / a
            b_over_a_square = b_over_a ** 2
            ccc = [1,
                    n[1] ** 2 + b_over_a_square * n[2] ** 2,
                    n[2] ** 2 + b_over_a_square * n[1] ** 2,
                    0,
                    -2 * n[1] * n[2] * (1 - b_over_a_square),
                    0,
                    0,
                    0,
                    2 * (n[2] * X[2] + b_over_a_square * n[1] * X[1]),
                    0]
            if verbose:
                txt = ""
                txt += "p=%f, q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (p, q, theta1, numpy.pi - theta1)
                txt += 'Ellipsoid of revolution a=%f \n'%a
                txt += 'Ellipsoid of revolution b=%f \n'%b
                txt += 'Ellipsoid of revolution c=%f \n'%c
                txt += 'Ellipsoid of revolution eccentricity: %f \n'%(c/a)
                txt += 'Optical element center at: [%f,%f,%f]\n'%(X[0], X[1], X[2])
                txt += 'Optical element normal: [%f,%f,%f]\n'%(n[0], n[1], n[2])
                print(txt)
        conic = S4Conic(ccc=numpy.array(ccc))
        return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)

    @classmethod
    def initialize_as_paraboloid_from_focal_distances(cls, p, q, theta1,
                                                      cylindrical=0, cylangle=0.0, switch_convexity=0,
                                                      method=1, verbose=1):
        """
        Creates an instance of S4Conic representing a paraboloid from factory parameters (p, q, theta).

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
        verbose : int, optional
            Flag for verbose output.
        method : int, optional
            See reference, 0: use Table 5, 1: use table 4.

        Returns
        -------
        Instance of S4Conic

        References
        ----------
        "Conic Surfaces and Transformations for X-Ray Beamline Optics Modeling", M. Sanchez del Rio and K. Goldberg,
        To be published (2024).

        """
        if method == 0:
            if p > q: # focusing
                ccc = [1,
                        numpy.sin(theta1) ** 2,
                        numpy.cos(theta1) ** 2,
                        0,
                        2 * numpy.cos(theta1) * numpy.sin(theta1),
                        0,
                        0,
                        0,
                        -4 * q * numpy.sin(theta1),
                        0,
                        ]
            else: # collimating
                ccc = [1,
                        numpy.sin(theta1) ** 2,
                        numpy.cos(theta1) ** 2,
                        0,
                        -2 * numpy.cos(theta1) * numpy.sin(theta1),
                        0,
                        0,
                        0,
                        -4 * p * numpy.sin(theta1),
                        0,
                        ]
        else:
            if p > q:
                a_p = q * numpy.sin(theta1) ** 2
                X = numpy.array([0, -q * numpy.sin(2 * theta1), q * numpy.cos(theta1) ** 2])
                N = numpy.array([0, 2 * q * numpy.sin(2 * theta1), 4 * a_p])
                n = N / numpy.sqrt(N[0] ** 2 + N[1] ** 2 + N[2] ** 2)
            else:
                a_p = p * numpy.sin(theta1) ** 2
                X = numpy.array([0, p * numpy.sin(2 * theta1), p * numpy.cos(theta1) ** 2])
                N = numpy.array([0, -2 * p * numpy.sin(2 * theta1), 4 * a_p])
                n = N / numpy.sqrt(N[0] ** 2 + N[1] ** 2 + N[2] ** 2)

            ccc = [1, n[2] ** 2, n[1] ** 2, 0, 2 * n[1] * n[2], 0, 0, 0, -4 * a_p / n[2], 0]

            if verbose:
                txt = ""
                if p > q:
                    txt += "Source is at infinity\n"
                    txt += "q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (q, theta1, numpy.pi - theta1)
                else:
                    txt += "Image is at infinity\n"
                    txt += "p=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (p, theta1, numpy.pi - theta1)
                txt += 'Parabloid of revolution [Y^2=2 PARAM Z] PARAM=%f \n' % (2 * a_p)
                txt += 'Parabloid of revolution [Y^2=4 a_p Z] a_p=%f \n' % (a_p)
                txt += 'Optical element center at: [%f,%f,%f]\n' % (X[0], X[1], X[2])
                txt += 'Optical element normal: [%f,%f,%f]\n' % (n[0], n[1], n[2])
                print(txt)

        conic = S4Conic(ccc=numpy.array(ccc))
        return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)

    @classmethod
    def initialize_as_hyperboloid_from_focal_distances(cls, p, q, theta1,
                                                       cylindrical=0, cylangle=0.0, switch_convexity=0,
                                                       method=0, verbose=1):
        """
        Creates an instance of S4Conic representing a hyperboloid from factory parameters (p, q, theta).

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
        verbose : int, optional
            Flag for verbose output.
        method : int, optional
            See reference, 0: use Table 5, 1: use table 4.

        Returns
        -------
        Instance of S4Conic

        References
        ----------
        "Conic Surfaces and Transformations for X-Ray Beamline Optics Modeling", M. Sanchez del Rio and K. Goldberg,
        To be published (2024).

        """
        if method == 0:
            ccc = [1,
                    numpy.sin(theta1) ** 2,
                    1 - (numpy.sin(theta1) * (p + q) / (p - q)) ** 2,
                    0,
                    -2 * numpy.sin(theta1) * numpy.cos(theta1) * (q + p) / (q - p),
                    0,
                    0,
                    0,
                    -4 * numpy.sin(theta1) * p * q / (q - p),
                    0,
                    ]
        else:
            a = numpy.abs(p - q) / 2
            c = 0.5 * numpy.sqrt(p ** 2 + q ** 2 - 2 * p * q * numpy.cos(2 * theta1))
            b = numpy.sqrt(c ** 2 - a ** 2)
            Yc = (p ** 2 - q ** 2) / (4 * c)
            X = numpy.array([0, Yc, b * numpy.sqrt((Yc / a) ** 2 - 1)])
            N = numpy.array([0, -2 * Yc / a ** 2, 2 * X[2] / b ** 2])
            if q > p:
                N *= -1
            n = N / numpy.sqrt(N[0] ** 2 + N[1] ** 2 + N[2] ** 2)

            b_over_a = b / a
            b_over_a_square = b_over_a ** 2

            ccc = [1,
                    n[1] ** 2 - b_over_a_square * n[2] ** 2,
                    n[2] ** 2 - b_over_a_square * n[1] ** 2,
                    0,
                    -2 * n[1] * n[2] * (1 + b_over_a_square),
                    0,
                    0,
                    0,
                    2 * (n[2] * X[2] - b_over_a_square * n[1] * X[1]),
                    0]
            if verbose:
                txt = ""
                txt += "p=%f, q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (p, q, theta1, numpy.pi - theta1)
                txt += 'Hyperboloid of revolution a=%f \n' % a
                txt += 'Hyperboloid of revolution b=%f \n' % b
                txt += 'Hyperboloid of revolution c=%f \n' % c
                txt += 'Hyperboloid of revolution eccentricity: %f \n' % (c / a)
                txt += 'Optical element center at: [%f,%f,%f]\n' % (X[0], X[1], X[2])
                txt += 'Optical element normal: [%f,%f,%f]\n' % (n[0], n[1], n[2])
                print(txt)
        conic = S4Conic(ccc=numpy.array(ccc))
        return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)


    #
    # required methods
    #

    def info(self):
        """
        Creates an info text.

        Returns
        -------
        str
        """
        txt = ""
        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        txt += "OE surface in form of conic equation: \n"
        txt += "  ccc[0]*X^2 + ccc[1]*Y^2 + ccc[2]*Z^2  \n"
        txt += "  ccc[3]*X*Y + ccc[4]*Y*Z + ccc[5]*X*Z  \n"
        txt += "  ccc[6]*X   + ccc[7]*Y   + ccc[8]*Z + ccc[9] = 0 \n"
        txt += " with \n"
        txt += " c[0] = %g \n "%self.ccc[0]
        txt += " c[1] = %g \n "%self.ccc[1]
        txt += " c[2] = %g \n "%self.ccc[2]
        txt += " c[3] = %g \n "%self.ccc[3]
        txt += " c[4] = %g \n "%self.ccc[4]
        txt += " c[5] = %g \n "%self.ccc[5]
        txt += " c[6] = %g \n "%self.ccc[6]
        txt += " c[7] = %g \n "%self.ccc[7]
        txt += " c[8] = %g \n "%self.ccc[8]
        txt += " c[9] = %g \n "%self.ccc[9]
        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'n"
        return txt

    def duplicate(self):
        """
        Duplicates an instance of S4Conic

        Returns
        -------
        instance of S4Conic.
        """
        return S4Conic.initialize_from_coefficients(self.ccc)

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
        # ; Calculates the normal at intercept points x2
        # ;
        normal = numpy.zeros_like(x2)

        normal[0,:] = 2 * self.ccc[1-1] * x2[0,:] + self.ccc[4-1] * x2[1,:] + self.ccc[6-1] * x2[2,:] + self.ccc[7-1]
        normal[1,:] = 2 * self.ccc[2-1] * x2[1,:] + self.ccc[4-1] * x2[0,:] + self.ccc[5-1] * x2[2,:] + self.ccc[8-1]
        normal[2,:] = 2 * self.ccc[3-1] * x2[2,:] + self.ccc[5-1] * x2[1,:] + self.ccc[6-1] * x2[0,:] + self.ccc[9-1]

        normalmod =  numpy.sqrt( normal[0, :]**2 + normal[1, :]**2 + normal[2, :]**2 )
        normal[0, :] /= normalmod
        normal[1, :] /= normalmod
        normal[2, :] /= normalmod

        return normal


    def calculate_intercept_and_choose_solution(self, x1, v1, reference_distance=10.0, method=0):
        """

        Calculates the intercept point (or stack of points) for a given ray or stack of rays,
        given a point XIN and director vector VIN.

        Parameters
        ----------
        XIN : numpy array
            The coordinates of a point of origin of the ray: shape [3, NRAYS].
        VIN : numpy array
            The coordinates of a director vector the ray: shape [3, NRAYS].
        reference_distance : float, optional
            A reference distance. The selected solution will be the closest to this refecrence_distance.
        method : int, optional
            0: automatic selection (essentially the same as in shadow3 but replacing TSOURCE (unavailable here) by
            reference_distance).
            1: use first solution.
            2: use second solution.

        Returns
        -------
        tuple
            The selected solution (time or flight path) as (t, iflag).

        """
        t1, t2, iflag = self.calculate_intercept(x1, v1)
        t = self.choose_solution(t1, t2, reference_distance=reference_distance, method=method)
        return t, iflag

    # todo: vectorize?
    def calculate_intercept(self, XIN, VIN):
        """
        Calculates the intercept point (or stack of points) for a given ray or stack of rays,
        given a point XIN and director vector VIN.

        Parameters
        ----------
        XIN : numpy array
            The coordinates of a point of origin of the ray: shape [3, NRAYS].
        VIN : numpy array
            The coordinates of a director vector the ray: shape [3, NRAYS].

        Returns
        -------
        tuple
            The two solutions (time or flight path) as (TPAR1.real, TPAR2.real, IFLAG).


        """
        CCC = self.ccc

        if XIN.shape == (3,):
            XIN.shape = (3, 1)
        if VIN.shape == (3,):
            VIN.shape = (3, 1)

        AA = CCC[1 - 1] * VIN[1 - 1, :] ** 2 \
             + CCC[2 - 1] * VIN[2 - 1, :] ** 2 \
             + CCC[3 - 1] * VIN[3 - 1, :] ** 2 \
             + CCC[4 - 1] * VIN[1 - 1, :] * VIN[2 - 1, :] \
             + CCC[5 - 1] * VIN[2 - 1, :] * VIN[3 - 1, :] \
             + CCC[6 - 1] * VIN[1 - 1, :] * VIN[3 - 1, :]

        BB = CCC[1 - 1] * XIN[1 - 1, :] * VIN[1 - 1, :] * 2 \
             + CCC[2 - 1] * XIN[2 - 1, :] * VIN[2 - 1, :] * 2 \
             + CCC[3 - 1] * XIN[3 - 1, :] * VIN[3 - 1, :] * 2 \
             + CCC[4 - 1] * (XIN[2 - 1, :] * VIN[1 - 1, :] \
                             + XIN[1 - 1, :] * VIN[2 - 1, :]) \
             + CCC[5 - 1] * (XIN[3 - 1, :] * VIN[2 - 1, :] \
                             + XIN[2 - 1, :] * VIN[3 - 1, :]) \
             + CCC[6 - 1] * (XIN[1 - 1, :] * VIN[3 - 1, :] \
                             + XIN[3 - 1, :] * VIN[1 - 1, :]) \
             + CCC[7 - 1] * VIN[1 - 1, :] \
             + CCC[8 - 1] * VIN[2 - 1, :] \
             + CCC[9 - 1] * VIN[3 - 1, :]

        CC = CCC[1 - 1] * XIN[1 - 1, :] ** 2 \
             + CCC[2 - 1] * XIN[2 - 1, :] ** 2 \
             + CCC[3 - 1] * XIN[3 - 1, :] ** 2 \
             + CCC[4 - 1] * XIN[2 - 1, :] * XIN[1 - 1, :] \
             + CCC[5 - 1] * XIN[2 - 1, :] * XIN[3 - 1, :] \
             + CCC[6 - 1] * XIN[1 - 1, :] * XIN[3 - 1, :] \
             + CCC[7 - 1] * XIN[1 - 1, :] \
             + CCC[8 - 1] * XIN[2 - 1, :] \
             + CCC[9 - 1] * XIN[3 - 1, :] \
             + CCC[10 - 1]

        # ;C
        # ;C Solve now the second deg. equation **
        # ;C
        TPAR1 = numpy.zeros_like(AA)
        TPAR2 = numpy.zeros_like(AA)
        IFLAG = numpy.ones_like(AA)

        # TODO: remove loop!
        for i in range(AA.size):
            if numpy.abs(AA[i]) < 1e-15:
                TPAR1[i] = - CC[i] / BB[i]
                TPAR2[i] = TPAR1[i]
            else:

                DENOM = 0.5 / AA[i]
                DETER = BB[i] ** 2 - CC[i] * AA[i] * 4

                if DETER < 0.0:
                    TPAR1[i] = 0.0
                    TPAR2[i] = 0.0
                    IFLAG[i] = -1
                else:
                    TPAR1[i] = -(BB[i] + numpy.sqrt(DETER)) * DENOM
                    TPAR2[i] = -(BB[i] - numpy.sqrt(DETER)) * DENOM

        if TPAR2.size == 1:
            TPAR2 = numpy.asscalar(TPAR2)

        # todo: return here the complex solutions and make the choice in choose_solution ?
        return TPAR1.real, TPAR2.real, IFLAG

    def choose_solution(self, TPAR1, TPAR2, reference_distance=10.0, method=0):
        # method = 0: new shadow4 way (essentially the same as in shadow3
        #             but replacing TSOURCE (unavailable here) by reference_distance
        # method = 1: use first solution
        # method = 2: use second solution
        TPAR = numpy.zeros(TPAR1.size)

        if method == 0:
            for i in range(TPAR1.size):
                if ( numpy.abs(TPAR1[i]-reference_distance) <= numpy.abs(TPAR2[i]-reference_distance)):
                   TPAR[i] = TPAR1[i]
                else:
                   TPAR[i] = TPAR2[i]
        elif method == 1:
            TPAR = TPAR1
        elif method == 2:
            TPAR = TPAR2

        return TPAR

    def surface_height(self, x, y, return_solution=0):
        """
         Calculates a 2D mesh array with the surface heights.

         Parameters
         ----------
         x : float (a scalar, vector or mesh)
             The x coordinate(s).
         y : float (a scalar, vector or mesh)
             The y coordinate(s).

         return_solution : int, optional
             Flag:
             0 = guess the solution with zero at pole,
             1 = get first solution,
             2 = get second solution.

         Returns
         -------
         2D numpy array
             the height scalar/vector/mesh depending on inputs.
         Notes
         -----
         x and y must be homogeneous, otherwise an error will occur:
              - both scalars
              - both mesh
              - one scalar and another vector
         """
        return self.height(y, x, return_solution=return_solution)

    def height(self, y=0, x=0, return_solution=0):
        """
         Calculates a 2D mesh array with the surface heights.
         (The same as surface_height but x,y interchanged!).

         Parameters
         ----------
         y : float (a scalar, vector or mesh), optional
             The y coordinate(s).
         x : float (a scalar, vector or mesh)
             The x coordinate(s).


         return_solution : int, optional
             Flag:
             0 = guess the solution with zero at pole,
             1 = get first solution,
             2 = get second solution.

         Returns
         -------
         2D numpy array
             the height scalar/vector/mesh depending on inputs.
         Notes
         -----
         y and x must be homogeneous, otherwise an error will occur:
              - both scalars
              - both mesh
              - one scalar and another vector
         """
        aa = self.ccc[2]
        bb = self.ccc[4] * y + self.ccc[5] * x + self.ccc[8]
        cc = self.ccc[0] * x**2 + self.ccc[1] * y**2 + self.ccc[3] * x * y + \
            self.ccc[6] * x + self.ccc[7] * y + self.ccc[9]

        if aa != 0:
            discr = bb**2 - 4 * aa * cc + 0j
            s1 = (-bb + numpy.sqrt(discr)) / 2 / aa
            s2 = (-bb - numpy.sqrt(discr)) / 2 / aa

            if return_solution == 0: # select the solution close to zero at pole
                if numpy.abs(s1).min() < numpy.abs(s2).min():
                    ss = s1
                else:
                    ss = s2
            elif return_solution == 1:
                ss = s1
            else:
                ss = s2
        else:
            ss = -cc / bb

        return numpy.real(ss)

    # todo: move the apply_* methods to the parent class
    def apply_specular_reflection_on_beam(self, beam):
        newbeam = beam.duplicate() # DONE!!!!! warning, input newbeam is changed... TODO: change this behaviour making a copy?

        # ;
        # ; TRACING...
        # ;

        x1 = newbeam.get_columns([1, 2, 3])  # numpy.array(a3.getshcol([1,2,3]))
        v1 = newbeam.get_columns([4, 5, 6])  # numpy.array(a3.getshcol([4,5,6]))
        flag = newbeam.get_column(10)  # numpy.array(a3.getshonecol(10))
        optical_path = newbeam.get_column(13)

        t1, t2, iflag = self.calculate_intercept(x1, v1)
        reference_distance = -newbeam.get_column(2).mean() + newbeam.get_column(3).mean()
        t = self.choose_solution(t1, t2, reference_distance=reference_distance)

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
        # v2 = self.vector_reflection(v1, normal)

        # shapes are here (3, NRAYS)
        # arrayofvectors.vector_reflection use shapes (NRAYS, 3) => Need transpose
        v2 = (vector_reflection(v1.T, normal.T)).T

        # ;
        # ; writes the mirr arrays
        # ;
        newbeam.set_column(1, x2[0])
        newbeam.set_column(2, x2[1])
        newbeam.set_column(3, x2[2])
        newbeam.set_column(4, v2[0])
        newbeam.set_column(5, v2[1])
        newbeam.set_column(6, v2[2])
        newbeam.set_column(10, flag)
        newbeam.set_column(13, optical_path + t)

        return newbeam, normal

    def apply_refraction_on_beam(self,
                                 beam,
                                 refraction_index_object,
                                 refraction_index_image,
                                 apply_attenuation=0,
                                 linear_attenuation_coefficient=0.0,  # in SI, i.e. m^-1
                                 ):

        # ;
        # ; TRACING...
        # ;
        newbeam = beam.duplicate()

        x1 = newbeam.get_columns([1, 2, 3])  # numpy.array(3, npoints)
        v1 = newbeam.get_columns([4, 5, 6])  # numpy.array(3, npoints)
        flag = newbeam.get_column(10)
        k_in_mod = newbeam.get_column(11)
        optical_path = newbeam.get_column(13)

        t1, t2, iflag = self.calculate_intercept(x1, v1)
        reference_distance = -newbeam.get_column(2).mean() + newbeam.get_column(3).mean()
        t = self.choose_solution(t1, t2, reference_distance=reference_distance)

        # for i in range(t.size):
        #     print(">>>> solutions: ",t1[i],t2[i],t[i])

        x2 = x1 + v1 * t
        for i in range(flag.size):
            if iflag[i] < 0: flag[i] = -100

        # ;
        # ; Calculates the normal at each intercept [see shadow's normal.F]
        # ;
        normal = self.get_normal(x2)

        # if surface is convex normal_z > 0;  if concave normal_z < 0
        # we always want normal_z > 0:
        # if normal[2,:].mean() < 0:
        #     normal *= (-1.0)
        #     print("Warning: o.e. NORMAL INVERTED")

        # ;
        # ; refraction
        # ;

        # note that sgn=None tells vector_refraction to compute the right sign of the sqrt.
        # This is equivalent to change the direction of the normal in the case that it is an inwards normal.
        v2t = vector_refraction(v1.T, normal.T, refraction_index_object, refraction_index_image, sgn=None)
        v2 = v2t.T

        # ;
        # ; writes the mirr.XX file
        # ;

        newbeam.set_column(1, x2[0])
        newbeam.set_column(2, x2[1])
        newbeam.set_column(3, x2[2])
        newbeam.set_column(4, v2[0])
        newbeam.set_column(5, v2[1])
        newbeam.set_column(6, v2[2])
        newbeam.set_column(10, flag)
        newbeam.set_column(11, k_in_mod * refraction_index_image / refraction_index_object)
        newbeam.set_column(13, optical_path + t * refraction_index_object)

        if apply_attenuation:
            att1 = numpy.sqrt(numpy.exp(-numpy.abs(t) * linear_attenuation_coefficient))
            print(">>> mu (object space): ", linear_attenuation_coefficient)
            print(">>> attenuation of amplitudes (object space): ", att1)
            newbeam.rays[:, 7 - 1 ] *= att1
            newbeam.rays[:, 8 - 1 ] *= att1
            newbeam.rays[:, 9 - 1 ] *= att1
            newbeam.rays[:, 16 - 1] *= att1
            newbeam.rays[:, 17 - 1] *= att1
            newbeam.rays[:, 18 - 1] *= att1

        return newbeam, normal

    def apply_grating_diffraction_on_beam(self, beam, ruling=[0.0], order=0, f_ruling=0):

        newbeam = beam.duplicate()

        x1 = newbeam.get_columns([1, 2, 3])  # numpy.array(a3.getshcol([1,2,3]))
        v1 = newbeam.get_columns([4, 5, 6])  # numpy.array(a3.getshcol([4,5,6]))
        flag = newbeam.get_column(10)  # numpy.array(a3.getshonecol(10))
        kin = newbeam.get_column(11) * 1e2 # in m^-1
        optical_path = newbeam.get_column(13)
        nrays = flag.size

        t1, t2, iflag = self.calculate_intercept(x1, v1)
        reference_distance = -newbeam.get_column(2).mean() + newbeam.get_column(3).mean()
        t = self.choose_solution(t1, t2, reference_distance=reference_distance)

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


    #
    # utilities
    #

    @classmethod
    def _transform_conic(cls, conic, cylindrical, cylangle, switch_convexity):
        if cylindrical:      conic._set_cylindrical(cylangle)
        if switch_convexity: conic._switch_convexity()
        return conic


    def _set_cylindrical(self, CIL_ANG):
        COS_CIL = numpy.cos(CIL_ANG)
        SIN_CIL = numpy.sin(CIL_ANG)

        A_1	 = self.ccc[0]
        A_2	 = self.ccc[1]
        A_3	 = self.ccc[2]
        A_4	 = self.ccc[3]
        A_5	 = self.ccc[4]
        A_6	 = self.ccc[5]
        A_7	 = self.ccc[6]
        A_8	 = self.ccc[7]
        A_9	 = self.ccc[8]
        A_10 = self.ccc[9]

        self.ccc[0] =  A_1 * SIN_CIL**4 + A_2 * COS_CIL**2 * SIN_CIL**2 - A_4 * COS_CIL * SIN_CIL**3
        self.ccc[1] =  A_2 * COS_CIL**4 + A_1 * COS_CIL**2 * SIN_CIL**2 - A_4 * COS_CIL**3 * SIN_CIL
        self.ccc[2] =  A_3						                     # Z^2
        self.ccc[3] =  - 2 * A_1 * COS_CIL * SIN_CIL**3 - 2 * A_2 * COS_CIL**3 * SIN_CIL + 2 * A_4 * COS_CIL**2 * SIN_CIL**2 # X Y
        self.ccc[4] =  A_5 * COS_CIL**2 - A_6 * COS_CIL * SIN_CIL	 # Y Z
        self.ccc[5] =  A_6 * SIN_CIL**2 - A_5 * COS_CIL * SIN_CIL	 # X Z
        self.ccc[6] =  A_7 * SIN_CIL**2 - A_8 * COS_CIL * SIN_CIL	 # X
        self.ccc[7] =  A_8 * COS_CIL**2 - A_7 * COS_CIL * SIN_CIL	 # Y
        self.ccc[8] =  A_9						                     # Z
        self.ccc[9] =  A_10

    def _switch_convexity(self):
        self.ccc[5-1]  = - self.ccc[5-1]
        self.ccc[6-1]  = - self.ccc[6-1]
        self.ccc[9-1]  = - self.ccc[9-1]

    def _rotation_surface_conic(self, alpha, axis ):
        if axis == 'x':
            self._rotation_surface_conic_x(alpha)
        elif axis == 'y':
            self._rotation_surface_conic_y(alpha)
        elif axis == 'z':
            self._rotation_surface_conic_z(alpha)

    def _rotation_surface_conic_x(self, alpha):
        a = numpy.cos(alpha)
        b = numpy.sin(alpha)

        c0 = self.ccc[0]
        c1 = self.ccc[1] * a ** 2 + self.ccc[2] * b ** 2 + self.ccc[4] * a * b
        c2 = self.ccc[1] * b ** 2 + self.ccc[2] * a ** 2 - self.ccc[4] * a * b
        c3 = self.ccc[3] * a + self.ccc[5] * b
        c4 = - 2 * self.ccc[1] * a * b + 2 * self.ccc[2] * a * b + self.ccc[4] * (a ** 2 - b ** 2)
        c5 = - self.ccc[3] * b + self.ccc[5] * a
        c6 = self.ccc[6]
        c7 = self.ccc[7] * a + self.ccc[8] * b
        c8 = - self.ccc[7] * b + self.ccc[8] * a
        c9 = self.ccc[9]

        self.ccc = numpy.array([c0, c1, c2, c3, c4, c5, c6, c7, c8, c9])

    def _rotation_surface_conic_y(self,alpha):

        a = numpy.cos(alpha)
        b = numpy.sin(alpha)


        c0 = self.ccc[0] * a**2 + self.ccc[2] * b**2 - self.ccc[5] * a * b
        c1 = self.ccc[1]
        c2 = self.ccc[0] * b**2 + self.ccc[2] * a**2 + self.ccc[5] * a * b
        c3 = self.ccc[3] * a - self.ccc[4] * b
        c4 = self.ccc[3] * b + self.ccc[4] * a
        c5 = 2 * self.ccc[0] * a * b - 2 * self.ccc[2] * a * b + self.ccc[5] * (a**2 - b**2)
        c6 = self.ccc[6] *a - self.ccc[8] * b
        c7 = self.ccc[7]
        c8 = self.ccc[6] * b + self.ccc[8] * a
        c9 = self.ccc[9]

        self.ccc = numpy.array([c0, c1, c2, c3, c4, c5, c6, c7, c8, c9])

    def _rotation_surface_conic_z(self, alpha):

        a = numpy.cos(alpha)
        b = numpy.sin(alpha)

        c0 = self.ccc[0] * a ** 2 + self.ccc[1] * b ** 2 + self.ccc[3] * a * b
        c1 = self.ccc[0] * b ** 2 + self.ccc[1] * a ** 2 - self.ccc[3] * a * b
        c2 = self.ccc[2]
        c3 = - 2 * self.ccc[0] * a * b + 2 * self.ccc[1] * a * b + self.ccc[3] * (a ** 2 - b ** 2)
        c4 = self.ccc[4] * a - self.ccc[5] * b
        c5 = self.ccc[4] * b + self.ccc[5] * a
        c6 = self.ccc[6] * a + self.ccc[7] * b
        c7 = - self.ccc[6] * b + self.ccc[7] * a
        c8 = self.ccc[8]
        c9 = self.ccc[9]

        self.ccc = numpy.array([c0, c1, c2, c3, c4, c5, c6, c7, c8, c9])

    def _translation_surface_conic(self, x0, axis = 'x'):

        if axis == 'x':
            self._translation_surface_conic_x(x0)
        elif axis == 'y':
            self._translation_surface_conic_y(x0)
        elif axis == 'z':
            self._translation_surface_conic_z(x0)

    def _translation_surface_conic_x(self, x0):

        c6 = - 2 * self.ccc[0] * x0 + self.ccc[6]
        c7 = - self.ccc[3] * x0 + self.ccc[7]
        c8 = - self.ccc[5] * x0 + self.ccc[8]
        c9 = self.ccc[0] * x0**2 + self.ccc[9] - self.ccc[6] * x0

        self.ccc = numpy.array([self.ccc[0], self.ccc[1], self.ccc[2], self.ccc[3], self.ccc[4], self.ccc[5], c6, c7, c8, c9])

    def _translation_surface_conic_y(self, y0):

        c6 = - self.ccc[3] * y0 + self.ccc[6]
        c7 = - 2 * self.ccc[1] * y0 + self.ccc[7]
        c8 = - self.ccc[4] * y0 + self.ccc[8]
        c9 = self.ccc[1] * y0**2 + self.ccc[9] - self.ccc[7] * y0

        self.ccc = numpy.array([self.ccc[0], self.ccc[1], self.ccc[2], self.ccc[3], self.ccc[4], self.ccc[5], c6, c7, c8, c9])

    def _translation_surface_conic_z(self, z0):

        c6 = - self.ccc[5] * z0 + self.ccc[6]
        c7 = - self.ccc[4] * z0 + self.ccc[7]
        c8 = - 2 * self.ccc[2] * z0 + self.ccc[8]
        c9 = self.ccc[2] * z0**2 + self.ccc[9] - self.ccc[8] * z0

        self.ccc = numpy.array([self.ccc[0], self.ccc[1], self.ccc[2], self.ccc[3], self.ccc[4], self.ccc[5], c6, c7, c8, c9])



    #
    # to delete
    #

    # # todo: delete? not used.
    # def z_vs_xy(self, x, y):
    #     if isinstance(x, numpy.ndarray):
    #         pass
    #     else:
    #         x = numpy.ndarray([x])
    #         y = numpy.ndarray([y])
    #
    #     ccc = self.ccc
    #
    #     AA = ccc[2] * numpy.ones_like(x)
    #     BB = ccc[4] * y + ccc[5] * x + ccc[8]
    #     CC = ccc[0] * x ** 2 + ccc[1] * y ** 2 + ccc[3] * x * y + ccc[6] * x + ccc[7] * y + ccc[9]
    #
    #     shape_x = x.shape
    #
    #     AAf = AA.flatten()
    #     BBf = BB.flatten()
    #     CCf = CC.flatten()
    #
    #     TPAR1 = numpy.zeros_like(CCf, dtype=complex)
    #     TPAR2 = numpy.zeros_like(CCf, dtype=complex)
    #
    #     for i in range(AAf.size):
    #         roots = numpy.roots([CCf[i], BBf[i], AAf[i]])
    #         TPAR1[i] = roots[0]
    #         TPAR2[i] = roots[1]
    #
    #     TPAR2.shape = shape_x
    #
    #     if TPAR2.size == 1:
    #         TPAR2 = numpy.asscalar(TPAR2)
    #
    #     return TPAR2.real
    #
    # def height(self, x, y, return_solution=0):
    #     return self.surface_height(y, x, return_solution=return_solution)
    #
    # # todo: delete
    # @classmethod
    # def initialize_as_ellipsoid_from_focal_distances_old(cls, p, q, theta1, cylindrical=0, cylangle=0.0,
    #                                                      switch_convexity=0):
    #     conic = S4Conic()
    #     conic.set_ellipsoid_from_focal_distances(p, q, theta1)
    #     return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)
    #
    # # todo: delete
    # @classmethod
    # def initialize_as_paraboloid_from_focal_distances_old(cls, p, q, theta1, cylindrical=0, cylangle=0.0,
    #                                                       switch_convexity=0):
    #     conic = S4Conic()
    #     conic.set_paraboloid_from_focal_distances(p, q, theta1)
    #     return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)
    #
    # # todo: delete
    # @classmethod
    # def initialize_as_hyperboloid_from_focal_distances_old(cls, p, q, theta1, cylindrical=0, cylangle=0.0,
    #                                                        switch_convexity=0):
    #     conic = S4Conic()
    #     conic.set_hyperboloid_from_focal_distances(p, q, theta1)
    #     return cls._transform_conic(conic, cylindrical, cylangle, switch_convexity)
    #
    # def set_sphere_from_focal_distances(self, ssour, simag, theta_grazing, verbose=True):
    #     # todo: implement also sagittal bending
    #     print("Theta grazing is: %f" % (theta_grazing))
    #     theta = (numpy.pi / 2) - theta_grazing
    #     rmirr = ssour * simag * 2 / numpy.cos(theta) / (ssour + simag)
    #
    #     if verbose:
    #         txt = ""
    #         txt += "p=%f, q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (ssour, simag, theta_grazing, theta)
    #         txt += "Radius= %f \n" % (rmirr)
    #         print(txt)
    #
    #     self.ccc[0] = 1.0  # X^2  # = 0 in cylinder case
    #     self.ccc[1] = 1.0  # Y^2
    #     self.ccc[2] = 1.0  # Z^2
    #     self.ccc[3] = .0  # X*Y   # = 0 in cylinder case
    #     self.ccc[4] = .0  # Y*Z
    #     self.ccc[5] = .0  # X*Z   # = 0 in cylinder case
    #     self.ccc[6] = .0  # X     # = 0 in cylinder case
    #     self.ccc[7] = .0  # Y
    #     self.ccc[8] = -2 * rmirr  # Z
    #     self.ccc[9] = .0  # G
    #
    # # todo: delete
    # def set_ellipsoid_from_focal_distances(self, ssour, simag, theta_grazing, verbose=True):
    #     tkt = self.calculate_ellipsoid_parameters_from_focal_distances(ssour, simag, theta_grazing, verbose=verbose)
    #
    #     AXMAJ = tkt["AXMAJ"]
    #     AXMIN = tkt["AXMIN"]
    #     RTCEN = tkt["RTCEN"]
    #     RNCEN = tkt["RNCEN"]
    #     YCEN = tkt["YCEN"]
    #     ZCEN = tkt["ZCEN"]
    #
    #     # ;C Computes now the quadric coefficient with the mirror center
    #     # ;C located at (0,0,0) and normal along (0,0,1)
    #
    #     A = 1 / AXMIN ** 2
    #     B = 1 / AXMAJ ** 2
    #     C = A
    #     self.ccc[0] = A
    #     self.ccc[1] = B * RTCEN[2 - 1] ** 2 + C * RTCEN[3 - 1] ** 2
    #     self.ccc[2] = B * RNCEN[2 - 1] ** 2 + C * RNCEN[3 - 1] ** 2
    #     self.ccc[3] = 0.0
    #     self.ccc[4] = 2 * (B * RNCEN[2 - 1] * RTCEN[2 - 1] + C * RNCEN[3 - 1] * RTCEN[3 - 1])
    #     self.ccc[5] = 0.0
    #     self.ccc[6] = 0.0
    #     self.ccc[7] = 0.0
    #     self.ccc[8] = 2 * (B * YCEN * RNCEN[2 - 1] + C * ZCEN * RNCEN[3 - 1])
    #     self.ccc[9] = 0.0
    #
    # # todo: delete
    # def set_paraboloid_from_focal_distances(self, SSOUR, SIMAG, theta_grazing, at_infinity=None, verbose=True):
    #     # ;C
    #     # ;C Computes the parabola
    #     # ;C
    #     theta = (numpy.pi / 2) - theta_grazing
    #     COSTHE = numpy.cos(theta)
    #     SINTHE = numpy.sin(theta)
    #
    #     if at_infinity is None:
    #         if SSOUR <= SIMAG:
    #             location = "q"  # q is infinite
    #         else:
    #             location = "p"
    #     else:
    #         if at_infinity == 0:
    #             location = "p"
    #         else:
    #             location = "q"
    #
    #     if location == "q":
    #         PARAM = 2 * SSOUR * COSTHE ** 2
    #         YCEN = -SSOUR * SINTHE ** 2
    #         ZCEN = -2 * SSOUR * SINTHE * COSTHE
    #         fact = -1.0
    #     elif location == "p":
    #         PARAM = 2 * SIMAG * COSTHE ** 2
    #         YCEN = - SIMAG * SINTHE ** 2
    #         ZCEN = -2 * SIMAG * SINTHE * COSTHE
    #         fact = 1.0
    #
    #     if verbose:
    #         txt = ""
    #         if location == "p":
    #             txt += "Source is at infinity\n"
    #             txt += "q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (SIMAG, theta_grazing, theta)
    #         else:
    #             txt += "Image is at infinity\n"
    #             txt += "p=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (SSOUR, theta_grazing, theta)
    #         txt += 'Parabloid of revolution PARAM=%f \n' % PARAM
    #         txt += 'Optical element center at: [0,%f,%f]\n' % (YCEN, ZCEN)
    #         print(txt)
    #
    #     self.ccc[0] = 1.0
    #     self.ccc[1] = COSTHE ** 2
    #     self.ccc[2] = SINTHE ** 2
    #     self.ccc[3] = 0.0
    #     self.ccc[4] = 2 * fact * COSTHE * SINTHE
    #     self.ccc[5] = 0.0
    #     self.ccc[6] = 0.0
    #     self.ccc[7] = 0.0
    #     self.ccc[8] = 2 * ZCEN * SINTHE - 2 * PARAM * COSTHE
    #     self.ccc[9] = 0.0
    #
    # def set_hyperboloid_from_focal_distances(self, SSOUR, SIMAG, theta_grazing, verbose=True):
    #
    #     theta = (numpy.pi / 2) - theta_grazing
    #     COSTHE = numpy.cos(theta)
    #     SINTHE = numpy.sin(theta)
    #
    #     AXMAJ = 0.5 * numpy.abs(SSOUR - SIMAG)
    #     AFOCI = 0.5 * numpy.sqrt(SSOUR ** 2 + SIMAG ** 2 - 2 * SSOUR * SIMAG * numpy.cos(2 * theta_grazing))
    #     AXMIN = numpy.sqrt(AFOCI ** 2 - AXMAJ ** 2)
    #
    #     ECCENT = AFOCI / numpy.abs(AXMAJ)
    #
    #     if SSOUR > SIMAG:
    #         YCEN = (SSOUR ** 2 - SIMAG ** 2) / 4 / AFOCI
    #         ZCEN = AXMIN * numpy.sqrt(YCEN ** 2 / AXMAJ ** 2 - 1.0)
    #
    #         RNCEN = numpy.array((0, -2 * YCEN / AXMAJ ** 2, 2 * ZCEN / AXMIN ** 2))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #         RNCEN /= numpy.sqrt(RNCEN[0] ** 2 + RNCEN[1] ** 2 + RNCEN[2] ** 2)
    #
    #         RTCEN = numpy.array((0, RNCEN[3 - 1], -RNCEN[2 - 1]))
    #     else:
    #         YCEN = (SSOUR ** 2 - SIMAG ** 2) / 4 / AFOCI
    #         ZCEN = AXMIN * numpy.sqrt(YCEN ** 2 / AXMAJ ** 2 - 1.0)
    #
    #         RNCEN = -numpy.array((0, -2 * YCEN / AXMAJ ** 2, 2 * ZCEN / AXMIN ** 2))  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    #         RNCEN /= numpy.sqrt(RNCEN[0] ** 2 + RNCEN[1] ** 2 + RNCEN[2] ** 2)
    #
    #         RTCEN = numpy.array((0, RNCEN[3 - 1], -RNCEN[2 - 1]))
    #
    #     if verbose:
    #         txt = ""
    #         txt += "p=%f, q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (SSOUR, SIMAG, theta_grazing, theta)
    #         txt += 'Hyperboloid of revolution a=%f \n' % AXMAJ
    #         txt += 'Hyperboloid of revolution b=%f \n' % AXMIN
    #         txt += 'Hyperboloid of revolution c=%f \n' % AFOCI
    #         txt += 'Hyperboloid of revolution focal distance c^2=%f \n' % (AFOCI ** 2)
    #         txt += 'Hyperboloid of revolution eccentricity: %f \n' % ECCENT
    #         txt += 'Optical element center at: [0,%f,%f]\n' % (YCEN, ZCEN)
    #         txt += 'Optical element normal: [%f,%f,%f]\n' % (RNCEN[0], RNCEN[1], RNCEN[2])
    #         txt += 'Optical element tangent: [%f,%f,%f]\n' % (RTCEN[0], RTCEN[1], RTCEN[2])
    #         print(txt)
    #
    #     # ;C
    #     # ;C Coefficients of the canonical form
    #     # ;C
    #     A = -1 / AXMIN ** 2
    #     B = 1 / AXMAJ ** 2
    #     C = A
    #     # ;C
    #     # ;C Rotate now in the mirror RF. The equations are the same as for the
    #     # ;C ellipse case.
    #     # ;C
    #     self.ccc[0] = A
    #     self.ccc[1] = B * RTCEN[2 - 1] ** 2 + C * RTCEN[3 - 1] ** 2
    #     self.ccc[2] = B * RNCEN[2 - 1] ** 2 + C * RNCEN[3 - 1] ** 2
    #     self.ccc[3] = 0.0
    #     self.ccc[4] = 2 * (B * RNCEN[2 - 1] * RTCEN[2 - 1] + C * RNCEN[3 - 1] * RTCEN[3 - 1])
    #     self.ccc[5] = 0.0
    #     self.ccc[6] = 0.0
    #     self.ccc[7] = 0.0
    #     self.ccc[8] = 2 * (B * YCEN * RNCEN[2 - 1] + C * ZCEN * RNCEN[3 - 1])
    #     self.ccc[9] = 0.0
    #
    # # todo: delete
    # @classmethod
    # def calculate_ellipsoid_parameters_from_focal_distances(cls, ssour, simag, theta_grazing, verbose=True):
    #     theta = (numpy.pi / 2) - theta_grazing
    #     COSTHE = numpy.cos(theta)
    #     SINTHE = numpy.sin(theta)
    #
    #     AXMAJ = (ssour + simag) / 2
    #     AXMIN = numpy.sqrt(simag * ssour) * COSTHE
    #
    #     AFOCI = numpy.sqrt(AXMAJ ** 2 - AXMIN ** 2)
    #     ECCENT = AFOCI / AXMAJ
    #     # ;C
    #     # ;C The center is computed on the basis of the object and image positions
    #     # ;C
    #     YCEN = (ssour - simag) * 0.5 / ECCENT
    #     ZCEN = -numpy.sqrt(1 - YCEN ** 2 / AXMAJ ** 2) * AXMIN
    #     # ;C
    #     # ;C Computes now the normal in the mirror center.
    #     # ;C
    #     RNCEN = numpy.zeros(3)
    #     RNCEN[1 - 1] = 0.0
    #     RNCEN[2 - 1] = -2 * YCEN / AXMAJ ** 2
    #     RNCEN[3 - 1] = -2 * ZCEN / AXMIN ** 2
    #     # ;CALL NORM(RNCEN,RNCEN)
    #     RNCEN = RNCEN / numpy.sqrt((RNCEN ** 2).sum())
    #     # ;C
    #     # ;C Computes the tangent versor in the mirror center.
    #     # ;C
    #     RTCEN = numpy.zeros(3)
    #     RTCEN[1 - 1] = 0.0
    #     RTCEN[2 - 1] = RNCEN[3 - 1]
    #     RTCEN[3 - 1] = -RNCEN[2 - 1]
    #
    #     # new
    #
    #     ELL_THE = numpy.arctan(ZCEN / YCEN)
    #
    #     YCEN2 = AXMAJ * AXMIN
    #     YCEN2 = YCEN2 / numpy.sqrt(AXMIN ** 2 + AXMAJ ** 2 * numpy.tan(ELL_THE) ** 2)
    #     ZCEN2 = YCEN2 * numpy.tan(ELL_THE)
    #     ZCEN2 = - numpy.abs(ZCEN2)
    #     if (numpy.cos(ELL_THE) < 0):
    #         YCEN2 = - numpy.abs(YCEN2)
    #     else:
    #         YCEN2 = numpy.abs(YCEN2)
    #
    #     if verbose:
    #         print("YCEN2,ZCEN2: ", YCEN2, ZCEN2)
    #         txt = ""
    #         txt += "p=%f, q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (ssour, simag, theta_grazing, theta)
    #         txt += 'Ellipsoid of revolution a=%f \n' % AXMAJ
    #         txt += 'Ellipsoid of revolution b=%f \n' % AXMIN
    #         txt += 'Ellipsoid of revolution c=sqrt(a^2-b^2)=%f \n' % AFOCI
    #         txt += 'Ellipsoid of revolution focal distance c^2=%f \n' % (AFOCI ** 2)
    #         txt += 'Ellipsoid of revolution eccentricity: %f \n' % ECCENT
    #         txt += 'Optical element center at: [0,%f,%f]\n' % (YCEN, ZCEN)
    #         txt += 'Optical element normal: [%f,%f,%f]\n' % (RNCEN[0], RNCEN[1], RNCEN[2])
    #         txt += 'Optical element tangent: [%f,%f,%f]\n' % (RTCEN[0], RTCEN[1], RTCEN[2])
    #         print(txt)
    #     return {
    #         "ssour": ssour, "simag": simag, "theta_grazing": theta_grazing, "theta": theta,
    #         "AXMAJ": AXMAJ, "AXMIN": AXMIN, "ELL_THE": ELL_THE,
    #         "AFOCI": AFOCI, "YCEN": YCEN, "ZCEN": ZCEN, "YCEN2": YCEN2, "ZCEN2": ZCEN2, "RNCEN": RNCEN, "RTCEN": RTCEN}

    # # todo: remove and use shadow4.tools.arrayofvectors.vector_reflection *** CHECK SHAPES BEFORE DOING IT, IT MAY NEED TRANSPOSE ****
    # def vector_reflection(self, v1, normal):
        #
        # # \vec{r} = \vec{i} - 2 (\vec{i} \vec{n}) \vec{n}
        # # \vec{r} = \vec{i} - 2 tmp3
        # tmp = v1 * normal
        # tmp2 = tmp[0, :] + tmp[1, :] + tmp[2, :]
        # tmp3 = normal.copy()
        #
        # for jj in (0, 1, 2):
        #     tmp3[jj, :] = tmp3[jj, :] * tmp2
        #
        # v2 = v1 - 2 * tmp3
        # v2mod = numpy.sqrt(v2[0, :] ** 2 + v2[1, :] ** 2 + v2[2, :] ** 2)
        # v2 /= v2mod
        #
        # return v2



if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()

    if True:
        ccc = S4Conic.initialize_as_plane()
        x2 = numpy.zeros((3,10))
        print("plane: ", ccc.get_normal(x2))

        ccc = S4Conic.initialize_as_sphere_from_external_parameters(1000.0)
        x2 = numpy.zeros((3,10))
        print("R=1000: ", ccc.get_normal(x2))

        ccc = S4Conic.initialize_as_sphere_from_focal_distances(1e15, 1000.0, 1e-3)
        x2 = numpy.zeros((3,10))
        print("R from (p,q): ", ccc.get_normal(x2))

        ccc = S4Conic.initialize_as_sphere_from_external_parameters(-1000.0)
        x2 = numpy.zeros((3,10))
        print("R=-1000: ", ccc.get_normal(x2))
    if True:
        p = 13.73 + 13.599
        q = 2.64
        theta1 = 0.02181
        # ccc = Conic.initialize_as_sphere_from_focal_distances(p, q, theta1, cylindrical=0, cylangle=0.0, switch_convexity=0)
        ccc = S4Conic.initialize_as_ellipsoid_from_focal_distances(p, q, theta1, cylindrical=0, cylangle=0.0, switch_convexity=0)
        # ccc = Conic.initialize_as_paraboloid_from_focal_distances(p, q, theta1, cylindrical=0, cylangle=0.0, switch_convexity=0)
        # ccc = Conic.initialize_as_hyperboloid_from_focal_distances(p, q, theta1, cylindrical=0, cylangle=0.0, switch_convexity=0)
        # print(ccc.info())

        y = numpy.linspace(-0.25, 0.25, 200)
        z = ccc.surface_height(y)
        from srxraylib.plot.gol import plot
        plot(y,z,xtitle="y",ytitle="z")

        #
        #
        #


        x = numpy.linspace(-0.15, 0.15, 100)
        Y = numpy.outer(numpy.ones_like(x),y)
        X = numpy.outer(x,numpy.ones_like(y))
        Z = ccc.surface_height(Y,X)

        from srxraylib.plot.gol import plot_image
        plot_image(Z,x,y)
        print(ccc.info())
        print("Ellipsoid parameters: ")
        ccc2 = S4Conic.initialize_as_ellipsoid_from_focal_distances(p, q, theta1)
        # tkt = S4Conic.calculate_ellipsoid_parameters_from_focal_distances(p, q, theta1)
        #
        # # using external parameters
        # ccc2 = S4Conic()
        # ccc2.set_ellipsoid_from_external_parameters(AXMAJ=tkt["AXMAJ"],AXMIN=tkt["AXMIN"],ELL_THE=tkt["ELL_THE"])
        # for key in tkt.keys():
        #     print(key,tkt[key])

        for i in range(10):
            print(ccc.get_coefficients()[i], ccc2.get_coefficients()[i])


    if True:
        p = 40.0
        q = 10.0
        theta1 = 0.003
        ccc = S4Conic.initialize_as_ellipsoid_from_focal_distances(p, q, theta1, cylindrical=0, cylangle=0.0, switch_convexity=0)

        y = numpy.linspace(-0.2, 0.2, 5000)
        x = numpy.linspace(-0.001, 0.001, 100)

        Y = numpy.outer(numpy.ones_like(x),y)
        X = numpy.outer(x,numpy.ones_like(y))
        Z = ccc.surface_height(X, Y)

        from srxraylib.plot.gol import plot_image
        plot_image(Z, x, y, aspect='auto')

        print(Z.shape, x.shape, y.shape)
        ccc.write_mesh_file(x, y,   filename="/tmp/mirror111.dat")
        ccc.write_mesh_h5file(x, y, filename="/tmp/mirror111.h5")

    if True:
        a = S4Conic.initialize_as_sphere_from_external_parameters(100)
        print(a.info())

    if True:
        ccc = S4Conic.initialize_as_hyperboloid_from_focal_distances(10.0, 3.0, 0.003, cylindrical=0, cylangle=0.0, switch_convexity=0)
        c = ccc.get_coefficients()
        print(c)
        s5 =  [-3703.714814855263, -0.03333333333333342, -3703.5998488688683, 0.0, -41.269717460357114, 0.0, 0.0, 0.0, -190.47647619130194, 0.0]
        s5 = numpy.array(s5)
        s5 = s5 / s5[0]
        for i in range(10):
            assert ( numpy.abs(s5[i] - c[i]) < 1e-3)

        ccc = S4Conic.initialize_as_hyperboloid_from_focal_distances(3, 10.0, 0.003, cylindrical=0, cylangle=0.0, switch_convexity=0)
        c = ccc.get_coefficients()
        print(c)
        s5 =  [-3703.714814855263, -0.033333333332666124, -3703.5998488688692, 0.0, 41.269717460237345, 0.0, 0.0, 3.0796371740537288e-12, 190.47647619130194, 0]
        s5 = numpy.array(s5)
        s5 = s5 / s5[0]
        for i in range(10):
            assert ( numpy.abs(s5[i] - c[i]) < 1e-3)


    # check method 0: Goldberg&Sanchez del Rio, 1: Sanchez del Rio&Goldberg.
    # shape = 2 # 0=parabola, 1=ellipse, 2=hyperbola
    if True:
        ntimes = 100
        for shape in [0,1,2]:
            for i in range(ntimes):
                p = 1000 * numpy.random.rand()
                q = 1000 * numpy.random.rand()
                theta1 = numpy.random.rand()

                # p = 200
                # q = 10
                # theta = 0.001

                if shape == 0:
                    conic = S4Conic.initialize_as_paraboloid_from_focal_distances(p, q, theta1,
                                                                               cylindrical=0, cylangle=0.0,
                                                                               switch_convexity=0,
                                                                               method=0, #0: Goldberg&Sanchez del Rio, 1: Sanchez del Rio&Goldberg.
                                                                               verbose=1)
                    print(conic.get_coefficients())

                    ccc = S4Conic.initialize_as_paraboloid_from_focal_distances(p, q, theta1,
                                                                               cylindrical=0, cylangle=0.0,
                                                                               switch_convexity=0,
                                                                               method=1, #0: Goldberg&Sanchez del Rio, 1: Sanchez del Rio&Goldberg.
                                                                               verbose=1)
                elif shape == 1:
                    conic = S4Conic.initialize_as_ellipsoid_from_focal_distances(p, q, theta1,
                                                                               cylindrical=0, cylangle=0.0, switch_convexity=0,
                                                                               method=0, verbose=1)
                    print(conic.get_coefficients())

                    ccc = S4Conic.initialize_as_ellipsoid_from_focal_distances(p, q, theta1,
                                                                               cylindrical=0, cylangle=0.0, switch_convexity=0,
                                                                               method=1, verbose=1)
                elif shape == 2:
                    conic = S4Conic.initialize_as_hyperboloid_from_focal_distances(p, q, theta1,
                                                                               cylindrical=0, cylangle=0.0, switch_convexity=0,
                                                                               method=0, verbose=1)
                    print(conic.get_coefficients())

                    ccc = S4Conic.initialize_as_hyperboloid_from_focal_distances(p, q, theta1,
                                                                               cylindrical=0, cylangle=0.0, switch_convexity=0,
                                                                               method=1, verbose=1)


                print(ccc.get_coefficients() * conic.get_coefficients()[0])

                c_p_1 = conic.get_coefficients()
                c_p_2 = ccc.get_coefficients() * c_p_1[0]
                for i in range(10):
                    diff = numpy.abs(c_p_1[i] - c_p_2[i])
                    if diff < 1e-5:
                        ss = ''
                    else:
                        ss = '<<<<<  PROBLEM >>>>>'
                    print(i, c_p_1[i], c_p_2[i], ss)
                    # assert (diff < 1e-5)


                print(conic.info())
