#
"""

Set of functions to make operations with an array of 3D vectors, v[n_vectors,3]

"""
import numpy

def vector_cross(u ,v):
    """

    Calculate the vector cross product.

    Parameters
    ----------
    u : numpy array shape (n_vectors,3)
        input vector 1.
    v : numpy array shape (n_vectors,3)
        input vector 2.

    Returns
    -------
     numpy array shape (n_vectors,3)
        Result vector.

    """
    # w = u X v
    # u = array (npoints,vector_index)

    w = numpy.zeros_like(u)
    w[: ,0] = u[: ,1] * v[: ,2] - u[: ,2] * v[: ,1]
    w[: ,1] = u[: ,2] * v[: ,0] - u[: ,0] * v[: ,2]
    w[: ,2] = u[: ,0] * v[: ,1] - u[: ,1] * v[: ,0]

    return w

def vector_modulus(u):
    """

    Calculates the modulus of a vector array.

    Parameters
    ----------
    u : numpy array shape (n_vectors,3)
        input vector 1.

    Returns
    -------
    numpy array shape (n_vectors)
        Result vector modulus.

    """
    return numpy.sqrt( u[: ,0 ]**2 + u[: ,1 ]**2 + u[: ,2 ]**2)

def vector_modulus_square(u):
    return ( u[: ,0 ]**2 + u[: ,1 ]**2 + u[: ,2 ]**2)

def vector_norm(u):
    """

    Calculate the normalized vector.

    Parameters
    ----------
    u : numpy array shape (n_vectors,3)
        input vector.

    Returns
    -------
     numpy array shape (n_vectors,3)
        Result vector.

    """
    # w = u / |u|
    u_norm = numpy.zeros_like(u)
    uu = numpy.sqrt( u[: ,0 ]**2 + u[: ,1 ]**2 + u[: ,2 ]**2)
    for i in range(3):
        u_norm[: ,i] = uu
    return u / u_norm

def vector_dot(u, v):
    """

    Calculate the dot product of two vectors.

    Parameters
    ----------
    u : numpy array shape (n_vectors,3)
        input vector 1.
    v : numpy array shape (n_vectors,3)
        input vector 2.

    Returns
    -------
     numpy array shape (n_vectors)
        Result dot product vector.

    """
    # w = u . v
    return u[: ,0] * v[: ,0] + u[: ,1] * v[: ,1] + u[: ,2] * v[: ,2]


def vector_sum(u, v):
    """

    Calculate the sum of two vectors.

    Parameters
    ----------
    u : numpy array shape (n_vectors,3)
        input vector 1.
    v : numpy array shape (n_vectors,3)
        input vector 2.

    Returns
    -------
     numpy array shape (n_vectors,3)
        Result vector.

    """
    # w = u + v
    w = numpy.zeros_like(u)
    for i in range(3):
        w[: ,i] = u[: ,i] + v[: ,i]
    return w

def vector_multiply_scalar(u, k):
    """

    Calculate the product of a vector by a scalar.

    Parameters
    ----------
    u : numpy array shape (n_vectors,3)
        input vector 1.
    k : numpy array shape (n_vectors)
        scalar values to multiply by.

    Returns
    -------
     numpy array shape (n_vectors)
        Result vector.

    """

    kk = numpy.array(k)

    if kk.size == 1:
        return u * kk
    else:
        w = numpy.zeros_like(u)
        for i in range(3):
            w[: ,i] = u[: ,i] * kk
        return w

def vector_add_scalar(u, k):
    """

    Calculate the sum of a vector and a scalar constant.

    Parameters
    ----------
    u : numpy array shape (n_vectors,3)
        input vector 1.
    k : numpy array shape (n_vectors)
        scalar values to be added.

    Returns
    -------
     numpy array shape (n_vectors)
        Result vector.

    """
    kk = numpy.array(k)

    if kk.size == 1:
        return u + kk
    else:
        w = numpy.zeros_like(u)
        for i in range(3):
            w[: ,i] = u[: ,i] + kk
        return w


def vector_diff(u, v):
    """

    Calculate the difference of two vectors u - v.

    Parameters
    ----------
    u : numpy array shape (n_vectors,3)
        input vector 1.
    v : numpy array shape (n_vectors,3)
        input vector 2.

    Returns
    -------
     numpy array shape (n_vectors,3)
        Result vector.

    """
    # w = u - v
    w = numpy.zeros_like(u)
    for i in range(3):
        w[: ,i] = u[: ,i] - v[: ,i]
    return w

def vector_reflection(v1, normal): # copied from s4_conic()
    """

    Calculate the reflection of a vector (ray) on a surface.

    Parameters
    ----------
    v1 : array of 3D vectors, shape: (n_vectors,3)
        Incident unit vector.
    normal : array of 3D vectors, shape: (n_vectors,3)
        Normal unit vector at the interface.

    Returns
    -------
    numpy array
        The unitary vector along the reflected direction with shape (n_vectors,3).
    """

    # \vec{r} = \vec{i} - 2 (\vec{i} \vec{n}) \vec{n}
    normal_norm = vector_norm(normal)
    return v1 - 2 * vector_multiply_scalar( normal_norm, vector_dot(v1, normal_norm))

def vector_refraction(vin, normal, n1, n2, sgn=1, do_check=0):
    """

    Calculate the refraction (transmission) vector using Snell's Law in vector form.

    Parameters
    ----------
    vin : array of 3D vectors, shape: (n_vectors,3)
        Incident unit vector.
    normal : array of 3D vectors, shape: (n_vectors,3)
        Normal unit vector at the interface.
    n1 : numpy array
        The refraction index of the incident mediuma (n_vectors)
    n2 : numpy array
        The refraction index of the object transmission (n_vectors).
    sgn : numpy array
        +1 or -1 to play and adjust the directions.
    do_check : int
        A flag to display incident and refracted angles.

    Returns
    -------
    numpy array
        The unitary vector along the refracted direction with shape (n_vectors,3).

    Reference
    ---------
    For a detailed explanation of the vector form of Snell's Law, see:
    https://www.starkeffects.com/snells-law-vector.shtml
    """

    if sgn is None: sgn = -numpy.sign(vector_dot(vin, normal))

    vin_norm = vector_norm(vin)
    normal_norm = vector_norm(normal)

    n_cross_vin = vector_cross(normal_norm, vin_norm)
    n_opp_cross_vin = vector_cross(normal_norm * (-1), vin_norm)

    sq2 = 1 - vector_dot(n_cross_vin, n_cross_vin) *  (n1 / n2) ** 2
    # vout = (n1/n2) * vector_cross(normal_norm,n_opp_cross_vin) - vector_multiply_scalar(normal_norm, numpy.sqrt(sq2))
    vout = vector_multiply_scalar(vector_cross(normal_norm, n_opp_cross_vin), (n1 / n2)) - vector_multiply_scalar(normal_norm, sgn * numpy.sqrt(sq2))

    if do_check:
        theta1 = numpy.arccos( vector_dot(vin_norm, normal_norm) * (-1))
        theta2 = numpy.arccos(vector_dot(vout, vector_multiply_scalar(normal_norm, -1)))
        print(">>>>> theta1: ", numpy.degrees(theta1))
        print(">>>>> theta2: ", numpy.degrees(theta2))
        print(">>>>> theta2 check: ", numpy.degrees(numpy.arcsin(n1/n2*numpy.sin(theta1))))
    return vout

def vector_scattering(K_IN, H, NORMAL):
    """

    Calculate the scattering of a wavevector.

    Method:

        - K_OUT_PARALLEL = K_IN_PARALLEL + H_PARALLEL

        - |K_OUT| = |K_IN|

    Parameters
    ----------
    K_IN  : numpy array of 3D vectors, shape: (n_vectors,3)
        Incident wavevector (modulus is 2 pi / wavelength) in cm**-1.
    H : numpy array of 3D vectors, shape: (n_vectors,3)
        The diffraction vector.
    NORMAL : numpy array of 3D vectors, shape: (n_vectors,3)
        The normal vector.

    Returns
    -------
    numpy array
        The diffracted wavevector with shape (n_vectors,3).
    """
    H_perp = vector_multiply_scalar(NORMAL, vector_dot(H, NORMAL))
    H_par = vector_diff(H, H_perp)

    K_IN_perp = vector_multiply_scalar(NORMAL, vector_dot(K_IN, NORMAL))
    K_IN_par = vector_diff(K_IN, K_IN_perp)

    K_OUT_par = vector_sum(K_IN_par, H_par)
    K_OUT_perp = vector_multiply_scalar(NORMAL, numpy.sqrt(vector_modulus_square(K_IN) - vector_modulus_square(K_OUT_par)))
    K_OUT = vector_sum(K_OUT_par, K_OUT_perp)
    return K_OUT

def vector_rotate_around_axis(u, rotation_axis, angle):
    """Rotates the vector around an axis. It uses the Rodrigues formula [rf]_

    Parameters
    ----------
    rotation_axis : numpy array of 3D vectors, shape: (n_vectors,3)
        Vector specifying the rotation axis (not necessarily unit vector).

    angle : float
        Rotation angle in radiants.

    Returns
    -------
    Vector instance
        Rotated vector as a new vector.

    References
    ----------
    .. [rf] http://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula

    """
    if isinstance(rotation_axis, numpy.ndarray):
        rotation_axis1 = rotation_axis
    else:
        rotation_axis1 = numpy.array(rotation_axis)

    if rotation_axis1.shape != u.shape:
        if len(rotation_axis1.shape) == 1: # e.g. rotation_axis=[1,0,0]
            rotation_axis2 = numpy.zeros_like(u)
            rotation_axis2[:, 0] = rotation_axis1[0]
            rotation_axis2[:, 1] = rotation_axis1[1]
            rotation_axis2[:, 2] = rotation_axis1[2]
            rotation_axis1 = rotation_axis2

    unit_rotation_axis = vector_norm(rotation_axis1)

    rotated_vector = vector_multiply_scalar(u, numpy.cos(angle))

    tmp_vector = vector_cross(unit_rotation_axis, u)
    tmp_vector = vector_multiply_scalar(tmp_vector, numpy.sin(angle))
    rotated_vector = vector_sum(rotated_vector, tmp_vector)

    scalar_factor = vector_dot(u, unit_rotation_axis) * (1.0 - numpy.cos(angle))
    tmp_vector = vector_multiply_scalar(unit_rotation_axis, scalar_factor)
    rotated_vector = vector_sum(rotated_vector, tmp_vector)

    return rotated_vector

def vector_default_efields(DIREC, pol_deg=1.0):
    """

    Creates two unitary vectors pointing on a sigma and pi-polarization directions. They are perpendicular to the
    incident direction, and the sum of the squarred moduli is one.

    Parameters
    ----------
    DIREC  : numpy array of 3D vectors, shape: (n_vectors,3)
        The vector with the direction.
    pol_deg : numpy array shape: (n_vectors)
        The polarization degree as defined in SHADOW.

    Returns
    -------
    tuple
        (E_S, E_P): the 3D vectors, shape: (n_vectors,3) for the S and P polarizations.
    """
    # Generates the normalized electric vectors perpendicular to DIREC
    # This is defined on the source plane, so that A_VEC is along the X-axis and AP_VEC is along Z-axis.
    # Then care must be taken so that A will be perpendicular to the ray direction.

    A_VEC = numpy.zeros_like(DIREC)
    A_VEC[:, 0] = 1.0

    # ! C   Rotate A_VEC so that it will be perpendicular to DIREC and with the
    # ! C   right components on the plane.
    A_TEMP = vector_cross(A_VEC, DIREC)
    A_VEC = vector_cross(DIREC, A_TEMP)
    A_VEC = vector_norm(A_VEC)
    AP_VEC = vector_cross(A_VEC, DIREC)
    AP_VEC = vector_norm(AP_VEC)

    #
    # obtain polarization for each ray (interpolation)
    #
    DENOM = numpy.sqrt(1.0 - 2.0 * pol_deg + 2.0 * pol_deg ** 2)
    AX = pol_deg / DENOM
    for i in range(3):
        A_VEC[:, i] *= AX

    AZ = (1.0 - pol_deg) / DENOM
    for i in range(3):
        AP_VEC[:, i] *= AZ

    return A_VEC, AP_VEC

if __name__ == "__main__":

    v1 = numpy.ones((200,3))
    assert (vector_dot(v1,v1)[0] == 3)

    v3 = numpy.ones((200,3)) * 3
    assert (vector_sum(v1, v3)[0,0] == 4)

    assert (vector_diff(v1, v3)[0,1] == -2)

    # reflection
    npoints = 1
    n = numpy.zeros((npoints,3))
    n[:,2] = 2.0

    v1 = numpy.zeros((npoints, 3))
    v1[:, 0] = numpy.sqrt(2) / 2
    v1[:, 2] = -numpy.sqrt(2) / 2
    v2 = vector_reflection(v1, n)
    print('v1: ', v1)
    print('n: ', n)
    print('v2: ', v2)

    # refraction
    # use sgn to make refraction independent of the direction of the normal
    npoints = 1
    n = numpy.zeros((npoints,3))
    n[:,2] = -1.0

    v1 = numpy.zeros((npoints, 3))
    v1[:, 0] = numpy.sqrt(2) / 2
    v1[:, 2] = numpy.sqrt(2) / 2
    v2 = vector_refraction(v1, n, n1=1.0, n2=1.5, sgn=1)
    print('\nv1: ', v1)
    print('n: ', n)
    print('v2: ', v2)
    print('sgn: ', vector_dot(v1, n))

    # http://www.starkeffects.com/snells-law-vector.shtml
    assert (numpy.abs(v2[0][0] - 0.471) < 1e-3)
    assert (numpy.abs(v2[0][1] - 0) < 1e-3)
    assert (numpy.abs(v2[0][2] - 0.882) < 1e-3)

    v2 = vector_refraction(v1, n, n1=1.0, n2=1.5, sgn=None)
    print('\nv1: ', v1)
    print('n: ', n)
    print('v2: ', v2)
    print('sgn: ', vector_dot(v1, n))


    n *= -1
    v2 = vector_refraction(v1, n, n1=1.0, n2=1.5, sgn=None)  # automatic
    print('\nv1: ', v1)
    print('n: ', n)
    print('v2: ', v2)
    print('sgn: ', vector_dot(v1, n))




    u = numpy.zeros((11,3))
    u[:, 0] = numpy.linspace(1, 2, 11) * 0
    u[:, 1] = numpy.linspace(2, 3, 11) * 0
    u[:, 2] = numpy.linspace(3, 4, 11) * 0 + 1
    print(vector_rotate_around_axis(u, [1,0,0], -numpy.radians(10)))


    axis = numpy.zeros((11, 3))
    axis[:,0] = 1
    print(axis.shape)
    print(vector_rotate_around_axis(u, axis, numpy.radians(10) ))