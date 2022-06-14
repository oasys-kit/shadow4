# operations with an array of 3D vectors, v[n_vectors,3]

import numpy

def vector_cross(u ,v):
    # w = u X v
    # u = array (npoints,vector_index)

    w = numpy.zeros_like(u)
    w[: ,0] = u[: ,1] * v[: ,2] - u[: ,2] * v[: ,1]
    w[: ,1] = u[: ,2] * v[: ,0] - u[: ,0] * v[: ,2]
    w[: ,2] = u[: ,0] * v[: ,1] - u[: ,1] * v[: ,0]

    return w

def vector_modulus(u):
    return numpy.sqrt( u[: ,0 ]**2 + u[: ,1 ]**2 + u[: ,2 ]**2)

def vector_modulus_square(u):
    return ( u[: ,0 ]**2 + u[: ,1 ]**2 + u[: ,2 ]**2)

def vector_norm(u):
    # w = u / |u|
    u_norm = numpy.zeros_like(u)
    uu = numpy.sqrt( u[: ,0 ]**2 + u[: ,1 ]**2 + u[: ,2 ]**2)
    for i in range(3):
        u_norm[: ,i] = uu
    return u / u_norm

def vector_dot(u, v):
    # w = u . v
    return u[: ,0] * v[: ,0] + u[: ,1] * v[: ,1] + u[: ,2] * v[: ,2]


def vector_sum(u, v):
    # w = u + v
    w = numpy.zeros_like(u)
    for i in range(3):
        w[: ,i] = u[: ,i] + v[: ,i]
    return w

def vector_multiply_scalar(u, k):
    kk = numpy.array(k)

    if kk.size == 1:
        return u * kk
    else:
        w = numpy.zeros_like(u)
        for i in range(3):
            w[: ,i] = u[: ,i] * kk
        return w

def vector_add_scalar(u, k):
    kk = numpy.array(k)

    if kk.size == 1:
        return u + kk
    else:
        w = numpy.zeros_like(u)
        for i in range(3):
            w[: ,i] = u[: ,i] + kk
        return w


def vector_diff(u, v):
    # w = u - v
    w = numpy.zeros_like(u)
    for i in range(3):
        w[: ,i] = u[: ,i] - v[: ,i]
    return w

def vector_reflection(v1,normal): # copied from s4_conic()
    # \vec{r} = \vec{i} - 2 (\vec{i} \vec{n}) \vec{n}
    normal_norm = vector_norm(normal)
    return v1 - 2 * vector_multiply_scalar( normal_norm, vector_dot(v1, normal_norm))

def vector_refraction(vin, normal, n1, n2):
    # http://www.starkeffects.com/snells-law-vector.shtml

    vin_norm = vector_norm(vin)
    normal_norm = vector_norm(normal)

    n_cross_vin = vector_cross(normal_norm, vin_norm)
    n_opp_cross_vin = vector_cross(normal_norm * (-1), vin_norm)

    sq2 = 1 - vector_multiply_scalar(vector_dot(n_cross_vin, n_cross_vin), (n1/n2)**2)
    vout = (n1/n2) * vector_cross(normal_norm,n_opp_cross_vin) - vector_multiply_scalar(normal_norm, numpy.sqrt(sq2))
    return vout

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
    npoints = 1
    n = numpy.zeros((npoints,3))
    n[:,2] = 1.0

    v1 = numpy.zeros((npoints, 3))
    v1[:, 0] = numpy.sqrt(2) / 2
    v1[:, 2] = -numpy.sqrt(2) / 2
    v2 = vector_refraction(v1, n, n1=1.0, n2=1.5)
    print('v1: ', v1)
    print('n: ', n)
    print('v2: ', v2)