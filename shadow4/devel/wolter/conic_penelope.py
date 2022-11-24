#
# implementation of generic equations for quadrics
# seee chapter 6 in https://www.oecd-nea.org/upload/docs/application/pdf/2020-10/penelope-2018__a_code_system_for_monte_carlo_simulation_of_electron_and_photon_transport.pdf
#

import numpy

from shadow4.optical_surfaces.s4_conic import S4Conic # for comparison
from shadow4.devel.wolter.conic_viewer import view_conic, compare_conics # for plot

#
# some tools
#


def reduced_quadric(kind='plane'):
    # I1 x^2 + I2 y^2 + I3 z^2 + I4 z + I5 = 0  (Eq 6.23 in Penelope Manual)

    if kind == 'plane':  # (Ibid, table 6.1
        return [0,0,0,1,-1]
    elif kind == 'pair of parallel planes':
        return [0,0,1,0,-1]
    elif kind == 'sphere':
        return [1,1,1,0,-1]
    elif kind == 'cylinder':
        return [1,1,0,0,-1]
    elif kind == 'hyperbolic cylinder':
        return [1,-1,0,0,-1]
    elif kind == 'hyperbolic cylinder 2':
        return [-1,1,0,0,-1]
    elif kind == 'cone':
        return [1,1,-1,0,0]
    elif kind == 'one sheet hyperboloid':
        return [1,1,-1,0,-1]
    elif kind == 'two sheet hyperboloid':
        return [1,1,-1,0,1]
    elif kind == 'paraboloid':
        return [1,1,0,-1,0]
    elif kind == 'parabolic cylinder':
        return [1,0,0,-1,0]
    elif kind == 'parabolic cylinder 2':
        return [0,1,0,-1,0]
    elif kind == 'hyperbolic paraboloid':
        return [1,-1,0,-1,0]
    elif kind == 'hyperbolic paraboloid 2':
        return [-1,1,0,-1,0]
    else:
        raise Exception('Invalid reduced-quadric name.')


def scale_reduced_quadric(reduced_quadric, xscale=1.0, yscale=1.0, zscale=1.0, return_list=True):
    # see ibid, eq 6.24
    out = reduced_quadric.copy()
    out[0] /= xscale**2
    out[1] /= yscale**2
    out[2] /= zscale**2
    out[3] /= zscale
    out[4] = float(out[4])

    if return_list:
        return out
    else: # return 3 matrices (ibid, eq 6.26)
        A33 = numpy.array( [[out[0],0,0],
                            [0,out[1],0],
                            [0,0,out[2]]])
        A03 = numpy.array( [[0,0,out[3]]])
        A00 = out[4]
        return A33,A03,A00

def expand_reduced_quadric(reduced_quadric_list):
    # returns the 10 coefficients of the generic quartic from the 5 coeffs of the reduced quartic

    QXX = reduced_quadric_list[0]
    QYY = reduced_quadric_list[1]
    QZZ = reduced_quadric_list[2]
    QXY = 0.0
    QYZ = 0.0
    QXZ = 0.0
    QX = 0.0
    QY = 0.0
    QZ = reduced_quadric_list[3]
    Q0 = reduced_quadric_list[4]

    return [QXX,QYY,QZZ,QXY,QYZ,QXZ,QX,QY,QZ,Q0]



def quartic_coefficients_matrices_to_list(A2,A1,A0, fix_zeros=False):
    AXX = A2[1-1,1-1]
    AXY = A2[1-1,2-1] + A2[2-1,1-1]
    AXZ = A2[1-1,3-1] + A2[3-1,1-1]
    AYY = A2[2-1,2-1]
    AYZ = A2[2-1,3-1] + A2[3-1,2-1]
    AZZ = A2[3-1,3-1]
    AX = A1[1-1]
    AY = A1[2-1]
    AZ = A1[3-1]
    out_list = [AXX, AYY, AZZ, AXY, AYZ, AXZ, AX, AY, AZ, A0]

    if fix_zeros:
        for i in range(len(out_list)):
            if numpy.abs(out_list[i]) < 1e-15:
                out_list[i] = 0.0

    return out_list

def quadric_coefficients_list_to_matrices(quartic_coefficients_list):

    AXX, AYY, AZZ, AXY, AYZ, AXZ, AX, AY, AZ, A0 = quartic_coefficients_list

    B2 = numpy.zeros((3, 3))
    # ibid eq. 6.18
    B2[1 - 1, 1 - 1] = AXX
    B2[1 - 1, 2 - 1] = 0.5 * AXY
    B2[1 - 1, 3 - 1] = 0.5 * AXZ
    B2[2 - 1, 1 - 1] = B2[1 - 1, 2 - 1]
    B2[2 - 1, 2 - 1] = AYY
    B2[2 - 1, 3 - 1] = 0.5 * AYZ
    B2[3 - 1, 1 - 1] = B2[1 - 1, 3 - 1]
    B2[3 - 1, 2 - 1] = B2[2 - 1, 3 - 1]
    B2[3 - 1, 3 - 1] = AZZ

    B1 = numpy.zeros((3))
    B1[1 - 1] = AX
    B1[2 - 1] = AY
    B1[3 - 1] = AZ

    B0 = A0

    return B2,B1,B0


def euler_rotation_matrix(omega,theta,phi, shortcut=False, fix_zeros=True):
    STHETA = numpy.sin(theta)
    CTHETA = numpy.cos(theta)
    SPHI = numpy.sin(phi)
    CPHI = numpy.cos(phi)
    SOMEGA = numpy.sin(omega)
    COMEGA = numpy.cos(omega)

    # ibidm eq, 6.9
    R = numpy.zeros((3,3))

    if shortcut:
        R[1-1,1-1] = 1
        R[1-1,2-1] = 0
        R[1-1,3-1] = 0
        R[2-1,1-1] = 0
        R[2-1,2-1] = CTHETA
        R[2-1,3-1] = -STHETA
        R[3-1,1-1] = 0
        R[3-1,2-1] = STHETA
        R[3-1,3-1] = CTHETA
    else:
        R[1-1,1-1] = CPHI*CTHETA*COMEGA-SPHI*SOMEGA
        R[1-1,2-1] = -CPHI*CTHETA*SOMEGA-SPHI*COMEGA
        R[1-1,3-1] = CPHI*STHETA
        R[2-1,1-1] = SPHI*CTHETA*COMEGA+CPHI*SOMEGA
        R[2-1,2-1] = -SPHI*CTHETA*SOMEGA+CPHI*COMEGA
        R[2-1,3-1] = SPHI*STHETA
        R[3-1,1-1] = -STHETA*COMEGA
        R[3-1,2-1] = STHETA*SOMEGA
        R[3-1,3-1] = CTHETA

    if fix_zeros:
        for i in range(3):
            for j in range(3):
                if numpy.abs(R[i,j]) < 1e-15:
                    R[i,j] = 0.0

    return R


#
# main routine for transforming the 10 conic coefficients
#


# translated from penelope cortran code
def rotate_and_shift_quartic(quartic_coefficients_list,
                             omega=0.0, theta=0.0, phi=0.0,
                             D=[0,0,0]):
#
# initial quartic
#
    B2 = numpy.zeros((3,3))

    AXX,AYY,AZZ,AXY,AYZ,AXZ,AX,AY,AZ,A0 = quartic_coefficients_list

    # ibid eq. 6.18
    B2[1-1,1-1] = AXX
    B2[1-1,2-1] = 0.5 * AXY
    B2[1-1,3-1] = 0.5 * AXZ
    B2[2-1,1-1] = B2[1-1,2-1]
    B2[2-1,2-1] = AYY
    B2[2-1,3-1] = 0.5 * AYZ
    B2[3-1,1-1] = B2[1-1,3-1]
    B2[3-1,2-1] = B2[2-1,3-1]
    B2[3-1,3-1] = AZZ

    B1 = numpy.zeros((3))
    B1[1-1] = AX
    B1[2-1] = AY
    B1[3-1] = AZ

    B0 = A0

    # so far, input [or scaled] F(x,y,z) can be written in matrix form (ibid, eq. 6.25)
    # notation: X=(x,y,z) ; B1=(AX,AY,AZ)
    # F = X' . (B2 X) + B1 . X + B0

    D1 = numpy.array(D)

    #
    #  ****  Rotation matrix.
    #
    R = euler_rotation_matrix(omega, theta, phi)

    A1 = numpy.zeros(3)
    A2 = numpy.zeros((3,3))
    for I in range(3):
        A1[I] = 0.0
        for J in range(3):
            A1[I] = A1[I] + R[I,J] * B1[J]
            A2[I,J] = 0.0
            for M in range(3):
                for K in range(3):
                    A2[I,J] = A2[I,J] + R[I,K] * B2[K,M] * R[J,M]


    #
    #  ****  Shifted-rotated quadric.
    #

    # 3rd equation 6.29 and complete the second equation
    for I in range(3):
        A2D = 0.0
        for J in range(3):
            A2D = A2D + A2[I,J] * D1[J]

        B1[I] = A1[I] - 2.0 * A2D
        B0 = B0 + D1[I] * (A2D - A1[I])


    AXX = A2[1-1,1-1]
    AXY = A2[1-1,2-1] + A2[2-1,1-1]
    AXZ = A2[1-1,3-1] + A2[3-1,1-1]
    AYY = A2[2-1,2-1]
    AYZ = A2[2-1,3-1] + A2[3-1,2-1]
    AZZ = A2[3-1,3-1]
    AX = B1[1-1]
    AY = B1[2-1]
    AZ = B1[3-1]
    A0 = B0

    transformed_coefficients = [AXX,AYY,AZZ,AXY,AYZ,AXZ,AX,AY,AZ,A0]
    for i in range(len(transformed_coefficients)):
        if numpy.abs(transformed_coefficients[i]) < 1e-15:
            transformed_coefficients[i] = 0.0

    return transformed_coefficients


# full matrix numpy implementation
def rotate_and_shift_quartic_NEW(quartic_coefficients_list,
                             omega=0.0, theta=0.0, phi=0.0,
                             D=[0.0,0.0,0.0]):
    #
    # initial quartic in matrix format
    #

    # so far, input [or scaled] F(x,y,z) can be written in matrix form (ibid, eq. 6.25)
    # notation: X=(x,y,z) ; B1=(AX,AY,AZ)
    # F = X' . (B2 X) + B1 . X + B0
    B2, B1, B0 = quadric_coefficients_list_to_matrices(quartic_coefficients_list)

    #
    #  ****  Rotation matrix.
    #
    R = euler_rotation_matrix(omega, theta, phi)

    #
    #  ****  Shifted-rotated quadric.
    #

    D1 = numpy.array(D)

    # first equation 6.29
    A2 = numpy.dot(R, numpy.dot(B2,R.T))

    # 2nd equation 6.29
    A1_first_term = numpy.dot(R, B1)
    A1 = A1_first_term - 2 * numpy.dot(A2,D1)

    # 3rd equation 6.29
    A0 = B0 + numpy.dot(D1.T,
                    (numpy.dot(A2, D1) - numpy.dot(R, B1)))

    transformed_coefficients = quartic_coefficients_matrices_to_list(A2,A1,A0,fix_zeros=True)
    return transformed_coefficients

def rotate_and_shift_quartic_MATHEMATICA(quartic_coefficients_list,
                             omega=0.0, theta=0.0, phi=0.0,
                             D=[0.0,0.0,0.0]):

    from numpy import sin as Sin
    from numpy import cos as Cos

    cxx, cyy, czz, cxy, cyz, cxz, cx, cy, cz, c0 = quartic_coefficients_list

    Amat = numpy.zeros((3,3))



    Amat[0,0] = cxx
    Amat[0,1] = (cxy*Cos(theta) - cxz*Sin(theta))/2.
    Amat[0,2] = (cxz*Cos(theta) + cxy*Sin(theta))/2.
    Amat[1,0] = (cxy*Cos(theta) - cxz*Sin(theta))/2.
    Amat[1,1] = cyy*Cos(theta)**2 - cyz*Cos(theta)*Sin(theta) + czz*Sin(theta)**2
    Amat[1,2] = (cyz*Cos(2*theta) + (cyy - czz)*Sin(2*theta))/2.
    Amat[2,0] = (cxz*Cos(theta) + cxy*Sin(theta))/2.
    Amat[2,1] = (cyz*Cos(2*theta) + (cyy - czz)*Sin(2*theta))/2.
    Amat[2,2] = czz*Cos(theta)**2 + cyz*Cos(theta)*Sin(theta) + cyy*Sin(theta)**2


    tx = D[0]
    ty = D[1]
    tz = D[2]

    Avec = numpy.zeros(3)


    Avec[0] = cx - 2*cxx*tx - (cxy*ty + cxz*tz)*Cos(theta) + (cxz*ty - cxy*tz)*Sin(theta)
    Avec[1] = -((cyy + czz)*ty) + (cy - cxy*tx)*Cos(theta) - (cyy*ty - czz*ty + cyz*tz)*Cos(2*theta) + (-cz + cxz*tx)*Sin(theta) + (cyz*ty - cyy*tz + czz*tz)*Sin(2*theta)
    Avec[2] = -((cyz*ty + 2*czz*tz)*Cos(theta)**2) + Sin(theta)*(cy - cxy*tx + (cyz*ty - 2*cyy*tz)*Sin(theta)) + Cos(theta)*(cz - cxz*tx - 2*(cyy*ty - czz*ty + cyz*tz)*Sin(theta))

    A0 = c0 + \
        tx*(-cx + cxx*tx + tz*((cxz*Cos(theta))/2. + (cxy*Sin(theta))/2.) + ty*((cxy*Cos(theta))/2. - (cxz*Sin(theta))/2.)) + \
        tz*(-(cz*Cos(theta)) - cy*Sin(theta) + tx*((cxz*Cos(theta))/2. + (cxy*Sin(theta))/2.) + \
        tz*(Sin(theta)*((cyz*Cos(theta))/2. + cyy*Sin(theta)) + Cos(theta)*(czz*Cos(theta) + (cyz*Sin(theta))/2.)) + \
        ty*(Sin(theta)*(cyy*Cos(theta) - (cyz*Sin(theta))/2.) + Cos(theta)*((cyz*Cos(theta))/2. - czz*Sin(theta)))) + \
        ty*(-(cy*Cos(theta)) + cz*Sin(theta) + tx*((cxy*Cos(theta))/2. - (cxz*Sin(theta))/2.) + \
        tz*(Cos(theta)*((cyz*Cos(theta))/2. + cyy*Sin(theta)) - Sin(theta)*(czz*Cos(theta) + (cyz*Sin(theta))/2.)) + \
        ty*(Cos(theta)*(cyy*Cos(theta) - (cyz*Sin(theta))/2.) - Sin(theta)*((cyz*Cos(theta))/2. - czz*Sin(theta))))

    print(">>>", Amat)
    print(">>>", Avec)
    print(">>>", A0)
    Alist = quartic_coefficients_matrices_to_list(Amat,Avec,A0, fix_zeros=True)
    print(">>>", Alist)
    return Alist

# full implementation from Mathematica code
def rotate_and_shift_quartic_MATHEMATICAFULLEULER(quartic_coefficients_list,
                             omega=0.0, theta=0.0, phi=0.0,
                             D=[0.0,0.0,0.0]):

    from numpy import sin as Sin
    from numpy import cos as Cos

    cxx, cyy, czz, cxy, cyz, cxz, cx, cy, cz, c0 = quartic_coefficients_list

    Amat = numpy.zeros((3,3))

    # this is for theta-rotation around X
    # Amat[0,0] = cxx
    # Amat[0,1] = (cxy*Cos(theta) - cxz*Sin(theta))/2.
    # Amat[0,2] = (cxz*Cos(theta) + cxy*Sin(theta))/2.
    # Amat[1,0] = (cxy*Cos(theta) - cxz*Sin(theta))/2.
    # Amat[1,1] = cyy*Cos(theta)**2 - cyz*Cos(theta)*Sin(theta) + czz*Sin(theta)**2
    # Amat[1,2] = (cyz*Cos(2*theta) + (cyy - czz)*Sin(2*theta))/2.
    # Amat[2,0] = (cxz*Cos(theta) + cxy*Sin(theta))/2.
    # Amat[2,1] = (cyz*Cos(2*theta) + (cyy - czz)*Sin(2*theta))/2.
    # Amat[2,2] = czz*Cos(theta)**2 + cyz*Cos(theta)*Sin(theta) + cyy*Sin(theta)**2

    # this includes the three euler angles
    Amat[0,0] = ((Cos(omega)*Cos(phi)*Cos(theta) - Sin(omega)*Sin(phi))*
        (Cos(phi)*Cos(theta)*(2*cxx*Cos(omega) - cxy*Sin(omega)) - (cxy*Cos(omega) + 2*cxx*Sin(omega))*Sin(phi) + cxz*Cos(phi)*Sin(theta)) +
        (Cos(phi)*Cos(theta)*Sin(omega) + Cos(omega)*Sin(phi))*(Cos(phi)*Cos(theta)*(-(cxy*Cos(omega)) + 2*cyy*Sin(omega)) + (2*cyy*Cos(omega) + cxy*Sin(omega))*Sin(phi) -
        cyz*Cos(phi)*Sin(theta)) + Cos(phi)*Sin(theta)*(Cos(phi)*Cos(theta)*(cxz*Cos(omega) - cyz*Sin(omega)) - (cyz*Cos(omega) + cxz*Sin(omega))*Sin(phi) + 2*czz*Cos(phi)*Sin(theta))
        )/2.

    Amat[0,1] = (Cos(2*omega)*(4*cxy*Cos(2*phi)*Cos(theta) + (cxx - cyy)*(3 + Cos(2*theta))*Sin(2*phi)) +
        4*Cos(2*phi)*((cxx - cyy)*Cos(theta)*Sin(2*omega) + (cyz*Cos(omega) + cxz*Sin(omega))*Sin(theta)) +
        Sin(2*phi)*(-(cxy*(3 + Cos(2*theta))*Sin(2*omega)) - 2*(cxx + cyy - 2*czz)*Sin(theta)**2 + 2*(cxz*Cos(omega) - cyz*Sin(omega))*Sin(2*theta)))/8.


    Amat[0,2] = (Cos(phi)*Cos(2*theta)*(cxz*Cos(omega) - cyz*Sin(omega)) - Cos(theta)*(cyz*Cos(omega) + cxz*Sin(omega))*Sin(phi) +
        (cxy*Cos(2*omega) + (cxx - cyy)*Sin(2*omega))*Sin(phi)*Sin(theta) - (Cos(phi)*(cxx + cyy - 2*czz + (cxx - cyy)*Cos(2*omega) - cxy*Sin(2*omega))*Sin(2*theta))/2.)/2.


    Amat[1,0] = (Cos(2*omega)*(4*cxy*Cos(2*phi)*Cos(theta) + (cxx - cyy)*(3 + Cos(2*theta))*Sin(2*phi)) +
        4*Cos(2*phi)*((cxx - cyy)*Cos(theta)*Sin(2*omega) + (cyz*Cos(omega) + cxz*Sin(omega))*Sin(theta)) +
        Sin(2*phi)*(-(cxy*(3 + Cos(2*theta))*Sin(2*omega)) - 2*(cxx + cyy - 2*czz)*Sin(theta)**2 + 2*(cxz*Cos(omega) - cyz*Sin(omega))*Sin(2*theta)))/8.


    Amat[1,1] = ((Cos(phi)*Sin(omega) + Cos(omega)*Cos(theta)*Sin(phi))*(Cos(phi)*(cxy*Cos(omega) + 2*cxx*Sin(omega)) + Cos(theta)*(2*cxx*Cos(omega) - cxy*Sin(omega))*Sin(phi) +
        cxz*Sin(phi)*Sin(theta)) + (Cos(omega)*Cos(phi) - Cos(theta)*Sin(omega)*Sin(phi))*
        (cxy*Cos(phi)*Sin(omega) + Cos(omega)*(2*cyy*Cos(phi) + cxy*Cos(theta)*Sin(phi)) + Sin(phi)*(-2*cyy*Cos(theta)*Sin(omega) + cyz*Sin(theta))) +
        Sin(phi)*Sin(theta)*(cxz*Cos(phi)*Sin(omega) + Cos(omega)*(cyz*Cos(phi) + cxz*Cos(theta)*Sin(phi)) + Sin(phi)*(-(cyz*Cos(theta)*Sin(omega)) + 2*czz*Sin(theta))))/2.


    Amat[1,2] = ((Cos(phi)*Sin(omega) + Cos(omega)*Cos(theta)*Sin(phi))*(cxz*Cos(theta) + (-2*cxx*Cos(omega) + cxy*Sin(omega))*Sin(theta)) +
        (Cos(omega)*Cos(phi) - Cos(theta)*Sin(omega)*Sin(phi))*(cyz*Cos(theta) + (-(cxy*Cos(omega)) + 2*cyy*Sin(omega))*Sin(theta)) +
        Sin(phi)*Sin(theta)*(2*czz*Cos(theta) + (-(cxz*Cos(omega)) + cyz*Sin(omega))*Sin(theta)))/2.


    Amat[2,0] = (Cos(phi)*Cos(2*theta)*(cxz*Cos(omega) - cyz*Sin(omega)) - Cos(theta)*(cyz*Cos(omega) + cxz*Sin(omega))*Sin(phi) +
        (cxy*Cos(2*omega) + (cxx - cyy)*Sin(2*omega))*Sin(phi)*Sin(theta) - (Cos(phi)*(cxx + cyy - 2*czz + (cxx - cyy)*Cos(2*omega) - cxy*Sin(2*omega))*Sin(2*theta))/2.)/2.


    Amat[2,1] = (Cos(phi)*Cos(theta)*(cyz*Cos(omega) + cxz*Sin(omega)) + Cos(2*theta)*(cxz*Cos(omega) - cyz*Sin(omega))*Sin(phi) -
        Cos(phi)*(cxy*Cos(2*omega) + (cxx - cyy)*Sin(2*omega))*Sin(theta) - ((cxx + cyy - 2*czz + (cxx - cyy)*Cos(2*omega) - cxy*Sin(2*omega))*Sin(phi)*Sin(2*theta))/2.)/2.


    Amat[2,2] = czz*Cos(theta)**2 + Cos(theta)*(-(cxz*Cos(omega)) + cyz*Sin(omega))*Sin(theta) + (cxx*Cos(omega)**2 - cxy*Cos(omega)*Sin(omega) + cyy*Sin(omega)**2)*Sin(theta)**2


    tx = D[0]
    ty = D[1]
    tz = D[2]

    Avec = numpy.zeros(3)


    Avec[0] = cx - 2*cxx*tx - (cxy*ty + cxz*tz)*Cos(theta) + (cxz*ty - cxy*tz)*Sin(theta)
    Avec[1] = -((cyy + czz)*ty) + (cy - cxy*tx)*Cos(theta) - (cyy*ty - czz*ty + cyz*tz)*Cos(2*theta) + (-cz + cxz*tx)*Sin(theta) + (cyz*ty - cyy*tz + czz*tz)*Sin(2*theta)
    Avec[2] = -((cyz*ty + 2*czz*tz)*Cos(theta)**2) + Sin(theta)*(cy - cxy*tx + (cyz*ty - 2*cyy*tz)*Sin(theta)) + Cos(theta)*(cz - cxz*tx - 2*(cyy*ty - czz*ty + cyz*tz)*Sin(theta))

    A0 = c0 + \
        tx*(-cx + cxx*tx + tz*((cxz*Cos(theta))/2. + (cxy*Sin(theta))/2.) + ty*((cxy*Cos(theta))/2. - (cxz*Sin(theta))/2.)) + \
        tz*(-(cz*Cos(theta)) - cy*Sin(theta) + tx*((cxz*Cos(theta))/2. + (cxy*Sin(theta))/2.) + \
        tz*(Sin(theta)*((cyz*Cos(theta))/2. + cyy*Sin(theta)) + Cos(theta)*(czz*Cos(theta) + (cyz*Sin(theta))/2.)) + \
        ty*(Sin(theta)*(cyy*Cos(theta) - (cyz*Sin(theta))/2.) + Cos(theta)*((cyz*Cos(theta))/2. - czz*Sin(theta)))) + \
        ty*(-(cy*Cos(theta)) + cz*Sin(theta) + tx*((cxy*Cos(theta))/2. - (cxz*Sin(theta))/2.) + \
        tz*(Cos(theta)*((cyz*Cos(theta))/2. + cyy*Sin(theta)) - Sin(theta)*(czz*Cos(theta) + (cyz*Sin(theta))/2.)) + \
        ty*(Cos(theta)*(cyy*Cos(theta) - (cyz*Sin(theta))/2.) - Sin(theta)*((cyz*Cos(theta))/2. - czz*Sin(theta))))

    Alist = quartic_coefficients_matrices_to_list(Amat,Avec,A0, fix_zeros=True)

    # print(">>>", Amat)
    # print(">>>", Avec)
    # print(">>>", A0)
    # print(">>>", Alist)
    return Alist

#
# specific conics
#

def sphere(ssour=10,simag=3,theta_grazing=3e-3):
    theta = (numpy.pi / 2) - theta_grazing
    rmirr = ssour * simag * 2 / numpy.cos(theta) / (ssour + simag)

    s1 = reduced_quadric('sphere')
    s2 = scale_reduced_quadric(s1, xscale=rmirr, yscale=rmirr, zscale=rmirr, return_list=True)
    s3 = expand_reduced_quadric(s2)
    s4 = rotate_and_shift_quartic(s3,
                             omega=0.0, theta=0.0, phi=0.0,
                             D=[0.0,0.0,rmirr])
    s5 = s4.copy()
    for i in range(10):
        s5[i] /= s4[0]

    print("Sphere: ")
    print("   R: ", rmirr)
    print("   reduced: ", s1)
    print("   scaled: ", s2)
    print("   expanded: ", s3)
    print("   rotated and shifted: ", s4)
    print("   normalized: ", s5)
    return {'p':ssour, 'q':simag, 'theta_grazing':theta_grazing, 'radius':rmirr, 'ccc':s5}

def paraboloid(ssour=10,simag=3,theta_grazing=3e-3, verbose=True):

    if ssour >= simag: # focusing
        a = simag * numpy.sin(theta_grazing)**2
        YCEN = - simag * numpy.sin(2 * theta_grazing)
        ZCEN = simag * numpy.cos(theta_grazing)**2

        CENTER = numpy.array([0, YCEN, ZCEN])

        NORMAL = numpy.array((0, -2 * YCEN, 4 * a))
        NORMAL_MOD = numpy.sqrt(NORMAL[0] ** 2 + NORMAL[1] ** 2 + NORMAL[2] ** 2)
        NORMAL /= NORMAL_MOD
        NORMAL_NEW = numpy.array((0, numpy.cos(theta_grazing), numpy.sin(theta_grazing)))
        # Euler angles
        omega = 1 / 2 * numpy.pi
        theta = numpy.pi / 2 - theta_grazing
        phi = 3 / 2 * numpy.pi
    else: # collimating
        a = ssour * numpy.sin(theta_grazing)**2
        YCEN = + ssour * numpy.sin(2 * theta_grazing)
        ZCEN = ssour * numpy.cos(theta_grazing)**2

        CENTER = numpy.array([0, YCEN, ZCEN])

        NORMAL = numpy.array((0, -2 * YCEN, 4 * a))
        NORMAL_MOD = numpy.sqrt(NORMAL[0] ** 2 + NORMAL[1] ** 2 + NORMAL[2] ** 2)
        NORMAL /= NORMAL_MOD
        NORMAL_NEW = numpy.array((0, -numpy.cos(theta_grazing), numpy.sin(theta_grazing)))
        # Euler angles
        omega = 1 / 2 * numpy.pi
        theta = -(numpy.pi / 2 - theta_grazing) #+ numpy.pi
        phi = 3 / 2 * numpy.pi


    PARAM = 2 * a
    if verbose:
        txt = ""
        if ssour >= simag:
            txt += "** Source is at infinity\n"
            txt += "** q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (simag, theta_grazing, theta)
        else:
            txt += "** Image is at infinity\n"
            txt += "** p=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (ssour, theta_grazing, theta)
        txt += '** Parabloid a=%f \n' % a
        txt += '** Parabloid of revolution PARAM=%f \n' % PARAM
        txt += '** 4 * a=%f \n' % (4 * a)
        txt += '** Optical element center at: (%f,%f,%f)\n' % (CENTER[0],CENTER[1],CENTER[2])
        txt += '** Normal: (%f,%f,%f) \n' % (NORMAL[0],NORMAL[1],NORMAL[2])
        txt += '** Normal NEW!!: (%f,%f,%f) \n' % (NORMAL_NEW[0],NORMAL_NEW[1],NORMAL_NEW[2])
        print(txt)



    ROTATED_CENTER = numpy.dot(euler_rotation_matrix(omega, theta, phi), CENTER)

    s1 = reduced_quadric('paraboloid')
    s2 = scale_reduced_quadric(s1, xscale=1.0, yscale=1.0, zscale=(1/4/a), return_list=True)
    s3 = expand_reduced_quadric(s2)

    s4 = rotate_and_shift_quartic_NEW(s3,
                                      omega=omega, theta=theta, phi=phi,
                                      D=-ROTATED_CENTER)

    s5 = s4.copy()
    # for i in range(10):
    #     s5[i] /= s4[0]

    print("**Paraboloid: ")
    print("**   a, theta[deg]: ", PARAM, theta*180/numpy.pi)
    print("**   reduced: ", s1)
    print("**   scaled: ", s2)
    print("**   expanded: ", s3)
    print("**   rotated and shifted: ", s4)
    print("**   normalized: ", s5)
    return {'p':ssour, 'q':simag, 'theta_grazing':theta_grazing, 'a':a, 'center':CENTER, 'normal':NORMAL,  'ccc':s5}


def ellipsoid(ssour=10, simag=3, theta_grazing=3e-3, verbose=True):

    a = 0.5 * (ssour + simag)
    b = numpy.sqrt(ssour * simag) * numpy.sin(theta_grazing)
    c = numpy.sqrt(a**2 - b**2)

    YCEN = (ssour**2 - simag**2) / 4 / c
    ZCEN = -b * numpy.sqrt(1 - YCEN**2 / a**2)
    CENTER = numpy.array([0,YCEN, ZCEN])

    NORMAL = numpy.array((0, -2 * YCEN / a**2, -2 * ZCEN / b**2))
    NORMAL_MOD = numpy.sqrt(NORMAL[0]**2 + NORMAL[1]**2 + NORMAL[2]**2)
    NORMAL /= NORMAL_MOD

    AXMAJ = a
    AXMIN = b
    AFOCI = c
    ECCENT = c / a

    #
    # Euler angles
    #
    omega = 1/2 * numpy.pi
    theta = numpy.arcsin(NORMAL[1]) # -theta_grazing # numpy.arccos(RNCEN[3-1])
    phi = 3/2 * numpy.pi

    if verbose:
        txt = ""
        txt += "** p=%f, q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (ssour, simag, theta_grazing, (numpy.pi / 2) - theta_grazing)
        txt += '** Ellipsoid of revolution a=%f \n' % AXMAJ
        txt += '** Ellipsoid of revolution b=%f \n' % AXMIN
        txt += '** Ellipsoid of revolution 1/a**2=%f \n' % (1/AXMAJ**2)
        txt += '** Ellipsoid of revolution 1/b**2=%f \n' % (1/AXMIN**2)
        txt += '** Ellipsoid of revolution c=sqrt(a^2-b^2)=%f \n' % AFOCI
        txt += '** Ellipsoid of revolution focal distance c^2=%f \n' % (AFOCI ** 2)
        txt += '** Ellipsoid of revolution excentricity: %f \n' % ECCENT
        txt += '** Optical element center at: [%f,%f,%f]\n' % (CENTER[0], CENTER[1], CENTER[2])
        txt += '** Optical element normal: [%f,%f,%f]\n' % (NORMAL[0], NORMAL[1], NORMAL[2])
        txt += '** Optical element tangent: [%f,%f,%f]\n' % (NORMAL[0], NORMAL[2], -NORMAL[1])
        txt += '** THETA from NORMAL %f deg\n' % (numpy.arcsin(NORMAL[1]) * 180 / numpy.pi)
        txt += '** THETA from NORMAL %f deg\n' % (numpy.arccos(NORMAL[2]) * 180 / numpy.pi)
        txt += '** THETA from EULER %f deg\n' % (theta * 180 / numpy.pi)
        txt += '** B nz ycen - A ny zcen: %f\n' % ((1/a**2) *NORMAL[2] * YCEN - (1/b**2) * NORMAL[1] * ZCEN)
        print(txt)

        print(txt)



    s1 = reduced_quadric('sphere')
    s2 = scale_reduced_quadric(s1, xscale=AXMIN, yscale=AXMAJ, zscale=AXMIN, return_list=True)
    s3 = expand_reduced_quadric(s2)
    D=-numpy.dot(euler_rotation_matrix(omega, theta, phi), CENTER)
    s4 = rotate_and_shift_quartic_NEW(s3,
                             omega=omega, theta=theta, phi=phi,
                             D=-numpy.dot(euler_rotation_matrix(omega, theta, phi), CENTER))
    s5 = s4.copy()
    # for i in range(10):
    #     s5[i] /= s4[0]

    print("**Ellipsoid: ")
    print("**   a,b, theta_grazing[rad]: ", AXMAJ, AXMIN, theta_grazing)
    print("**   euler [deg]: ", omega * 180 / numpy.pi, theta * 180 / numpy.pi, phi * 180 / numpy.pi)
    print("**   D: ", D[0],D[1],D[2])
    print("**   reduced: ", s1)
    print("**   scaled: ", s2)
    print("**   expanded: ", s3)
    print("**   rotated and shifted: ", s4)
    print("**   normalized: ", s5)

    return {'p':ssour, 'q':simag, 'theta_grazing':theta_grazing,
                   'a':a, 'b':b, 'c':c,
                   'center':CENTER, 'normal':NORMAL,  'ccc':s5}


def hyperboloid(ssour=10, simag=3, theta_grazing=3e-3, verbose=True):

    a = 0.5 * numpy.abs(ssour - simag)
    c = 0.5 * numpy.sqrt(ssour**2 + simag**2 - 2 * ssour * simag * numpy.cos(2 * theta_grazing))
    b = numpy.sqrt(c**2 - a**2)

    # Large p (p>q)
    if ssour > simag: # select center in first quadrant
        YCEN = (ssour**2 - simag**2) / 4 / c
        ZCEN = b * numpy.sqrt(YCEN**2 / a**2 - 1)
        NORMAL = numpy.array((0, -2 * YCEN / a ** 2, 2 * ZCEN / b ** 2)) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        NORMAL_MOD = numpy.sqrt(NORMAL[0] ** 2 + NORMAL[1] ** 2 + NORMAL[2] ** 2)
        NORMAL /= NORMAL_MOD
        #
        # Euler angles
        #
        omega = 1 / 2 * numpy.pi
        theta = numpy.arcsin(NORMAL[1]) #<<<<<<<<<<<<<<<<<<<<<
        phi = 3 / 2 * numpy.pi
    # Large q
    else: # center in 2nd quadrant
        YCEN = (ssour**2 - simag**2) / 4 / c
        ZCEN = b * numpy.sqrt(YCEN**2 / a**2 - 1)
        NORMAL = -numpy.array((0, -2 * YCEN / a ** 2, 2 * ZCEN / b ** 2)) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        NORMAL_MOD = numpy.sqrt(NORMAL[0] ** 2 + NORMAL[1] ** 2 + NORMAL[2] ** 2)
        NORMAL /= NORMAL_MOD
        #
        # Euler angles
        #
        omega = 1 / 2 * numpy.pi
        theta = -numpy.arccos(NORMAL[2]) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        phi = 3 / 2 * numpy.pi

    CENTER = numpy.array([0,YCEN, ZCEN])

    ecc = c / a


    # #
    # # Euler angles
    # #
    # omega = 1/2 * numpy.pi
    # theta = numpy.arcsin(NORMAL[1])# numpy.pi - numpy.abs(numpy.arcsin(NORMAL[1])) # numpy.arcsin(NORMAL[1]) # -theta_grazing # numpy.arccos(RNCEN[3-1])
    # phi = 3/2 * numpy.pi

    if verbose:
        txt = ""
        txt += "** p=%f, q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (ssour, simag, theta_grazing, (numpy.pi / 2) - theta_grazing)
        txt += '** Hyperboloid of revolution a=%f \n' % a
        txt += '** Hyperboloid of revolution b=%f \n' % b
        txt += '** Hyperboloid of revolution 1/a**2=%f \n' % (1/a**2)
        txt += '** Hyperboloid of revolution 1/b**2=%f \n' % (1/b**2)
        txt += '** Hyperboloid of revolution c=sqrt(a^2+b^2)=%f \n' % c
        txt += '** Hyperboloid of revolution focal distance c^2=%f \n' % (c ** 2)
        txt += '** Hyperboloid of revolution excentricity: %f \n' % ecc
        txt += '** Optical element center at: [%f,%f,%f]\n' % (CENTER[0], CENTER[1], CENTER[2])
        txt += '** Optical element normal: [%f,%f,%f]\n' % (NORMAL[0], NORMAL[1], NORMAL[2])
        txt += '** Optical element tangent: [%f,%f,%f]\n' % (NORMAL[0], NORMAL[2], -NORMAL[1])
        txt += '** THETA from NORMAL %f rad\n' % (numpy.arcsin(NORMAL[1]) )
        txt += '** THETA from NORMAL %f rad\n' % (numpy.arccos(NORMAL[2]) )
        txt += '** THETA from euler %f rad\n' % (theta )
        # txt += '** B nz ycen - A ny zcen: %f\n' % ((1/a**2) *NORMAL[2] * YCEN - (1/b**2) * NORMAL[1] * ZCEN)
        print(txt)

        print(txt)


    s1 = [-1,1,-1,0,-1] # reduced_quadric('one sheet hyperboloid')
    s2 = scale_reduced_quadric(s1, xscale=b, yscale=a, zscale=b, return_list=True)
    s3 = expand_reduced_quadric(s2)

    D = -numpy.dot(euler_rotation_matrix(omega, theta, phi), CENTER)
    # D = numpy.dot(euler_rotation_matrix(0, numpy.pi, 0.0), D) # make the incident beam in the negative part

    s4 = rotate_and_shift_quartic_NEW(s3,
                             omega=omega, theta=theta, phi=phi,
                             D=D)
    # if ssour < simag:
    #     # make the incident beam in the negative part
    #     s4 = rotate_and_shift_quartic_NEW(s4,
    #                              omega=numpy.pi/2, theta=numpy.pi, phi=numpy.pi*3/2,
    #                              D=[0,0,0])
    s5 = s4.copy()
    # for i in range(10):
    #     s5[i] /= s4[0]

    print("**Hyperboloid: ")
    print("**   a,b, theta_grazing[rad]: ", a, b, theta_grazing)
    print("**   euler [deg]: ", omega * 180 / numpy.pi, theta * 180 / numpy.pi, phi * 180 / numpy.pi)
    print("**   D: ", D[0],D[1],D[2])
    print("**   rotated N: ", numpy.dot(euler_rotation_matrix(omega, theta, phi), NORMAL))
    print("**   reduced: ", s1)
    print("**   scaled: ", s2)
    print("**   expanded: ", s3)
    print("**   rotated and shifted: ", s4)
    print("**   normalized: ", s5)
    A = -1/b**2
    B = 1/a**2
    ny = NORMAL[1]
    nz = NORMAL[2]
    ccc = [A, A*ny**2+B*nz**2,A*nz**2+B*ny**2,0,2*(B-A)*ny*nz,0.,0.,0.,2*(B*ny*YCEN+A*nz*ZCEN),0.]
    print("**   using SHADOW way: ", ccc)

    return {'p':ssour, 'q':simag, 'theta_grazing':theta_grazing,
                   'a':a, 'b':b, 'c':c,
                   'center':CENTER, 'normal':NORMAL,  'ccc':s5}
#
# TESTING ROUTINES
#
def sphere_check():
    s5 = sphere(ssour=10,simag=3,theta_grazing=3e-3)
    c=numpy.zeros(10)
    c[1-1]= 1
    c[2-1]= 1
    c[3-1]= 1
    c[4-1]= 0
    c[5-1]= 0
    c[6-1]= 0
    c[7-1]= 0
    c[8-1]= 0
    c[9-1]= -3076.93
    c[10-1]= 0

    for i in range(10):
        print(s5['ccc'][i] , c[i])
        assert(numpy.abs(s5['ccc'][i] - c[i]) < 1e-2)

def ellipsoid_check(ssour=10,simag=3,theta_grazing=3e-3, do_plot=False):


    ccc = S4Conic.initialize_as_ellipsoid_from_focal_distances(ssour, simag, theta_grazing,
                                        cylindrical=0, cylangle=0.0, switch_convexity=0)
    print("ccc: ", ccc.get_coefficients())

    s5 = ellipsoid(ssour=ssour,simag=simag,theta_grazing=theta_grazing)

    c = ccc.get_coefficients()
    print("ccc: ", c)
    print("s5: ", s5['ccc'])


    for i in range(10):
        print(i, c[i], s5['ccc'][i])
        assert(numpy.abs(s5['ccc'][i] - c[i]) < 1e-2)

    # view_conic(s5, x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1)
    if do_plot:
        compare_conics(s5['ccc'], ccc.get_coefficients(), x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1,
                       titles=['s5','ccc'])

def hyperboloid_check(ssour=10,simag=3,theta_grazing=3e-3, do_plot=False):


    # ccc = S4Conic.initialize_as_hyperboloid_from_focal_distances(ssour, simag, theta_grazing,
    #                                     cylindrical=0, cylangle=0.0, switch_convexity=0)
    # print("ccc: ", ccc.get_coefficients())

    if ssour < simag:
        c = [-3703.714814855263, -0.033333333332666124, -3703.5998488688692, 0.0, 41.269717460237345, 0.0, 0.0, 3.0795921368564905e-12, 190.47647619130197, 0.0]
        branch = 0
    else:
        c = [-3703.714814855263, -0.03333333333333342, -3703.5998488688683, 0.0, -41.269717460357114, 0.0, 0.0, 0.0, -190.47647619130197, 0.0]
        branch = 0

    s5 = hyperboloid(ssour=ssour,simag=simag,theta_grazing=theta_grazing)

    print("ccc: ", c)
    print("s5: ", s5['ccc'])


    for i in range(10):
        print(i, c[i], s5['ccc'][i])
        assert(numpy.abs(s5['ccc'][i] - c[i]) < 1e-2)

    # view_conic(s5, x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1)
    if do_plot:
        compare_conics(s5['ccc'], c, x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1,
                       titles=['s5','ccc'], branch=branch)

    return s5

def parabola_check(ssour=10,simag=10,theta_grazing=3e-3, do_plot=False):


    ccc = S4Conic.initialize_as_paraboloid_from_focal_distances(ssour, simag, theta_grazing,
                                        cylindrical=0, cylangle=0.0, switch_convexity=0)
    print("ccc: ", ccc.get_coefficients())

    s5 = paraboloid(ssour=ssour,simag=simag,theta_grazing=theta_grazing)

    print("ccc: ", ccc.get_coefficients())
    print("s5: ", s5['ccc'])

    c = ccc.get_coefficients()
    for i in range(10):
        print(i, c[i], s5['ccc'][i])

    for i in range(10):
        print(s5['ccc'][i] , c[i])
        assert(numpy.abs(s5['ccc'][i] - c[i]) < 1e-2)

    # view_conic(s5, x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1)
    if do_plot:
        compare_conics(s5, ccc.get_coefficients(), x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1,
                       titles=['s5','ccc'])


def height(ccc,y=0,x=0,return_solution=0):
    """

    :param y: a scalar, vector or mesh
    :param x: a scalar, vector or mesh
        y and x must be homogeneous, otherwise an error will occur:
         both scalars
         both mesh
         one scalar and another vector
    :param return_solution: 0 = guess the solution with zero at pole,
                            1 = get first solution
                            2 = get second solution
    :return: the height scalar/vector/mesh depending on inputs
    """
    aa = ccc[2]
    bb = ccc[4] * y + ccc[5] * x + ccc[8]
    cc = ccc[0] * x**2 + ccc[1] * y**2 + ccc[3] * x * y + \
        ccc[6] * x + ccc[7] * y + ccc[9]

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

if __name__ == "__main__":



    sphere_check()

    parabola_check(ssour=1e8,simag=10,theta_grazing=3e-3, do_plot=0)
    parabola_check(ssour=10,simag=1e9,theta_grazing=3e-3, do_plot=0)
    #
    ellipsoid_check(ssour=10,simag=3,theta_grazing=3e-3, do_plot=0)
    #
    #
    hyperboloid_check(ssour=10,simag=3,theta_grazing=3e-3, do_plot=0)
    hyperboloid_check(ssour=3, simag=10, theta_grazing=3e-3, do_plot=0)


