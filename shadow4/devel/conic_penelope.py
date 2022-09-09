import numpy

from shadow4.optical_surfaces.s4_conic import S4Conic
from conic_viewer import view_conic, compare_conics

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


def rotate_and_shift_quartic_NEW(quartic_coefficients_list,
                             omega=0.0, theta=0.0, phi=0.0,
                             D=[0,0,0]):
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



"""
C  *********************************************************************
C                       SUBROUTINE ROTSHF
C  *********************************************************************
      SUBROUTINE ROTSHF(OMEGA,THETA,PHI,DX,DY,DZ,
     1                  AXX,AXY,AXZ,AYY,AYZ,AZZ,AX,AY,AZ,A0)
C
C     This subroutine rotates and shifts a quadric surface.
CC  *********************************************************************
C                       SUBROUTINE ROTSHF
C  *********************************************************************
      SUBROUTINE ROTSHF(OMEGA,THETA,PHI,DX,DY,DZ,
     1                  AXX,AXY,AXZ,AYY,AYZ,AZZ,AX,AY,AZ,A0)
C
C     This subroutine rotates and shifts a quadric surface.
C
C  Input parameters:
C     OMEGA, THETA, PHI ... Euler rotation angles,
C     DX, DY, DZ .......... components of the displacement vector,
C     AXX, ..., A0 ........ coefficients of the initial quadric.
C
C  Output parameters:
C     AXX, ..., A0 ........ coefficients of the transformed quadric.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION R(3,3),A2(3,3),B2(3,3),A1(3),B1(3),D1(3)
C
C  ****  Initial quadric.
C
      B2(1,1)=AXX
      B2(1,2)=0.5D0*AXY
      B2(1,3)=0.5D0*AXZ
      B2(2,1)=B2(1,2)
      B2(2,2)=AYY
      B2(2,3)=0.5D0*AYZ
      B2(3,1)=B2(1,3)
      B2(3,2)=B2(2,3)
      B2(3,3)=AZZ
      B1(1)=AX
      B1(2)=AY
      B1(3)=AZ
      B0=A0
      D1(1)=DX
      D1(2)=DY

C  Input parameters:
C     OMEGA, THETA, PHI ... Euler rotation angles,
C     DX, DY, DZ .......... components of the displacement vector,
C     AXX, ..., A0 ........ coefficients of the initial quadric.
C
C  Output parameters:
C     AXX, ..., A0 ........ coefficients of the transformed quadric.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z), INTEGER*4 (I-N)
      DIMENSION R(3,3),A2(3,3),B2(3,3),A1(3),B1(3),D1(3)
C
C  ****  Initial quadric.
C
      B2(1,1)=AXX
      B2(1,2)=0.5D0*AXY
      B2(1,3)=0.5D0*AXZ
      B2(2,1)=B2(1,2)
      B2(2,2)=AYY
      B2(2,3)=0.5D0*AYZ
      B2(3,1)=B2(1,3)
      B2(3,2)=B2(2,3)
      B2(3,3)=AZZ
      B1(1)=AX
      B1(2)=AY
      B1(3)=AZ
      B0=A0
      D1(1)=DX
      D1(2)=DY
      D1(3)=DZ
C
C  ****  Rotation matrix.
C
      STHETA=SIN(THETA)
      CTHETA=COS(THETA)
      SPHI=SIN(PHI)
      CPHI=COS(PHI)
      SOMEGA=SIN(OMEGA)
      COMEGA=COS(OMEGA)
C
      R(1,1)=CPHI*CTHETA*COMEGA-SPHI*SOMEGA
      R(1,2)=-CPHI*CTHETA*SOMEGA-SPHI*COMEGA
      R(1,3)=CPHI*STHETA
      R(2,1)=SPHI*CTHETA*COMEGA+CPHI*SOMEGA
      R(2,2)=-SPHI*CTHETA*SOMEGA+CPHI*COMEGA
      R(2,3)=SPHI*STHETA
      R(3,1)=-STHETA*COMEGA
      R(3,2)=STHETA*SOMEGA
      R(3,3)=CTHETA
C
C  ****  Rotated quadric.
C
      DO I=1,3
        A1(I)=0.0D0
        DO J=1,3
          A1(I)=A1(I)+R(I,J)*B1(J)
          A2(I,J)=0.0D0
          DO M=1,3
            DO K=1,3
              A2(I,J)=A2(I,J)+R(I,K)*B2(K,M)*R(J,M)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
C  ****  Shifted-rotated quadric.
C
      DO I=1,3
        A2D=0.0D0
        DO J=1,3
          A2D=A2D+A2(I,J)*D1(J)
        ENDDO
        B1(I)=A1(I)-2.0D0*A2D
        B0=B0+D1(I)*(A2D-A1(I))
      ENDDO
C
      AXX=A2(1,1)
      IF(ABS(AXX).LT.1.0D-15) AXX=0.0D0
      AXY=A2(1,2)+A2(2,1)
      IF(ABS(AXY).LT.1.0D-15) AXY=0.0D0
      AXZ=A2(1,3)+A2(3,1)
      IF(ABS(AXZ).LT.1.0D-15) AXZ=0.0D0
      AYY=A2(2,2)
      IF(ABS(AYY).LT.1.0D-15) AYY=0.0D0
      AYZ=A2(2,3)+A2(3,2)
      IF(ABS(AYZ).LT.1.0D-15) AYZ=0.0D0
      AZZ=A2(3,3)
      IF(ABS(AZZ).LT.1.0D-15) AZZ=0.0D0
      AX=B1(1)
      IF(ABS(AX).LT.1.0D-15) AX=0.0D0
      AY=B1(2)
      IF(ABS(AY).LT.1.0D-15) AY=0.0D0
      AZ=B1(3)
      IF(ABS(AZ).LT.1.0D-15) AZ=0.0D0
      A0=B0
      IF(ABS(A0).LT.1.0D-15) A0=0.0D0
      RETURN
      END

"""

def euler_rotation_matrix(omega,theta,phi, shortcut=False):

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

    return R

def rotate_and_shift_quartic(quartic_coefficients_list,
                             omega=0.0, theta=0.0, phi=0.0,
                             D=[0,0,0]):
#
# initial quartic
#
    B2 = numpy.zeros((3,3))

    AXX,AYY,AZZ,AXY,AYZ,AXZ,AX,AY,AZ,A0 = quartic_coefficients_list

      # B2(1,1)=AXX
      # B2(1,2)=0.5D0*AXY
      # B2(1,3)=0.5D0*AXZ
      # B2(2,1)=B2(1,2)
      # B2(2,2)=AYY
      # B2(2,3)=0.5D0*AYZ
      # B2(3,1)=B2(1,3)
      # B2(3,2)=B2(2,3)
      # B2(3,3)=AZZ
      # B1(1)=AX
      # B1(2)=AY
      # B1(3)=AZ
      # B0=A0
      # D1(1)=DX
      # D1(2)=DY

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

    # so far, input [or acaled] F(x,y,z) can be written in matrix form (ibid, eq. 6.25)
    # notation: X=(x,y,z) ; B1=(AX,AY,AZ)
    # F = X' . (B2 X) + B1 . X + B0

    D1 = numpy.array(D)

    #
    #  ****  Rotation matrix.
    #
    R = euler_rotation_matrix(omega, theta, phi)

    #
    #  ****  Rotated quadric.
    #
      # DO I=1,3
      #   A1(I)=0.0D0
      #   DO J=1,3
      #     A1(I)=A1(I)+R(I,J)*B1(J)
      #     A2(I,J)=0.0D0
      #     DO M=1,3
      #       DO K=1,3
      #         A2(I,J)=A2(I,J)+R(I,K)*B2(K,M)*R(J,M)
      #       ENDDO
      #     ENDDO
      #   ENDDO
      # ENDDO
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

    # The loop before can be splited in
    # first equation 6.29
    # for I in range(3):
    #     for J in range(3):
    #         A2[I,J] = 0.0
    #         for M in range(3):
    #             for K in range(3):
    #                 A2[I,J] = A2[I,J] + R[I,K] * B2[K,M] * R[J,M]
    # First term (R A) of the 2nd equation 6.29
    # for I in range(3):
    #     A1[I] = 0.0
    #     for J in range(3):
    #         A1[I] = A1[I] + R[I,J] * B1[J]


    #
    #  ****  Shifted-rotated quadric.
    #


#       DO I=1,3
#         A2D=0.0D0
#         DO J=1,3
#           A2D=A2D+A2(I,J)*D1(J)
#         ENDDO
#         B1(I)=A1(I)-2.0D0*A2D
#         B0=B0+D1(I)*(A2D-A1(I))
#       ENDDO
    # 3rd equation 6.29 and complete the second equation
    for I in range(3):
        A2D = 0.0
        for J in range(3):
            A2D = A2D + A2[I,J] * D1[J]

        B1[I] = A1[I] - 2.0 * A2D
        B0 = B0 + D1[I] * (A2D - A1[I])


    #       AXX=A2(1,1)
    #       AXY=A2(1,2)+A2(2,1)
    #       AXZ=A2(1,3)+A2(3,1)
    #       AYY=A2(2,2)
    #       AYZ=A2(2,3)+A2(3,2)
    #       AZZ=A2(3,3)
    #       AX=B1(1)
    #       AY=B1(2)
    #       AZ=B1(3)
    #       A0=B0


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
    """
    C  ****  Expanded quadric.
 102  CONTINUE
      QXX=KQ(1)/XSCALE**2
      QXY=0.0D0
      QXZ=0.0D0
      QYY=KQ(2)/YSCALE**2
      QYZ=0.0D0
      QZZ=KQ(3)/ZSCALE**2
      QX=0.0D0
      QY=0.0D0
      QZ=KQ(4)/ZSCALE
      Q0=KQ(5)
    """

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
    return s5

def paraboloid(ssour=10,simag=3,theta_grazing=3e-3, verbose=True):

    theta = (numpy.pi / 2) - theta_grazing


    if ssour >= simag: # focusing
        a = simag * numpy.sin(theta_grazing)**2
        YCEN = - simag * numpy.sin(2 * theta_grazing)
        ZCEN = simag * numpy.cos(theta_grazing)**2
        # Euler angles
        omega = 1 / 2 * numpy.pi
        theta = numpy.pi / 2 - theta_grazing
        phi = 3 / 2 * numpy.pi
    else: # collimating
        a = ssour * numpy.sin(theta_grazing)**2
        YCEN = - ssour * numpy.sin(2 * theta_grazing)
        ZCEN = ssour * numpy.cos(theta_grazing)**2
        # Euler angles
        omega = 1 / 2 * numpy.pi
        theta = numpy.pi / 2 - theta_grazing #+ numpy.pi
        phi = 3 / 2 * numpy.pi


    CENTER = numpy.array([0,YCEN, ZCEN])


    NORMAL = numpy.array((0, -2 * YCEN, 4 * a))
    NORMAL_MOD = numpy.sqrt(NORMAL[0]**2 + NORMAL[1]**2 + NORMAL[2]**2)
    NORMAL /= NORMAL_MOD


    PARAM = 2 * a
    if verbose:
        txt = ""
        if ssour >= simag:
            txt += "Source is at infinity\n"
            txt += "q=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (simag, theta_grazing, theta)
        else:
            txt += "Image is at infinity\n"
            txt += "p=%f, theta_grazing=%f rad, theta_normal=%f rad\n" % (ssour, theta_grazing, theta)
        txt += 'Parabloid of revolution PARAM=%f \n' % PARAM
        txt += '4 * a=%f \n' % (4 * a)
        txt += 'Optical element center at: (%f,%f,%f)\n' % (CENTER[0],CENTER[1],CENTER[2])

        txt += 'Normal: (%f,%f,%f) \n' % (NORMAL[0],NORMAL[1],NORMAL[2])
        print(txt)


    print("Euler-rotated 010: ", numpy.dot(euler_rotation_matrix(omega,theta,phi),[0,1,0]))


    ROTATED_CENTER = numpy.dot(euler_rotation_matrix(omega, theta, phi), CENTER)
    print("Euler-rotated Center: ", ROTATED_CENTER)

    s1 = reduced_quadric('paraboloid')
    s2 = scale_reduced_quadric(s1, xscale=1.0, yscale=1.0, zscale=(1/4/a), return_list=True)
    s3 = expand_reduced_quadric(s2)

    if ssour >= simag: # focusing
        s4 = rotate_and_shift_quartic_NEW(s3,
                                          omega=omega, theta=theta, phi=phi,
                                          D=-ROTATED_CENTER)
    else: # collimating
        s4 = rotate_and_shift_quartic_NEW(s3,
                                          omega=omega, theta=theta, phi=phi,
                                          D=-ROTATED_CENTER)
        s4 = rotate_and_shift_quartic_NEW(s4,
                                          omega=numpy.pi, theta=0.0, phi=0.0,
                                          D=[0,0,0])

    s5 = s4.copy()
    # for i in range(10):
    #     s5[i] /= s4[0]

    print("Paraboloid: ")
    print("   a, theta[deg]: ", PARAM, theta*180/numpy.pi)
    print("   reduced: ", s1)
    print("   scaled: ", s2)
    print("   expanded: ", s3)
    print("   rotated and shifted: ", s4)
    print("   normalized: ", s5)
    return s5


def ellipsoid(ssour=10,simag=3,theta_grazing=3e-3):
    theta = (numpy.pi / 2) - theta_grazing
    COSTHE = numpy.cos(theta)
    SINTHE = numpy.sin(theta)

    AXMAJ = (ssour + simag) / 2
    AXMIN = numpy.sqrt(simag * ssour) * COSTHE

    AFOCI = numpy.sqrt(AXMAJ ** 2 - AXMIN ** 2)
    ECCENT = AFOCI / AXMAJ
    # ;C
    # ;C The center is computed on the basis of the object and image positions
    # ;C
    YCEN = (ssour - simag) * 0.5 / ECCENT
    ZCEN = -numpy.sqrt(1 - YCEN ** 2 / AXMAJ ** 2) * AXMIN

    RNCEN = numpy.zeros(3)
    RNCEN[1 - 1] = 0.0
    RNCEN[2 - 1] = -2 * YCEN / AXMAJ ** 2
    RNCEN[3 - 1] = -2 * ZCEN / AXMIN ** 2
    # ;CALL NORM(RNCEN,RNCEN)
    RNCEN = RNCEN / numpy.sqrt((RNCEN ** 2).sum())


    # Euler angles
    #
    omega = 1/2 * numpy.pi
    theta = numpy.arccos(RNCEN[3-1])
    phi = 3/2 * numpy.pi

    # Theta = numpy.arccos(RNCEN[2])

    print("N: ", RNCEN)
    print("Euler-rotated 001: ", numpy.dot(euler_rotation_matrix(omega,theta,phi),[0,0,1]))
    print( " or: ",
          numpy.cos(phi)*numpy.sin(theta),
          numpy.sin(phi) * numpy.sin(theta),
          numpy.cos(theta))


    print("Center: ", 0,YCEN, ZCEN)
    CENTER = numpy.array([0, YCEN, ZCEN])
    ROTATED_CENTER = numpy.dot(euler_rotation_matrix(omega, theta, phi), CENTER)
    # ROTATED_CENTER = numpy.array([0, YCEN, ZCEN])
    # ROTATED_CENTER = numpy.dot(euler_rotation_matrix(-phi, theta, -omega), [0, YCEN, ZCEN])
    print("Euler-rotated Center: ", ROTATED_CENTER)

    # omega = numpy.arccos(-RNCEN[2-1]/numpy.sqrt(1-RNCEN[3-1]**2))
    # theta = numpy.arccos(RNCEN[3-1])
    # Y3 = numpy.sin(theta) * numpy.cos()
    # phi = numpy.arccos(RNCEN[3-1]/numpy.sqrt(1-RNCEN[3-1]**2))

    # print("Euler omega, theta, phi: ", omega*180/numpy.pi, theta*180/numpy.pi, phi*180/numpy.pi)
    s1 = reduced_quadric('sphere')
    s2 = scale_reduced_quadric(s1, xscale=AXMIN, yscale=AXMAJ, zscale=AXMIN, return_list=True)
    s3 = expand_reduced_quadric(s2)
    s4 = rotate_and_shift_quartic_NEW(s3,
                             omega=omega, theta=theta, phi=phi,
                             D=-ROTATED_CENTER)
    s5 = s4.copy()
    # for i in range(10):
    #     s5[i] /= s4[0]

    print("Ellipsoid: ")
    print("   a,b, theta[deg]: ", AXMAJ, AXMIN, theta*180/numpy.pi)
    print("   reduced: ", s1)
    print("   scaled: ", s2)
    print("   expanded: ", s3)
    print("   rotated and shifted: ", s4)
    print("   normalized: ", s5)
    return s5

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
        print(s5[i] , c[i])
        assert(numpy.abs(s5[i] - c[i]) < 1e-2)

def ellipsoid_check(ssour=10,simag=3,theta_grazing=3e-3, do_plot=False):


    ccc = S4Conic.initialize_as_ellipsoid_from_focal_distances(ssour, simag, theta_grazing,
                                        cylindrical=0, cylangle=0.0, switch_convexity=0)
    print("ccc: ", ccc.get_coefficients())

    s5 = ellipsoid(ssour=ssour,simag=simag,theta_grazing=theta_grazing)

    print("ccc: ", ccc.get_coefficients())
    print("s5: ", s5)

    for i in range(10):
        print(i, ccc.get_coefficients()[i], s5[i])

    # view_conic(s5, x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1)
    if do_plot:
        compare_conics(s5, ccc.get_coefficients(), x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1,
                       titles=['s5','ccc'])

def parabola_check(ssour=10,simag=10,theta_grazing=3e-3, do_plot=False):


    ccc = S4Conic.initialize_as_paraboloid_from_focal_distances(ssour, simag, theta_grazing,
                                        cylindrical=0, cylangle=0.0, switch_convexity=0)
    print("ccc: ", ccc.get_coefficients())

    s5 = paraboloid(ssour=ssour,simag=simag,theta_grazing=theta_grazing)

    print("ccc: ", ccc.get_coefficients())
    print("s5: ", s5)

    c = ccc.get_coefficients()
    for i in range(10):
        print(i, c[i], s5[i])

    for i in range(10):
        print(s5[i] , c[i])
        assert(numpy.abs(s5[i] - c[i]) < 1e-2)

    # view_conic(s5, x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1)
    if do_plot:
        compare_conics(s5, ccc.get_coefficients(), x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1,
                       titles=['s5','ccc'])

if __name__ == "__main__":

    # sphere_check()


    # sphere_check()
    # ellipsoid_check(do_plot=True)


    parabola_check(ssour=1e8,simag=10,theta_grazing=3e-3, do_plot=0)
    parabola_check(ssour=10,simag=1e9,theta_grazing=3e-3, do_plot=0)