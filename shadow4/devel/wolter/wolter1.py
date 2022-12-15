import numpy
from shadow4.devel.wolter.conic_penelope import ellipsoid, hyperboloid, paraboloid, rotate_and_shift_quartic, euler_rotation_matrix

from shadow4.devel.wolter.conics_from_factory_parameters import ellipsoid as fac_ellipsoid
from shadow4.devel.wolter.conics_from_factory_parameters import paraboloid as fac_paraboloid
from shadow4.devel.wolter.conics_from_factory_parameters import hyperboloid as fac_hyperboloid
from shadow4.devel.wolter.conics_from_factory_parameters import ken_ellipsoid, ken_hyperboloid, ken_paraboloid

def cyl(ccc):
    ccc1 = ccc.copy()
    ccc1[0] = 0
    ccc1[3] = 0
    ccc1[5] = 0
    ccc1[6] = 0
    return ccc1


def recipe1(
    p_ell = 10.0,
    q_ell = 3.0,
    distance = 0.3,
    theta = 0.003,
    ratio_hyp = 3.0, # ratio_hyp = q_hyp / p_ell > 1.0
    verbose = 1,
                ):


    q_hyp = q_ell - distance
    a_hyp = 0.5 * q_hyp * (1 - 1/ratio_hyp)
    p_hyp = q_hyp - 2 * a_hyp

    tkt_ell = ellipsoid(p_ell,q_ell,theta)
    tkt_hyp = hyperboloid(p_hyp,q_hyp,theta)

    if verbose:
        print("ell p,q", p_ell, q_ell)
        print("hyp p,q", p_hyp, q_hyp)

        print("position ell p,q", p_ell, distance / 2)
        print("position hyp p,q", distance / 2, p_hyp)

        print("length=", p_ell + distance + p_hyp)
        m_ell = (q_ell / p_ell)
        m_hyp = (p_hyp / q_hyp)
        print("magnification ell: ", m_ell)
        print("magnification hyp: ", m_hyp)
        print("magnification: ", m_ell * m_hyp)
        print("demagnification: ", 1 / (m_ell * m_hyp))

    return tkt_ell,tkt_hyp


def recipe2(
    p_ell = 10.0,
    distance = 0.3,
    p_hyp = 0.9,
    theta = 0.003,
    # m_ell = 0.3,
    m_hyp = 1/3,
    verbose = 1,
            ):
    q_hyp = p_hyp / m_hyp
    # a_hyp = 0.5 * (q_hyp - p_hyp)


    q_ell = q_hyp + distance

    m_ell = q_ell / p_ell


    if q_ell >= p_ell:
        raise Exception("must be: q_ell < p_ell")

    if q_ell >= p_ell:
        raise Exception("must be: q_ell < p_ell")

    if q_hyp >= q_ell:
        raise Exception("must be: q_hyp < q_ell")

    if q_hyp <= p_hyp:
        raise Exception("must be: q_hyp > p_hyp")

    tkt_ell = ellipsoid(p_ell,q_ell,theta)
    tkt_hyp = hyperboloid(p_hyp,q_hyp,theta)

    if verbose:
        print("distance:", distance)
        print("distance 2 (check):", q_ell - q_hyp)

        print("ell p,q", p_ell, q_ell)
        print("hyp p,q", p_hyp, q_hyp)

        print("position ell p,q", p_ell, distance / 2)
        print("position hyp p,q", distance / 2, p_hyp)

        print("length=", p_ell + distance + p_hyp)
        print("magnification ell: ", m_ell)
        print("magnification hyp: ", m_hyp)
        print("magnification: ", m_ell * m_hyp)
        print("demagnification: ", 1 / (m_ell * m_hyp))

    return tkt_ell,tkt_hyp

def recipe3(    #  common kick point
    p_ell = 100.0,
    q_ell = 10.0,
    p_hyp = 0.9,
    theta = 0.003,
    verbose = 1,
    method=0, # 0=Penelope, 1=Factory, 2=Ken
            ):

    q_hyp = q_ell

    if p_ell > 1e10: # use parabola
        if method == 0:
            tkt_ell = paraboloid(p_ell, q_ell, theta)
            ccc_ell = tkt_ell['ccc']
        elif method == 1:
            ccc_ell = fac_paraboloid(p_ell, q_ell, theta)
        elif method == 2:
            ccc_ell = ken_paraboloid(p_ell, q_ell, theta)

        f_ell = q_ell
    else:
        if method == 0:
            tkt_ell = ellipsoid(p_ell, q_ell, theta)
            ccc_ell = tkt_ell['ccc']
        elif method == 1:
            ccc_ell = fac_ellipsoid(p_ell, q_ell, theta)
        elif method == 2:
            ccc_ell = ken_ellipsoid(p_ell, q_ell, theta)

        f_ell = 1.0 / (1.0/p_ell + 1.0/q_ell)



    if method == 0:
        tkt_hyp = hyperboloid(p_hyp, q_hyp, theta)
        ccc_hyp = tkt_hyp['ccc']
    elif method == 1:
        ccc_hyp = fac_hyperboloid(p_hyp, q_hyp, theta)
    elif method == 2:
        ccc_hyp = ken_hyperboloid(p_hyp, q_hyp, theta)

    f_hyp = q_hyp / p_hyp

    # if p_hyp <= q_hyp:
    #     pass
    # else:
    #     raise Exception("must be: p_hyp <= q_hyp")


    if verbose:


        print("ell p,q", p_ell, q_ell)
        print("hyp p,q", p_hyp, q_hyp)

        print("length=", p_ell + p_hyp)
        print("f ell: ", f_ell)
        print("f hyp: ", f_hyp)
        print("f: ", f_ell * f_hyp)

    ccc_ell, ccc_hyp = numpy.array(ccc_ell), numpy.array(ccc_hyp)

    ccc_ell[numpy.argwhere(numpy.abs(ccc_ell) < 1e-14)] = 0
    ccc_hyp[numpy.argwhere(numpy.abs(ccc_hyp) < 1e-14)] = 0

    tkt_ell = {}
    tkt_ell['ccc'] = ccc_ell

    tkt_hyp = {}
    tkt_hyp['ccc'] = ccc_hyp

    return tkt_ell, tkt_hyp

def recipe4(    #  centered system parabola-hyperbola
    f11 = -0.00194644,
    f12 = 0.0,
    f21 = 1.905,
    f22 = 0.0,
    theta = 0.0159872,
    verbose = 1,
    method=0, # 0=Penelope, 1=Factory, 2=Ken
            ):

    if f12 != 0.0 or f22 != 0:
        raise Exception("Is your origin at the common focus?")

    p = numpy.abs(f11-f12)

    # intersection point at the parabola matching angle (https://doi.org/10.1107/S1600577522004593)
    c0 = (p / 2) / (numpy.tan(theta)) ** 2 - (p/2)
    # y^2 = 2px + p^2
    c1 = numpy.sqrt( 2 * p * c0 + p**2)

    ccc_centered_parabola = [1, 1, 0, 0, 0, 0, 0, 0, -2 * p, -p ** 2]  # normal incidence

    #
    # hyperbola
    #

    c = f21/2

    # get a from the hyperbola
    # (x-c)/a)^2 - (y/b)^2 = 1
    # b^2 = c^2 - a^2

    A = 1
    B = -c0**2 + 2 * c * c0 - 2 * c**2 - c1**2
    C = c**4 + c**2 * c0**2 - 2 * c**3 * c0
    S1 = numpy.sqrt( (-B + numpy.sqrt(B**2 - 4*A*C)) / (2 * A) )
    S2 = numpy.sqrt( (-B - numpy.sqrt(B ** 2 - 4 * A * C)) / (2 * A) )
    a = numpy.min((S1,S2))
    b = numpy.sqrt(c**2 - a**2)

    #
    # coeffs

    # Parabola
    # y^2 = p(2 * x + p) (Underwood)
    # x^2 + z^2 = 2 p y + p^2 (Shadow)
    ccc_centered_parabola = [1,1,0, 0,0,0,     0,0,-2*p, -p**2] # normal incidence
    ccc_centered_parabola = numpy.array(ccc_centered_parabola)

    # Hyperbola
    # (x-c)^2/a^2 - y^2/b^2 = 1 (Underwood)
    # (z-c)^2/a^2 - (y^2+x*2)/b^2 = 1 (Shadow)
    # # normal incidence (Underwood x->z, y->y  ->x)
    ccc_centered_hyperbola = [-1/b**2, -1/b**2, 1/a**2, \
                              0, 0, 0, \
                              0, 0, -2*c/a**2, \
                              (c/a)**2 - 1]
    ccc_centered_hyperbola = numpy.array(ccc_centered_hyperbola)


    if verbose:
        print("f11, f12", f11, f12)
        print("f21, f22", f21, f22)

        print("c0, c1: ", c0, c1)
        print("S1, S2: ", S1, S2)


        print("   theta grazing [deg]: ", theta * 180 / numpy.pi)
        print("   Parabola p [m]:", p)
        print("   Hyperbola a, b, c [m]: ", a, b, c)

        print("\n\nCalculated parameters: ")
        x_pmin = c0
        y_pmin = c1
        print("   ** Origin is at parabola focus (=far hyperbola focus)**")
        print("   Common point: x_pmin, y_pmin: ", x_pmin, y_pmin)
        print(" ")


        print("\n   normalized ccc_centered_parabola", ccc_centered_parabola / ccc_centered_parabola[0])
        print("   normalized ccc_centered_hyperbola", ccc_centered_hyperbola / ccc_centered_hyperbola[0])

        print("\n   ccc_centered_parabola", ccc_centered_parabola)
        print("   ccc_centered_hyperbola", ccc_centered_hyperbola)


    return  {'ccc':ccc_centered_parabola}, {'ccc':ccc_centered_hyperbola}

    # ccc_ell, ccc_hyp = numpy.array(ccc_ell), numpy.array(ccc_hyp)
    #
    # ccc_ell[numpy.argwhere(numpy.abs(ccc_ell) < 1e-14)] = 0
    # ccc_hyp[numpy.argwhere(numpy.abs(ccc_hyp) < 1e-14)] = 0
    #
    # tkt_ell = {}
    # tkt_ell['ccc'] = ccc_ell
    #
    # tkt_hyp = {}
    # tkt_hyp['ccc'] = ccc_hyp
    #
    # return tkt_ell, tkt_hyp
if __name__ == "__main__":

    if False:
        tkt_ell, tkt_hyp = recipe1(
            p_ell = 10.0,
            q_ell = 3.0,
            distance = 0.3,
            theta = 0.003,
            ratio_hyp = 3.0,
        )

        # print(tkt_ell)
        # print(tkt_hyp)

        print("\n\n>>>>>\n\n")
        # correct for incidence in the negative Y
        ccc1 = tkt_hyp['ccc']
        ccc2 = rotate_and_shift_quartic(ccc1, omega=0.0, theta=0.0, phi=numpy.pi, )
        print(ccc2)

        # """
        # [3703.7148148348147, 0.03333333333333333, 3703.7051501405344, 0.0, 11.965776068354531, 0.0, 0.0, 0.0, -102.56425641041795, 0.0]
        # [-45724.874256651594, -0.411522633747064, -45723.22816611661, 0.0, -548.6951988985766, 0.0, 0.0, 2.9816149549333204e-12, 740.7418518484149, -1.3322676295501878e-15]
        # """


    if False:
        # tkt_ell, tkt_hyp = recipe2(
        #     p_ell=10.0,
        #     distance=0.3,
        #     p_hyp=0.9,
        #     theta=0.003,
        #     # m_ell=0.3,
        #     m_hyp=1/3,
        #     verbose=1,
        # )

        tkt_ell, tkt_hyp = recipe2(
            p_ell=10.0,
            distance=0.5,
            p_hyp=2.5,
            theta=0.003,
            # m_ell=0.3,
            m_hyp=1/3,
            verbose=1,
        )

        print("\n\n>>>>>\n\n")
        # correct for incidence in the negative Y
        ccc1 = tkt_hyp['ccc']
        ccc2 = rotate_and_shift_quartic(ccc1, omega=0.0, theta=0.0, phi=numpy.pi, )
        print(ccc2)

        # """
        # [3703.7148148348147, 0.03333333333333333, 3703.7051501405344, 0.0, 11.965776068354531, 0.0, 0.0, 0.0, -102.56425641041795, 0.0]
        # [-45724.874256651594, -0.411522633747064, -45723.22816611661, 0.0, -548.6951988985766, 0.0, 0.0, 2.9816149549333204e-12, 740.7418518484149, -1.3322676295501878e-15]
        # """

    if False: # Underwood off-axis
        tkt_ell, tkt_hyp = recipe3( p_ell = 1e11,
                                    q_ell = 3.808047,
                                    p_hyp = 1.904995,
                                    theta = 0.0159872,
                                    verbose = 1,
                                    method = 0)
        ccc_ell = tkt_ell['ccc']
        ccc_hyp = tkt_hyp['ccc']

        print("\n\nNormalized Parabola: ",  ccc_ell/ccc_ell[0])
        print("Normalized Hyperbola: ", ccc_hyp/ccc_hyp[0])

        print("\n\nParabola: ",  ccc_ell)
        print("Hyperbola: ", ccc_hyp)


    if True: # Underwood centered
        tkt_ell, tkt_hyp = recipe4(verbose=1)
        ccc_ell = tkt_ell['ccc']
        ccc_hyp = tkt_hyp['ccc']

        # print("\n\nNormalized Parabola: ",  ccc_ell/ccc_ell[0])
        # print("Normalized Hyperbola: ", ccc_hyp/ccc_hyp[0])
        #
        # print("\n\nParabola: ",  ccc_ell)
        # print("Hyperbola: ", ccc_hyp)



