import numpy
from conic_penelope import ellipsoid, hyperboloid, rotate_and_shift_quartic, euler_rotation_matrix

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
        # print(cyl(tkt_ell['ccc']))
        print(tkt_ell['ccc'])
        # correct for incidence in the negative Y
        ccc1 = tkt_hyp['ccc']
        ccc2 = rotate_and_shift_quartic(ccc1, omega=0.0, theta=0.0, phi=numpy.pi, )
        print(ccc2)

        # """
        # [3703.7148148348147, 0.03333333333333333, 3703.7051501405344, 0.0, 11.965776068354531, 0.0, 0.0, 0.0, -102.56425641041795, 0.0]
        # [-45724.874256651594, -0.411522633747064, -45723.22816611661, 0.0, -548.6951988985766, 0.0, 0.0, 2.9816149549333204e-12, 740.7418518484149, -1.3322676295501878e-15]
        # """


    if True:
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
        # print(tkt_ell)
        # print(tkt_hyp)

        print("\n\n>>>>>\n\n")
        # print(cyl(tkt_ell['ccc']))
        print(tkt_ell['ccc'])
        # correct for incidence in the negative Y
        ccc1 = tkt_hyp['ccc']
        ccc2 = rotate_and_shift_quartic(ccc1, omega=0.0, theta=0.0, phi=numpy.pi, )
        print(ccc2)

        # """
        # [3703.7148148348147, 0.03333333333333333, 3703.7051501405344, 0.0, 11.965776068354531, 0.0, 0.0, 0.0, -102.56425641041795, 0.0]
        # [-45724.874256651594, -0.411522633747064, -45723.22816611661, 0.0, -548.6951988985766, 0.0, 0.0, 2.9816149549333204e-12, 740.7418518484149, -1.3322676295501878e-15]
        # """


