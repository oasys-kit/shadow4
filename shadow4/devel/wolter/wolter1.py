import numpy
from conic_penelope import ellipsoid, hyperboloid, rotate_and_shift_quartic, euler_rotation_matrix

def cyl(ccc):
    ccc1 = ccc.copy()
    ccc1[0] = 0
    ccc1[3] = 0
    ccc1[5] = 0
    ccc1[6] = 0
    return ccc1


if __name__ == "__main__":
    p_ell = 10.0
    q_ell = 3.0

    distance = 0.3
    theta = 0.003
    ratio_hyp = 3.0 # q_hyp / p_ell > 1.0



    q_hyp = q_ell - distance
    a_hyp = 0.5 * q_hyp * (1 - 1/ratio_hyp)
    p_hyp = q_hyp - 2 * a_hyp

    tkt_ell = ellipsoid(p_ell,q_ell,theta)
    tkt_hyp = hyperboloid(p_hyp,q_hyp,theta)


    print(tkt_ell)
    print(tkt_hyp)

    print(">>>>>\n\n")
    # print(cyl(tkt_ell['ccc']))
    print(tkt_ell['ccc'])


    ccc1 = tkt_hyp['ccc']
    ccc2 = rotate_and_shift_quartic(ccc1, omega=0.0, theta=0.0, phi=numpy.pi,)
    print(ccc2)

    # ccc3 = cyl(ccc1)
    # print(ccc3)

    print("ell p,q", p_ell,q_ell)
    print("hyp p,q", p_hyp, q_hyp)

    print("position ell p,q", p_ell,distance/2)
    print("position hyp p,q", distance/2, p_hyp)

    print("length=", p_ell + distance + p_hyp)
    m_ell = (q_ell/p_ell)
    m_hyp = (p_hyp/q_hyp)
    print("magnification ell: ", m_ell )
    print("magnification hyp: ", m_hyp)
    print("magnification: ", m_ell * m_hyp)
    print("demagnification: ", 1/(m_ell * m_hyp))