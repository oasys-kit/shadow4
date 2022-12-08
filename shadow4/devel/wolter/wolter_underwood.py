import numpy
from shadow4.devel.wolter.conic_penelope import ellipsoid, hyperboloid, rotate_and_shift_quartic, euler_rotation_matrix
from conic_penelope import ellipsoid, hyperboloid, paraboloid, rotate_and_shift_quartic, euler_rotation_matrix

def cyl(ccc):
    ccc1 = ccc.copy()
    ccc1[0] = 0
    ccc1[3] = 0
    ccc1[5] = 0
    ccc1[6] = 0
    return ccc1

def get_x_parabola(y, p_u=1):
    return (y**2 - p_u**2) / 2 / p_u

if __name__ == "__main__":

    theta = 0.916 * numpy.pi / 180

    #
    # parabola
    #
    p = 0.194644e-2
    # x^2 + z^2 = 2 p y + p^2
    # f = a = p/2
    f_parabola = p/2


    y_pmin = p / numpy.tan(theta)
    print("theta, p, y_pmin: ", theta, p, y_pmin)

    f_eff = 190.5e-2
    y_pmin = f_eff * numpy.sin(4*theta)
    x_pmin = get_x_parabola(y_pmin, p_u=p)
    r_pmin = numpy.sqrt(x_pmin**2 + y_pmin**2)
    print("x_pmin, y_pmin, r_pmin: ", x_pmin, y_pmin, r_pmin)

    ccc_p = paraboloid(ssour=1e10, simag=r_pmin, theta_grazing=theta, verbose=True)


    #
    # hyperbola
    #
    print("\n\n\n")
    p_hyp = r_pmin
    c = f_eff / 2
    r = numpy.sqrt(y_pmin**2 + (x_pmin - 2*c)**2)
    q_hyp = r
    theta_new = 0.5 * numpy.arccos( ((2 * c)**2 - r_pmin**2 - r**2) / (-2 * r * r_pmin))
    print("theta, theta_new: ", theta, theta_new)
    # ccc_h = hyperboloid(p_hyp, q_hyp, theta_new)
    ccc_h = hyperboloid(q_hyp, p_hyp, theta_new)
    print("theta_new in deg: ", theta_new * 180 / numpy.pi)
    print("r: %20.15f" % r)
    print("f_eff", f_eff, r_pmin, 0.907256, 1./(1/r_pmin+1/0.907256))

    a_hyp = 0.951526
    b_hyp = 0.043058
    c_hyp = numpy.sqrt(a_hyp**2 + b_hyp**2)


    ccc_centered_parabola = [1,0,1, 0,-2*p,0,  0,0,0, -p**2]
    ccc_centered_parabola = [1,1,0, 0,0,0,  0,0,-2*p, -p**2] # normal incidence
    print("ccc_centered_parabola", ccc_centered_parabola)

    ccc_centered_hyperbola = [-1/b_hyp**2, 1/a_hyp**2, -1/b_hyp**2, \
                              0, 0, 0, \
                              0, -2*c_hyp/a_hyp, (c_hyp/a_hyp)**2 - 1]
    print("ccc_centered_hyperbola", ccc_centered_hyperbola)