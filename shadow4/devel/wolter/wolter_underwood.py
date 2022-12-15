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

def draft1():
    theta = 0.916 * numpy.pi / 180

    #
    # parabola
    #
    p = 0.194644e-2
    # x^2 + z^2 = 2 p y + p^2
    # f = a = p/2
    y_pmin = p / numpy.tan(theta)
    print("theta, p, y_pmin: ", theta, p, y_pmin)

    f_eff = 190.5e-2
    y_pmin = f_eff * numpy.sin(4*theta)
    x_pmin = get_x_parabola(y_pmin, p_u=p)
    r_pmin = numpy.sqrt(x_pmin**2 + y_pmin**2)
    print("x_pmin, y_pmin, r_pmin: ", x_pmin, y_pmin, r_pmin)



    A = 14.84e-4
    y_pmax = numpy.sqrt(y_pmin**2 + A/numpy.pi)
    x_pmax = get_x_parabola(y_pmax, p_u=p)
    r_pmax = numpy.sqrt(x_pmax**2 + x_pmin**2)
    print("x_pmax, y_pmax, r_pmax: ", x_pmax, y_pmax, r_pmax)

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
    ccc_h = hyperboloid(q_hyp, p_hyp, theta)
    print(ccc_h)
    print("theta_new in deg: ", theta_new * 180 / numpy.pi)
    print("r: %20.15f" % r)
    print("f_eff", f_eff, r_pmin, 0.907256, 1./(1/r_pmin+1/0.907256))

    a_hyp = ccc_h['a'] # 0.951526
    b_hyp = ccc_h['b'] # 0.043058
    c_hyp = ccc_h['c']  # numpy.sqrt(a_hyp**2 + b_hyp**2)


    ccc_centered_parabola = [1,0,1, 0,-2*p,0,  0,0,0,    -p**2]
    ccc_centered_parabola = [1,1,0, 0,0,0,     0,0,-2*p, -p**2] # normal incidence
    print("ccc_centered_parabola", ccc_centered_parabola)

    # ccc_centered_hyperbola = [-1/b_hyp**2, 1/a_hyp**2, -1/b_hyp**2, \
    #                           0, 0, 0, \
    #                           0, -2*c_hyp/a_hyp, (c_hyp/a_hyp)**2 - 1]
    # normal incidence (Underwood x->z, y->y  ->x)
    ccc_centered_hyperbola = [-1/b_hyp**2, -1/b_hyp**2, 1/a_hyp**2, \
                              0, 0, 0, \
                              0, 0, -2*c_hyp/a_hyp**2, \
                              (c_hyp/a_hyp)**2 - 1]
    print("ccc_centered_hyperbola", ccc_centered_hyperbola)
    print("2*y_pmin=", 2*y_pmin)
    print("2*y_pmax=", 2 * y_pmax)

    print("x_pmin-2*c_hyp=", x_pmin - 2 * c_hyp, x_pmin + 2 * c_hyp)
    print("2*y_pmin=", 2 * y_pmin)
    print("2*c_hyp=", 2 * c_hyp)
    print("a_hyp+c_hyp=", a_hyp + c_hyp)

    # from shadow4.devel.wolter.conic_viewer import view_conic
    # view_conic(ccc_centered_hyperbola, x_min=-0.14, x_max=0.14, y_min=-0.14, y_max=0.14,)
    from srxraylib.plot.gol import plot

def centered_system(
    theta = 0.916 * numpy.pi / 180,  # grazing angle for both mirrors at kick point
    p = 0.194644e-2,                 # parabola p
    c = (190.5e-2 / 2),              # hyperbola c (=f_eff/2)
    ):

    #
    # parabola
    #

    # y^2 = p(2 * x + p) (Underwood)
    # x^2 + z^2 = 2 p y + p^2 (Shadow)
    ccc_centered_parabola = [1,1,0, 0,0,0,     0,0,-2*p, -p**2] # normal incidence


    y_pmin = p / numpy.tan(theta)
    x_pmin = get_x_parabola(y_pmin, p_u=p)
    r_pmin = numpy.sqrt(x_pmin**2 + y_pmin**2)

    #
    # hyperbola
    #

    p_hyp = r_pmin
    q_hyp = numpy.sqrt(y_pmin**2 + (x_pmin - 2*c)**2)

    ccc_h = hyperboloid(q_hyp, p_hyp, theta, verbose=False)
    #
    a_hyp = ccc_h['a'] # 0.951526
    b_hyp = ccc_h['b'] # 0.043058
    c_hyp = ccc_h['c']  # numpy.sqrt(a_hyp**2 + b_hyp**2)

    # (x-c)^2/a^2 - y^2/b^2 = 1 (Underwood)
    # (z-c)^2/a^2 - (y^2+x*2)/b^2 = 1 (Shadow)
    # # normal incidence (Underwood x->z, y->y  ->x)
    ccc_centered_hyperbola = [-1/b_hyp**2, -1/b_hyp**2, 1/a_hyp**2, \
                              0, 0, 0, \
                              0, 0, -2*c_hyp/a_hyp**2, \
                              (c_hyp/a_hyp)**2 - 1]

    print("\n\nInputs: ")
    print("   theta grazing [deg]: ", theta*180/numpy.pi)
    print("   Parabola p [m]:", p)
    print("   Hyperbola c [m]: ", c)


    print("\n\nCalculated parameters: ")
    print("   ** Origin is at parabola focus (=far hyperbola focus)**")
    print("   Kick point yp_min=y_hmax=p/tan(theta)", y_pmin)
    print("   Kick point: x_pmin, y_pmin, r_pmin: ", x_pmin, y_pmin, r_pmin)
    print("   hyperbola (input) p,q, theta, q/p: ", p_hyp, q_hyp, theta, q_hyp/p_hyp)
    print("   hyperbola (output) a,b,c: ", a_hyp, b_hyp, c_hyp)
    print(" ")

    ccc_centered_parabola = numpy.array(ccc_centered_parabola)
    ccc_centered_hyperbola = numpy.array(ccc_centered_hyperbola)
    print("\n   normalized ccc_centered_parabola", ccc_centered_parabola/ccc_centered_parabola[0])
    print("   normalized ccc_centered_hyperbola", ccc_centered_hyperbola/ccc_centered_hyperbola[0])

    print("\n   ccc_centered_parabola", ccc_centered_parabola)
    print("   ccc_centered_hyperbola", ccc_centered_hyperbola)


def noncentered_system(
    theta = 0.916 * numpy.pi / 180,  # grazing angle for both mirrors at kick point
    p = 0.194644e-2,                 # parabola p
    c = (190.5e-2 / 2),              # hyperbola c (=f_eff/2)):
    ):

    #
    # parabola
    #

    y_pmin = p / numpy.tan(theta)
    x_pmin = get_x_parabola(y_pmin, p_u=p)
    r_pmin = numpy.sqrt(x_pmin**2 + y_pmin**2)


    ccc_p = paraboloid(ssour=1e10, simag=r_pmin, theta_grazing=theta, verbose=True)
    #
    # hyperbola
    #
    p_hyp = r_pmin
    q_hyp = numpy.sqrt(y_pmin**2 + (x_pmin - 2*c)**2)
    ccc_h = hyperboloid(q_hyp, p_hyp, theta)
    print(ccc_h)

    print("\n\nInputs: ")
    print("   theta grazing [deg]: ", theta*180/numpy.pi)
    print("   Parabola p [m]:", p)
    print("   Hyperbola c [m]: ", c)


    print("\n\nCalculated parameters: ")
    print("\n\n   Normalized ccc_parabola: ", ccc_p['ccc']/ccc_p['ccc'][0])
    print("   Normalized ccc_hyperbola (for shadow, change sign of yz term!!!!!)): ", ccc_h['ccc']/ccc_h['ccc'][0])
    print("\n   ccc_parabola: ", ccc_p['ccc'])
    print("   ccc_hyperbola (for shadow, change sign of yz term!!!!!)): ", ccc_h['ccc'])



if __name__ == "__main__":
    # draft1()

    centered_system()

    # noncentered_system()


