
from shadow4.optical_surfaces.s4_conic import S4Conic # for comparison
from shadow4.devel.wolter.conic_viewer import view_conic, compare_conics # for plot

import numpy

from numpy import sin as Sin
from numpy import cos as Cos
from numpy import tan as Tan
from numpy import sqrt as Sqrt

def Cot(x):
    return 1/Tan(x)

def Sec(x):
    return 1/Cos(x)

def Csc(x):
    return 1/Sin(x)

# see conics_penelope_paraboloid_focusing.nb
def paraboloid(p=1e10, q=10,theta=3e-3):
    if p > q:
        return paraboloid_focusing(q=q, theta=theta)
    else:
        return paraboloid_collimating(p=p, theta=theta)

def paraboloid_focusing(q=10,theta=3e-3):
    return [1, Sin(theta)**2, Cos(theta)**2, 0, 2*Cos(theta)*Sin(theta), 0, 0, 0, -4*q*Sin(theta),0 ]

# see conics_penelope_paraboloid_collimating.nb
def paraboloid_collimating(p=10,theta=3e-3):
    return [1, Sin(theta) ** 2, Cos(theta) ** 2, 0, -2 * Cos(theta) * Sin(theta), 0, 0, 0, -4 * p * Sin(theta), 0]

# see conics_penelope_ellipsoid.nb
def ellipsoid(p=10,q=3,theta=3e-3):

    return [Csc(theta)**2/(p*q),
            1/(p*q),
            (-(p - q)**2 + (p + q)**2*Csc(theta)**2)/(p*q*(p + q)**2),
            0,
    (2*(p - q)*Sqrt(((p + q)**2*Cos(theta)**2)/(p**2 + q**2 + 2*p*q*Cos(2*theta)))*Sqrt(p**2 + q**2 + 2*p*q*Cos(2*theta))*
        Sqrt(Csc(theta)**2/(p + q)**2))/(p*q*(p + q)),
            0,
            0,
    (4*(p - q)*(-(Sqrt(p*q)*Sqrt((p*q*Cos(theta)**2)/(p**2 + q**2 + 2*p*q*Cos(2*theta)))*Csc(theta)) +
        p*q*Sqrt(((p + q)**2*Cos(theta)**2)/(p**2 + q**2 + 2*p*q*Cos(2*theta)))*Sqrt(Csc(theta)**2/(p + q)**2)))/
        (p*q*(p + q)*Sqrt(p**2 + q**2 + 2*p*q*Cos(2*theta))*Sqrt(Csc(theta)**2/(p + q)**2)),
    -(Sqrt(Csc(theta)**2/(p + q)**2)*(-2*(p**2 - q**2)**2*Cot(theta)**2 +
    Csc(theta)**2*((p - q)**2*(p**2 + 6*p*q + q**2) + (p - q)**4*Cos(2*theta) +
       (8*(p*q)**1.5*Cos(theta)*Sqrt(((p + q)**2*Cos(theta)**2)/(p**2 + q**2 + 2*p*q*Cos(2*theta)))*Cot(theta))/
        (Sqrt((p*q*Cos(theta)**2)/(p**2 + q**2 + 2*p*q*Cos(2*theta)))*Sqrt(Csc(theta)**2/(p + q)**2))))*Sin(theta)**2)/
    (2.*p*q*(p**2 + q**2 + 2*p*q*Cos(2*theta))),
            0]

def hyperboloid(p=10,q=3,theta=3e-3):
    if p >= q:
        return hyperboloid_large_p(p=p,q=q,theta=theta)
    else:
        return hyperboloid_large_q(p=p, q=q, theta=theta)

def hyperboloid_large_p(p=10,q=3,theta=3e-3):
    if p < q:
        raise Exception("p<q")
    return [
        -(Csc(theta)**2/(p*q)),-(((p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))*Csc(theta)**2)/
        (p*q*(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2))),
        (4*(p + q)**2 - ((p - q)**4*Csc(theta)**2*(1 + Csc(theta)**2))/(p*q))/((p - q)**2*(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)),0,
        (-2*(p + q)*(p**2 + q**2 - 2*p*q*Cos(2*theta))*Csc(theta)**2*Sqrt(1/(1 + (p + q)**2/((p - q)**2*(1 + Csc(theta)**2)))))/
        (p*(p - q)**3*q*Sqrt(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))*
        Sqrt((2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)/((p - q)**2*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))),0,0,
        -2*((-2*(p + q)*Sqrt(1/(1 + (p + q)**2/((p - q)**2*(1 + Csc(theta)**2)))))/
        ((p - q)*Sqrt(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))) +
        (Sqrt(2)*(-p**2 + q**2)*Sqrt(-((p*q*(-3 + Cos(2*theta)))/(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))*
        Sqrt(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))*
        Sqrt((2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)/((p - q)**2*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))*Sin(theta))/
        (Sqrt(p*q)*(-2*(p**2 - p*q + q**2) + (p**2 + q**2)*Cos(2*theta)))),
        -((((p + q)**2*(p**2 + q**2 - 2*p*q*Cos(2*theta))*Csc(theta)*Sqrt(1/(2 + (2*(p + q)**2)/((p - q)**2*(1 + Csc(theta)**2))))*
        (4*Sqrt(p*q)*Sqrt(-((p*q*(-3 + Cos(2*theta)))/(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta)))) +
        (p - q)**2*Sqrt(((p - q)**2*(-3 + Cos(2*theta)))/(-2*(p**2 - p*q + q**2) + (p**2 + q**2)*Cos(2*theta)))*Csc(theta)*
        Sqrt((2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)/((p - q)**2*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))))/(p*q) +
        2*(4*(p + q)**2 - ((p - q)**4*Csc(theta)**2*(1 + Csc(theta)**2))/(p*q))*
        ((p + q)**2/(2.*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))*
        Sqrt((2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)/((p - q)**2*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))) -
        Sqrt(2)*Sqrt(p*q)*Sqrt(-((p*q*(-3 + Cos(2*theta)))/(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))*
        Sqrt(1/(1 + (p + q)**2/((p - q)**2*(1 + Csc(theta)**2))))*Sin(theta)))/
        ((p - q)**2*(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2))),0]

def hyperboloid_large_q(p=3,q=10,theta=3e-3):
    if p > q:
        raise Exception("p>q")
    return [
        -(Csc(theta)**2/(p*q)),-(((p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))*Csc(theta)**2)/
        (p*q*(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2))),
        (4*(p + q)**2 - ((p - q)**4*Csc(theta)**2*(1 + Csc(theta)**2))/(p*q))/((p - q)**2*(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)),0,
        (2*(p**2 + q**2 - 2*p*q*Cos(2*theta))*Sqrt(-((p*q*(-3 + Cos(2*theta)))/(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))*Csc(theta)**3*
        Sqrt((p + q)**2/(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)))/
        ((p - q)**2*(p*q)**1.5*Sqrt((4*(p**2 + q**2) + 2*(p - q)**2*Csc(theta)**2)/
        ((p - q)**2*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))),0,0,
        -2*(-(Csc(theta)*(-4*p*q*Sqrt((p + q)**2/(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)) +
        ((p - q)*(p + q)*Csc(theta)**2)/
        (Sqrt(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))*
        Sqrt((2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)/((p - q)**2*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))))*
        (1 + Sin(theta)**2))/
        (2.*Sqrt(p*q)*(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)*
        Sqrt((p*q*(1 + Sin(theta)**2))/(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta)))) +
        ((p**2 + q**2 - 2*p*q*Cos(2*theta))*Sqrt(-((p*q*(-3 + Cos(2*theta)))/(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))*Csc(theta)**3*
        Sqrt((p + q)**2/(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2))*
        ((-2*p*q*(-3 + Cos(2*theta)))/
        ((p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))*
        Sqrt((2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)/((p - q)**2*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))) +
        ((p - q)*(p + q)*Sqrt(((p + q)**2*Sin(theta)**2)/((p - q)**2 + 2*(p**2 + q**2)*Sin(theta)**2)))/
        Sqrt(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))/
        (2.*(p - q)**2*(p*q)**1.5*Sqrt((4*(p**2 + q**2) + 2*(p - q)**2*Csc(theta)**2)/
        ((p - q)**2*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta)))))),
        -2*((2*(p + q)*((p + q)**2/(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2))**1.5)/
        ((p - q)*Sqrt(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))) +
        ((p - q)*(-3 + Cos(2*theta))*(-2*(p**2 - p*q + q**2) + (p**2 + q**2)*Cos(2*theta))*Csc(theta)**4*
        (p*(Sqrt((p + q)**2/(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)) -
        Sqrt(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))*
        Sqrt((2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)/((p - q)**2*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))))) +
        q*(Sqrt((p + q)**2/(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)) +
        Sqrt(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))*
        Sqrt((2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)/((p - q)**2*(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta)))))))/
        (Sqrt(p**2 - 4*p*q + q**2 + 2*p*q*Cos(2*theta))*(2*(p**2 + q**2) + (p - q)**2*Csc(theta)**2)**2)),0]


#
# Ken
#


def ken_paraboloid_focusing(q=10,theta=3e-3):
    return [1, Sin(theta)**2, Cos(theta)**2, 0, 2*Cos(theta)*Sin(theta), 0, 0, 0, -4*q*Sin(theta),0 ]

# see conics_penelope_paraboloid_collimating.nb
def ken_paraboloid_collimating(p=10,theta=3e-3):
    ccc = ken_paraboloid_focusing(q=p,theta=theta)
    ccc[4] = -ccc[4]
    return ccc

def ken_hyperboloid_large_q_old(p=3,q=10,theta=3e-3):
    c = Cos(theta)
    s = Sin(theta)
    return [
        0,
        -0.25 * s**2 * (q-p)**2,
        p*q - 0.25 * c**2 * (p+q)**2,
        0,
        0.5 * s * (q-p) * c * (p+q),
        0,
        0,
        0,
        s * (q-p) * p * q,
        0
    ]

def ken_ellipsoid(p=3,q=10,theta=3e-3):
    c = Cos(theta)
    s = Sin(theta)
    h = (p - q) * c
    return [
        (p+q)**2,
        (p+q)**2 * s**2,
        h**2 + 4 * p * q,
        0,
        2 * s * (p + q) * h,
        0,
        0,
        0,
        -4 * s * (p + q)  * p * q,
        0
    ]

def ken_hyperboloid(p=3,q=10,theta=3e-3):
    c = Cos(theta)
    s = Sin(theta)
    return [
        1,
        s**2,
        c**2 - 4 * p *q * s**2 / (q - p)**2,
        0,
        - 2 * s * c * (p + q) / (q - p),
        0,
        0,
        0,
        -4 * s * p * q  / (q - p),
        0
    ]

def ken_hyperboloid_large_q(p=3,q=10,theta=3e-3):
    return ken_hyperboloid(p,q,theta)

def ken_hyperboloid_large_p(p=3,q=10,theta=3e-3):
    return ken_hyperboloid(p,q,theta)

def ken_hyperboloid_large_p_old(p=3,q=10,theta=3e-3):
    return ken_hyperboloid_large_q_old(p,q,theta)

#
# tools
#
def cylinder(c_in):
    c_out = c_in.copy()
    c_out[0] = 0.0
    c_out[3] = 0.0
    c_out[5] = 0.0
    c_out[6] = 0.0
    return c_out

def normalize(c_in, index=0, clean=True):
    c_out = [0] * 10
    for i in range(10):
        c_out[i] = c_in[i] / c_in[index]
        if clean:
            if numpy.abs(c_out[i]) < 1e-15:
                c_out[i] = 0.0
    return c_out


#
# checks
#
#
#
#
def ellipsoid_check(ssour=10,simag=3,theta_grazing=3e-3, do_plot=False):


    ccc = S4Conic.initialize_as_ellipsoid_from_focal_distances(ssour, simag, theta_grazing,
                                        cylindrical=0, cylangle=0.0, switch_convexity=0)
    print("ccc: ", ccc.get_coefficients())

    s5 = ellipsoid(ssour,simag,theta_grazing)

    c = ccc.get_coefficients()
    print("ccc: ", c)
    print("s5: ", s5)


    for i in range(10):
        print(i, c[i], s5[i])
        assert(numpy.abs(s5[i] - c[i]) < 1e-2)

    # view_conic(s5, x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1)
    if do_plot:
        compare_conics(s5, ccc.get_coefficients(), x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1,
                       titles=['s5','ccc'])

def parabola_check(ssour=10,simag=10,theta_grazing=3e-3, do_plot=False):


    ccc = S4Conic.initialize_as_paraboloid_from_focal_distances(ssour, simag, theta_grazing,
                                        cylindrical=0, cylangle=0.0, switch_convexity=0)
    print("ccc: ", ccc.get_coefficients())

    if simag < ssour:
        s5 = paraboloid_focusing(simag,theta_grazing)
    else:
        s5 = paraboloid_collimating(ssour, theta_grazing)


    print("ccc: ", ccc.get_coefficients())
    print("s5: ", s5)

    c = ccc.get_coefficients()
    for i in range(10):
        print(i, c[i], s5[i])

    for i in range(10):
        print(i, s5[i] , c[i])
        assert(numpy.abs(s5[i] - c[i]) < 1e-2)

    # view_conic(s5, x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1)
    if do_plot:
        compare_conics(s5, ccc.get_coefficients(), x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1,
                       titles=['s5','ccc'])

def hyperbola_check(ssour=10,simag=3,theta_grazing=3e-3, do_plot=False):

    if ssour > simag:
        s5 = hyperboloid_large_p(ssour, simag,theta_grazing)
        c = [-3703.714814834814, -0.033331264034642545, -3703.599850917718, 0, -41.26934605862947, 0, 0, 2.220446049250313e-16, -190.48238875690424, 0]

    else:
        s5 = hyperboloid_large_q(ssour, simag,theta_grazing)
        c = [-3703.714814834814, -0.033331264034642545, -3703.599850917718, 0, 41.26934605862948, 0, 0, 2.220446049250313e-16, 190.48238875690424, 0]


    print("ccc: ", c)
    print("s5: ", s5)

    for i in range(10):
        print(i, c[i], s5[i])

    for i in range(10):
        print(i, s5[i] , c[i])
        assert(numpy.abs(s5[i] - c[i]) < 1e-2)

    # view_conic(s5, x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1)
    if do_plot:
        compare_conics(s5, ccc.get_coefficients(), x_min=-0.01, x_max=0.01, y_min=-0.1, y_max=0.1,
                       titles=['s5','ccc'])


if __name__ == "__main__":
    # print(ellipsoid(p=10,q=3,theta=3e-3))
    # print(paraboloid_focusing(q=10,theta=3e-3))
    # print(paraboloid_collimating(p=10, theta=3e-3))
    # print(hyperboloid_large_p(p=10, q=3, theta=3e-3))
    # print(hyperboloid_large_q(p=3, q=10, theta=3e-3))
    #
    # ellipsoid_check(ssour=10,simag=3,theta_grazing=3e-3, do_plot=False)
    # parabola_check(ssour=10e10,simag=10,theta_grazing=3e-3, do_plot=False)
    # parabola_check(ssour=10, simag=10e10, theta_grazing=3e-3, do_plot=False)
    # hyperbola_check(ssour=10, simag=3, theta_grazing=3e-3, do_plot=False)
    # hyperbola_check(ssour=3, simag=10, theta_grazing=3e-3, do_plot=False)

    # print("cylinder, p<q:")
    # print(normalize(cylinder(hyperboloid_large_q(p=3, q=10, theta=3e-3)), index=2))
    # print(normalize(ken_hyperboloid_large_q_old(p=3, q=10, theta=3e-3),index=2))
    #
    # print("cylinder, p>q:")
    # print(normalize(cylinder(hyperboloid_large_p(p=10, q=3, theta=3e-3)), index=2))
    # print(normalize(ken_hyperboloid_large_p_old(p=10, q=3, theta=3e-3),index=2))

    print("hyperboloid, p<q:")
    print(normalize(hyperboloid_large_q(p=3, q=10, theta=3e-3), index=0))
    print(ken_hyperboloid_large_q(p=3, q=10, theta=3e-3))

    print("hyperboloid, p>q:")
    print(normalize(hyperboloid_large_p(p=10, q=3, theta=3e-3), index=0))
    print(ken_hyperboloid_large_p(p=10, q=3, theta=3e-3))

    print("ellipsoid")
    print(normalize(ellipsoid(p=10, q=3, theta=3e-3), index=0))
    print(normalize(ken_ellipsoid(p=10, q=3, theta=3e-3), index=0))

    print("parabola")
    print(normalize(paraboloid_focusing(q=10, theta=3e-3), index=0))
    print(normalize(ken_paraboloid_focusing(q=10, theta=3e-3), index=0))