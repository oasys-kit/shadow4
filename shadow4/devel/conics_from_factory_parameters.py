
from shadow4.optical_surfaces.s4_conic import S4Conic # for comparison
from conic_viewer import view_conic, compare_conics # for plot

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
def paraboloid_focusing(q=10,theta=3e-3):
    return [1, Sin(theta)**2, Cos(theta)**2, 0, 2*Cos(theta)*Sin(theta), 0, 0, 0, -4*q*Sin(theta),0 ]

# see conics_penelope_paraboloid_collimating.nb
def paraboloid_collimating(p=10,theta=3e-3):
    return [1, Sin(theta) ** 2, Cos(theta) ** 2, 0, -2 * Cos(theta) * Sin(theta), 0, 0, 0, -4 * p * Sin(theta), 0]

# see conics_penelope_ellipsoid.nb
def ellipsoid(p=10,q=3,theta=3e-3):

    return [Csc(theta)**2/(p*q),1/(p*q),
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


if __name__ == "__main__":
    print(ellipsoid(p=10,q=3,theta=3e-3))
    print(paraboloid_focusing(q=10,theta=3e-3))
    print(paraboloid_collimating(p=10, theta=3e-3))

    ellipsoid_check(ssour=10,simag=3,theta_grazing=3e-3, do_plot=False)
    parabola_check(ssour=10e10,simag=10,theta_grazing=3e-3, do_plot=False)
    parabola_check(ssour=10, simag=10e10, theta_grazing=3e-3, do_plot=False)