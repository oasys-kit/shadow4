import numpy
from minishadow.optical_surfaces.conic import Conic
from minishadow.optical_surfaces.conic import Conic
from srxraylib.plot.gol import plot_surface

from mesh import Mesh

if __name__ == "__main__":



    fmirr         = 1
    p             = 1000.0       # source-mirror
    q             = 300.0        # mirror-image
    alpha         = 0.0      # mirror orientation angle
    theta_grazing = 5e-3     # grazing angle, rad
    fcyl          = 0 # oe1.FCYL
    f_convex      = 0 # oe1.F_CONVEX



    print("fmirr = %s, p=%f, q=%f, alpha=%f, theta_grazing=%f rad, fcyl=%d"%\
              (fmirr,p,q,alpha,theta_grazing,fcyl))

    ccc = Conic()

    ccc.set_sphere_from_focal_distances(p,q,theta_grazing) # ,itype=fmirr,cylindrical=fcyl)
    print(ccc.info())

    x0 = numpy.array([0.2, -0.50, 0.5])
    v0 = numpy.array([0.0, 1.0, -1.0])

    t_exact = ccc.calculate_intercept(x0,v0,keep=0)

    print("exact solution: ",t_exact[0][0])

    #
    # mesh
    #
    mm = Mesh()
    mm.set_ray(x0,v0)
    mm.set_surface(ccc.z_vs_xy)

    x_start = 0.4

    t_mesh = mm.solve(x_start)

    print("mesh solution: ",t_mesh)
    print("   line: ",mm.line(t_mesh))
    print("   surface: ",mm.surface_vs_t(t_mesh))


    #
    # plot surface
    #


    # x = numpy.linspace(-0.25,0.25,20)
    # y = numpy.linspace(-0.5,0.5,100)
    # X = numpy.outer(x,numpy.ones_like(y))
    # Y = numpy.outer(numpy.ones_like(x),y)
    #
    # Z = ccc.z_vs_xy(X,Y)
    # print(">>>>>",Z.shape)
    # plot_surface(Z.real,x,y)