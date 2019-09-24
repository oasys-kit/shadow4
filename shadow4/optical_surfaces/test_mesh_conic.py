import numpy

from shadow4.optical_surfaces.conic import Conic
from srxraylib.plot.gol import plot_surface


from mesh import Mesh


def write_shadow_surface(s,xx,yy,filename='presurface.dat'):
    """
      write_shadowSurface: writes a mesh in the SHADOW/presurface format
      SYNTAX:
           out = write_shadowSurface(z,x,y,filename=filename)
      INPUTS:
           z - 2D array of heights z(x,y)
           x - 1D array of spatial coordinates along mirror width.
           y - 1D array of spatial coordinates along mirror length.

      OUTPUTS:
           filename - output file in SHADOW format. If undefined, the
                     file is names "presurface.dat"

    """

    try:
       fs = open(filename, 'w')
    except IOError:
       out = 0
       print ("Error: can\'t open file: "+filename)
       return
    else:
        # dimensions
        fs.write( "%d  %d \n"%(xx.size,yy.size))
        # y array
        for i in range(yy.size):
            fs.write("%20.18g  "%(yy[i]))
        fs.write("\n")
        # for each x element, the x value followed by the corresponding z(y) profile
        for i in range(xx.size):
            tmps = ""
            for j in range(yy.size):
                tmps += "%20.18g  "%(s[i,j])
            fs.write("%20.18g    %s \n"%(xx[i],tmps))
        fs.close()
        print ("write_shadow_surface: File for SHADOW "+filename+" written to disk.")




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

    print("exact solution: ",t_exact)#[0][0])

    #
    # mesh object using exact surface
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


    x = numpy.linspace(-0.25,0.25,20)
    y = numpy.linspace(-0.5,0.5,100)
    X = numpy.outer(x,numpy.ones_like(y))
    Y = numpy.outer(numpy.ones_like(x),y)

    Z = ccc.z_vs_xy(X,Y)

    # plot_surface(Z.real,x,y)


    #
    # create surface mesh and load it
    #
    write_shadow_surface(Z,x,y,filename="test_mesh_conic.dat")

    mm.load_file("test_mesh_conic.dat")

    x_start = 0.4

    t_mesh_file = mm.solve(x_start)

    print("mesh file solution: ",t_mesh_file)
    print("   line: ",mm.line(t_mesh_file))
    print("   surface: ",mm.surface_vs_t(t_mesh_file))

    assert (numpy.abs(t_mesh - t_exact) < 1e-4)
    assert (numpy.abs(t_mesh_file - t_mesh) < 1e-4)