
# https://stackoverflow.com/questions/13360062/python-curves-intersection-with-fsolve-and-function-arguments-using-numpy

# https://www.reddit.com/r/matlab/comments/pd7rr/finding_the_point_of_intersection_between_a_line/

from scipy.optimize import fsolve
from srxraylib.plot.gol import plot,plot_image, plot_surface


import numpy



def plot_surface_and_line(Z,x,y,zz,xx,yy,show=True):

    import matplotlib.pylab as plt
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter

    #
    # plot
    #

    X = numpy.outer(x,numpy.ones_like(y))
    Y = numpy.outer(numpy.ones_like(x),y)

    fig = plt.figure(figsize=None)
    ax = fig.gca(projection='3d')
    #

    cmap = cm.coolwarm

    ax = fig.gca(projection='3d')
    ax.plot(xx, yy, zz, label='parametric curve')
    ax.set_xlim(y.min(),y.max())
    ax.set_ylim(y.min(),y.max())

    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cmap, linewidth=0, antialiased=False)

    ax.set_zlim(zz.min(),zz.max())
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.title("")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ax.legend()

    if show:
        plt.show()

    return fig

class Mesh(object):
    def __init__(self,surface=None):
        self.__x0 = [0.0,0.0,0.0]
        self.__v0 = [0.0,0.0,0.0]
        self.surface = surface


    def line(self,t):
        return (self.__x0[0] + self.__v0[0] * t,
                self.__x0[1] + self.__v0[1] * t,
                self.__x0[2] + self.__v0[2] * t)

    def line_z(self,t):
        return self.__x0[2] + self.__v0[2] * t

    def surface_z_vs_t(self,t):
        return self.surface_vs_t(t)[2]

    def surface_vs_t(self,t):
        x1 = self.__x0[0] + self.__v0[0] * t
        y1 = self.__x0[1] + self.__v0[1] * t
        return x1,y1,self.surface(x1,y1)

    def solve(self,x_start):
        t_solution = fsolve(lambda t:self.surface_z_vs_t(t)-self.line_z(t), x_start)
        return t_solution[0]

    def set_ray(self,x0,v0):
        self.__x0 = x0
        self.__v0 = v0

    def set_surface(self,surface):
        self.surface = surface

def sphere(x, y, radius=5.0):
        return radius - numpy.sqrt(radius**2 - x**2 - y**2)

if __name__ == "__main__":
    x = numpy.linspace(-0.25,0.25,20)
    y = numpy.linspace(-0.5,0.5,100)
    X = numpy.outer(x,numpy.ones_like(y))
    Y = numpy.outer(numpy.ones_like(x),y)



    Z = sphere(X,Y)

    x0 = [0.2, y.min(), numpy.abs(y.min())]
    v0 = [0.0, 1.0, -1.0]

    t = numpy.linspace(0,y.max()-y.min(),100)
    xx = x0[0] + v0[0] * t
    yy = x0[1] + v0[1] * t
    zz = x0[2] + v0[2] * t

    # plot_surface_and_line(Z,x,y,zz,xx,yy)


    mm = Mesh()
    mm.set_ray(x0,v0)
    mm.set_surface(sphere)

    x_start = 0.4

    t_solution = mm.solve(x_start)

    print("t_solution: ",t_solution)
    print("line: ",mm.line(t_solution))
    print("surface: ",mm.surface_vs_t(t_solution))