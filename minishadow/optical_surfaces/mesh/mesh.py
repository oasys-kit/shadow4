
# https://stackoverflow.com/questions/13360062/python-curves-intersection-with-fsolve-and-function-arguments-using-numpy

# https://www.reddit.com/r/matlab/comments/pd7rr/finding_the_point_of_intersection_between_a_line/

from scipy.optimize import fsolve
from scipy import interpolate
from srxraylib.plot.gol import plot,plot_image, plot_surface, plot_scatter
import sys

import numpy


#
# rays[index,column]
# vector[xyz,index]
#


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

    def load_file(self,filename,kind='cubic'):
        x,y,z = self.read_surface_error_file(filename)
        self.surface = interpolate.interp2d(x,y,z.T, kind=kind)

    #copied from shadowOui util/shadow_util
    @classmethod
    def read_surface_error_file(cls, filename):

        file = open(filename, "r")

        rows = file.readlines()

        dimensions = rows[0].split()
        n_x = int(dimensions[0])
        n_y = int(dimensions[1])

        if n_x > 500:
            raise Exception("Malformed file: maximum allowed point in X direction is 500")

        x_coords = numpy.zeros(0)
        y_coords = numpy.zeros(0)
        z_values = numpy.zeros((n_x, n_y))


        index = 1
        dim_y_row = len(rows[index].split())
        is_ycoord = True
        first_x_row_index = 0

        while(is_ycoord):
            y_values = rows[index].split()

            if len(y_values) == dim_y_row:
                for y_value in y_values:
                    y_coords = numpy.append(y_coords, float(y_value))
            else:
                first_x_row_index = index
                is_ycoord = False

            index +=1

        first_x_row = rows[first_x_row_index].split()

        if len(first_x_row) == 2:
            x_index = 0
            z_index = 0

            for index in range(first_x_row_index, len(rows)):
                if z_index == 0:
                    values = rows[index].split()
                    x_coords = numpy.append(x_coords, float(values[0]))
                    z_value = float(values[1])
                else:
                    z_value = float(rows[index])

                z_values[x_index, z_index] = z_value
                z_index += 1

                if z_index == n_y:
                    x_index += 1
                    z_index = 0
        else:
            x_rows = []

            for index in range(2, len(rows)):

                x_row = rows[index].split("\t")

                if len(x_row) != 1 + n_y:
                    x_row = rows[index].split()

                if len(x_row) != 1 + n_y:
                    raise Exception("Malformed file: check format")

                x_rows.append(x_row)

            for x_index in range(0, len(x_rows)):
                x_coords = numpy.append(x_coords, float(x_rows[x_index][0]))

                for z_index in range(0, len(x_rows[x_index]) - 1):
                    z_value = float(x_rows[x_index][z_index + 1])

                    z_values[x_index, z_index] = z_value

        return x_coords, y_coords, z_values

    def get_normal(self,x2):
        # ;
        # ; Calculates the normal at intercept points x2 [see shadow's normal.F]
        # ;

        normal = numpy.zeros_like(x2)

        eps = 100 * sys.float_info.epsilon

        X_0 = x2[0,:]
        Y_0 = x2[1,:]
        Z_0 = x2[2,:]

        N_0 = numpy.zeros_like(X_0)
        N_1 = numpy.zeros_like(X_0)
        N_2 = numpy.ones_like(X_0)
        for i in range(X_0.size):
            z00 = self.surface(X_0[i],Y_0[i])
            N_0[i] = -1.0 * (self.surface(X_0[i]+eps,Y_0[i]) - z00) / eps
            N_1[i] = -1.0 * (self.surface(X_0[i],Y_0[i]+eps) - z00) / eps



        n2 = numpy.sqrt(N_0**2 + N_1**2 + N_2**2)
        #
        normal[0,:] = N_0 / n2
        normal[1,:] = N_1 / n2
        normal[2,:] = N_2 / n2

        return normal

    def calculate_intercept(self,XIN,VIN,keep=0):

        npoints = XIN.shape[1]
        answer = numpy.zeros(npoints)
        i_flag = numpy.ones(npoints)


        print(">>>>> main loop to find solutions (slow...)")
        for i in range(npoints):
            self.__x0 = XIN[:,i]
            self.__v0 = VIN[:,i]
            try:
                t_solution = self.solve(0.0)
                answer[i] = t_solution
            except:
                i_flag[i] = -1
        print(">>>>> done main loop to find solutions (Thanks for waiting!)")
        return answer,i_flag

    def apply_specular_reflection_on_beam(self,newbeam):
        # ;
        # ; TRACING...
        # ;

        x1 =   newbeam.get_columns([1,2,3])
        v1 =   newbeam.get_columns([4,5,6])
        flag = newbeam.get_column(10)

        t,iflag = self.calculate_intercept(x1,v1)

        x2 = numpy.zeros_like(x1)
        x2[0,:] = x1[0,:] + v1[0,:] * t
        x2[1,:] = x1[1,:] + v1[1,:] * t
        x2[2,:] = x1[2,:] + v1[2,:] * t

        for i in range(flag.size):
            if iflag[i] < 0: flag[i] = -100

        # # ;
        # # ; Calculates the normal at each intercept
        # # ;

        normal = self.get_normal(x2)

        # ;
        # ; reflection
        # ;

        v2 = self.vector_reflection(v1,normal)

        # # ;
        # # ; writes the mirr.XX file
        # # ;

        newbeam.rays[:,1-1] = x2[0,:]
        newbeam.rays[:,2-1] = x2[1,:]
        newbeam.rays[:,3-1] = x2[2,:]
        newbeam.rays[:,4-1] = v2[0,:]
        newbeam.rays[:,5-1] = v2[1,:]
        newbeam.rays[:,6-1] = v2[2,:]
        newbeam.rays[:,10-1] = flag
        #
        return newbeam,t,x1,v1,x2,v2

    # todo: move to superclass or Vector class
    def vector_reflection(self,v1,normal):
        tmp = v1 * normal
        tmp2 = tmp[0,:] + tmp[1,:] + tmp[2,:]
        tmp3 = normal.copy()

        for jj in (0,1,2):
            tmp3[jj,:] = tmp3[jj,:] * tmp2

        v2 = v1 - 2 * tmp3
        v2mod = numpy.sqrt(v2[0,:]**2 + v2[1,:]**2 + v2[2,:]**2)
        v2 /= v2mod

        return v2

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

    #
    # now real surface from file
    #
    # x,y,z = mm.read_surface_error_file("test_mesh_conic.dat")
    #
    # print(z.min(),z.max())
    # plot_surface(z,x,y)

    mm.load_file("test_mesh_conic.dat")