from shadow4.optical_surfaces.s4_optical_surface import S4OpticalSurface

# https://stackoverflow.com/questions/13360062/python-curves-intersection-with-fsolve-and-function-arguments-using-numpy

# https://www.reddit.com/r/matlab/comments/pd7rr/finding_the_point_of_intersection_between_a_line/

# from scipy.optimize import fsolve
from scipy.optimize import root, root_scalar, brentq
from scipy import interpolate

from srxraylib.plot.gol import plot,plot_image, plot_surface, plot_scatter
import sys
import time
import numpy

from shadow4.optical_surfaces.s4_optical_surface import S4OpticalSurface
from shadow4.tools.arrayofvectors import vector_refraction, vector_scattering
from shadow4.tools.arrayofvectors import vector_cross, vector_dot, vector_multiply_scalar, vector_sum, vector_diff
from shadow4.tools.arrayofvectors import vector_modulus_square, vector_modulus, vector_norm, vector_rotate_around_axis



# TODO: interp2d is deprecated in SciPy 1.10 and will be removed in SciPy 1.13.0.

class S4Mesh(S4OpticalSurface):
    def __init__(self, surface=None, mesh_x=None, mesh_y=None, mesh_z=None, kind='linear'):
        self.__x0 = [0.0,0.0,0.0]
        self.__v0 = [0.0,0.0,0.0]
        self.surface = surface # Surface must be the function defining height(x,y) or a scipy.interpolate.interp2d instance
        self.mesh_x = mesh_x # not used if surface is defined
        self.mesh_y = mesh_y # not used if surface is defined
        self.mesh_z = mesh_z # not used if surface is defined
        self.kind = kind

        if (surface is None) and (mesh_z is not None) and (mesh_x is not None) and (mesh_x is not None):
            self.calculate_surface_from_mesh()

    def set_ray(self,x0,v0):
        self.__x0 = x0
        self.__v0 = v0

    def set_surface(self, surface):
        if isinstance(surface, interpolate.interp2d):
            self.surface = surface
        else:
            try:
                z = surface(0, 0)
                self.surface = surface
            except:
                raise Exception("Surface must be the function defining height(x,y) or a scipy.interpolate.interp2d instance")

    def set_mesh(self, z, x, y):
        if (z.shape[0] != x.size) or (z.shape[1] != y.size):
            raise Exception("Bad input. It must be z[n,m], x[n], y[m]")
        self.mesh_x = x
        self.mesh_y = y
        self.mesh_z = z
        self.surface = None # reset

    def get_mesh_x_y(self):
        return self.mesh_x, self.mesh_y

    def get_mesh_z(self):
        return self.mesh_z

    def calculate_surface_from_mesh(self, alg=1):
        if alg == 0:
            self.surface = interpolate.interp2d(self.mesh_x, self.mesh_y, self.mesh_z.T, kind=self.kind)
        elif alg == 1:
            self.__TMP = interpolate.RectBivariateSpline(self.mesh_x, self.mesh_y, self.mesh_z, kx=1, ky=3)
            self.surface = lambda x, y: self.__TMP.ev(x, y)

    def surface_height(self, x, y):
        return self.surface(x, y)

    def add_to_mesh(self, z1):
        print(">>>>>>>>>>>>>>ADDING TO MESH", z1.shape, self.mesh_z.shape, z1, self.mesh_z)
        if self.mesh_z is None:
            raise Exception("Cannot add to None")

        if isinstance(z1, float):
            self.mesh_z += z1
        elif isinstance(z1, int):
            self.mesh_z += z1
        elif isinstance(z1, numpy.ndarray):
            if z1.shape != self.mesh_z.shape:
                raise Exception("Cannot add array [%,%] to mesh_z[%d,%d]" % (z1.mesh[0],
                                                                             z1.mesh[1],
                                                                             self.mesh_z.shape[0],
                                                                             self.mesh_z.shape[1]))
            self.mesh_z += z1
        else:
            print(">>>>Entered data type: ", type(z1) )
            raise Exception("Entry type not supported")


        self.calculate_surface_from_mesh()

    def load_h5file(self,filename):
        x,y,z = self.read_surface_error_h5file(filename)
        self.mesh_x = x
        self.mesh_y = y
        self.mesh_z = z
        self.calculate_surface_from_mesh()

    def load_surface_data(self, surface_data_object):
        self.mesh_x = surface_data_object._xx.copy()
        self.mesh_y = surface_data_object._yy.copy()
        self.mesh_z = surface_data_object._zz.copy().T
        self.calculate_surface_from_mesh()

    def load_surface_data_arrays(self,x,y,Z):
        self.mesh_x = x
        self.mesh_y = y
        self.mesh_z = Z
        self.calculate_surface_from_mesh()

    def load_file(self,filename):
        x,y,z = self.read_surface_error_file(filename)
        self.mesh_x = x
        self.mesh_y = y
        self.mesh_z = z
        self.calculate_surface_from_mesh()

        # self.surface = interpolate.interp2d(x,y,z.T, kind=kind)

    #
    #
    #
    def line(self,t):
        return (self.__x0[0] + self.__v0[0] * t,
                self.__x0[1] + self.__v0[1] * t,
                self.__x0[2] + self.__v0[2] * t)

    def line_z(self,t):
        return self.__x0[2] + self.__v0[2] * t

    def surface_z_vs_t(self,t):
        return self.surface_vs_t(t)[2]

    def surface_vs_t(self, t):
        x1 = self.__x0[0] + self.__v0[0] * t
        y1 = self.__x0[1] + self.__v0[1] * t
        return x1,y1,self.surface(x1,y1)

    def equation_to_solve(self, t):
        t = numpy.array(t)
        # surface height
        x1 = self.__x0[0] + self.__v0[0] * t
        y1 = self.__x0[1] + self.__v0[1] * t

        if t.size == 1:
            z1 =  self.surface(x1,y1)
        else:
            z1 = numpy.zeros_like(x1)
            for i in range(z1.size):
                tmp = self.surface(x1[i],y1[i])
                z1[i] = tmp

        # line vs t
        l1 = self.__x0[2] + self.__v0[2] * t

        return z1 - l1

    def solve(self, x_start):
        # t_solution = fsolve(lambda t:self.surface_z_vs_t(t)-self.line_z(t), x_start)
        # return t_solution[0]

        # t_solution = fsolve(self.equation_to_solve, x_start)
        # return t_solution[0]

        # t_solution = brentq(self.equation_to_solve, 0, 200)
        # return t_solution

        # t_solution = root(self.equation_to_solve, x_start, method='hybr')                 # 6.2955/5.742
        # t_solution = root(self.equation_to_solve, x_start, method='hybr', tol=1e-11)      # 5.795/5.742
        # t_solution = root(self.equation_to_solve, x_start, method='broyden2', tol=1e-11)  # 6.084/5.742
        # t_solution = root(self.equation_to_solve, x_start, method='anderson')             # 5.82/5.742
        # t_solution = root(self.equation_to_solve, x_start, method='Krylov')             # 8.12/5.742
        # t_solution = root(self.equation_to_solve, x_start, method='Krylov', tol=1e-11)  # 6.284/5.742
        # t_solution = root(self.equation_to_solve, x_start, method='diagbroyden', tol=1e-11)  # 6.08/5.742


        # t_solution = root(self.equation_to_solve, x_start, method='hybr', tol=1e-11)  # 6.17/5.27
        # t_solution = root(self.equation_to_solve, x_start, method='hybr')  # 6.17/5.27
        # t_solution = root(self.equation_to_solve, x_start, method='anderson') # 6.36/5.27

        t_solution = root(self.equation_to_solve, x_start, method='hybr', tol=1e-13)  # 6.07/5.27

        # print(t_solution['message'])
        return t_solution['x']


    #
    #
    #
    @classmethod
    def read_surface_error_h5file(cls, filename):
        import h5py
        f = h5py.File(filename, 'r')
        x = f["/surface_file/X"][:]
        y = f["/surface_file/Y"][:]
        Z = f["/surface_file/Z"][:]
        f.close()
        return x, y, Z.T.copy()


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

        z00 = self.surface(X_0, Y_0)

        N_0 = -1.0 * (self.surface(X_0 + eps, Y_0) - z00) / eps
        N_1 = -1.0 * (self.surface(X_0, Y_0 + eps) - z00) / eps
        N_2 = numpy.ones_like(X_0)

        n2 = numpy.sqrt(N_0**2 + N_1**2 + N_2**2)

        normal[0,:] = N_0 / n2
        normal[1,:] = N_1 / n2
        normal[2,:] = N_2 / n2

        return normal

    def calculate_intercept(self,XIN,VIN,keep=0):
        npoints = XIN.shape[1]
        answer = numpy.zeros(npoints)
        i_flag = numpy.ones(npoints)

        print("\n\n>>>>> main loop to find solutions (slow...)")
        t0 = time.time()
        for i in range(npoints):
            self.__x0 = XIN[:,i]
            self.__v0 = VIN[:,i]
            try:
                t_solution = self.solve(0.0)
                answer[i] = t_solution
            except:
                i_flag[i] = -1

        t1 = time.time()
        print(">>>>", answer)
        print(">>>>> done main loop to find solutions (Thanks for waiting!) Spent: %g s for %d rays (%g ms/ray)\n\n" % \
              (t1-t0, npoints, 1000 * (t1-t0) / npoints))

        return answer,i_flag

    def calculate_intercept_and_choose_solution(self, x1, v1, reference_distance=10.0, method=0):
        return self.calculate_intercept(x1, v1)

    # todo: move the apply_* methods to the parent class
    def apply_specular_reflection_on_beam(self,newbeam):
        # ;
        # ; TRACING...
        # ;

        import time
        t0 = time.time()

        x1 =   newbeam.get_columns([1,2,3])
        v1 =   newbeam.get_columns([4,5,6])
        flag = newbeam.get_column(10)
        optical_path = newbeam.get_column(13)

        t,iflag = self.calculate_intercept(x1,v1)

        print("check point 1", time.time() - t0, "s")

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

        print("check point 2", time.time() - t0, "s")

        # ;
        # ; reflection
        # ;

        v2 = self.vector_reflection(v1,normal)

        print("check point 3", time.time() - t0, "s")

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
        newbeam.rays[:,13-1] = optical_path + t
        #
        return newbeam,normal,t,x1,v1,x2,v2

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

    #
    # grating routines
    #
    def apply_grating_diffraction_on_beam(self, beam, ruling=[0.0], order=0, f_ruling=0):

        newbeam = beam.duplicate()

        x1 = newbeam.get_columns([1, 2, 3])  # numpy.array(a3.getshcol([1,2,3]))
        v1 = newbeam.get_columns([4, 5, 6])  # numpy.array(a3.getshcol([4,5,6]))
        flag = newbeam.get_column(10)  # numpy.array(a3.getshonecol(10))
        kin = newbeam.get_column(11) * 1e2 # in m^-1
        optical_path = newbeam.get_column(13)
        nrays = flag.size

        # t1, t2, iflag = self.calculate_intercept(x1, v1)
        reference_distance = -newbeam.get_column(2).mean() + newbeam.get_column(3).mean()
        # t = self.choose_solution(t1, t2, reference_distance=reference_distance)
        t, iflag = self.calculate_intercept_and_choose_solution(x1, v1, reference_distance=reference_distance)

        x2 = x1 + v1 * t
        for i in range(flag.size):
            if iflag[i] < 0: flag[i] = -100

        # ;
        # ; Calculates the normal at each intercept [see shadow's normal.F]
        # ;

        normal = self.get_normal(x2)

        # ;
        # ; reflection
        # ;
        # v2 =  v1.T - 2 * vector_multiply_scalar(normal.T, vector_dot(v1.T, normal.T))
        # V_OUT = v2.copy()
        # v2 = v2.T

        # ;
        # ; grating scattering
        # ;
        if True:
            DIST = x2[1]
            RDENS = 0.0
            for n in range(len(ruling)):
                RDENS += ruling[n] * DIST**n

            PHASE = optical_path + 2 * numpy.pi * order * DIST * RDENS / kin
            G_MOD = 2 * numpy.pi * RDENS * order


            # capilatized vectors are [:,3] as required for vector_* operations
            VNOR = normal.T
            VNOR = vector_multiply_scalar(VNOR, -1.0) # outward normal


            # print(">>>> VNOR: (%20.18g,%20.18g,%20.18f) mod: %20.18f" % (VNOR[-1, 0], VNOR[-1, 1], VNOR[-1, 2],
            #                                          (VNOR[-1, 0]**2 + VNOR[-1, 1]**2 + VNOR[-1, 2]**2)))

            # versors
            X_VRS = numpy.zeros((nrays,3))
            X_VRS[:,0] = 1
            Y_VRS = numpy.zeros((nrays, 3))
            Y_VRS[:,1] = 1

            if f_ruling == 0:
                G_FAC = vector_dot(VNOR, Y_VRS)
                G_FAC = numpy.sqrt(1 - G_FAC**2)
            elif f_ruling == 1:
                G_FAC = 1.0
            elif f_ruling == 5:
                G_FAC = vector_dot(VNOR, Y_VRS)
                G_FAC = numpy.sqrt(1 - G_FAC**2)

            G_MODR = G_MOD * G_FAC


            K_IN = vector_multiply_scalar(v1.T, kin)
            K_IN_NOR = vector_multiply_scalar(VNOR, vector_dot(K_IN, VNOR) )
            K_IN_PAR = vector_diff(K_IN, K_IN_NOR)


            VTAN = vector_cross(VNOR, X_VRS)
            GSCATTER = vector_multiply_scalar(VTAN, G_MODR)


            K_OUT_PAR = vector_sum(K_IN_PAR, GSCATTER)
            K_OUT_NOR = vector_multiply_scalar(VNOR,  numpy.sqrt(kin**2 - vector_modulus_square(K_OUT_PAR)))
            K_OUT = vector_sum(K_OUT_PAR, K_OUT_NOR)
            V_OUT = vector_norm(K_OUT)

        # ;
        # ; writes the mirr.XX file
        # ;

        newbeam.set_column(1, x2[0])
        newbeam.set_column(2, x2[1])
        newbeam.set_column(3, x2[2])
        newbeam.set_column(4, V_OUT.T[0])
        newbeam.set_column(5, V_OUT.T[1])
        newbeam.set_column(6, V_OUT.T[2])
        newbeam.set_column(10, flag)
        newbeam.set_column(13, optical_path + t)

        return newbeam, normal
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




if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()


    def sphere(x, y, radius=5.0):
        return radius - numpy.sqrt(radius ** 2 - x ** 2 - y ** 2)

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

    plot_surface_and_line(Z,x,y,zz,xx,yy)

    #
    # mesh object
    #
    mm = S4Mesh()
    mm.set_ray(x0,v0)
    mm.set_surface(sphere)

    zz = mm.surface_height(X, Y)
    plot_surface(zz, x, y, xtitle="X")

    #
    # intercept
    #
    x_start = 0

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

    # mm.load_file("test_mesh_conic.dat")
    x2 = numpy.zeros((3,10))
    x2 = numpy.random.rand(30).reshape((3,10))
    print("normal: ", mm.get_normal(x2))
