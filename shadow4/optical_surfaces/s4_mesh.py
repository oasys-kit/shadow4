from typing import Callable, Tuple, List, Union
from scipy.optimize import root
from scipy import interpolate

from srxraylib.plot.gol import plot_surface
import sys
import time
import numpy

from shadow4.optical_surfaces.s4_optical_surface import S4OpticalSurface
from shadow4.tools.arrayofvectors import vector_cross, vector_dot, vector_multiply_scalar, vector_sum, vector_diff
from shadow4.tools.arrayofvectors import vector_modulus_square, vector_norm
from shadow4.beam.s4_beam import S4Beam


class S4Mesh(S4OpticalSurface):
    def __init__(self,
                 surface: Callable = None,
                 mesh_x: numpy.ndarray = None,
                 mesh_y: numpy.ndarray = None,
                 mesh_z: numpy.ndarray = None):
        self.__x0 = None
        self.__v0 = None
        self._surface = surface # Surface must be the function defining height(x,y)
        self._mesh_x  = mesh_x  # not used if surface is defined
        self._mesh_y  = mesh_y  # not used if surface is defined
        self._mesh_z  = mesh_z  # not used if surface is defined

        if (surface is None) and not (mesh_z is None or mesh_x is None or mesh_x is None):
            self._calculate_surface_from_mesh()

    def set_rays(self, x0: numpy.ndarray, v0: numpy.ndarray):
        self.__x0 = x0
        self.__v0 = v0

    def set_surface(self, surface: Callable):
        try:
            _ = surface(0, 0)
            self._surface = surface
        except:
            raise Exception("Surface must be the function defining height(x,y) or a scipy.interpolate.interp2d instance")

    def set_mesh(self, z: numpy.ndarray, x: numpy.ndarray, y: numpy.ndarray):
        if z is None or x is None or y is None: raise ValueError("Bad input. Arrays cannot be None")
        if (z.shape[0] != x.size) or (z.shape[1] != y.size): raise ValueError("Bad input. It must be z[n,m], x[n], y[m]")
        self._mesh_x  = x
        self._mesh_y  = y
        self._mesh_z  = z
        self._surface = None # reset

    def get_mesh_x_y(self) -> Tuple[numpy.ndarray, numpy.ndarray]:
        return self._mesh_x, self._mesh_y

    def get_mesh_z(self) -> numpy.ndarray:
        return self._mesh_z

    def surface_height(self, x: Union[float, numpy.ndarray], y: Union[float, numpy.ndarray]):
        return self._surface(x, y)

    def add_to_mesh(self, z1: Union[int, float, numpy.ndarray]):
        print(">>>>>>>>>>>>>>ADDING TO MESH", z1.shape, self._mesh_z.shape, z1, self._mesh_z)
        if self._mesh_z is None: raise ValueError("Cannot add to None")

        if isinstance(z1, float) or isinstance(z1, int):  self._mesh_z += z1
        elif isinstance(z1, numpy.ndarray):
            if z1.shape != self._mesh_z.shape:
                raise ValueError("Cannot add array [%,%] to mesh_z[%d,%d]" % (z1.mesh[0],
                                                                              z1.mesh[1],
                                                                              self._mesh_z.shape[0],
                                                                              self._mesh_z.shape[1]))
            self._mesh_z += z1
        else:
            print(">>>>Entered data type: ", type(z1) )
            raise ValueError("Entry type not supported")

        self._calculate_surface_from_mesh()

    def load_file(self, filename: str):
        x, y, z = self._read_surface_error_file(filename)

        self._mesh_x = x
        self._mesh_y = y
        self._mesh_z = z
        self._calculate_surface_from_mesh()

    def load_h5file(self, filename: str):
        x, y, z = self._read_surface_error_h5file(filename)

        self._mesh_x = x
        self._mesh_y = y
        self._mesh_z = z
        self._calculate_surface_from_mesh()

    def load_surface_data(self, surface_data_object):
        self._mesh_x = surface_data_object._xx.copy()
        self._mesh_y = surface_data_object._yy.copy()
        self._mesh_z = surface_data_object._zz.copy().T

        self._calculate_surface_from_mesh()

    def load_surface_data_arrays(self, x: numpy.ndarray, y: numpy.ndarray, Z: numpy.ndarray):
        self._mesh_x = x
        self._mesh_y = y
        self._mesh_z = Z
        self._calculate_surface_from_mesh()

    # todo: move the apply_* methods to the parent class
    def apply_specular_reflection_on_beam(self, newbeam: S4Beam) -> Tuple[S4Beam,
                                                                          numpy.ndarray,
                                                                          numpy.ndarray,
                                                                          numpy.ndarray,
                                                                          numpy.ndarray,
                                                                          numpy.ndarray,
                                                                          numpy.ndarray ]:
        # ;
        # ; TRACING...
        # ;
        x1           = newbeam.get_columns([1,2,3])
        v1           = newbeam.get_columns([4,5,6])
        flag         = newbeam.get_column(10)
        optical_path = newbeam.get_column(13)

        t, iflag = self.calculate_intercept(x1, v1)

        x2 = numpy.zeros_like(x1)
        x2[0, :] = x1[0, :] + numpy.multiply(v1[0, :], t)
        x2[1, :] = x1[1, :] + numpy.multiply(v1[1, :], t)
        x2[2, :] = x1[2, :] + numpy.multiply(v1[2, :], t)

        flag[numpy.where(iflag < 0)] = -100

        # # ;
        # # ; Calculates the normal at each intercept
        # # ;
        normal = self.get_normal(x2)

        # ;
        # ; reflection
        # ;
        v2 = self.vector_reflection(v1, normal)

        newbeam.rays[:,2-1] = x2[1,:]
        newbeam.rays[:,3-1] = x2[2,:]
        newbeam.rays[:,4-1] = v2[0,:]
        newbeam.rays[:,5-1] = v2[1,:]
        newbeam.rays[:,6-1] = v2[2,:]
        newbeam.rays[:,10-1] = flag
        newbeam.rays[:,13-1] = optical_path + t
        #
        return newbeam, normal, t, x1, v1, x2, v2

    def get_normal(self, x2: numpy.ndarray):
        # ;
        # ; Calculates the normal at intercept points x2 [see shadow's normal.F]
        # ;

        normal = numpy.zeros_like(x2)

        eps = 100 * sys.float_info.epsilon

        X_0 = x2[0, :]
        Y_0 = x2[1, :]

        z00 = self._surface(X_0, Y_0)

        N_0 = -1.0 * (self._surface(X_0 + eps, Y_0) - z00) / eps
        N_1 = -1.0 * (self._surface(X_0, Y_0 + eps) - z00) / eps
        N_2 = numpy.ones_like(X_0)

        n2 = numpy.sqrt(N_0 ** 2 + N_1 ** 2 + N_2 ** 2)

        normal[0, :] = N_0 / n2
        normal[1, :] = N_1 / n2
        normal[2, :] = N_2 / n2

        return normal

    def calculate_intercept(self, XIN: numpy.ndarray, VIN: numpy.ndarray, keep=0):
        npoints = XIN.shape[1]

        print("\n\n>>>>> main loop to find solutions (slow...)")
        t0 = time.time()

        self.__x0 = XIN
        self.__v0 = VIN

        i_flag = numpy.ones(npoints)
        answer, success = self._solve(x_start=numpy.zeros(npoints))
        i_flag[numpy.where(success == False)] = -1

        t1 = time.time()
        print(">>>>> done main loop to find solutions, spent: %g s for %d rays (%g ms/ray)\n\n" % (t1 - t0, npoints, 1000 * (t1 - t0) / npoints))

        return answer, i_flag

    def calculate_intercept_and_choose_solution(self, x1: numpy.ndarray, v1: numpy.ndarray, reference_distance: float = 10.0, method: int = 0):
        return self.calculate_intercept(x1, v1)

    # todo: move to superclass or Vector class
    def vector_reflection(self, v1, normal):
        tmp = v1 * normal
        tmp2 = tmp[0, :] + tmp[1, :] + tmp[2, :]
        tmp3 = normal.copy()

        for jj in (0, 1, 2): tmp3[jj, :] = tmp3[jj, :] * tmp2

        v2 = v1 - 2 * tmp3
        v2mod = numpy.sqrt(v2[0, :] ** 2 + v2[1, :] ** 2 + v2[2, :] ** 2)
        v2 /= v2mod

        return v2

    #
    # grating routines
    #
    def apply_grating_diffraction_on_beam(self, beam, ruling=[0.0], order=0, f_ruling=0):
        newbeam = beam.duplicate()

        x1 = newbeam.get_columns([1, 2, 3])
        v1 = newbeam.get_columns([4, 5, 6])
        flag = newbeam.get_column(10)
        kin = newbeam.get_column(11) * 1e2  # in m^-1
        optical_path = newbeam.get_column(13)
        nrays = flag.size

        reference_distance = -newbeam.get_column(2).mean() + newbeam.get_column(3).mean()
        t, iflag = self.calculate_intercept_and_choose_solution(x1, v1, reference_distance=reference_distance)

        x2 = x1 + numpy.multiply(v1, t)
        flag[numpy.where(iflag < 0)] = -100

        # ;
        # ; Calculates the normal at each intercept [see shadow's normal.F]
        # ;
        normal = self.get_normal(x2)

        # ;
        # ; grating scattering
        # ;
        if True:
            DIST = x2[1]
            RDENS = 0.0
            for n in range(len(ruling)): RDENS += ruling[n] * DIST ** n

            PHASE = optical_path + 2 * numpy.pi * order * DIST * RDENS / kin
            G_MOD = 2 * numpy.pi * RDENS * order

            # capilatized vectors are [:,3] as required for vector_* operations
            VNOR = normal.T
            VNOR = vector_multiply_scalar(VNOR, -1.0)  # outward normal

            # versors
            X_VRS = numpy.zeros((nrays, 3))
            X_VRS[:, 0] = 1
            Y_VRS = numpy.zeros((nrays, 3))
            Y_VRS[:, 1] = 1

            if f_ruling == 0:
                G_FAC = vector_dot(VNOR, Y_VRS)
                G_FAC = numpy.sqrt(1 - G_FAC ** 2)
            elif f_ruling == 1:
                G_FAC = 1.0
            elif f_ruling == 5:
                G_FAC = vector_dot(VNOR, Y_VRS)
                G_FAC = numpy.sqrt(1 - G_FAC ** 2)

            G_MODR = G_MOD * G_FAC

            K_IN = vector_multiply_scalar(v1.T, kin)
            K_IN_NOR = vector_multiply_scalar(VNOR, vector_dot(K_IN, VNOR))
            K_IN_PAR = vector_diff(K_IN, K_IN_NOR)

            VTAN = vector_cross(VNOR, X_VRS)
            GSCATTER = vector_multiply_scalar(VTAN, G_MODR)

            K_OUT_PAR = vector_sum(K_IN_PAR, GSCATTER)
            K_OUT_NOR = vector_multiply_scalar(VNOR, numpy.sqrt(kin ** 2 - vector_modulus_square(K_OUT_PAR)))
            K_OUT = vector_sum(K_OUT_PAR, K_OUT_NOR)
            V_OUT = vector_norm(K_OUT)

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
    # PROTECTED METHODS
    #

    def _calculate_surface_from_mesh(self):
        self.__interpolating_agent = interpolate.RectBivariateSpline(self._mesh_x, self._mesh_y, self._mesh_z, kx=1, ky=3)

        self._surface = lambda x, y: self.__interpolating_agent.ev(x, y)

    def _line(self, t: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
        return self.__x0[0, :] + numpy.multiply(self.__v0[0, :], t), \
               self.__x0[1, :] + numpy.multiply(self.__v0[1, :], t), \
               self.__x0[2, :] + numpy.multiply(self.__v0[2, :], t)

    def _line_x_y(self, t: numpy.ndarray) ->  numpy.ndarray:
        return self.__x0[0, :] + numpy.multiply(self.__v0[0, :], t), \
               self.__x0[1, :] + numpy.multiply(self.__v0[1, :], t)

    def _line_z(self, t: numpy.ndarray) ->  numpy.ndarray:
        return self.__x0[2, :] + numpy.multiply(self.__v0[2, :], t)

    def _surface_z_vs_t(self, t: numpy.ndarray) -> numpy.ndarray:
        _, _, surface_z = self._surface_vs_t(t)

        return surface_z

    def _surface_vs_t(self, t: numpy.ndarray) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]:
        x1, y1 = self._line_x_y(t)

        return x1, y1, self._surface(x1, y1)

    def _equation_to_solve(self, t: numpy.ndarray) ->  numpy.ndarray:
        _, _, z1 = self._surface_vs_t(t)

        return z1 - self._line_z(t)

    def _solve(self, x_start: numpy.ndarray):
        t_solution = root(self._equation_to_solve, x_start, method='df-sane', tol=None)

        return t_solution['x'], t_solution['success']


    #
    # CLASS METHODS
    #
    @classmethod
    def _read_surface_error_h5file(cls, filename):
        import h5py
        f = h5py.File(filename, 'r')
        x = f["/surface_file/X"][:]
        y = f["/surface_file/Y"][:]
        Z = f["/surface_file/Z"][:]
        f.close()
        return x, y, Z.T.copy()


    #copied from shadowOui util/shadow_util
    @classmethod
    def _read_surface_error_file(cls, filename: str):
        file = open(filename, "r")
        rows = file.readlines()

        dimensions = rows[0].split()
        n_x = int(dimensions[0])
        n_y = int(dimensions[1])

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
                for y_value in y_values: y_coords = numpy.append(y_coords, float(y_value))
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

                if len(x_row) != 1 + n_y: x_row = rows[index].split()
                if len(x_row) != 1 + n_y: raise Exception("Malformed file: check format")

                x_rows.append(x_row)

            for x_index in range(0, len(x_rows)):
                x_coords = numpy.append(x_coords, float(x_rows[x_index][0]))

                for z_index in range(0, len(x_rows[x_index]) - 1):
                    z_values[x_index, z_index] = float(x_rows[x_index][z_index + 1])

        return x_coords, y_coords, z_values

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

    t_solution = mm._solve(x_start)

    print("t_solution: ",t_solution)
    print("line: ", mm._line(t_solution))
    print("surface: ", mm._surface_vs_t(t_solution))

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
