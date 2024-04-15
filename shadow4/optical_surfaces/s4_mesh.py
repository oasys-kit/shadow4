"""

Defines the shadow4 Mesh class to deal with a numerical surfaces (defined by an array of points).

"""
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

from shadow4.tools.arrayofvectors import vector_reflection
from shadow4.tools.logger import is_verbose, is_debug

class S4Mesh(S4OpticalSurface):
    """
    Class to manage optical surfaces defined by numeric arrays.

    Parameters
    ----------
    surface : callable, optional
        A function to return the surface height z from coordinates x, y as arguments. Nor needed if mesh_* are given.
    mesh_x : numpy array, optional
        The array with the X coordinate in m. Nor used if surface is given.
    mesh_y : numpy array, optional
        The array with the X coordinate in m. Nor used if surface is given.
    mesh_z : numpy array, optional
        The array with the X coordinate in m. Nor used if surface is given.
    """
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

    #
    # setters + getters
    #

    def get_mesh_x_y(self) -> Tuple[numpy.ndarray, numpy.ndarray]:
        """
        Returns the x and y arrays.

        Returns
        -------
        tuple
        (x, y) numpy 1D arrays.
        """
        return self._mesh_x, self._mesh_y

    def get_mesh_z(self) -> numpy.ndarray:
        """
        returns the z 2D array.

        Returns
        -------
        numpy 2D array
        """
        return self._mesh_z

    #
    # overloaded methods
    #

    def info(self):
        """
        Creates an info text.

        Returns
        -------
        str
        """
        txt = ""


        txt += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
        if self.__x0     is not None: txt += "\nLoaded rays. Shape:" + repr(self.__x0.shape)
        if self._mesh_x  is not None: txt += "\nX shape:"            + repr(self._mesh_x.shape)
        if self._mesh_y  is not None: txt += "\nY shape:"            + repr(self._mesh_y.shape)
        if self._mesh_z  is not None: txt += "\nZ shape:"            + repr(self._mesh_z.shape)
        if self._surface is not None: txt += "\nfunction:"           + repr(self._surface)
        txt += "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

        return txt

    def duplicate(self):
        """
        Duplicates an instance of S4Toroid

        Returns
        -------
        instance of S4Toroid.
        """
        return S4Mesh(
                      surface = self._surface,
                      mesh_x = self._mesh_x,
                      mesh_y = self._mesh_y,
                      mesh_z = self._mesh_z,
                      )

    def get_normal(self, x2: numpy.ndarray):
        """
        Calculates the normal vector (or stack of vectors) at a point on the surface.

        Parameters
        ----------
        x2 : numpy array
            The coordinates vector(s) of shape [3, NRAYS].

        Returns
        -------
        numpy array
            The normal vector(s) of shape [3, NRAYS].

        """
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

    def calculate_intercept_and_choose_solution(self, x1: numpy.ndarray, v1: numpy.ndarray,
                                                reference_distance: float = 10.0, method: int = 0, # required by upper class
                                                ):
        """

        Calculates the intercept point (or stack of points) for a given ray or stack of rays,
        given a point XIN and director vector VIN.

        Parameters
        ----------
        XIN : numpy array
            The coordinates of a point of origin of the ray: shape [3, NRAYS].
        VIN : numpy array
            The coordinates of a director vector the ray: shape [3, NRAYS].
        reference_distance : float, optional
            Not used in S4Mesh.
        method : int, optional
            Not used in S4Mesh.

        Returns
        -------
        tuple
            (answer, i_flag) The selected solution (time or flight path numpy array) and the flag numpy array.

        """
        return self.calculate_intercept(x1, v1)

    def calculate_intercept(self, XIN: numpy.ndarray, VIN: numpy.ndarray, keep=0):
        """

        Calculates the intercept point (or stack of points) for a given ray or stack of rays,
        given a point XIN and director vector VIN.

        Parameters
        ----------
        XIN : numpy array
            The coordinates of a point of origin of the ray: shape [3, NRAYS].
        VIN : numpy array
            The coordinates of a director vector the ray: shape [3, NRAYS].

        Returns
        -------
        tuple
            (answer, i_flag) The selected solution (time or flight path numpy array) and the flag numpy array.

        """
        npoints = XIN.shape[1]

        if is_debug(): print("\n\n>>>>> main loop to find solutions")
        t0 = time.time()

        self.__x0 = XIN
        self.__v0 = VIN

        i_flag = numpy.ones(npoints)
        answer, success = self._solve(x_start=numpy.zeros(npoints))

        i_flag[numpy.where(success == False)] = -1

        t1 = time.time()
        if is_debug(): print(">>>>> done main loop to find solutions, spent: %g s for %d rays (%g ms/ray)\n\n" % (t1 - t0, npoints, 1000 * (t1 - t0) / npoints))

        return answer, i_flag


    def surface_height(self, x: Union[float, numpy.ndarray], y: Union[float, numpy.ndarray]):
        """
        Calculates a 2D mesh array with the surface heights.

        Parameters
        ----------
        x : float or numpy 2D array (mesh)
            The x coordinate(s).
        y : float or numpy 2D array
            The y coordinate(s).

        Returns
        -------
        2D numpy array
            the height mesh.
        """
        return self._surface(x, y)

    #
    # Other calculations
    #

    def add_to_mesh(self, z1: Union[int, float, numpy.ndarray]):
        """
        Add a new numerical mesh to the existing one.

        Parameters
        ----------
        z1 : int, float, or numpy array
            The 2D Mesh (or the scalar value).

        """
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
            if is_debug(): print(">>>>Entered data type: ", type(z1) )
            raise ValueError("Entry type not supported")

        self._calculate_surface_from_mesh()

    def load_file(self, filename: str):
        """
        Loads a numeric mesh from an text file (with the SHADOW3 presurface preprocessor format).

        Parameters
        ----------
        filename : str
            The file name.
        """
        x, y, z = self._read_surface_error_file(filename)

        self._mesh_x = x
        self._mesh_y = y
        self._mesh_z = z
        self._calculate_surface_from_mesh()

    def load_h5file(self, filename: str):
        """
        Loads a numeric mesh from an hdf5 file (with the standard OASYS surface description).

        Parameters
        ----------
        filename : str
            The file name.
        """
        x, y, z = self._read_surface_error_h5file(filename)

        self._mesh_x = x
        self._mesh_y = y
        self._mesh_z = z
        self._calculate_surface_from_mesh()

    def load_surface_data(self, surface_data_object):
        """
        Loads a numeric mesh from an OASYS surface data object.

        Parameters
        ----------
        filename : instance of OasysSurfaceData
            The data object.
        """
        self._mesh_x = surface_data_object._xx.copy()
        self._mesh_y = surface_data_object._yy.copy()
        self._mesh_z = surface_data_object._zz.copy().T

        self._calculate_surface_from_mesh()

    def load_surface_data_arrays(self, x: numpy.ndarray, y: numpy.ndarray, Z: numpy.ndarray):
        """
        Loads a numeric mesh from arrays.

        Parameters
        ----------
        x : numpy array
            1D array with the x coordinates in m.
        y : numpy array
            1D array with the y coordinates in m.
        Z : numpy array
            2D array with the z(x,y) coordinates in m.
        """
        self._mesh_x = x
        self._mesh_y = y
        self._mesh_z = Z
        self._calculate_surface_from_mesh()


    #
    # PROTECTED METHODS
    #
    def _set_rays(self, x0: numpy.ndarray, v0: numpy.ndarray):
        self.__x0 = x0
        self.__v0 = v0

    def _set_surface(self, surface: Callable):
        try:
            _ = surface(0, 0)
            self._surface = surface
        except:
            raise Exception("Surface must be the function defining height(x,y) or a scipy.interpolate.interp2d instance")

    def _calculate_surface_from_mesh(self):
        self.__interpolating_agent = interpolate.RectBivariateSpline(self._mesh_x, self._mesh_y, self._mesh_z, kx=3, ky=3)

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

    @classmethod
    def _read_surface_error_h5file(cls, filename):
        import h5py
        f = h5py.File(filename, 'r')
        x = f["/surface_file/X"][:]
        y = f["/surface_file/Y"][:]
        Z = f["/surface_file/Z"][:]
        f.close()
        return x, y, Z.T.copy()

    @classmethod
    def _read_surface_error_file(cls, filename: str): #copied from shadowOui util/shadow_util
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



if __name__ == "__main__":
    from srxraylib.plot.gol import set_qt
    set_qt()


    def sphere(x, y, radius=5.0):
        return radius - numpy.sqrt(radius ** 2 - x ** 2 - y ** 2)


    def plot_surface_and_line(Z, x, y, zz, xx, yy, show=True):

        try:
            import matplotlib.pylab as plt
            from mpl_toolkits.mplot3d import Axes3D
            from matplotlib import cm
            from matplotlib.ticker import LinearLocator, FormatStrFormatter
        except:
            print("Please install matplotlib to allow graphics")
        #
        # plot
        #

        X = numpy.outer(x, numpy.ones_like(y))
        Y = numpy.outer(numpy.ones_like(x), y)

        fig = plt.figure(figsize=None)
        ax = fig.add_subplot(projection='3d')
        #

        cmap = cm.coolwarm

        ax = fig.add_subplot(projection='3d')
        ax.plot(xx, yy, zz, label='parametric curve')
        ax.set_xlim(y.min(), y.max())
        ax.set_ylim(y.min(), y.max())

        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cmap, linewidth=0, antialiased=False)

        ax.set_zlim(zz.min(), zz.max())
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

