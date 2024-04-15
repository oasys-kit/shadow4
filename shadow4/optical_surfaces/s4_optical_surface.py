"""

Abstract class defining the interfaces of the optical surfaces to be implemented in the subclasses.

It also defines common utilities to write the surfaces to files.

"""
import numpy
import os
import h5py
import time

from shadow4.tools.arrayofvectors import vector_reflection, vector_refraction, vector_scattering
from shadow4.tools.arrayofvectors import vector_cross, vector_dot, vector_multiply_scalar, vector_sum, vector_diff
from shadow4.tools.arrayofvectors import vector_modulus_square, vector_modulus, vector_norm, vector_rotate_around_axis
from shadow4.tools.logger import is_verbose, is_debug

class S4OpticalSurface(object):

    def info(self):
        raise NotImplementedError("Subclasses should implement this!")

    def duplicate(self):
        raise NotImplementedError("Subclasses should implement this!")

    def surface_height(self, x, y, **kwargs):
        raise NotImplementedError("Subclasses should implement this!")

    def get_normal(self, x, **kwargs):
        raise NotImplementedError("Subclasses should implement this!")

    def calculate_intercept(self, XIN, VIN, **kwargs): # todo: remove?
        raise NotImplementedError("Subclasses should implement this!")

    def choose_solution(self, TPAR1, TPAR2, **kwargs):
        raise NotImplementedError("Subclasses should implement this!")

    def calculate_intercept_and_choose_solution(self, XIN, VIN, **kwargs): # todo: common implementation it here
        raise NotImplementedError("Subclasses should implement this!")

    def apply_specular_reflection_on_beam(self, beam):
        newbeam = beam.duplicate()

        # ;
        # ; TRACING...
        # ;

        x1 = newbeam.get_columns([1, 2, 3])
        v1 = newbeam.get_columns([4, 5, 6])
        flag = newbeam.get_column(10)
        optical_path = newbeam.get_column(13)

        reference_distance = -newbeam.get_column(2).mean() + newbeam.get_column(3).mean()
        t, iflag = self.calculate_intercept_and_choose_solution(x1, v1, reference_distance=reference_distance, method=0)

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
        v2 = (vector_reflection(v1.T, normal.T)).T

        # ;
        # ; writes the mirr arrays
        # ;
        newbeam.set_column(1, x2[0])
        newbeam.set_column(2, x2[1])
        newbeam.set_column(3, x2[2])
        newbeam.set_column(4, v2[0])
        newbeam.set_column(5, v2[1])
        newbeam.set_column(6, v2[2])
        newbeam.set_column(10, flag)
        newbeam.set_column(13, optical_path + t)

        return newbeam, normal, t, x1, v1, x2, v2

    def apply_refraction_on_beam(self,
                                 beam,
                                 refraction_index_object,
                                 refraction_index_image,
                                 apply_attenuation=0,
                                 linear_attenuation_coefficient=0.0,  # in SI, i.e. m^-1
                                 ):

        # ;
        # ; TRACING...
        # ;
        newbeam = beam.duplicate()

        x1 = newbeam.get_columns([1, 2, 3])
        v1 = newbeam.get_columns([4, 5, 6])
        flag = newbeam.get_column(10)
        k_in_mod = newbeam.get_column(11)
        optical_path = newbeam.get_column(13)

        reference_distance = -newbeam.get_column(2).mean() + newbeam.get_column(3).mean()
        t, iflag = self.calculate_intercept_and_choose_solution(x1, v1, reference_distance=reference_distance, method=0)

        x2 = x1 + v1 * t
        for i in range(flag.size):
            if iflag[i] < 0: flag[i] = -100

        # ;
        # ; Calculates the normal at each intercept [see shadow's normal.F]
        # ;
        normal = self.get_normal(x2)

        # ;
        # ; refraction
        # ;

        # note that sgn=None tells vector_refraction to compute the right sign of the sqrt.
        # This is equivalent to change the direction of the normal in the case that it is an inwards normal.
        v2t = vector_refraction(v1.T, normal.T, refraction_index_object, refraction_index_image, sgn=None)
        v2 = v2t.T

        # ;
        # ; writes the beam arrays
        # ;

        newbeam.set_column(1, x2[0])
        newbeam.set_column(2, x2[1])
        newbeam.set_column(3, x2[2])
        newbeam.set_column(4, v2[0])
        newbeam.set_column(5, v2[1])
        newbeam.set_column(6, v2[2])
        newbeam.set_column(10, flag)
        newbeam.set_column(11, k_in_mod * refraction_index_image / refraction_index_object)
        newbeam.set_column(13, optical_path + t * refraction_index_object)

        if apply_attenuation:
            att1 = numpy.sqrt(numpy.exp(-numpy.abs(t) * linear_attenuation_coefficient))
            if is_debug(): print(">>> mu (object space): ", linear_attenuation_coefficient)
            if is_debug(): print(">>> attenuation of amplitudes (object space): ", att1)
            newbeam.rays[:, 7 - 1 ] *= att1
            newbeam.rays[:, 8 - 1 ] *= att1
            newbeam.rays[:, 9 - 1 ] *= att1
            newbeam.rays[:, 16 - 1] *= att1
            newbeam.rays[:, 17 - 1] *= att1
            newbeam.rays[:, 18 - 1] *= att1

        return newbeam, normal


    def apply_grating_diffraction_on_beam(self, beam, ruling=[0.0], order=0, f_ruling=0):

        newbeam = beam.duplicate()

        x1 = newbeam.get_columns([1, 2, 3])
        v1 = newbeam.get_columns([4, 5, 6])
        flag = newbeam.get_column(10)
        kin = newbeam.get_column(11) * 1e2 # in m^-1
        optical_path = newbeam.get_column(13)
        nrays = flag.size

        reference_distance = -newbeam.get_column(2).mean() + newbeam.get_column(3).mean()
        t, iflag = self.calculate_intercept_and_choose_solution(x1, v1, reference_distance=reference_distance, method=0)

        x2 = x1 + v1 * t
        for i in range(flag.size):
            if iflag[i] < 0: flag[i] = -100

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
            for n in range(len(ruling)):
                RDENS += ruling[n] * DIST**n

            PHASE = optical_path + 2 * numpy.pi * order * DIST * RDENS / kin
            G_MOD = 2 * numpy.pi * RDENS * order


            # capilatized vectors are [:,3] as required for vector_* operations
            VNOR = normal.T
            VNOR = vector_multiply_scalar(VNOR, -1.0) # outward normal


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
        # ; writes the beam arrays
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
    # common utilities
    #
    def write_mesh_file(self, x, y, filename="surface.dat"):
        """
        Writes the optical surface in a file with SHADOW3 format (presurface preprocessor).
        The arrays for the x and y coordinates should be provided to evaluate the height z(x, y).

        Parameters
        ----------
        x : numpy array
            The array with the X coordinates.
        y : numpy array
            The array with the Y coordinates.
        filename : str, optional
            The file name (for output).

        Returns
        -------
        int
            1=Succes, 0=Error
        """
        X = numpy.outer(x, numpy.ones_like(y))
        Y = numpy.outer(numpy.ones_like(x), y)
        Z = self.surface_height(X, Y)
        return self.write_shadow_surface(Z.T, x, y, outFile=filename)

    def write_mesh_h5file(self, x, y, filename="surface.h5",  subgroup_name="surface_file", overwrite=True):
        """
        Writes the optical surface in a hdf5 file.
        The arrays for the x and y coordinates should be provided to evaluate the height z(x, y).

        Parameters
        ----------
        x : numpy array
            The array with the X coordinates.
        y : numpy array
            The array with the Y coordinates.
        filename : str, optional
            The file name (for output).
        subgroup_name : str, optional
            The h5 subgroup name (surface_file for OASYS compatibility).
        overwrite : int, optional
            Fleg to overwritte file: 0=Append to existing, 1=Overwrite existing file.

        """
        X = numpy.outer(x, numpy.ones_like(y))
        Y = numpy.outer(numpy.ones_like(x), y)
        Z = self.surface_height(X, Y)

        if (os.path.isfile(filename)) and (overwrite == True): os.remove(filename)

        if not os.path.isfile(filename):  # if file doesn't exist, create it.
            file = h5py.File(filename, 'w')
            # points to the default data to be plotted
            file.attrs['default'] = subgroup_name + '/Z'
            # give the HDF5 root some more attributes
            file.attrs['file_name'] = filename
            file.attrs['file_time'] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            file.attrs['creator'] = 'write_surface_file'
            file.attrs['code'] = 'Oasys'
            file.attrs['HDF5_Version'] = h5py.version.hdf5_version
            file.attrs['h5py_version'] = h5py.version.version
            file.close()

        file = h5py.File(filename, 'a')

        try:
            f1 = file.create_group(subgroup_name)
        except:
            f1 = file[subgroup_name]

        f1z = f1.create_dataset("Z", data=Z.T)
        f1x = f1.create_dataset("X", data=x)
        f1y = f1.create_dataset("Y", data=y)

        # NEXUS attributes for automatic plot
        f1.attrs['NX_class'] = 'NXdata'
        f1.attrs['signal'] = "Z"
        f1.attrs['axes'] = [b"Y", b"X"]

        f1z.attrs['interpretation'] = 'image'
        f1x.attrs['long_name'] = "X [m]"
        f1y.attrs['long_name'] = "Y [m]"

        file.close()
        print("File %s written to disk." % filename)

    @classmethod
    def write_shadow_surface(cls, s, xx, yy, outFile='presurface.dat'):
        """
        Writes a mesh in the SHADOW3/presurface format

        Parameters
        ----------
        s : numpy array
            The array with the heights s(xx, yy) in m.
        xx : numpy array
            The array with the X spatial coordinates (sagittal, along mirror width) in m.
        yy : numpy array
            The array with the Y spatial coordinates (tangential, along mirror length.) in m.
        outFile : str, optional
            The file name (for output).

        Returns
        -------
        int
            1=Succes, 0=Error

        """
        # modified from shadow3 ShadowTools
        out = 1

        try:
            fs = open(outFile, 'w')
        except IOError:
            out = 0
            print("Error: can\'t open file: " + outFile)
            return
        else:
            # dimensions
            fs.write(repr(xx.size) + " " + repr(yy.size) + " \n")
            # y array
            for i in range(yy.size):
                fs.write(' ' + repr(yy[i]))
            fs.write("\n")
            # for each x element, the x value and the corresponding z(y) profile
            for i in range(xx.size):
                tmps = ""
                for j in range(yy.size):
                    tmps = tmps + "  " + repr(s[j, i])
                fs.write(' ' + repr(xx[i]) + " " + tmps)
                fs.write("\n")
            fs.close()
            if is_verbose(): print("write_shadow_surface: File for SHADOW " + outFile + " written to disk.")
        return out

if __name__ == "__main__":
    a = S4OpticalSurface()
    a.info()