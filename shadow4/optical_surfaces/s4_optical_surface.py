"""

Abstract class defining the interfaces of the optical surfaces to be implemented in the subclasses.

It also defines common utilities to write the surfaces to files.

"""
import numpy
import os
import h5py
import time



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

    def apply_specular_reflection_on_beam(self, beam, **kwargs):  # todo: common implementation it here?
        raise NotImplementedError("Subclasses should implement this!")

    def apply_refraction_on_beam(self, beam, **kwargs):  # todo: common implementation it here?
        raise NotImplementedError("Subclasses should implement this!")

    def apply_grating_diffraction_on_beam(self, beam, **kwargs):  # todo: common implementation it here?
        raise NotImplementedError("Subclasses should implement this!")

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
    def write_shadow_surface(cls, s, xx, yy, outFile='presurface.dat', verbose=1):
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
        verbose : int, optional
            Set to 1 for verbose output.

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
            if verbose: print("write_shadow_surface: File for SHADOW " + outFile + " written to disk.")
        return out

if __name__ == "__main__":
    a = S4OpticalSurface()
    a.info()