
# abstract class defining the interfaces of the optical surfaces to be implemented in the children classes
# it also defines common utilities
import numpy
import os
import h5py
import time



class S4OpticalSurface(object):

    def info(self):
        raise NotImplementedError

    def duplicate(self):
        raise NotImplementedError

    def surface_height(self, x, y, **kwargs):
        raise NotImplementedError

    def get_normal(self, x, **kwargs):
        raise NotImplementedError

    def calculate_intercept(self, XIN, VIN, **kwargs):
        raise NotImplementedError

    def apply_specular_reflection_on_beam(self, beam, **kwargs):
        raise NotImplementedError

    def apply_refraction_on_beam(self, beam, **kwargs):
        raise NotImplementedError

    def apply_crystal_diffraction_bragg_symmetric_on_beam(self, beam, **kwargs):
        raise NotImplementedError

    #
    # common utilities
    #
    def write_mesh_file(self, x, y, filename="surface.dat"):
        X = numpy.outer(x, numpy.ones_like(y))
        Y = numpy.outer(numpy.ones_like(x), y)
        Z = self.surface_height(X, Y)
        write_shadow_surface(Z.T, x, y, outFile=filename)

    def write_mesh_h5file(self, x, y, filename="surface.h5",  subgroup_name="surface_file", overwrite=True):

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

# copied from shadow3 ShadowTools
def write_shadow_surface(s, xx, yy, outFile='presurface.dat'):
    """
      write_shadowSurface: writes a mesh in the SHADOW/presurface format
      SYNTAX:
           out = write_shadowSurface(z,x,y,outFile=outFile)
      INPUTS:
           z - 2D array of heights
           x - 1D array of spatial coordinates along mirror width.
           y - 1D array of spatial coordinates along mirror length.

      OUTPUTS:
           out - 1=Success, 0=Failure
           outFile - output file in SHADOW format. If undefined, the
                     file is names "presurface.dat"

    """
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
        # for each x element, the x value and the corresponding z(y)
        # profile
        for i in range(xx.size):
            tmps = ""
            for j in range(yy.size):
                tmps = tmps + "  " + repr(s[j, i])
            fs.write(' ' + repr(xx[i]) + " " + tmps)
            fs.write("\n")
        fs.close()
        print("write_shadow_surface: File for SHADOW " + outFile + " written to disk.")

if __name__ == "__main__":
    a = S4OpticalSurface()
    a.info()