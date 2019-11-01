import numpy
from syned.storage_ring.magnetic_structure import MagneticStructure
from scipy.ndimage import gaussian_filter1d

class MagneticStructure1DField(MagneticStructure):
    def __init__(self,y=None,B=None):
            self._y = y
            self._B = B

    @classmethod
    def initialize_from_file(cls,filename,skiprows=0):
        """
        Initializes an instance of the class from a file
        :param filename: String with file name
        :param skiprows: Number of rows to skip (for jumping comments)
        :return: An instance of MagneticStructure1DField
        """
        m = MagneticStructure1DField()
        m.load_from_file(filename,skiprows=skiprows)
        return m

    @classmethod
    def initialize_from_arrays(cls,B,y):
        """
        Initializes an instance of the class from arrays
        :param B: array with magnetic field
        :param y: array with abscissas
        :return: An instance of MagneticStructure1DField
        """
        m = MagneticStructure1DField()
        m.set_from_arrays(B,y)
        return m

    @classmethod
    def initialize_bending_magnet(cls,B0,length,npoints=2000,overdimension_factor=1.2,smooth_sigma=None):
        """
        Initializes an instance of the class containing a centered bending magnet
        :param B0: scalar with magnetic field
        :param length: bending magnet length
        :param npoints: number of points for the arrays
        :param overdimension_factor: a factor (larger than one, default 1.2) to add zeros at the bending magnet edges
        :param smooth_sigma: (default = None, not applied) a value with the sigma for Gaussian smooth (in pixels)
        :return: An instance of MagneticStructure1DField
        """

        m = MagneticStructure1DField()
        m.set_interval_and_zero_field(-0.5*length*overdimension_factor,
                                      0.5*length*overdimension_factor,
                                      npoints)
        m.add_bending_magnet(B0,length,0.0)
        if smooth_sigma is not None:
            m.smooth_edges(smooth_sigma)

        return m

    def get_magnetic_field(self):
        """
        Returns array with magnetic field
        :return: array with magnetic field
        """
        return self._B.copy()

    def get_abscissas(self):
        """
        Returns arrays with abscissas
        :return: array with abscissas
        """
        return self._y.copy()

    def reset(self):
        """
        Reset internal arrays
        :return:
        """
        self._y = None
        self._B = None

    def load_from_file(self,filename,skiprows=0):
        """
        Load a magnetic field from a file

        :param filename: String with file name
        :param skiprows: Number of rows to skip (for jumping comments)
        :return:
        """
        try:
            a = numpy.loadtxt(filename,skiprows=skiprows)
            if a.shape[0] > a.shape[1]:
                self._y = a[:,0]
                self._B = a[:,1]
            else:
                self._y = a[:,1]
                self._B = a[:,0]
        except:
            raise Exception("Failed to load file: %s%"%filename)

    def set_from_arrays(self,B,y):
        """
        Sets the abscissas axis and magnetic field from external arrays.

        :param B: Array with magnetic field in T
        :param y: Array with abscissas in meters
        :return:
        """
        if B.size != y.size:
            raise Exception("B and y must have the same dimension")
        self._B = B.copy()
        self._y = y.copy()

    def set_interval_and_zero_field(self,y_from,y_to,npoints):
        """
        Initializes internal arrays with abscissas axis and zero magnetic field
        :param y_from: start abscissa (in meters)
        :param y_to: end abscissa (in meters)
        :param npoints: number of points
        :return:
        """
        self._y = numpy.linspace(y_from,y_to,npoints)
        self._B = numpy.zeros_like(self._y)

    def add_bending_magnet(self,B0,bending_magnet_length,bending_magnet_center=0.0):
        """
        Defines a constant field (bending magnet) over an interval.
        Requires that the magnetic field had been initialized (with set_interval_and_zero_field() )
        :param B0: (scalar) The magnetic field value
        :param bending_magnet_length: Length of bending magnet in meters
        :param bending_magnet_center: Center of bending magnet in meters
        :return:
        """
        igood = numpy.argwhere( numpy.abs(self._y - bending_magnet_center) <= bending_magnet_length / 2 )
        if len(igood) > 0:
            self._B[igood] = B0


    def info(self):
        """
        Return a string some information
        :return: string
        """
        txt = ""
        txt += "Info on 1D profile: \n"
        txt += "  number of points: %d: \n"%self._y.size
        txt += "  y min: %f  max: %f: \n" % (self._y.min(), self._y.max())
        txt += "  B min: %f  max: %f: \n" % (self._B.min(), self._B.max())
        txt += "  Integral of B: %f \n" % numpy.trapz(self._B,self._y)

        return (txt)

    def flip_B(self):
        """
        Multiplies by -1 the magnetic field

        :return:
        """
        self._B *= -1.0

    def add_spatial_shift(self,shift):
        """
        Adds a shift to the ordinate axis (y)
        :param shift: the shift value to be added to y
        :return:
        """
        self._y += shift

    def smooth_edges(self,sigma=2.5):
        """
        Makes a Gaussian smoothing of the magnetic field.

        :param sigma: The sigma for the Gaussian smooth in steps (pixels)
        :return:
        """
        self._B = gaussian_filter1d(self._B, sigma)


if __name__ == "__main__":

    from srxraylib.plot.gol import plot, set_qt

    L = 1.605  # m


    o = MagneticStructure1DField()
    o.set_interval_and_zero_field(-L/2,L/2,2000)
    o.add_bending_magnet(-0.876,0.5,-0.5)
    o.add_bending_magnet(0.16, 0.35, 0.)
    o.add_bending_magnet(-0.8497, 0.5, +0.5)

    print(o.info())
    print(type(o))

    assert(isinstance(o,MagneticStructure1DField))
    assert(isinstance(o,MagneticStructure))


    #
    # plot
    #
    B0 = o.get_magnetic_field()
    o.smooth_edges(sigma=2.5)
    B1 = o.get_magnetic_field()
    y = o.get_abscissas()
    plot(y,B0,y,B1)


    # another example
    o2 = MagneticStructure1DField.initialize_bending_magnet(0.8,0.5,npoints=1000,smooth_sigma=10)
    plot(o2.get_abscissas(), o2.get_magnetic_field())