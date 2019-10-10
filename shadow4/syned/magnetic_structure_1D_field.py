import numpy
from syned.storage_ring.magnetic_structure import MagneticStructure

class MagneticStructure1DField(MagneticStructure):
    def __init__(self,y=None,B=None):
            self.y = y
            self.B = B

    @classmethod
    def initialize_from_file(cls,filename,skiprows=0):
        m = MagneticStructure1DField()
        m.load_from_file(filename,skiprows=skiprows)
        return m

    def reset(self):
        self.y = None
        self.B = None

    def load_from_file(self,filename,skiprows=0):
        try:
            a = numpy.loadtxt(filename,skiprows=skiprows)
            if a.shape[0] > a.shape[1]:
                self.y = a[:,0]
                self.B = a[:,1]
            else:
                self.y = a[:,1]
                self.B = a[:,0]
        except:
            raise Exception("Failed to load file: %s%"%filename)

    def info(self):
        txt = ""
        txt += "Info on 1D profile: \n"
        txt += "  number of points: %d: \n"%self.y.size
        txt += "  y min: %f  max: %f: \n" % (self.y.min(), self.y.max())
        txt += "  B min: %f  max: %f: \n" % (self.B.min(), self.B.max())
        txt += "  Integral of B: %f \n" % numpy.trapz(self.B,self.y)

        return (txt)

    def flip_B(self):
        self.B *= -1.0

    def add_spatial_shift(self,shift):
        self.y += shift


if __name__ == "__main__":

    from srxraylib.plot.gol import plot, set_qt
    from scipy.ndimage import gaussian_filter1d


    set_qt()
    #
    # create a test file with magnetic field
    #

    L = 1605.0 #mm

    y = numpy.linspace(0,L, 2000)

    B = y * 0.0


    for i in range(y.size):
        if y[i] > 75 and y[i] < 575: B[i] = -0.876
        if y[i] > 650 and y[i] < 975: B[i] = 0.16
        if y[i] > 1030 and y[i] < 1530: B[i] = -0.8497

    # plot(y, B)



    B2 = gaussian_filter1d(B, 2.5)

    y -= y[y.size//2]
    y *= 1e-3

    plot(y, B, y, B2, legend=["original","smoothed"],xtitle="y / m",ytitle="B / T")

    filename = "BM_smooth.b"
    f = open(filename, "w")
    for i in range(y.size):
        f.write("%f  %f\n" % (y[i], B2[i]))
    f.close()
    print("File written to disk: %s"%filename)


    #
    #
    #

    o = MagneticStructure1DField.initialize_from_file(filename)

    print(o.info())
    print(type(o))
    print(isinstance(o,MagneticStructure1DField))
    print(isinstance(o,MagneticStructure))

    plot(o.y,o.B)

    a = numpy.vstack((o.y,o.B)).T

    print(a)
    print(">>>>>>>>>>>",a.shape)


