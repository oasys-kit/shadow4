
import numpy
from shadow4.beam.beam import Beam

class SourceGridCartesian(object):

    def __init__(self,
                 real_space,
                 direction_space,
                 real_space_points,
                 direction_space_points,
                 real_space_center,
                 direction_space_center,
                 ):
        """
        This defines a grid source, so points starting in a cube-like volume in real space
        and directions gridded in X,Z

        :param real_space: the widths of the real_space volume (parallellepipedal) [Dx,Dy,Dz]
        :param direction_space: The "angular" aperture [Dx',Dz']
        :param real_space_points: Number of points [Nx,Ny,Nz]
        :param direction_space_points: Number of points [Nx',Nz']
        :param real_space_center: Center in real space [Cx,Cy,Cz]
        :param direction_space_center:  Center in diraction space [Cx',Cz'] (note that (Cx')^2+(Cz')^2 < 1)
        :return:
        """
        self._real_space = real_space
        self._direction_space = direction_space
        self._real_space_points = real_space_points
        self._direction_space_points = direction_space_points
        self._real_space_center = real_space_center
        self._direction_space_center = direction_space_center



    @classmethod
    def initialize_point_source(cls,
                 direction_space=[1e-6,1e-6],
                 direction_space_points=[100,100],
                 direction_space_center=[0.0,0.0] ):
        """
        Initializes a point source

        :param direction_space: Default: [1e-6,1e-6],
        :param direction_space_points: Default: [100,100],
        :param direction_space_center: Default: [0.0,0.0] ):
        :return:
        """
        return SourceGridCartesian(real_space=[0,0,0],
                 direction_space=direction_space,
                 real_space_points=[1,1,1],
                 direction_space_points=direction_space_points,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=direction_space_center,)

    @classmethod
    def initialize_collimated_source(cls,
                 real_space=[1e-6,0.0,1e-6],
                 real_space_points=[100,1,100],
                 real_space_center=[0.0,0.0,0.0] ):
        """

        :param real_space: Default: [1e-6,0,1e-6]
        :param real_space_points: Default: [100,1,100]
        :param real_space_center: Default: [0.0,0.0,0.0]
        :return:
        """
        return SourceGridCartesian(real_space=real_space,
                 direction_space=[0.0,0.0],
                 real_space_points=real_space_points,
                 direction_space_points=[1,1],
                 real_space_center=real_space_center,
                 direction_space_center=[0.0,0.0],)

    #
    # getters
    #

    def get_number_of_points(self):
        """
        Returns the total number of points
        :return:
        """
        return self._real_space_points[0] * self._real_space_points[1] * self._real_space_points[2]* \
            self._direction_space_points[0] * self._direction_space_points[1]

    def get_arrays_real_space(self):
        """
        Returns three arrays with the spatial coordinates 1D arrays
        :return: x,y,z
        """
        if self._real_space_points[0] <= 1:
            x = numpy.array([self._real_space_center[0]])
        else:
            x = numpy.linspace(-0.5*self._real_space[0],
                                0.5*self._real_space[0],
                               self._real_space_points[0]) + self._real_space_center[0]

        if self._real_space_points[1] <= 1:
            y = numpy.array([self._real_space_center[1]])
        else:
            y = numpy.linspace(-0.5*self._real_space[1],
                                0.5*self._real_space[1],
                               self._real_space_points[1]) + self._real_space_center[1]

        if self._real_space_points[2] <= 1:
            z = numpy.array([self._real_space_center[2]])
        else:
            z = numpy.linspace(-0.5*self._real_space[2],
                                0.5*self._real_space[2],
                               self._real_space_points[2]) + self._real_space_center[2]

        return x,y,z

    def get_arrays_direction_space(self):
        """
        Returns two arrays with the direction angles (in fact components of the direction vector)
        :return: xp,zp
        """
        if self._direction_space_points[0] <= 1:
            x = numpy.array([self._direction_space_center[0]])
        else:
            hdiv1 = 0.5*self._direction_space[0]
            hdiv2 = -0.5*self._direction_space[0]
            xmax1 = numpy.tan(hdiv1)
            xmax2 = numpy.tan(hdiv2)

            # x = numpy.linspace(hdiv1,hdiv2,self._direction_space_points[0]) + self._direction_space_center[0]
            x = numpy.linspace(0,1,self._direction_space_points[0])
            x = x * (xmax1 - xmax2) + xmax2 + self._direction_space_center[0]

        if self._direction_space_points[1] <= 1:
            y = numpy.array([self._direction_space_center[1]])
        else:
            vdiv1 =  0.5*self._direction_space[1]
            vdiv2 = -0.5*self._direction_space[1]
            ymax1 = numpy.tan(vdiv1)
            ymax2 = numpy.tan(vdiv2)

            #y = numpy.linspace(vdiv1,vdiv2,self._direction_space_points[1]) + self._direction_space_center[1]
            y = numpy.linspace(0,1,self._direction_space_points[1])
            y = y * (ymax1 - ymax2) + ymax2 + self._direction_space_center[1]

        return x,y


    def get_mesh_divergences(self):
        """
        Returns two mesh arrays (Nx,Nz) with the Xp and Zp velues
        :return: Xp,Zp
        """
        xp,zp = self.get_arrays_direction_space()

        XP = numpy.array(numpy.outer(xp,numpy.ones_like(zp)))
        YP = numpy.array(numpy.outer(numpy.ones_like(xp),zp))

        thetar = numpy.arctan(numpy.sqrt(XP*XP + YP*YP))
        phir = numpy.arctan2(YP,XP)


        return numpy.cos(phir) * numpy.sin(thetar), numpy.sin(phir) * numpy.sin(thetar)


    def get_mesh_real_space(self):
        """
        Returns two mesh arrays with the spatial cross section coordinates X,Z
        :return: X,Z
        """
        x,y,z =  self.get_arrays_real_space()
        return numpy.array(numpy.outer(x,numpy.ones_like(z))), \
               numpy.array(numpy.outer(numpy.ones_like(x),z))


    def get_volume_divergences(self):
        """
        Returns an array (3,npoints) with xp,yp,zp (first index 0,1,2, respectively) with the
        direction vectors
        :return: xpypzp array
        """
        XP,ZP = self.get_mesh_divergences()
        YP = numpy.sqrt(numpy.ones_like(XP) -XP**2 -ZP**2 )
        tmp = numpy.vstack((XP.flatten(),YP.flatten(),ZP.flatten()))
        return tmp

    def get_volume_real_space(self):
        """
        Returns an array (3,npoints) with x,y,z (first index 0,1,2, respectively) with the
        spatial coordinates
        :return: xyz
        """
        x,y,z = self.get_arrays_real_space()
        x.flatten()
        y.flatten()
        z.flatten()
        X = numpy.outer(x,numpy.ones_like(y))
        Y = numpy.outer(numpy.ones_like(x),y)
        X.flatten()
        Y.flatten()
        XX = numpy.outer(X,numpy.ones_like(z))
        YY = numpy.outer(Y,numpy.ones_like(z))
        ZZ = numpy.outer(numpy.ones_like(X),z)
        return numpy.vstack((XX.flatten(),YY.flatten(),ZZ.flatten()))

    def get_volume(self):
        """
        Returns an array (6,npoints) with x,y,z,xp,yp,zp (first index 0,1,2,3,4,5 respectively) with the
        spatial and direction coordinates
        :return: xyzxpypzp array
        """
        v1 = self.get_volume_real_space()
        v2 = self.get_volume_divergences()

        v1x = v1[0,:].copy().flatten()
        v1y = v1[1,:].copy().flatten()
        v1z = v1[2,:].copy().flatten()
        v2x = v2[0,:].copy().flatten()
        v2y = v2[1,:].copy().flatten()
        v2z = v2[2,:].copy().flatten()

        V1x = numpy.outer(v1x,numpy.ones_like(v2x)).flatten()
        V1y = numpy.outer(v1y,numpy.ones_like(v2x)).flatten()
        V1z = numpy.outer(v1z,numpy.ones_like(v2x)).flatten()

        V2x = numpy.outer(numpy.ones_like(v1x),v2x).flatten()
        V2y = numpy.outer(numpy.ones_like(v1x),v2y).flatten()
        V2z = numpy.outer(numpy.ones_like(v1x),v2z).flatten()

        return numpy.vstack((V1x,V1y,V1z,V2x,V2y,V2z))


    #
    # info
    #
    def info(self):
        """
        Returns an array of strings with info.
        :return:
        """
        txt = ""

        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
        txt += "Grid source: \n"
        txt += "Number of points: %d \n"%self.get_number_of_points()
        txt += "Gridding in real space:      %d, %d, %d \n"%(self._real_space_points[0],
                                                                self._real_space_points[1],
                                                                self._real_space_points[2])
        txt += "Gridding in direction space: %d, %d \n"%(self._direction_space_points[0],
                                                                self._direction_space_points[1])
        txt += "\n"
        txt += "real_space "+repr(self._real_space) + "\n"
        txt += "direction_space "+repr(self._direction_space) + "\n"
        txt += "real_space_points "+repr(self._real_space_points) + "\n"
        txt += "direction_space_points "+repr(self._direction_space_points) + "\n"
        txt += "real_space_center "+repr(self._real_space_center) + "\n"
        txt += "direction_space_center "+repr(self._direction_space_center) + "\n"

        txt += "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'n"

        return txt


    def get_beam(self,wavelength=1e-10):
        """
        Returns a Beam
        :param wavelength: the photon wavelength in m
        :return:
        """

        rays = numpy.zeros((self.get_number_of_points(),18))
        rays[:, 0] = self.get_volume()[0,:]
        rays[:, 1] = self.get_volume()[1,:]
        rays[:, 2] = self.get_volume()[2,:]
        rays[:, 3] = self.get_volume()[3,:]
        rays[:, 4] = self.get_volume()[4,:]
        rays[:, 5] = self.get_volume()[5,:]
        rays[:,6] = 1.0 # Es
        rays[:,9] = 1   # flag
        rays[:,10] = 2 * numpy.pi / (wavelength * 1e2) # wavenumber
        rays[:,11] = numpy.arange(self.get_number_of_points(),dtype=float) # index

        return Beam.initialize_from_array(rays)


if __name__ == "__main__":
    from srxraylib.plot.gol import plot_scatter, set_qt
    from shadow4.compatibility.beam3 import Beam3
    set_qt()

    a = SourceGridCartesian.initialize_point_source(
                direction_space        = [2e-3,2e-3],
                direction_space_points = [20,  20],
                direction_space_center = [0.0, 0.0] )
    print(a.info())

    x,y,z = a.get_arrays_real_space()
    print("x:",x)
    print("y:",y)
    print("z:",z)

    xp,zp = a.get_arrays_direction_space()
    # print("xp:",xp)
    # print("zp:",zp)

    XP,ZP = a.get_mesh_divergences()
    print("XP ZP.shape",XP.shape,ZP.shape)

    VP = a.get_volume_divergences()
    print("VP",VP.shape,VP.size)


    Vx = a.get_volume_real_space()
    print("Vx: ",Vx.shape)

    V = a.get_volume()
    print("V: ",V.shape)



    beam = a.get_beam()
    beam_shadow3 = Beam3.initialize_from_shadow4_beam( beam )

    beam_shadow3.write("begin.dat")

    import Shadow
    Shadow.ShadowTools.plotxy(beam_shadow3,4,6)


    #
    #
    a = SourceGridCartesian.initialize_collimated_source(real_space=[10.,0.0,10.0],real_space_points=[100,1,100])
    print(a.info())


    beam = a.get_beam()
    plot_scatter(beam.get_column(1), beam.get_column(3))

    beam_shadow3 = Beam3.initialize_from_shadow4_beam( beam )

    import Shadow
    Shadow.ShadowTools.plotxy(beam_shadow3,1,3)
    #

    beam_shadow3.write("begin.dat")