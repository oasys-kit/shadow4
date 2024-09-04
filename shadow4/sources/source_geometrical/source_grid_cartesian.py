"""
Grid source defined in cartesian coordinates.
"""
import numpy
from shadow4.beam.s4_beam import S4Beam
from shadow4.sources.s4_light_source_base import S4LightSourceBase

from shadow4.tools.arrayofvectors import vector_default_efields

class SourceGridCartesian(S4LightSourceBase):
    """
    Defines a grid source, so points starting in a cube-like volume in real space and directions gridded in X,Z

    Parameters
    ----------
    real_space_width : list, optional
        the widths of the real_space volume (parallellepipedal) [Dx,Dy,Dz].
    direction_space_width : list, optional
        The "angular" aperture [Dx',Dz'].
    real_space_points : list, optional
        Number of points [Nx,Ny,Nz].
    direction_space_points : list, optional
        Number of points [Nx',Nz']
    real_space_center : list, optional
        Center coordinates in real space [Cx,Cy,Cz].
    direction_space_center : list, optional
        Center coordinates in divergence space [Cx',Cz']. Note that (Cx')^2+(Cz')^2 < 1.
    name : str, optional
        A name.
    nrays : int, optional
        Number of rays generated using SourceGaussian.get_beam()
    seed : int, optional
        Seed for the Monte Carlo generator.
    wavelength : float, optional
        The photon wavelength in m.
    polarization_degree : float
        The polarization degree (cos_s / (cos_s + cos_p).
    polarization_phase_deg : float, optional
        The polarization phase in degrees (0=linear).
    coherent_beam : int, optional
        (0) random (incoherent), or (1) constant (coherent) s-phases.
    """
    def __init__(self,
                 real_space_width=[1e-3,1e-3,1e-3],
                 direction_space_width=[0,0],
                 real_space_points=[10,10,10],
                 direction_space_points=[1,1],
                 real_space_center=[0,0,0],
                 direction_space_center=[0,0],
                 name="Undefined",
                 nrays=0, # not used
                 seed=0, # not used
                 wavelength=1e-10,
                 polarization_degree=1.0,
                 polarization_phase_deg=0.0,
                 coherent_beam=1,
                 ):
        super().__init__(name=name, nrays=nrays, seed=seed)
        self._real_space_width = real_space_width
        self._direction_space_width = direction_space_width
        self._real_space_points = real_space_points
        self._direction_space_points = direction_space_points
        self._real_space_center = real_space_center
        self._direction_space_center = direction_space_center
        self._wavelength = wavelength
        self._polarization_degree = polarization_degree
        self._polarization_phase_deg = polarization_phase_deg
        self._coherent_beam = coherent_beam

        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._set_support_text([
            ("real_space_width", "the widths of the real_space volume (parallellepipedal) [Dx,Dy,Dz].", ""),
            ("direction_space_width", "The angular aperture [Dx',Dz'].", ""),
            ("real_space_points", "Number of points [Nx,Ny,Nz].", ""),
            ("direction_space_points", "Number of points [Nx',Nz']", ""),
            ("real_space_center", "Center coordinates in real space [Cx,Cy,Cz].", ""),
            ("direction_space_center",
             "Center coordinates in divergence space [Cx',Cz']. Note that (Cx')^2+(Cz')^2 < 1.", ""),
            ("name", "the name", ""),
            ("nrays", "Number of rays generated using SourceGaussian.get_beam()", ""),
            ("seed", "Seed for the Monte Carlo generator.", ""),
            ("wavelength", "The photon wavelength in m.", ""),
            ("polarization_degree", "The polarization degree (cos_s / (cos_s + cos_p).", ""),
            ("polarization_phase_deg", "The polarization phase in degrees (0=linear).", ""),
            ("coherent_beam", "Random(incoherent) (0) or constant (coherent) s-phases.", ""),
        ] )


    @classmethod
    def initialize_point_source(cls,
                 direction_space_width=[1e-6,1e-6],
                 direction_space_points=[100,100],
                 direction_space_center=[0.0,0.0] ):
        """
        Initializes a point source (zero size).

        Parameters
        ----------
        direction_space_width : list, optional
            Interval for the direction space [Xwidth, Zwidth].
        direction_space_points : list, optional
            Number of points for the direction space [Nx, Nz].
        direction_space_center : list, optional
            Center for the direction space [Xcenter, Zcenter].

        Returns
        -------
        instance of SourceGridCartesian.

        """
        return SourceGridCartesian(real_space_width=[0,0,0],
                 direction_space_width=direction_space_width,
                 real_space_points=[1,1,1],
                 direction_space_points=direction_space_points,
                 real_space_center=[0.0,0.0,0.0],
                 direction_space_center=direction_space_center,)

    @classmethod
    def initialize_collimated_source(cls,
                 real_space_width=[1e-6,0.0,1e-6],
                 real_space_points=[100,1,100],
                 real_space_center=[0.0,0.0,0.0] ):
        """
        Initializes a collimated source (zero divergence).
        Parameters
        ----------
        real_space_width : list, optional
            real_space_width [Xwidth, Ywidth, Zwidth].
        real_space_points : list, optional
            real_space_points [Nx, Ny, Nz].
        real_space_center : list, optional
            real_space_points [Xcenter, Ycenter, Zcenter].

        Returns
        -------

        """
        return SourceGridCartesian(real_space_width=real_space_width,
                 direction_space_width=[0.0,0.0],
                 real_space_points=real_space_points,
                 direction_space_points=[1,1],
                 real_space_center=real_space_center,
                 direction_space_center=[0.0,0.0],)

    #
    # getters
    #

    def get_number_of_points(self):
        """
        Returns the total number of points or rays.

        Returns
        -------
        int
        """
        return self._real_space_points[0] * self._real_space_points[1] * self._real_space_points[2]* \
            self._direction_space_points[0] * self._direction_space_points[1]

    def get_arrays_real_space(self):
        """
        Returns three arrays with the sampled spatial coordinates.

        Returns
        -------
        tuple
            (x, y, z).

        """
        if self._real_space_points[0] <= 1:
            x = numpy.array([self._real_space_center[0]])
        else:
            x = numpy.linspace(-0.5*self._real_space_width[0],
                                0.5*self._real_space_width[0],
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
            z = numpy.linspace(-0.5*self._real_space_width[2],
                                0.5*self._real_space_width[2],
                               self._real_space_points[2]) + self._real_space_center[2]

        return x,y,z

    def get_arrays_direction_space(self):
        """
        Returns two arrays with the sampled angles (in fact, the components of the direction vector).

        Returns
        -------
        tuple
            (x', z')

        """
        if self._direction_space_points[0] <= 1:
            x = numpy.array([self._direction_space_center[0]])
        else:
            hdiv1 = 0.5*self._direction_space_width[0]
            hdiv2 = -0.5*self._direction_space_width[0]
            xmax1 = numpy.tan(hdiv1)
            xmax2 = numpy.tan(hdiv2)

            x = numpy.linspace(0,1,self._direction_space_points[0])
            x = x * (xmax1 - xmax2) + xmax2 + self._direction_space_center[0]

        if self._direction_space_points[1] <= 1:
            y = numpy.array([self._direction_space_center[1]])
        else:
            vdiv1 =  0.5*self._direction_space_width[1]
            vdiv2 = -0.5*self._direction_space_width[1]
            ymax1 = numpy.tan(vdiv1)
            ymax2 = numpy.tan(vdiv2)

            y = numpy.linspace(0,1,self._direction_space_points[1])
            y = y * (ymax1 - ymax2) + ymax2 + self._direction_space_center[1]

        return x,y


    def get_mesh_divergences(self):
        """
        Returns two mesh arrays (Nx, Nz) with the Xp and Zp values.

        Returns
        -------
        tuple
            (X',Z')
        """
        xp,zp = self.get_arrays_direction_space()

        XP = numpy.array(numpy.outer(xp,numpy.ones_like(zp)))
        YP = numpy.array(numpy.outer(numpy.ones_like(xp),zp))

        thetar = numpy.arctan(numpy.sqrt(XP*XP + YP*YP))
        phir = numpy.arctan2(YP,XP)


        return numpy.cos(phir) * numpy.sin(thetar), numpy.sin(phir) * numpy.sin(thetar)


    def get_mesh_real_space(self):
        """
        Returns two mesh arrays with the spatial cross section coordinates X,Z.

        Returns
        -------
        tuple
            (x,z)

        """
        x,y,z =  self.get_arrays_real_space()
        return numpy.array(numpy.outer(x,numpy.ones_like(z))), \
               numpy.array(numpy.outer(numpy.ones_like(x),z))


    def get_volume_divergences(self):
        """
        Returns an array (3,npoints) with xp,yp,zp (first index 0,1,2, respectively) with the direction vectors.

        Returns
        -------
        numpy array

        """
        XP,ZP = self.get_mesh_divergences()
        YP = numpy.sqrt(numpy.ones_like(XP) -XP**2 -ZP**2 )
        tmp = numpy.vstack((XP.flatten(),YP.flatten(),ZP.flatten()))
        return tmp

    def get_volume_real_space(self):
        """
        Returns an array (3,npoints) with x,y,z (first index 0,1,2, respectively) with the spatial coordinates.

        Returns
        -------
        numpy array
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
        Returns an array (6, npoints) with x,y,z,xp,yp,zp (first index 0,1,2,3,4,5 respectively) with the
        spatial and direction coordinates.

        Returns
        -------
        numpy array

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



    # get_info

    def get_info(self):
        """
        Returns an array of strings with info.

        Returns
        -------
        str
        """
        txt = ""

        txt += "Gridding in real space:      %d, %d, %d \n"%(self._real_space_points[0],
                                                                self._real_space_points[1],
                                                                self._real_space_points[2])
        txt += "Gridding in direction space: %d, %d \n"%(self._direction_space_points[0],
                                                                self._direction_space_points[1])
        txt += "\n"
        txt += "real_space_width "+repr(self._real_space_width) + "\n"
        txt += "direction_space_width "+repr(self._direction_space_width) + "\n"
        txt += "real_space_points "+repr(self._real_space_points) + "\n"
        txt += "direction_space_points "+repr(self._direction_space_points) + "\n"
        txt += "real_space_center "+repr(self._real_space_center) + "\n"
        txt += "direction_space_center "+repr(self._direction_space_center) + "\n"

        txt += "\nWavelength = %g m" % self._wavelength
        txt += "Degree of polarization is %f. Angular difference in phase is %f\n" % \
               (self._polarization_degree, self._polarization_phase_deg)


        return txt


    def get_beam(self):
        """
        Returns an instance of S4Beam with the sampled rays.

        Returns
        -------
        instance of S4Beam

        """
        N = self.get_number_of_points()
        rays = numpy.zeros((N, 18))
        rays[:, 0] = self.get_volume()[0, :]
        rays[:, 1] = self.get_volume()[1, :]
        rays[:, 2] = self.get_volume()[2, :]
        rays[:, 3] = self.get_volume()[3, :]
        rays[:, 4] = self.get_volume()[4, :]
        rays[:, 5] = self.get_volume()[5, :]
        rays[:,9] = 1   # flag
        rays[:,10] = 2 * numpy.pi / (self._wavelength * 1e2) # wavenumber in cm**-1
        rays[:,11] = numpy.arange(self.get_number_of_points(),dtype=float) # index
        if not self._coherent_beam:
            rays[:, 13] = numpy.random.random(N) * 2 * numpy.pi # Phase s
        rays[:, 14] = rays[:, 13] + numpy.radians(self._polarization_phase_deg) # Phase p

        DIREC = rays[:,3:6]
        A_VEC, AP_VEC = vector_default_efields(DIREC, pol_deg=self._polarization_degree)
        rays[:, 6:9] = A_VEC
        rays[:, 15:18] = AP_VEC

        return S4Beam.initialize_from_array(rays)

    def to_python_code(self):
        """
        Returns the python code to recreate the grid source.

        Returns
        -------
        str
            The python code.
        """

        txt = ""
        txt += "\n#\n#\n#"

        txt += "\nfrom shadow4.sources.source_geometrical.source_grid_cartesian import SourceGridCartesian"
        txt += "\nlight_source = SourceGridCartesian(name='%s', " % (self.get_name())
        txt += "\n   real_space_width = [%f, %f, %f]," % (tuple(self._real_space_width))
        txt += "\n   real_space_center = [%f, %f, %f]," % (tuple(self._real_space_center))
        txt += "\n   real_space_points = [%d, %d, %d]," % (tuple(self._real_space_points))
        txt += "\n   direction_space_width = [%f, %f]," % (tuple(self._direction_space_width))
        txt += "\n   direction_space_center = [%f, %f]," % (tuple(self._direction_space_center))
        txt += "\n   direction_space_points = [%d, %d]," % (tuple(self._direction_space_points))
        txt += "\n   wavelength = %g," % self._wavelength
        txt += "\n   polarization_degree = %g," % self._polarization_degree
        txt += "\n   polarization_phase_deg = %g," % self._polarization_phase_deg
        txt += "\n   coherent_beam = %d)" % self._coherent_beam
        txt += "\nbeam = light_source.get_beam()"

        return txt


if __name__ == "__main__":
    from srxraylib.plot.gol import plot_scatter

    a = SourceGridCartesian.initialize_point_source(
                direction_space_width  = [2e-3,2e-3],
                direction_space_points = [20,  20],
                direction_space_center = [0.0, 0.0] )
    print(a.get_info())

    x,y,z = a.get_arrays_real_space()
    print("x:",x)
    print("y:",y)
    print("z:",z)

    xp,zp = a.get_arrays_direction_space()

    XP,ZP = a.get_mesh_divergences()
    print("XP ZP.shape",XP.shape, ZP.shape)

    VP = a.get_volume_divergences()
    print("VP",VP.shape,VP.size)


    Vx = a.get_volume_real_space()
    print("Vx: ",Vx.shape)

    V = a.get_volume()
    print("V: ",V.shape)

    beam = a.get_beam()
    plot_scatter(beam.get_column(4), beam.get_column(6), plot_histograms=0, title="Point source. Cols 4,6")

    print("check orthogonality", beam.efields_orthogonal())


    #
    #
    #

    a = SourceGridCartesian.initialize_collimated_source(real_space_width=[10.,0.0,10.0],real_space_points=[100,1,100])
    print(a.info()) # syned
    print(a.get_info()) # local

    beam = a.get_beam()
    plot_scatter(beam.get_column(1), beam.get_column(3), plot_histograms=0, title="Collimated source. Cols 1,3")


    print(a.to_python_code())




