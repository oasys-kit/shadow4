"""
Grid source defined in polar coordinates.
"""
import numpy
from shadow4.beam.s4_beam import S4Beam
from shadow4.sources.s4_light_source_base import S4LightSourceBase

from shadow4.tools.arrayofvectors import vector_default_efields

class SourceGridPolar(S4LightSourceBase):

    def __init__(self,
                 real_space_width=[1e-6,0,1e-6],
                 direction_space_width=[1e-6,1e-6],
                 real_space_points=[100,36],
                 direction_space_points=[1,1],
                 real_space_center=[0,0,0],
                 direction_space_center=[0,1,0],
                 name="Undefined",
                 nrays=0, # not used
                 seed=0, # not used
                 wavelength=1e-10,
                 polarization_degree=1.0,
                 polarization_phase_deg=0.0,
                 coherent_beam=1,
                 ):
        """
        Defines a grid source, so points starting in a ellipsoid-like volume in real space and angularly gridded in X,Z.

        Parameters
        ----------
        real_space_width : list, optional
            The widths of the real_space_width volume [2a,2b,2c] of the ellipsoid.
        direction_space_width : list, optional
            The "angular" aperture [Dx',Dz'].
        real_space_points : list, optional
            Number of points [Nradial,Nangular].
        direction_space_points : list, optional
            Number of points [Nradial',Nangular'].
        real_space_center : list, optional
            Center cartesian coordinates in real space [Xc,Yc,Zc].
        direction_space_center : list, optional
            Center coordinates in divergence space [X'_c,Z'_c]. Note that (X')^2+(Z')^2 < 1.
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
        super().__init__(name=name, nrays=nrays, seed=seed)

        if real_space_width[1] != 0:
            raise Exception("Finite source depth not yet implemented!") # TODO: implement it!

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


    @classmethod
    def initialize_point_source(cls,
                real_space_center=[0.0, 0.0, 0.0],
                direction_space_width=[1e-6,1e-6],
                direction_space_points=[5,36],
                direction_space_center=[0.0,0.0],
                                ):
        """
        Initializes a point source.

        Parameters
        ----------
        direction_space_width : list, optional
            The "angular" aperture [Dx',Dz'].
        direction_space_points : list, optional
            Number of points [Nradial', Nangular'].
        real_space_center : list, optional
            Center cartesian coordinates in real space [Xc,Yc,Zc].
        direction_space_center : list, optional
            Center coordinates in divergence space [X'_c,Z'_c]. Note that (X')^2+(Z')^2 < 1.

        Returns
        -------
        instance of SourceGridPolar.

        """
        return SourceGridPolar(real_space_width=[0,0,0],
                 direction_space_width=direction_space_width,
                 real_space_points=[1,1,1],
                 direction_space_points=direction_space_points,
                 real_space_center=real_space_center,
                 direction_space_center=direction_space_center,
                )

    @classmethod
    def initialize_collimated_source(cls,
                real_space_width=[1e-6,0.0,1e-6],
                real_space_points=[100,36],
                real_space_center=[0.0,0.0,0.0],
                direction_space_center=[0.0,0.0],
                                ):
        """

        Parameters
        ----------
        real_space_width : list, optional
            The widths of the real_space_width volume [2a, 2b, 2c] of the ellipsoid.
        real_space_points : list, optional
            Number of points [Nradial, Nangular].
        real_space_center : list, optional
            Center cartesian coordinates in real space [Xc, Yc, Zc].
        direction_space_center : list, optional
            Center coordinates in divergence space [X'_c,Z'_c]. Note that (X')^2+(Z')^2 < 1.

        Returns
        -------
        instance of SourceGridPolar.
        """
        return SourceGridPolar(real_space_width=real_space_width,
                 direction_space_width=[0.0,0.0],
                 real_space_points=real_space_points,
                 direction_space_points=[1,1],
                 real_space_center=real_space_center,
                 direction_space_center=direction_space_center,)

    #
    # getters
    #

    def get_number_of_points_real_space(self):
        """
        Returns the number of points in real space.

        Returns
        -------
        int
        """
        return self._real_space_points[0] * self._real_space_points[1]

    def get_number_of_points_direction_space(self):
        """
        Returns the number of points in direction space.

        Returns
        -------
        int
        """
        return self._direction_space_points[0] * self._direction_space_points[1]

    def get_number_of_points(self):
        """
        Returns the total number of points.

        Returns
        -------
        int
        """
        return self.get_number_of_points_real_space() * self.get_number_of_points_direction_space()


    def _get_arrays_real_space(self):
        # Returns three arrays with the spatial coordinates 1D arrays
        radial_ratio = numpy.linspace(0, 1, self._real_space_points[0] )

        angular = numpy.arange(self._real_space_points[1]) / self._real_space_points[1]
        angular *= 2 * numpy.pi


        npoints = self.get_number_of_points_real_space()

        x = numpy.zeros(npoints)
        z = numpy.zeros(npoints)

        i = -1
        for radius_ratio in radial_ratio:
            for angle in angular:
                i += 1
                x[i] = self._real_space_width[0]/2 * radius_ratio * numpy.cos(angle)
                z[i] = self._real_space_width[2]/2 * radius_ratio * numpy.sin(angle)
        y = numpy.zeros_like(x)
        return x,y,z

    def _get_arrays_direction_space(self):
        # Returns two arrays with the direction angles (in fact components of the direction vector)
        radial_ratio = numpy.linspace(0, 1, self._direction_space_points[0] )

        angular = numpy.arange(self._direction_space_points[1]) / self._direction_space_points[1]
        angular *= 2 * numpy.pi

        npoints = self.get_number_of_points_direction_space() # self._real_space_points[0] * self._real_space_points[1]

        vx = numpy.zeros(npoints)
        vz = numpy.zeros(npoints)

        i = -1
        for radius_ratio in radial_ratio:
            for angle in angular:
                i += 1
                vx[i] = self._direction_space_width[0]/2 * radius_ratio * numpy.cos(angle)
                vz[i] = self._direction_space_width[0]/2 * radius_ratio * numpy.sin(angle)

        return vx,vz

        return numpy.array([0]), numpy.array([0])

    def _get_volume(self):
        # Returns an array (6,npoints) with x,y,z,xp,yp,zp (first index 0,1,2,3,4,5 respectively) with the
        # spatial and direction coordinates
        X, Y, Z = self._get_arrays_real_space()
        VX, VZ = self._get_arrays_direction_space()

        npoint = self.get_number_of_points()
        V1x = numpy.zeros(npoint)
        V1y = numpy.zeros(npoint)
        V1z = numpy.zeros(npoint)
        V2x = numpy.zeros(npoint)
        V2z = numpy.zeros(npoint)
        V2y = numpy.zeros(npoint)

        ij = -1
        for i in range(self.get_number_of_points_real_space()):
            for j in range(self.get_number_of_points_direction_space()):
                ij += 1
                V1x[ij] = X[i]
                V1y[ij] = Y[i]
                V1z[ij] = Z[i]
                V2x[ij] = VX[j]
                V2z[ij] = VZ[j]

        V1x += self._real_space_center[0]
        V1y += self._real_space_center[1]
        V1z += self._real_space_center[2]
        V2x += self._direction_space_center[0]
        V2z += self._direction_space_center[1]

        try:
            V2y = numpy.sqrt(1 - V2x ** 2 - V2z ** 2)
        except:
            raise Exception('Failed to normalize directions. Try smalled angular width and/or angular center.')

        return numpy.vstack((V1x,V1y,V1z,V2x,V2y,V2z))


    #
    # info
    #
    def get_info(self):
        """
        Returns an array of strings with info.

        Returns
        -------
        str
        """
        txt = ""
        txt += "Number of points: Nreal_space: %d, Ndirection_space: %d, Total: %d \n"%(
                                                self.get_number_of_points_real_space(),
                                                self.get_number_of_points_direction_space(),
                                                self.get_number_of_points())
        txt += "Gridding in real space:      %d, %d\n"%(self._real_space_points[0],
                                                                self._real_space_points[1])
        txt += "Gridding in direction space: %d, %d \n"%(self._direction_space_points[0],
                                                                self._direction_space_points[1])
        txt += "\n"
        txt += "real_space widths"+repr(self._real_space_width) + "\n"
        txt += "direction_space widths "+repr(self._direction_space_width) + "\n"
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
        volume = self._get_volume()
        rays[:, 0] = volume[0, :]
        rays[:, 1] = volume[1, :]
        rays[:, 2] = volume[2, :]
        rays[:, 3] = volume[3, :]
        rays[:, 4] = volume[4, :]
        rays[:, 5] = volume[5, :]
        rays[:,9] = 1   # flag
        rays[:,10] = 2 * numpy.pi / (self._wavelength * 1e2) # wavenumber in cm**-1
        rays[:,11] = numpy.arange(self.get_number_of_points(),dtype=float) # index

        if not self._coherent_beam:
            rays[:, 13] = numpy.random.random(N) * 2 * numpy.pi # Phase s
        rays[:, 14] = rays[:, 13] + numpy.radians(self._polarization_phase_deg) # Phase p

        DIREC = rays[:, 3:6]
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

        txt += "\nfrom shadow4.sources.source_geometrical.source_grid_polar import SourceGridPolar"
        txt += "\nlight_source = SourceGridPolar(name='%s', " % (self.get_name())
        txt += "\n   real_space_width = [%f, %f, %f]," % (tuple(self._real_space_width))
        txt += "\n   real_space_center = [%f, %f, %f]," % (tuple(self._real_space_center))
        txt += "\n   real_space_points = [%d, %d]," % (tuple(self._real_space_points))
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

    if True:
        a = SourceGridPolar.initialize_point_source(
                    direction_space_width  = [2e-3,2e-3],
                    direction_space_points = [20,  5],
                    direction_space_center = [0.0, 0.0] )
        print(a.get_info())

        beam = a.get_beam()
        plot_scatter(beam.get_column(4) * 1e6, beam.get_column(6) * 1e6, title="Polar grid in direction space")
        print("check orthogonality", beam.efields_orthogonal())
    #
    #
    if False:
        a = SourceGridPolar.initialize_collimated_source(
            real_space_width=[2e-6, 0.0, 1e-6],
            real_space_points=[10, 4],
            real_space_center=[0.0, 0.0, 0.0]
        )
        print(a.get_info())


        beam = a.get_beam()
        plot_scatter(beam.get_column(1)*1e6, beam.get_column(3)*1e6, xrange=[-1,1], yrange=[-1,1],
                     title="Polar grid in real space")

    if False:
        a = SourceGridPolar(
            real_space_width=[1e-6, 0.0, 1e-6],
            real_space_points=[2, 4],
            real_space_center=[0.0, 0.0, 0.0],
            direction_space_width=[2e-6, 2e-6],
            direction_space_points=[5, 5],
            direction_space_center=[0.0, 0.0])

        print(a.get_info())


        beam = a.get_beam()
        plot_scatter(beam.get_column(1)*1e6, beam.get_column(3)*1e6, xrange=[-1,1], yrange=[-1,1],
                     title="Real space")
        plot_scatter(beam.get_column(4)*1e6, beam.get_column(6)*1e6, xrange=[-1,1], yrange=[-1,1],
                     title="Directions space")
        plot_scatter(beam.get_column(1)*1e6, beam.get_column(4)*1e6, xrange=[-1,1], yrange=[-1,1],
                     title="Phase space X")
        plot_scatter(beam.get_column(3)*1e6, beam.get_column(6)*1e6, xrange=[-1,1], yrange=[-1,1],
                     title="Phase space Z")

    if False:
        a = SourceGridPolar(
            real_space_width=[2e-3, 0.0, 2e-3],
            real_space_points=[2, 8],
            real_space_center=[0.0, 0.0, 0.0],
            direction_space_width=[20e-3, 20e-3],
            direction_space_points=[3, 359],
            direction_space_center=[0.0, 0.0])

        print(a.get_info())


        beam = a.get_beam()
        plot_scatter(beam.get_column(1)*1e6, beam.get_column(3)*1e6, xrange=[-1.1e3,1.1e3], yrange=[-1.1e3,1.1e3], title="Real space")
        plot_scatter(beam.get_column(4)*1e6, beam.get_column(6)*1e6, xrange=[-11e3,11e3], yrange=[-11e3,11e3], title="Directions space")
        plot_scatter(beam.get_column(1)*1e6, beam.get_column(4)*1e6, xrange=[-1.1e3,1.1e3], yrange=[-11e3,11e3], title="Phase space X")
        plot_scatter(beam.get_column(3)*1e6, beam.get_column(6)*1e6, xrange=[-1.1e3,1.1e3], yrange=[-11e3,11e3], title="Phase space Z")

    if False:
        a = SourceGridPolar(
            real_space_width=[2e-3, 0.0, 2e-3],
            real_space_points=[2, 8],
            real_space_center=[0.0, 0.0, 0.0],
            direction_space_width=[20e-3, 20e-3],
            direction_space_points=[3, 359],
            direction_space_center=[0.0, 0.0])

        print(a.to_python_code())