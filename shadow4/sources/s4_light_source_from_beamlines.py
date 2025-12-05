"""
Defines the a LightSource merging several beamlines.
"""
import numpy
from syned.storage_ring.empty_light_source import EmptyLightSource

class S4LightSourceFromBeamlines(EmptyLightSource):
    """
    Class to create a light source from a S4Beam in an H5 file.

    Parameters
    ----------
    name : str, optional
        A name.
    beamlines : list, optiona
        A list with the beamlines (instances of S4Beamline).
    ids : list, optional
        A list with the ids or names (str).
    weights : list, optional
        A list with the weights (float).

    """
    def __init__(self, name="Undefined", beamlines=None, ids=None, weights=None):
        super().__init__(name=name)

        if  beamlines is None:
            self._beamlines = []
        else:
            if not isinstance(beamlines, list): raise Exception("Invalid beamlines. Must be a list")
            self._beamlines = beamlines

        if ids is None:
            self._set_default_ids()
        else:
            if not isinstance(ids, list): raise Exception("Invalid ids. Must be a list")
            if len(ids) != self.number_of_beamlines(): raise Exception("Invalid ids.")
            self._ids = ids

        if weights is None:
            self._set_default_weights()
        else:
            if not isinstance(weights, list): raise Exception("Invalid weights. Must be a list")
            if len(weights) != self.number_of_beamlines(): raise Exception("Invalid weights.")
            self._weights = weights

        self._set_support_text([
            ("name","Name",""),
            ("beamlines","Beamlines",""),
            ("weights", "Weights", ""),
            ("ids", "ID names", ""),
            ] )

    def number_of_beamlines(self):
        """
        returns the number of beamlines.

        Returns
        -------
        int
            The number of stored beamlines.
        """
        return len(self._beamlines)

    def append_beamline(self, beamline, id="appended beamline", weight=1.0):
        """
        Appends a beamline

        Parameters
        ----------
        beamline : instance of S4Beamline
            the beamline to be appended.
        id : str, optional
            the beamline id.
        weight : float, optional
            the intensity weigh for the beamline.

        """

        self._beamlines.append(beamline)
        self._ids.append(id)
        self._weights.append(weight)

    def _set_default_ids(self):
        self._ids = []
        n = self.number_of_beamlines()
        for i in range(n):
            self._ids.append("beamline_%d" % (i + 1))

    def _set_default_weights(self):
        self._weights = []
        n = self.number_of_beamlines()
        for i in range(n):
            self._weights.append(1.0)

    def _get_beams(self):
        out = []
        n = self.number_of_beamlines()
        for i in range(n):

            bl = self._beamlines[i]
            print(">>>>>>>>>>>>>>!!!!", bl)
            output_beam, output_mirr = bl.run_beamline()
            out.append(output_beam)
        return out

    def get_beam(self):
        """
        Retirns the S4 beam.

        Parameters
        ----------
        copy : int
            Returns the beam stored in the class (0) or a copy of it (1).

        Returns
        -------
            S4Beam instance
            The S4 beam.

        """
        beam = None
        beams = self._get_beams()
        for i in range(len(beams)):
            output_beam, output_mirr = self._beamlines[i].run_beamline()
            if i == 0:
                beam = output_beam
            else:
                beam.append_beam(output_beam, update_column_index=True)
                if self._weights[i] != 1.0:
                    beam.apply_attenuation(numpy.sqrt(self._weights[i])) # weight is intensity, attenuator is amplitude!
        return beam


    def to_python_code(self, **kwargs):
        """
        Returns the python code.

        Returns
        -------
        str
            The python code.
        """

        n = self.number_of_beamlines()
        txt = ""

        for i in range(n):
            script = self._beamlines[i].to_python_code()

            indented_script = '\n'.join('    ' + line for line in script.splitlines())

            txt += "def run_beamline_%d():\n" % (i + 1)
            txt += indented_script
            txt += "\n    return beamline"
            txt += "\n\n"



        txt += "\n#\n#\n#"
        txt += "\nfrom shadow4.sources.s4_light_source_from_beamlines import S4LightSourceFromBeamlines"
        txt += "\nlight_source = S4LightSourceFromBeamlines(name='%s')" % (self.get_name())
        txt += "\n"
        for i in range(n):
            txt += "\nlight_source.append_beamline(run_beamline_%d(), id='%s', weight=%f)" % \
                   (i + 1, self._ids[i], self._weights[i])

        txt += "\nbeam = light_source.get_beam()"
        return txt

    def get_info(self):
        """
        Returns specific information.

        Returns
        -------
        str
        """
        n = self.number_of_beamlines()
        txt = "LightSource created by merging %d beamlines:" % n

        for i in range(n):
            txt += "\n beamline %d, id: '%s', weight: %f" % (i + 1, self._ids[i], self._weights[i])

        return txt

if __name__ == "__main__":
    def get_beamline():
        from shadow4.beamline.s4_beamline import S4Beamline

        beamline = S4Beamline()
        #
        #
        #
        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
        light_source = SourceGeometrical(name='Geometrical Source', nrays=10000, seed=12345)
        light_source.set_spatial_type_gaussian(sigma_h=0.000553, sigma_v=0.000029)
        light_source.set_depth_distribution_off()
        light_source.set_angular_distribution_uniform(hdiv1=-0.000000, hdiv2=0.000000, vdiv1=-0.000030, vdiv2=0.000030)
        light_source.set_energy_distribution_singleline(10850.000000, unit='eV')
        light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
        beam = light_source.get_beam()

        beamline.set_light_source(light_source)

        # optical element number XX
        boundary_shape = None

        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
        optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator (1)',
                                   boundary_shape=boundary_shape,
                                   i_abs=0,  # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
                                   i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=31.15, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
        beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # test plot
        if 0:
            from srxraylib.plot.gol import plot_scatter
            plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
                         title='(Intensity,Photon Energy)', plot_histograms=0)
            plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1),
                         title='(X,Z) in microns')

        return beamline

    a = S4LightSourceFromBeamlines(beamlines=[get_beamline(), get_beamline()], ids=['11', '22'], weights=None)
    a.append_beamline(get_beamline())
    print(a.info())
    print(a.get_info())
    beam = a.get_beam()
    print(beam, beam.N)
    print('script: ', a.to_python_code())

