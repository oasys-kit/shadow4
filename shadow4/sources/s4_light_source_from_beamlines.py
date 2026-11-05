"""
Defines the a LightSource merging several beamlines.
"""
import numpy
from syned.storage_ring.empty_light_source import EmptyLightSource

class S4LightSourceFromBeamlines(EmptyLightSource):
    """
    Class to create a light source merging several beamlines.

    Parameters
    ----------
    name : str, optional
        A name.
    beamlines : list, optional
        A list with the beamlines (instances of S4Beamline).
    ids : list, optional
        A list with the ids or names (str).
    weights : list, optional
        A list with the weights (float).
    flag_homogeneous : int, optional
        1 if the beamlines have identical light sources (e.g., a parallel run where only the seed changes). In this
        case, calculate_spectrum() returns the spectrum of the light source of the first beamline (the weights are
        not applied). 0 (default) for heterogeneous beamlines (the spectrum cannot be calculated).

    """
    def __init__(self, name="Undefined", beamlines=None, ids=None, weights=None, flag_homogeneous=0):
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

        self._flag_homogeneous = flag_homogeneous

        self._set_support_text([
            ("name","Name",""),
            ("beamlines","Beamlines",""),
            ("weights", "Weights", ""),
            ("ids", "ID names", ""),
            ("flag_homogeneous", "Flag homogeneous beamlines", ""),
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
            output_beam, output_mirr = bl.run_beamline()
            out.append(output_beam)
        return out

    def get_beam(self):
        """
        Returns the S4 beam.

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
            script_i = self._beamlines[i].to_python_code_packed(
                add_head="def dry_run_beamline_%d():" % (i + 1),
                add_in_lightsource=None,
                add_return="return beamline",
                add_main="",
                filename="",
                dry_run=True)

            txt += script_i
            txt += "\n\n"


        txt += "\n#\n#\n#"
        txt += "\nfrom shadow4.sources.s4_light_source_from_beamlines import S4LightSourceFromBeamlines"
        txt += "\nlight_source = S4LightSourceFromBeamlines(name='%s', flag_homogeneous=%d)" % (self.get_name(), self._flag_homogeneous)
        txt += "\n"
        for i in range(n):
            txt += "\nlight_source.append_beamline(dry_run_beamline_%d(), id='%s', weight=%f)" % \
                   (i + 1, self._ids[i], self._weights[i])

        txt += "\nbeam = light_source.get_beam()"
        return txt

    def calculate_spectrum(self, **params):
        """
        Calculates the spectrum of the light source, for homogeneous beamlines (flag_homogeneous=1).

        Parameters
        ----------
        **params
            Parameters passed to the calculate_spectrum() method of the light source of the first beamline.

        Returns
        -------
        tuple
            The spectrum, as returned by the calculate_spectrum() method of the light source of the first beamline
            (the weights are not applied).

        Raises
        ------
        Exception
            If the beamlines are heterogeneous (flag_homogeneous=0), or if there are no beamlines.

        """
        if self.number_of_beamlines() == 0:
            raise Exception("Spectrum cannot be calculated: no beamlines")
        if self._flag_homogeneous:
            return self._beamlines[0].get_light_source().calculate_spectrum(**params)
        else:
            raise Exception("Spectrum cannot be calculated for heterogeneous beamlines")

    def get_info(self):
        """
        Returns specific information.

        Returns
        -------
        str
        """
        n = self.number_of_beamlines()
        if self._flag_homogeneous:
            yy = "homogeneous"
        else:
            yy = "heterogeneous"
        txt = "LightSource created by merging %d %s beamlines:" % (n, yy)

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
        # electron beam
        from shadow4.sources.s4_electron_beam import S4ElectronBeam
        electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
        electron_beam.set_sigmas_all(sigma_x=3.34603e-05, sigma_y=7.2815e-06, sigma_xp=4.5174e-06, sigma_yp=1.9407e-06)
        electron_beam.set_dispersion_all(0, 0, 0, 0)

        # magnetic structure
        from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian
        source = S4UndulatorGaussian(
            period_length=0.017,  # syned Undulator parameter (length in m)
            number_of_periods=117.647,  # syned Undulator parameter
            photon_energy=18000.0,  # Photon energy (in eV)
            delta_e=20.0,  # Photon energy width (in eV)
            ng_e=100,  # Photon energy scan number of points
            flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
            flag_energy_spread=0,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
            harmonic_number=1,  # harmonic number
            flag_autoset_flux_central_cone=1,  # value to set the flux peak
            flux_central_cone=668967299633458.4,  # value to set the flux peak
        )

        # light source
        from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource
        light_source = S4UndulatorGaussianLightSource(name='Undulator Gaussian', electron_beam=electron_beam,
                                                      magnetic_structure=source, nrays=10000, seed=5676561)
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

    a = S4LightSourceFromBeamlines(beamlines=[get_beamline(), get_beamline()], ids=['11', '22'], weights=None, flag_homogeneous=1)
    a.append_beamline(get_beamline())
    print(a.info())
    print(a.get_info())
    beam = a.get_beam()
    print(beam, beam.N)
    print('script: ', a.to_python_code())
    spectrum = a.calculate_spectrum()
    print('spectrum: %d arrays of %d points' % (len(spectrum), spectrum[0].size))

