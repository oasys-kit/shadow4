import numpy

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.optical_elements.ideal_elements.ideal_fzp import IdealFZP

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_optical_element_decorators import S4OpticalElementDecorator
from shadow4.beamline.s4_beamline_element import S4BeamlineElement


class S4IdealFZP(IdealFZP, S4OpticalElementDecorator):
    """
    Defines an ideal Fresnel Zone Plate.

    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the optical element.
    focusing_direction : int
        0=None, 1=x (sagittal), 2=z (meridional), 3=2D focusing.
    focal : float
        The focal length in meters.
    nominal_wavelength : float
        The nominal wavelength in m for where the focal length is defined.
    diameter : float
        The FZP diameter in m.
    """
    def __init__(self,
                 name="Ideal FZP",
                 focusing_direction=3,      # 0=None, 1=x (sagittal), 2=z (meridional), 3=2D focusing.
                 focal=1.0,                 #  focal distance (m)
                 nominal_wavelength=1e-10,  # nominal wavelength in m
                 diameter=0.001,            # FZP diameter in m
                 ):
        super().__init__(name=name,
                         focusing_direction=focusing_direction,
                         focal=focal,
                         nominal_wavelength=nominal_wavelength,
                         # r0=r0,
                         diameter=diameter,
                         )

    def to_python_code(self, **kwargs):
        """
        Creates the python code for defining the element.

        Parameters
        ----------
        **kwargs

        Returns
        -------
        str
            Python code.
        """
        txt_pre = """

from shadow4.beamline.optical_elements.ideal_elements.s4_ideal_fzp import S4IdealFZP
optical_element = S4IdealFZP(name='{name:s}',
    focusing_direction={focusing_direction:d}, # 0=None, 1=x (sagittal), 2=z (meridional), 3=2D focusing.
    focal={focal:g}, # focal distance (m)
    nominal_wavelength={nominal_wavelength:g}, # nominal wavelength in m
    diameter={diameter:g}, # FZP diameter in m
    )
"""
        txt = txt_pre.format(**{ 'name': self.get_name(),
                                 'focusing_direction': self.focusing_direction(),
                                 'focal': self.focal(),
                                 'nominal_wavelength': self.nominal_wavelength(),
                                 'diameter': self.diameter(),
                                })
        return txt

class S4IdealFZPElement(S4BeamlineElement):
    """
    Constructor.

    Parameters
    ----------
    optical_element : instance of OpticalElement, optional
        The syned optical element.
    coordinates : instance of ElementCoordinates, optional
        The syned element coordinates.
    input_beam : instance of S4Beam, optional
        The S4 incident beam.

    Returns
    -------
    instance of S4IdealFZPElement.
    """
    def __init__(self,
                 optical_element : S4IdealFZP = None,
                 coordinates : ElementCoordinates = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4IdealFZP(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         input_beam=input_beam)

    def trace_beam(self, **params):
        """
        Runs (ray tracing) the input beam through the element.

        Parameters
        ----------
        **params

        Returns
        -------
        tuple
            (output_beam, footprint) instances of S4Beam.
        """
        input_beam = self.get_input_beam().duplicate()

        p, q = self.get_coordinates().get_p_and_q()
        if p != 0.0: input_beam.retrace(p, resetY=True)

        #
        ##############################################################################
        #

        output_beam = input_beam.duplicate()

        focusing_direction = self.get_optical_element().focusing_direction()

        if focusing_direction > 0:
            #
            lambda1 = input_beam.get_column(19)  # lambda in Angstroms
            lambda1m = lambda1 * 1e-10
            x = input_beam.get_column(1)
            z = input_beam.get_column(3)
            xpin = input_beam.get_column(4)
            zpin = input_beam.get_column(6)

            if focusing_direction == 1:
                z *= 0
            elif focusing_direction == 2:
                x *= 0

            #
            Kmod = 2 * numpy.pi / lambda1m  # wavevector modulus in m-1
            r = numpy.sqrt(x ** 2. + z ** 2.)  # distance to center
            Kxin = Kmod * xpin
            Kzin = Kmod * zpin

            nrays = input_beam.N

            n = numpy.zeros(nrays)
            d = numpy.zeros(nrays)

            #
            # calculate n (index of n-th zone) and d (radius if the nth zone minus
            # radius of the n-a zone)
            #
            # Rays that arrive onto the inner zone
            # IN are the indices of rays that arrive inside the inner zone

            R0 = self.get_optical_element().r0_exact()
            focal = self.get_optical_element().focal()
            nomlambda = self.get_optical_element().nominal_wavelength()
            DDm = self.get_optical_element().diameter()


            IN = numpy.where(r <= R0)
            IN = numpy.array(IN)
            if IN.size > 0:
                n[IN] = 0.0
                d[IN] = 0.0

            # Rays that arrive outside the inner zone
            # (see formulas in A.G. Michette, "X-ray science and technology"
            #  Institute of Physics Publishing (1993))

            OUT = numpy.where(r >= R0)
            OUT = numpy.array(OUT)
            if OUT.size > 0:
                n[OUT] = (r[OUT] ** 2 - R0 ** 2) / (nomlambda * focal)  # eq 8.56
                # d-spacing: we suppose the ray is in the middle of two r's with fractional n
                d[OUT] = numpy.sqrt((n[OUT] + .5) * nomlambda * focal + R0 ** 2) - \
                         numpy.sqrt((n[OUT] - .5) * nomlambda * focal + R0 ** 2)

            # computing G (the "grating" wavevector)

            Gx = -numpy.pi / d * (x / r)
            Gz = -numpy.pi / d * (z / r)
            # capture infinities
            Gx[d == 0] = 0.0
            Gz[d == 0] = 0.0

            # computing kout

            Kxout = Kxin + Gx
            Kzout = Kzin + Gz
            xpout = Kxout / Kmod
            zpout = Kzout / Kmod

            # Handle rays that arrive outside the FZP
            # flag for lost rays

            LOST = numpy.where(r > DDm / 2)
            LOST = numpy.array(LOST)
            if LOST.size > 0:
                output_beam.rays[LOST, 9] = -100.0

            output_beam.rays[:, 3] = xpout
            output_beam.rays[:, 4] = numpy.sqrt(1 - xpout ** 2 - zpout ** 2)
            output_beam.rays[:, 5] = zpout


        footprint = output_beam.duplicate()
        footprint.rotate(numpy.pi / 2, axis=1)

        return output_beam, footprint

    def to_python_code(self, **kwargs):
        """
        Creates the python code for defining the element.

        Parameters
        ----------
        **kwargs

        Returns
        -------
        str
            Python code.
        """
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        txt += self.to_python_code_coordinates()
        txt += "\nfrom shadow4.beamline.optical_elements.ideal_elements.s4_ideal_fzp import S4IdealFZPElement"
        txt += "\nbeamline_element = S4IdealFZPElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt

if __name__ == "__main__":

    if 0:
        # a = S4IdealFZP()
        # print(a.to_python_code())
        # a = S4SuperIdealFZP()
        # print(a.to_python_code())
        b = S4IdealFZPElement()
        print(b.to_python_code())
        # b = S4SuperIdealFZPElement()
        # print(b.to_python_code())


    #==============================

    if 1:
        from shadow4.beamline.s4_beamline import S4Beamline

        beamline = S4Beamline()

        #
        #
        #
        from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical
        light_source = SourceGeometrical(name='SourceGeometrical', nrays=5000, seed=5676561)
        light_source.set_spatial_type_rectangle(width=0.000800,height=0.000800)
        light_source.set_depth_distribution_off()
        light_source.set_angular_distribution_flat(hdiv1=-0.000000,hdiv2=0.000000,vdiv1=-0.000000,vdiv2=0.000000)
        light_source.set_energy_distribution_singleline(1.540000, unit='A')
        light_source.set_polarization(polarization_degree=1.000000, phase_diff=0.000000, coherent_beam=0)
        beam = light_source.get_beam()

        beamline.set_light_source(light_source)

        # optical element number XX
        boundary_shape = None

        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
        optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator', boundary_shape=boundary_shape,
            i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
            i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=1, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
        beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX

        from shadow4.beamline.optical_elements.ideal_elements.s4_ideal_fzp import S4IdealFZP
        optical_element = S4IdealFZP(name='Ideal FZP',
            focusing_direction=3, # 0=None, 1=x (sagittal), 2=z (meridional), 3=2D focusing.
            focal=1, # focal distance (m)
            nominal_wavelength=1.54e-10, # nominal wavelength in m
            diameter=0.000618, # FZP diameter in m
            )

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=0, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
        from shadow4.beamline.optical_elements.ideal_elements.s4_ideal_fzp import S4IdealFZPElement
        beamline_element = S4IdealFZPElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        # optical element number XX
        boundary_shape = None

        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4Screen
        optical_element = S4Screen(name='Generic Beam Screen/Slit/Stopper/Attenuator', boundary_shape=boundary_shape,
            i_abs=0, # 0=No, 1=prerefl file_abs, 2=xraylib, 3=dabax
            i_stop=0, thick=0, file_abs='<specify file name>', material='Au', density=19.3)

        from syned.beamline.element_coordinates import ElementCoordinates
        coordinates = ElementCoordinates(p=1, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
        from shadow4.beamline.optical_elements.absorbers.s4_screen import S4ScreenElement
        beamline_element = S4ScreenElement(optical_element=optical_element, coordinates=coordinates, input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)


        # test plot
        if True:
           from srxraylib.plot.gol import plot_scatter
           # plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1), title='(Intensity,Photon Energy)', plot_histograms=0)
           plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns')


