import numpy

from syned.beamline.optical_element import OpticalElement

from shadow4.beamline.s4_optical_element_decorators import S4OpticalElementDecorator

from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements

from shadow4.tools.logger import is_verbose

from shadow4.beamline.optical_elements.crystals.s4_conic_crystal import S4ConicCrystal, S4ConicCrystalElement
from shadow4.beamline.optical_elements.mirrors.s4_conic_mirror import S4ConicMirror, S4ConicMirrorElement


class S4Compound(OpticalElement, S4OpticalElementDecorator):
    """
    Shadow4 Mirror Class
    This is a base class for mirrors.
    Use derived classes for plane or other curved mirror surfaces.

    Constructor.

    Parameters
    ----------
    name : str, optional
        The name of the mirror.
    boundary_shape : instance of BoundaryShape, optional
        The boundary shape of the mirror.
    surface_shape : instance of SurfaceShape, optional
        The surface shape of the mirror.
    f_reflec : int, optional
         the reflectivity of surface:
            - 0=no reflectivity,
            - 1=full polarization.
    f_refl : int, optional
        A flag to indicate the source of reflectivities:
            * 0=prerefl file,
            * 1=electric susceptibility,
            * 2=user defined file (1D angle in mrad, reflectivity),
            * 3=user defined file (1D energy in eV, reflectivity),
            * 4=user defined file (2D energy in eV, angle in mrad, reflectivity),
            * 5=direct calculation using xraylib,
            * 6=direct calculation using dabax.
    file_refl : str, optional
            name of user defined file (for f_refl=0).
    refraction_index : complex, optional
            complex scalar with refraction index n (for f_refl=1).
    material : str, optional
            string with material formula (for f_refl=5,6)
    density : float, optional
            material density in g/cm^3 (for f_refl=5,6)

    Returns
    -------
    instance of S4Mirror.
    """
    def __init__(self,
                 name="Undefined",
                 oe_list=None,
                 ):

        OpticalElement.__init__(self,
                        name=name,
                        )

        if oe_list is None:
            oe_list = []
        else:
            self._oe_list = oe_list


        # support text containg name of variable, help text and unit. Will be stored in self._support_dictionary
        self._add_support_text([
            ("oe_list",            "S4: list of optical elements",                 ""),
        ] )

    def get_info(self):
        """
        Returns the specific information of the S4 mirror optical element.

        Returns
        -------
        str
        """
        txt = "\n\n"
        txt += "COMPOUND OPTICAL ELEMENT\n"

        for i, element in enumerate(self._oe_list):
            txt += "   ELEMENT number %d\n" % (i + 1)
            txt += "\n" + element.get_info()

        return txt

    def to_python_code_boundary_shape(self):
        """
        Creates a code block with information of boundary shape.

        Returns
        -------
        str
            The text with the code.
        """
        txt = "" # "\nfrom shadow4.beamline.optical_elements.mirrors.s4_plane_mirror import S4PlaneMirror"
        # bs = self._boundary_shape
        # if bs is None:
        #     txt += "\nboundary_shape = None"
        # elif isinstance(bs, Rectangle):
        #     txt += "\nfrom syned.beamline.shape import Rectangle"
        #     txt += "\nboundary_shape = Rectangle(x_left=%g, x_right=%g, y_bottom=%g, y_top=%g)" % bs.get_boundaries()
        # elif isinstance(bs, Ellipse):
        #     txt += "\nfrom syned.beamline.shape import Ellipse"
        #     txt += "\nboundary_shape = Ellipse(a_axis_min=%g, a_axis_max=%g, b_axis_min=%g, b_axis_max=%g)" % bs.get_boundaries()
        return txt


class S4CompoundElement(S4BeamlineElement):
    """
    The base class for Shadow4 mirror element.
    It is made of a S4Mirror and an ElementCoordinates instance. It also includes the input beam.

    Use derived classes for plane or other curved crystal surfaces.

    Constructor.

    Parameters
    ----------
    optical_element : instance of OpticalElement, optional
        The syned optical element.
    coordinates : instance of ElementCoordinates, optional
        The syned element coordinates.
    movements : instance of S4BeamlineElementMovements, optional
        The S4 element movements.
    input_beam : instance of S4Beam, optional
        The S4 incident beam.

    Returns
    -------
    instance of S4MirrorElement.
    """
    def __init__(self,
                 optical_element : S4Compound = None,
                 coordinates : ElementCoordinates = None,
                 movements : S4BeamlineElementMovements = None,
                 input_beam : S4Beam = None):
        super().__init__(optical_element=optical_element if optical_element is not None else S4Compound(),
                         coordinates=coordinates if coordinates is not None else ElementCoordinates(),
                         movements=movements,
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
        flag_lost_value = params.get("flag_lost_value", -1)
        change_reference_system_in = params.get("change_reference_system_in", True)
        change_reference_system_out = params.get("change_reference_system_out", True)

        if is_verbose():
            if not change_reference_system_in:
                print("change_reference_system_in = False: skipping reference change to o.e.")
            if not change_reference_system_out:
                print("change_reference_system_out = False: skipping reference change from o.e. to image")

        p = self.get_coordinates().p()
        q = self.get_coordinates().q()
        theta_grazing1 = numpy.pi / 2 - self.get_coordinates().angle_radial()
        theta_grazing2 = numpy.pi / 2 - self.get_coordinates().angle_radial_out()
        alpha1 = self.get_coordinates().angle_azimuthal()

        print(">>>> main coord: ", p, q, numpy.degrees(theta_grazing1), numpy.degrees(theta_grazing2), numpy.degrees(alpha1))
        #
        input_beam = self.get_input_beam().duplicate()

        #
        # put beam in mirror reference system
        #
        if change_reference_system_in:
            input_beam.rotate(alpha1, axis=2)
            input_beam.rotate(theta_grazing1, axis=1)
            input_beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])
            print(">>> compound changes to in")

        # mirror movement:
        movements = self.get_movements()
        if movements is not None:
            if movements.f_move:
                input_beam.rot_for(OFFX=movements.offset_x,
                                   OFFY=movements.offset_y,
                                   OFFZ=movements.offset_z,
                                   X_ROT=movements.rotation_x,
                                   Y_ROT=movements.rotation_y,
                                   Z_ROT=movements.rotation_z)


        #
        # reflect beam in the mirror surface
        #
        elements = self.get_optical_element()._oe_list

        print(">>> elements: ", type(elements), elements)

        beam = input_beam.duplicate()
        coordinates = ElementCoordinates(p=0, q=0,
                                         angle_radial=0,
                                         angle_azimuthal=0,
                                         angle_radial_out=0)

        for i, element in enumerate(elements):
            print(">>> element number %d : " % (i + 1), element)

            print("    >>>> vin: ", beam.get_columns([4, 5, 6])[:, 0])
            if isinstance(element, S4ConicMirror):
                print("    >>> Conic mirror")
                be = S4ConicMirrorElement(optical_element=element,
                                           coordinates=coordinates,
                                           movements=None,
                                           input_beam=beam,
                                           )
                beam, footprint = be.trace_beam(flag_lost_value=flag_lost_value,
                                                change_reference_system_in=False,
                                                change_reference_system_out=False)
            elif isinstance(element, S4ConicCrystal):
                print("    >>> Conic crystal")
                be = S4ConicCrystalElement(optical_element=element,
                                           coordinates=coordinates,
                                           movements=None,
                                           input_beam=beam,
                                           )
                beam, footprint = be.trace_beam(flag_lost_value=flag_lost_value,
                                                change_reference_system_in=False,
                                                change_reference_system_out=False)

            else:
                raise Exception("Not implemented %s in CompoundElement", type(element))

            print("    >>>> intercept: ", footprint.get_columns([1, 2, 3])[:, 0])
            print("    >>>> vout: ", footprint.get_columns([4, 5, 6])[:, 0])
        #
        # from mirror reference system to image plane
        #

        output_beam = footprint.duplicate()
        if change_reference_system_out:
            output_beam.change_to_image_reference_system(theta_grazing2, q)
            print(">>> compound changes to out")

        return output_beam, footprint

def get_optical_element_instance_channel_cut():
    try:    name = self.getNode().title
    except: name = "Channel Cut Crystal Monochromator"

    boundary_shape = None

    from shadow4.beamline.optical_elements.crystals.s4_conic_crystal import S4ConicCrystal

    crystal_separation = 0.0005

    ccc1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.5 * crystal_separation]
    ccc2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.5 * crystal_separation]
    pitch = 17e-6
    roll = 0
    yaw = 0

    R_I = [[1,0,0],
         [0,1,0],
         [0,0,1]]

    T = [0,0,0]

    #Rx, beta
    R_pitch =  [[1, 0,                0],
                [0, numpy.cos(pitch), -numpy.sin(pitch)],
                [0, numpy.sin(pitch), numpy.cos(pitch)]]

    # Ry, gamma
    R_roll = [  [numpy.cos(roll), 0, numpy.sin(roll)],
                [0,               1, 0              ],
                [-numpy.sin(roll),0, numpy.cos(roll)]   ]

    # Rz, alpha
    R_yaw = [  [numpy.cos(yaw), -numpy.sin(yaw), 0],
               [numpy.sin(yaw),  numpy.cos(yaw), 0],
               [0,               0,              1] ]


    # R = Rz Rx Ry
    R = numpy.array(R_yaw) @ numpy.array(R_pitch) @ numpy.array(R_roll)

    print(">>> R: ", R)

    conic_coefficients1 = rotate_and_translate_coefficients(ccc1, R, T)

    conic_coefficients2 = rotate_and_translate_coefficients(ccc2, R, T)

    optical_element1 = S4ConicCrystal(name='Generic Crystal',
                                      boundary_shape=boundary_shape,
                                      conic_coefficients=conic_coefficients1,
                                      material='Si', miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                      f_bragg_a=False, asymmetry_angle=0.0,
                                      is_thick=1, thickness=0.001,
                                      f_central=1, f_phot_cent=0, phot_cent=5000.0,
                                      file_refl='bragg.dat',
                                      f_ext=0,
                                      material_constants_library_flag=0,
                                      # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                      )

    optical_element2 = S4ConicCrystal(name='Generic Crystal',
                                      boundary_shape=boundary_shape,
                                      conic_coefficients=conic_coefficients2,
                                      material='Si', miller_index_h=1, miller_index_k=1, miller_index_l=1,
                                      f_bragg_a=False, asymmetry_angle=0.0,
                                      is_thick=1, thickness=0.001,
                                      f_central=0, f_phot_cent=0, phot_cent=5000.0,
                                      file_refl='bragg.dat',
                                      f_ext=0,
                                      material_constants_library_flag=0,
                                      # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
                                      )
    #
    #
    #
    return S4Compound(name=name, oe_list=[optical_element1, optical_element2])

import numpy as np
def rotate_and_translate_coefficients(coe_list, R_M, T):
    R_M = numpy.array(R_M, dtype=float)
    T = numpy.array(T, dtype=float)

    axx, ayy, azz, axy, ayz, axz, ax, ay, az, a0 = coe_list
    A2 = np.array([[axx,axy/2,axz/2],
    [axy/2,ayy,ayz/2],
    [axz/2,ayz/2,azz]])
    A1 = np.array([ax,ay,az])
    A0 = a0
    B2 = np.dot(R_M, np.dot(A2,R_M.T)) # first equation 6.29
    B1 = np.dot(R_M, A1) - 2 * np.dot(B2, T) # 2nd equation 6.29
    B0 = A0 + np.dot(T.T, (np.dot(B2, T) - \
                           np.dot(R_M, A1))) # 3rd equation 6.29
    return [B2[0,0], B2[1,1], B2[2,2],
            B2[0,1] + B2[1,0], B2[1,2] + B2[2,1], B2[0,2] + B2[2,0],
            B1[0], B1[1], B1[2], B0]

# # factory parameters (input)
# p = 10.0
# q = 3.0
# theta = 0.003
# # ellipse parameters
# a = (p + q) / 2
# b = np.sqrt(p * q) * np.sin(theta)
# c = np.sqrt(a**2 - b**2)
# # mirror center
# yc = (p ** 2 - q ** 2) / 4 / c
# zc = -b * np.sqrt(1 - yc ** 2 / a ** 2)
# X = np.array([0, yc, zc])
# # normal to the mirror at center
# N = np.array((0, -2 * yc / a ** 2, -2 * zc / b ** 2))
# n = N / np.sqrt((N**2).sum())
# # angle between N and Z
# Theta = np.arcsin(n[1])
# # rotation matrix
# R_M = np.array([[1,0,0],
# [0,np.cos(Theta),-np.sin(Theta)],
# [0,np.sin(Theta),np.cos(Theta)]])
# # translation vector
# T = -np.dot(R_M,X)
# # coefficients of the ellipsoid at the centered system
# c_in = [1/b**2,1/a**2,1/b**2,0,0,0,0,0,0,-1]
# # transformed coeffcients
# c_out = rotate_and_translate_coefficients(c_in,R_M,T)
# print("coeffs in centered frame: ", c_in)
# print("coeffs in local frame: ", c_out)
# # results of run
# # coeffs in centered frame:
# # [3703.715,0.02367,3703.715,0,0,0,0,0,0,-1]
# # coeffs in local frame:
# # [3703.715,0.0333353,3703.705,0,11.966,0,0,0,-102.564,0]

if __name__ == "__main__":
    def plot_2d(beam, footprint, xrange=None, yrange=None, show=1):
        from srxraylib.plot.gol import plot, plot_image, plot_image_with_histograms, plot_show

        #
        ticket = beam.histo2(1, 3, nbins_h=100, nbins_v=100, xrange=[-0.001789094021090414, 0.0017219218954817816],
                             yrange=[-0.001789094021090414, 0.0017219218954817816], nolost=1, ref=23)

        title = "BEAM I: %.1f " % ticket['intensity']
        if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
        if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

        plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                                   title=title, xtitle="column 1", ytitle="column 3",
                                   cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=0)
        #
        ticket = footprint.histo2(2, 1, nbins_h=100, nbins_v=100, xrange=[-0.00367611372770682, 0.003385538247030695],
                                  yrange=[-0.0017426214951288058, 0.0016776050111402614], nolost=1, ref=23)

        title = "FOOTPRINT I: %.1f " % ticket['intensity']
        if ticket['fwhm_h'] is not None: title += "FWHM H: %f " % ticket['fwhm_h']
        if ticket['fwhm_v'] is not None: title += "FWHM V: %f " % ticket['fwhm_v']

        plot_image_with_histograms(ticket['histogram'], ticket['bin_h_center'], ticket['bin_v_center'],
                                   title=title, xtitle="column 2", ytitle="column 1",
                                   cmap='jet', add_colorbar=True, figsize=(8, 8), histo_path_flag=1, show=1)

    # two parallel mirrors i) -z=0, ii) z=0
    if 0:
        from shadow4.beamline.s4_beamline import S4Beamline
        beamline = S4Beamline()

        # electron beam
        from shadow4.sources.s4_electron_beam import S4ElectronBeam

        electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
        electron_beam.set_sigmas_all(sigma_x=3.01836e-05, sigma_y=3.63641e-06, sigma_xp=4.36821e-06,
                                     sigma_yp=1.37498e-06)
        electron_beam.set_dispersion_all(0, 0, 0, 0)

        # magnetic structure
        from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian

        source = S4UndulatorGaussian(
            period_length=0.042,  # syned Undulator parameter (length in m)
            number_of_periods=38.571,  # syned Undulator parameter
            photon_energy=5000.0,  # Photon energy (in eV)
            delta_e=2.0,  # Photon energy width (in eV)
            ng_e=100,  # Photon energy scan number of points
            flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
            flag_energy_spread=0,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
            harmonic_number=1,  # harmonic number
            flag_autoset_flux_central_cone=0,  # value to set the flux peak
            flux_central_cone=10000000000.0,  # value to set the flux peak
        )

        # light source
        from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource

        light_source = S4UndulatorGaussianLightSource(name='Undulator Gaussian', electron_beam=electron_beam,
                                                      magnetic_structure=source, nrays=50000, seed=5676561)
        beam = light_source.get_beam()

        beamline.set_light_source(light_source)

        # optical element number XX
        boundary_shape = None

        from shadow4.beamline.optical_elements.mirrors.s4_conic_mirror import S4ConicMirror

        optical_element1 = S4ConicMirror(name='Generic Mirror', boundary_shape=boundary_shape,
                                        conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
                                        f_reflec=0, f_refl=5, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Ni', coating_density=8.902, coating_roughness=0)

        optical_element2 = S4ConicMirror(name='Generic Mirror', boundary_shape=boundary_shape,
                                        conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.0005],
                                        f_reflec=0, f_refl=5, file_refl='<none>', refraction_index=0.99999 + 0.001j,
                                        coating_material='Ni', coating_density=8.902, coating_roughness=0)

        from syned.beamline.element_coordinates import ElementCoordinates

        coordinates = ElementCoordinates(p=38.1, q=1.012,
                                         angle_radial=1.164204424,
                                         angle_azimuthal=0,
                                         angle_radial_out=numpy.radians(180) - 1.164204424)

        #
        #
        #

        ce = S4CompoundElement(optical_element=S4Compound(name='Compound',
                                                          oe_list=[optical_element1, optical_element2]),
                               coordinates=coordinates,
                               movements=None,
                               input_beam=beam)
        beam, footprint = ce.trace_beam()

        beamline.append_beamline_element(ce)

        plot_2d(beam, footprint)


    # channel-cut
    if 1:
        from shadow4.beamline.s4_beamline import S4Beamline

        beamline = S4Beamline()

        # electron beam
        from shadow4.sources.s4_electron_beam import S4ElectronBeam

        electron_beam = S4ElectronBeam(energy_in_GeV=6, energy_spread=0.001, current=0.2)
        electron_beam.set_sigmas_all(sigma_x=3.01836e-05, sigma_y=3.63641e-06, sigma_xp=4.36821e-06,
                                     sigma_yp=1.37498e-06)
        electron_beam.set_dispersion_all(0, 0, 0, 0)

        # magnetic structure
        from shadow4.sources.undulator.s4_undulator_gaussian import S4UndulatorGaussian

        source = S4UndulatorGaussian(
            period_length=0.042,  # syned Undulator parameter (length in m)
            number_of_periods=38.571,  # syned Undulator parameter
            photon_energy=5000.0,  # Photon energy (in eV)
            delta_e=2.0,  # Photon energy width (in eV)
            ng_e=100,  # Photon energy scan number of points
            flag_emittance=1,  # when sampling rays: Use emittance (0=No, 1=Yes)
            flag_energy_spread=0,  # when sampling rays: Use e- energy spread (0=No, 1=Yes)
            harmonic_number=1,  # harmonic number
            flag_autoset_flux_central_cone=0,  # value to set the flux peak
            flux_central_cone=10000000000.0,  # value to set the flux peak
        )

        # light source
        from shadow4.sources.undulator.s4_undulator_gaussian_light_source import S4UndulatorGaussianLightSource

        light_source = S4UndulatorGaussianLightSource(name='Undulator Gaussian', electron_beam=electron_beam,
                                                      magnetic_structure=source, nrays=50000, seed=5676561)
        beam = light_source.get_beam()

        beamline.set_light_source(light_source)

        # # optical element number XX
        # boundary_shape = None
        #
        # from shadow4.beamline.optical_elements.crystals.s4_conic_crystal import S4ConicCrystal
        #
        # optical_element1 = S4ConicCrystal(name='Generic Crystal',
        #                                  boundary_shape=boundary_shape,
        #                                  conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
        #                                  material='Si', miller_index_h=1, miller_index_k=1, miller_index_l=1,
        #                                  f_bragg_a=False, asymmetry_angle=0.0,
        #                                  is_thick=1, thickness=0.001,
        #                                  f_central=1, f_phot_cent=0, phot_cent=5000.0,
        #                                  file_refl='bragg.dat',
        #                                  f_ext=0,
        #                                  material_constants_library_flag=0,
        #                                  # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
        #                                  )
        #
        # optical_element2 = S4ConicCrystal(name='Generic Crystal',
        #                                  boundary_shape=boundary_shape,
        #                                  conic_coefficients=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.0005],
        #                                  material='Si', miller_index_h=1, miller_index_k=1, miller_index_l=1,
        #                                  f_bragg_a=False, asymmetry_angle=0.0,
        #                                  is_thick=1, thickness=0.001,
        #                                  f_central=0, f_phot_cent=0, phot_cent=5000.0,
        #                                  file_refl='bragg.dat',
        #                                  f_ext=0,
        #                                  material_constants_library_flag=0,
        #                                  # 0=xraylib,1=dabax,2=preprocessor v1,3=preprocessor v2
        #                                  )

        optical_element = get_optical_element_instance_channel_cut()

        from syned.beamline.element_coordinates import ElementCoordinates

        coordinates = ElementCoordinates(p=38.1, q=1,
                                         angle_radial=numpy.radians(66.703995),
                                         angle_azimuthal=0,
                                         angle_radial_out=numpy.radians(180 - 66.703995))

        beamline_element = S4CompoundElement(
            optical_element=optical_element,
            coordinates=coordinates,
            movements=None,
            input_beam=beam)

        beam, footprint = beamline_element.trace_beam()

        beamline.append_beamline_element(beamline_element)

        plot_2d(beam, footprint)

