import numpy

from shadow4.syned.shape import Rectangle, Ellipse, TwoEllipses # TODO from syned.beamline.shape
from shadow4.syned.shape import Toroidal, Conic, NumericalMesh, Plane, Sphere, Ellipsoid, Paraboloid, Hyperboloid # TODO from syned.beamline.shape
from shadow4.syned.shape import SphericalCylinder # TODO from syned.beamline.shape


from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.beamline import BeamlineElement

from syned.beamline.optical_element_with_surface_shape import OpticalElementsWithSurfaceShape
from syned.beamline.optical_elements.mirrors.mirror import Mirror as SyMirror

from shadow4.optical_surfaces.conic import Conic as S4Conic
from shadow4.optical_surfaces.toroid import Toroid as S4Toroid
from shadow4.optical_surfaces.mesh import Mesh as S4Mesh
from shadow4.physical_models.prerefl.prerefl import PreRefl


class Mirror(object):

    def __init__(self, beamline_element_syned=None):
        if beamline_element_syned is None:

            self._beamline_element_syned = BeamlineElement(
                SyMirror(name="Undefined",
                 surface_shape=Conic(conic_coefficients=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0, 0.0]), #plane
                 boundary_shape=None,
                 coating=None,
                 coating_thickness=None),
                ElementCoordinates())

        else:

            # if not isinstance(beamline_element_syned._optical_element,SyMirror):
            #     raise Exception("Please initialize shadow4 Mirror with syned Mirror")
            #
            # if (isinstance(beamline_element_syned._optical_element._surface_shape, Conic)) or \
            #     (isinstance(beamline_element_syned._optical_element._surface_shape, Toroidal)) or \
            #     (isinstance(beamline_element_syned._optical_element._surface_shape, NumericalMesh)):
            #     pass
            # else:
            #     print(">>>",beamline_element_syned._optical_element._surface_shape)
            #     raise Exception("Only Conic, Toroid and NumericalMesh syned-surface-shapes are accepted by shadow4")
            #
            #
            # self._beamline_element_syned = beamline_element_syned
            if isinstance(beamline_element_syned._optical_element, SyMirror):
                pass
            else:
                raise Exception("Please initialize shadow4 Mirror with syned Mirror")

            ok = False
            for obj in [Conic, Toroidal, NumericalMesh, Plane, Sphere, Ellipsoid, Paraboloid, Hyperboloid, SphericalCylinder]:
                if isinstance(beamline_element_syned._optical_element._surface_shape, obj): ok = True
            if ok:
                self._beamline_element_syned = beamline_element_syned
            else:
                raise Exception(
                    "Please initialize shadow4 Mirror with syned Mirror surface shape as Conic, Toroidal, NumericalMesh, Plane, Sphere, Ellipsoid, Paraboloid, Hyperboloid, SphericalCylinder")



    def info(self):
        if self._beamline_element_syned is not None:
            return (self._beamline_element_syned.info())

    def trace_beam(self,beam_in,undo_shadow_orientation_angle_rotation=False):

        p = self._beamline_element_syned.get_coordinates().p()
        q = self._beamline_element_syned.get_coordinates().q()
        theta_grazing1 = numpy.pi / 2 - self._beamline_element_syned.get_coordinates().angle_radial()
        alpha1 = self._beamline_element_syned.get_coordinates().angle_azimuthal()

        #
        beam = beam_in.duplicate()

        #
        # put beam in mirror reference system
        #
        beam.rotate(alpha1, axis=2)
        beam.rotate(theta_grazing1, axis=1)
        beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        #
        # reflect beam in the mirror surface
        #
        soe = self._beamline_element_syned.get_optical_element()

        v_in = beam.get_columns([4,5,6])
        if not isinstance(soe, OpticalElementsWithSurfaceShape): # undefined
            raise Exception("Undefined mirror")
        else:
            surshape = self._beamline_element_syned.get_optical_element().get_surface_shape()
            if isinstance(surshape,Conic):
                print(">>>>> Conic mirror")
                conic = self._beamline_element_syned.get_optical_element().get_surface_shape()
                ccc = S4Conic.initialize_from_coefficients(conic._conic_coefficients)
                mirr, normal = ccc.apply_specular_reflection_on_beam(beam)
            elif isinstance(surshape, Toroidal):
                print(">>>>> Toroidal mirror",self._beamline_element_syned.get_optical_element().get_surface_shape()._min_radius,
                      self._beamline_element_syned.get_optical_element().get_surface_shape()._maj_radius)
                toroid = S4Toroid()
                toroid.set_toroid_radii( \
                    self._beamline_element_syned.get_optical_element().get_surface_shape()._maj_radius,
                    self._beamline_element_syned.get_optical_element().get_surface_shape()._min_radius,)
                mirr, normal = toroid.apply_specular_reflection_on_beam(beam)
            elif isinstance(surshape, NumericalMesh):
                print(">>>>> NumericalMesh mirror")
                num_mesh = S4Mesh()
                num_mesh.load_h5file(self._beamline_element_syned.get_optical_element().get_surface_shape()._h5file)
                mirr,normal,t,x1,v1,x2,v2 = num_mesh.apply_specular_reflection_on_beam(beam)
            elif isinstance(surshape, Plane):
                print(">>>>> Plane mirror")
                ccc = S4Conic.initialize_as_plane()
                mirr, normal = ccc.apply_specular_reflection_on_beam(beam)
            elif isinstance(surshape, SphericalCylinder):  # Note this check must come before Sphere as SphericalCylinder is Sphere
                print(">>>>> SphericalCylinder mirror")
                if surshape.get_direction() == 0:
                    cylangle = 0.0
                elif surshape.get_direction() == 1:
                    cylangle = numpy.pi / 2
                else:
                    raise Exception("Undefined cylinder direction")

                ccc = S4Conic.initialize_as_sphere_from_curvature_radius(surshape.get_radius(),cylindrical=True,cylangle=cylangle)
                mirr, normal = ccc.apply_specular_reflection_on_beam(beam)

            elif isinstance(surshape, Sphere):
                print(">>>>> Sphere mirror")
                ccc = S4Conic.initialize_as_sphere_from_curvature_radius(surshape.get_radius(),cylindrical=False,cylangle=0.0)
                mirr, normal = ccc.apply_specular_reflection_on_beam(beam)

            elif isinstance(surshape, Ellipsoid):
                print(">>>>> Ellipsoid mirror",surshape)
                raise Exception("Not enough information in Syned Ellipsoid object!!!")
            elif isinstance(surshape, Paraboloid):
                print(">>>>> Paraboloid mirror")
                raise Exception("Not enough information in Syned Paraboloid object!!!")
            elif isinstance(surshape, Hyperboloid):
                print(">>>>> Hyperboloid mirror")
                raise Exception("Not enough information in Syned Hyperboloid object!!!")
            else:
                raise Exception("cannot trace this surface shape")

        #
        # apply mirror boundaries
        #
        mirr.apply_boundaries_syned(self._beamline_element_syned.get_optical_element().get_boundary_shape())

        #
        # apply mirror reflectivity
        # TODO: add phase
        #

        coating = self._beamline_element_syned.get_optical_element()._coating
        print(">>>>>>>>>>>>>>>>>> COATING: ", coating)
        if coating is not None:
            try:
                prerefl_file = coating
                print(">>>>>>>>>>> PREREFL FILE", prerefl_file)
                pr = PreRefl()
                pr.read_preprocessor_file(prerefl_file)
                print(pr.info())
            except:
                raise Exception("the syned coating in mirror definition must contain the prerefl preprocessor file")

            v_out = beam.get_columns([4, 5, 6])
            angle_in = numpy.arccos( v_in[0,:] * normal[0,:] +
                                     v_in[1,:] * normal[1,:] +
                                     v_in[2,:] * normal[2,:])

            angle_out = numpy.arccos( v_out[0,:] * normal[0,:] +
                                     v_out[1,:] * normal[1,:] +
                                     v_out[2,:] * normal[2,:])

            grazing_angle_mrad = 1e3 * (numpy.pi / 2 - angle_in)

            Rs, Rp, Ru = pr.reflectivity_fresnel(grazing_angle_mrad=grazing_angle_mrad,
                                                 photon_energy_ev=beam.get_column(-11),
                                                 roughness_rms_A=0.0)

            beam.apply_reflectivities(numpy.sqrt(Rs), numpy.sqrt(Rp))


        #
        # TODO: write angle.xx for comparison
        #


        #
        # from mirror reference system to image plane
        #

        beam_out = mirr.duplicate()
        beam_out.rotate(theta_grazing1, axis=1)
        # do not undo alpha rotation: newbeam.rotate(-alpha, axis=2)
        if undo_shadow_orientation_angle_rotation:
            beam_out.rotate(-alpha1, axis=2)

        beam_out.retrace(q, resetY=True)

        return beam_out, mirr

    #
    # i/o utilities
    #
    def set_positions(self, p, q, theta_grazing, theta_azimuthal=None):
        self._beamline_element_syned.get_coordinates()._p = p
        self._beamline_element_syned.get_coordinates()._q = q
        self._beamline_element_syned.get_coordinates()._angle_radial = numpy.pi / 2 - theta_grazing
        if theta_azimuthal is not None:
            self._beamline_element_syned.get_coordinates()._angle_azimuthal = theta_azimuthal

    def set_surface_plane(self):
        self._beamline_element_syned._optical_element._surface_shape = \
            Conic(conic_coefficients=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0])

    def set_surface_conic(self, conic_coefficients=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0]):
        self._beamline_element_syned._optical_element._surface_shape = Conic(conic_coefficients=conic_coefficients)

    def set_surface_toroid(self, min_radius=0.0, maj_radius=0.0):
        self._beamline_element_syned._optical_element = Toroidal(min_radius=min_radius, maj_radius=maj_radius)

    def set_boundaries_rectangle(self, x_left=-1e3, x_right=1e3, y_bottom=-1e3, y_top=1e3):
        self._beamline_element_syned._optical_element._boundary_shape = \
            Rectangle(x_left=x_left, x_right=x_right, y_bottom=y_bottom, y_top=y_top)



if __name__ == "__main__":
    pass