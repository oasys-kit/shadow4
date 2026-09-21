# SPDX-License-Identifier: CECILL-2.1
# Copyright (c) 2026 ESRF - the European Synchrotron

"""
Tests for shadow4.beamline.optical_elements.phase_deflectors.

Covers both the generic base (S4PhaseDeflector / S4PhaseDeflectorElement,
and the module-level kernel) and the concrete numerical-mesh element
(S4NumericalMeshPhaseDeflector / S4NumericalMeshPhaseDeflectorElement), which
replaces the old refractors-based S4ThinTransmission.
"""

import h5py
import numpy
import pytest

from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.shape import Rectangle

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
from shadow4.beamline.optical_elements.phase_deflectors.s4_phase_deflector import (
    S4PhaseDeflector,
    S4PhaseDeflectorElement,
    apply_phase_deflection,
    _validate_thickness_map,
    _as_ray_array,
    _interp2d_regular,
)
from shadow4.beamline.optical_elements.phase_deflectors.s4_numerical_mesh_phase_deflector import (
    S4NumericalMeshPhaseDeflector,
    S4NumericalMeshPhaseDeflectorElement,
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _write_h5_mesh(path, x_axis, y_axis, z_of_x_y):
    """
    Write a synthetic OASYS-style surface h5 file.

    z_of_x_y(X, Y) is evaluated on a meshgrid with indexing="ij" (shape
    (nx, ny)); the on-disk convention stores Z with shape (ny, nx) (verified
    empirically against shadow4.optical_surfaces.s4_mesh.S4Mesh).
    """
    X, Y = numpy.meshgrid(x_axis, y_axis, indexing="ij")  # shape (nx, ny)
    z_nx_ny = z_of_x_y(X, Y)
    with h5py.File(path, "w") as f:
        f["/surface_file/X"] = numpy.asarray(x_axis, dtype=float)
        f["/surface_file/Y"] = numpy.asarray(y_axis, dtype=float)
        f["/surface_file/Z"] = z_nx_ny.T.copy()  # (ny, nx) on-disk convention
    return z_nx_ny  # (nx, ny); z_nx_ny.T is the (ny, nx) "profile" convention


def _make_phase_deflector(tmp_path, z_of_x_y, nx=21, ny=15, half_width=1e-3,
                           half_height=2e-3, **kwargs):
    x_axis = numpy.linspace(-half_width, half_width, nx)
    y_axis = numpy.linspace(-half_height, half_height, ny)
    z_nx_ny = _write_h5_mesh(tmp_path / "mesh.h5", x_axis, y_axis, z_of_x_y)
    oe = S4NumericalMeshPhaseDeflector(name="pd", surface_data_file=str(tmp_path / "mesh.h5"), **kwargs)
    return oe, x_axis, y_axis, z_nx_ny


def _pencil_beam(n=1, k=1.0e11):
    """
    A beam of n rays at the origin, in the element-local reference frame used
    by ``apply_phase_deflection`` (transverse = X, Y; forward = -Z, matching
    what S4PhaseDeflectorElement.trace_beam's rotate/translate produces for
    angle_radial=0). Es=(1,0,0), Ep=0, alive.
    """
    beam = S4Beam.initialize_as_pencil(n)
    beam.set_column(5, 0.0)   # Vy = 0 (initialize_as_pencil() defaults to the +Y convention)
    beam.set_column(6, -1.0)  # Vz = -1 (element-local forward direction)
    beam.set_column(11, numpy.full(n, k))
    return beam


# ---------------------------------------------------------------------------
# S4PhaseDeflector: generic base class
# ---------------------------------------------------------------------------

class TestS4PhaseDeflectorBase:
    def test_default_construction_does_not_raise(self):
        S4PhaseDeflector()

    def test_rejects_invalid_f_r_ind(self):
        with pytest.raises(ValueError):
            S4PhaseDeflector(f_r_ind=4)

    def test_get_thickness_map_is_abstract(self):
        with pytest.raises(NotImplementedError):
            S4PhaseDeflector().get_thickness_map()

    def test_get_optical_surface_instance_is_abstract(self):
        with pytest.raises(NotImplementedError):
            S4PhaseDeflector().get_optical_surface_instance()

    def test_get_info_does_not_crash_without_a_mesh(self):
        txt = S4PhaseDeflector().get_info()
        assert "PHASE DEFLECTOR" in txt

    def test_get_optical_constants_user_mode_matches_manual_formula(self):
        oe = S4PhaseDeflector(f_r_ind=0, r_ind=1.0 - 3e-6, r_attenuation=12.0)
        k = numpy.array([1e11, 2e11])
        delta, beta = oe.get_optical_constants(photon_energy_eV=numpy.array([10000.0, 20000.0]), k=k)
        assert numpy.allclose(delta, 3e-6)
        assert numpy.allclose(beta, 12.0 / (2.0 * k))

    def test_support_dictionary_registers_new_fields(self):
        oe = S4PhaseDeflector(material="Be", density=1.85)
        for key in ("surface_shape", "material", "density", "f_r_ind", "r_ind", "r_attenuation", "file_r_ind"):
            assert key in oe._support_dictionary
        # regression: a DabaxXraylib instance caches a non-JSON-serializable
        # SpecFile once used; "dabax" must stay out of the support dictionary
        # (matching S4Interface/S4CRL/S4Lens), or to_json()/.info() can crash
        # after tracing in dabax mode.
        assert "dabax" not in oe._support_dictionary
        # base attrs from OpticalElement.__init__ must still be present too
        assert "name" in oe._support_dictionary
        assert "boundary_shape" in oe._support_dictionary


# ---------------------------------------------------------------------------
# low-level helper functions (generic array math, unaffected by the
# refractors -> phase_deflectors move)
# ---------------------------------------------------------------------------

class TestValidateThicknessMap:
    def _valid_dict(self):
        return {
            "profile": numpy.ones((3, 4)),
            "x_axis": numpy.linspace(-1, 1, 4),
            "y_axis": numpy.linspace(-1, 1, 3),
        }

    def test_accepts_valid_map(self):
        profile, x_axis, y_axis = _validate_thickness_map(self._valid_dict())
        assert profile.shape == (3, 4)

    def test_missing_key_raises_keyerror(self):
        d = self._valid_dict()
        del d["x_axis"]
        with pytest.raises(KeyError):
            _validate_thickness_map(d)

    def test_shape_mismatch_raises(self):
        d = self._valid_dict()
        d["profile"] = numpy.ones((3, 5))
        with pytest.raises(ValueError):
            _validate_thickness_map(d)

    def test_non_monotonic_axis_raises(self):
        d = self._valid_dict()
        d["x_axis"] = numpy.array([0.0, 1.0, 0.5, -1.0])
        with pytest.raises(ValueError):
            _validate_thickness_map(d)

    def test_negative_profile_raises(self):
        d = self._valid_dict()
        d["profile"][0, 0] = -1.0
        with pytest.raises(ValueError):
            _validate_thickness_map(d)

    def test_decreasing_axis_is_accepted(self):
        d = self._valid_dict()
        d["x_axis"] = d["x_axis"][::-1]
        profile, x_axis, y_axis = _validate_thickness_map(d)
        assert x_axis[0] > x_axis[-1]


class TestAsRayArray:
    def test_scalar_broadcast(self):
        arr = _as_ray_array("k", 5.0, 4)
        assert arr.shape == (4,)
        assert numpy.all(arr == 5.0)

    def test_correct_length_array_passthrough(self):
        values = numpy.array([1.0, 2.0, 3.0])
        arr = _as_ray_array("k", values, 3)
        assert numpy.allclose(arr, values)

    def test_wrong_length_raises(self):
        with pytest.raises(ValueError):
            _as_ray_array("k", numpy.array([1.0, 2.0]), 3)

    def test_non_1d_raises(self):
        with pytest.raises(ValueError):
            _as_ray_array("k", numpy.ones((2, 2)), 4)


class TestInterp2dRegular:
    def test_exact_on_linear_function(self):
        x_axis = numpy.linspace(0, 3, 4)
        y_axis = numpy.linspace(0, 2, 3)
        X, Y = numpy.meshgrid(x_axis, y_axis)  # shape (ny, nx)
        values = 2.0 * X + 3.0 * Y + 1.0

        xs = numpy.array([0.5, 1.5, 2.5])
        ys = numpy.array([0.5, 1.0, 1.5])
        out = _interp2d_regular(x_axis, y_axis, values, xs, ys)
        expected = 2.0 * xs + 3.0 * ys + 1.0
        assert numpy.allclose(out, expected)

    def test_out_of_bounds_is_nan(self):
        x_axis = numpy.linspace(0, 1, 3)
        y_axis = numpy.linspace(0, 1, 3)
        values = numpy.zeros((3, 3))
        out = _interp2d_regular(x_axis, y_axis, values, numpy.array([-1.0, 2.0]), numpy.array([0.5, 0.5]))
        assert numpy.all(numpy.isnan(out))

    def test_decreasing_axes_supported(self):
        x_axis = numpy.linspace(3, 0, 4)  # decreasing
        y_axis = numpy.linspace(0, 2, 3)
        X, Y = numpy.meshgrid(x_axis, y_axis)
        values = 2.0 * X + 1.0
        out = _interp2d_regular(x_axis, y_axis, values, numpy.array([1.5]), numpy.array([1.0]))
        assert numpy.allclose(out, [2.0 * 1.5 + 1.0])


# ---------------------------------------------------------------------------
# apply_phase_deflection: physics kernel
# ---------------------------------------------------------------------------

class TestApplyPhaseDeflection:
    def _flat_thickness(self, t, half=1e-3):
        x_axis = numpy.linspace(-half, half, 5)
        y_axis = numpy.linspace(-half, half, 5)
        profile = numpy.full((5, 5), t)
        return {"profile": profile, "x_axis": x_axis, "y_axis": y_axis}

    def test_flag_lost_value_must_be_negative(self):
        beam = _pencil_beam(1)
        thickness = self._flat_thickness(1e-6)
        with pytest.raises(ValueError):
            apply_phase_deflection(beam, thickness, k=1e11, delta=1e-6, beta=1e-8, flag_lost_value=0.0)

    def test_uniform_slab_beer_lambert_attenuation(self):
        n = 5
        beam = _pencil_beam(n, k=1e11)
        t = 2e-6
        thickness = self._flat_thickness(t)
        k = 1e11
        beta = 3e-9
        out = apply_phase_deflection(beam, thickness, k=k, delta=1e-6, beta=beta)

        expected_amplitude = numpy.exp(-k * beta * t)
        assert numpy.allclose(out.rays[:, 6], expected_amplitude)  # Es_x scaled
        assert numpy.allclose(out.rays[:, 7], 0.0)
        assert numpy.allclose(out.rays[:, 15:18], 0.0)  # Ep untouched (was zero)

    def test_uniform_slab_opd_and_phase(self):
        beam = _pencil_beam(3, k=1e11)
        t = 1.5e-6
        thickness = self._flat_thickness(t)
        k = 1e11
        delta = 4e-6
        out = apply_phase_deflection(beam, thickness, k=k, delta=delta, beta=0.0)

        expected_opd = -delta * t
        expected_phase = -k * delta * t
        assert numpy.allclose(out.rays[:, 12], expected_opd)
        assert numpy.allclose(out.rays[:, 13], expected_phase)
        assert numpy.allclose(out.rays[:, 14], expected_phase)

    def test_uniform_slab_no_angular_kick(self):
        beam = _pencil_beam(2, k=1e11)
        thickness = self._flat_thickness(1e-6)
        out = apply_phase_deflection(beam, thickness, k=1e11, delta=5e-6, beta=1e-8)
        assert numpy.allclose(out.rays[:, 3], 0.0)   # Vx unchanged
        assert numpy.allclose(out.rays[:, 4], 0.0)   # Vy unchanged
        assert numpy.allclose(out.rays[:, 5], -1.0)  # Vz unchanged (unit vector already)

    def test_prism_ramp_deflects_towards_thinner_side(self):
        # t(x,y) = t0 + a*x, so dt/dx = a everywhere (exact for numpy.gradient on linear data).
        half = 1e-3
        x_axis = numpy.linspace(-half, half, 7)
        y_axis = numpy.linspace(-half, half, 7)
        t0 = 5e-3
        a = 2.0
        X, _ = numpy.meshgrid(x_axis, y_axis)  # shape (ny, nx)
        profile = t0 + a * X
        thickness = {"profile": profile, "x_axis": x_axis, "y_axis": y_axis}

        beam = _pencil_beam(1, k=1e11)
        beam.set_column(1, numpy.array([0.3e-3]))  # interior point, X
        beam.set_column(2, numpy.array([0.0]))  # Y (local-frame transverse axis)

        delta = 7e-6  # typical x-ray delta > 0 (n < 1)
        out = apply_phase_deflection(beam, thickness, k=1e11, delta=delta, beta=0.0)

        expected_kick_x = -delta * a
        assert numpy.isclose(out.rays[0, 3], expected_kick_x, rtol=1e-6)
        assert numpy.isclose(out.rays[0, 4], 0.0, atol=1e-12)  # no y-gradient
        # sanity: for delta>0 and a>0 (thickness increasing with +x) the kick bends toward -x
        assert expected_kick_x < 0.0

    def test_rays_outside_support_are_marked_lost_and_zeroed(self):
        thickness = self._flat_thickness(1e-6, half=1e-3)
        beam = _pencil_beam(2, k=1e11)
        beam.set_column(1, numpy.array([0.0, 5e-3]))  # second ray far outside the grid
        out = apply_phase_deflection(beam, thickness, k=1e11, delta=1e-6, beta=1e-8, flag_lost_value=-7.0)
        assert out.rays[0, 9] > 0.0  # still alive
        assert out.rays[1, 9] == -7.0  # newly lost
        assert numpy.allclose(out.rays[1, [6, 7, 8, 15, 16, 17]], 0.0)

    def test_apply_to_lost_false_leaves_already_lost_rays_untouched(self):
        thickness = self._flat_thickness(1e-6)
        beam = _pencil_beam(1, k=1e11)
        beam.set_column(10, numpy.array([-3.0]))  # already lost before this element
        original_es = beam.rays[0, 6]
        out = apply_phase_deflection(beam, thickness, k=1e11, delta=1e-6, beta=1e-8, apply_to_lost=False)
        assert out.rays[0, 6] == original_es  # untouched
        assert out.rays[0, 9] == -3.0  # flag preserved as-is

    def test_apply_to_lost_true_applies_kicks_but_preserves_flag(self):
        thickness = self._flat_thickness(1e-6)
        beam = _pencil_beam(1, k=1e11)
        beam.set_column(10, numpy.array([-3.0]))
        k, beta = 1e11, 5e-9
        out = apply_phase_deflection(beam, thickness, k=k, delta=1e-6, beta=beta, apply_to_lost=True)
        assert numpy.isclose(out.rays[0, 6], numpy.exp(-k * beta * 1e-6))
        assert out.rays[0, 9] == -3.0  # flag preserved, not reset to alive

    def test_duplicate_true_does_not_modify_input(self):
        thickness = self._flat_thickness(1e-6)
        beam = _pencil_beam(1, k=1e11)
        original = beam.rays.copy()
        out = apply_phase_deflection(beam, thickness, k=1e11, delta=1e-6, beta=1e-8, duplicate=True)
        assert out is not beam
        assert numpy.allclose(beam.rays, original)
        assert not numpy.allclose(out.rays, original)

    def test_duplicate_false_modifies_in_place(self):
        thickness = self._flat_thickness(1e-6)
        beam = _pencil_beam(1, k=1e11)
        out = apply_phase_deflection(beam, thickness, k=1e11, delta=1e-6, beta=1e-8, duplicate=False)
        assert out is beam


# ---------------------------------------------------------------------------
# S4NumericalMeshPhaseDeflector: construction, thickness map, code generation
# ---------------------------------------------------------------------------

class TestS4NumericalMeshPhaseDeflectorConstruction:
    def test_rejects_non_finite_thickness_scaling(self):
        with pytest.raises(ValueError):
            S4NumericalMeshPhaseDeflector(thickness_scaling=float("nan"))

    def test_rejects_non_positive_coordinate_scaling(self):
        with pytest.raises(ValueError):
            S4NumericalMeshPhaseDeflector(coordinate_scaling=0.0)
        with pytest.raises(ValueError):
            S4NumericalMeshPhaseDeflector(coordinate_scaling=-1.0)

    def test_rejects_invalid_f_r_ind(self):
        with pytest.raises(ValueError):
            S4NumericalMeshPhaseDeflector(f_r_ind=4)

    def test_default_construction_does_not_raise(self):
        S4NumericalMeshPhaseDeflector()

    def test_inherits_phase_deflector_and_numerical_mesh_decorator(self):
        from shadow4.beamline.s4_optical_element_decorators import S4NumericalMeshOpticalElementDecorator
        oe = S4NumericalMeshPhaseDeflector()
        assert isinstance(oe, S4PhaseDeflector)
        assert isinstance(oe, S4NumericalMeshOpticalElementDecorator)

    def test_support_dictionary_registers_mesh_specific_fields(self):
        oe = S4NumericalMeshPhaseDeflector()
        for key in ("surface_data_file", "apply_to_lost", "shift_thickness_to_zero",
                    "thickness_scaling", "coordinate_scaling", "invert_surface"):
            assert key in oe._support_dictionary
        assert "dabax" not in oe._support_dictionary


class TestGetThicknessMap:
    def test_requires_file(self):
        oe = S4NumericalMeshPhaseDeflector()
        with pytest.raises(ValueError):
            oe.get_thickness_map()

    def test_nonsquare_mesh_orientation_is_correct(self, tmp_path):
        # regression test: S4Mesh.get_mesh_z() is (nx, ny); get_thickness_map must
        # transpose it to the (ny, nx) convention used by apply_phase_deflection.
        oe, x_axis, y_axis, z_nx_ny = _make_phase_deflector(
            tmp_path, lambda X, Y: numpy.abs(X) + 0.0 * Y, nx=21, ny=9,
            material="Be", density=1.85, f_r_ind=0,
        )
        tm = oe.get_thickness_map()
        assert tm["profile"].shape == (9, 21)
        assert numpy.allclose(tm["x_axis"], x_axis)
        assert numpy.allclose(tm["y_axis"], y_axis)
        assert numpy.allclose(tm["profile"], z_nx_ny.T)

    def test_thickness_scaling_multiplies_profile_only(self, tmp_path):
        oe, x_axis, y_axis, z_nx_ny = _make_phase_deflector(
            tmp_path, lambda X, Y: numpy.abs(X) + numpy.abs(Y), thickness_scaling=3.0,
        )
        tm = oe.get_thickness_map()
        assert numpy.allclose(tm["profile"], 3.0 * z_nx_ny.T)
        assert numpy.allclose(tm["x_axis"], x_axis)

    def test_coordinate_scaling_multiplies_axes_only(self, tmp_path):
        oe, x_axis, y_axis, z_nx_ny = _make_phase_deflector(
            tmp_path, lambda X, Y: numpy.abs(X) + numpy.abs(Y), coordinate_scaling=2.0,
        )
        tm = oe.get_thickness_map()
        assert numpy.allclose(tm["x_axis"], 2.0 * x_axis)
        assert numpy.allclose(tm["y_axis"], 2.0 * y_axis)
        assert numpy.allclose(tm["profile"], z_nx_ny.T)

    def test_negative_profile_without_shift_raises(self, tmp_path):
        oe, *_ = _make_phase_deflector(tmp_path, lambda X, Y: X + 0.0 * Y)  # symmetric -> negative half
        with pytest.raises(ValueError):
            oe.get_thickness_map()

    def test_negative_profile_validate_false_is_allowed(self, tmp_path):
        oe, *_ = _make_phase_deflector(tmp_path, lambda X, Y: X + 0.0 * Y)
        tm = oe.get_thickness_map(validate=False)
        assert numpy.any(tm["profile"] < 0.0)

    def test_shift_thickness_to_zero(self, tmp_path):
        oe, *_ = _make_phase_deflector(
            tmp_path, lambda X, Y: X + 0.0 * Y, shift_thickness_to_zero=True,
        )
        tm = oe.get_thickness_map()
        assert numpy.isclose(numpy.nanmin(tm["profile"]), 0.0)
        assert numpy.all(tm["profile"] >= -1e-15)

    def test_invert_surface_with_shift(self, tmp_path):
        oe, x_axis, y_axis, z_nx_ny = _make_phase_deflector(
            tmp_path, lambda X, Y: numpy.abs(X) + numpy.abs(Y),
            invert_surface=True, shift_thickness_to_zero=True,
        )
        tm = oe.get_thickness_map()
        profile_no_invert = z_nx_ny.T
        expected = -profile_no_invert - numpy.nanmin(-profile_no_invert)
        assert numpy.allclose(tm["profile"], expected)


class TestToPythonCode:
    def test_none_fields_are_not_quoted(self):
        oe = S4NumericalMeshPhaseDeflector(name="pd2")
        code = oe.to_python_code()
        assert "surface_data_file=None," in code
        assert "'None'" not in code

    def test_string_fields_are_quoted(self, tmp_path):
        oe, *_ = _make_phase_deflector(tmp_path, lambda X, Y: numpy.abs(X) + numpy.abs(Y),
                                        f_r_ind=1, file_r_ind="material.dat")
        code = oe.to_python_code()
        assert "file_r_ind='material.dat'," in code
        assert "surface_data_file='%s'," % str(tmp_path / "mesh.h5") in code

    def test_generated_code_reproduces_thickness_map(self, tmp_path):
        oe, x_axis, y_axis, z_nx_ny = _make_phase_deflector(
            tmp_path, lambda X, Y: numpy.abs(X) + numpy.abs(Y), nx=17, ny=11,
            material="Be", density=1.85, f_r_ind=0,
            r_ind=0.999999, r_attenuation=10.0,
        )
        tm_before = oe.get_thickness_map()

        code = oe.to_python_code()
        namespace = {}
        exec(code, namespace)
        regenerated = namespace["optical_element"]

        tm_after = regenerated.get_thickness_map()
        assert numpy.allclose(tm_after["profile"], tm_before["profile"])
        assert numpy.allclose(tm_after["x_axis"], tm_before["x_axis"])
        assert regenerated._material == "Be"
        assert regenerated._r_attenuation == 10.0


# ---------------------------------------------------------------------------
# S4NumericalMeshPhaseDeflectorElement.trace_beam: plumbing / integration
# ---------------------------------------------------------------------------

class TestTraceBeam:
    def _element(self, tmp_path, **kwargs):
        oe, x_axis, y_axis, z_nx_ny = _make_phase_deflector(
            tmp_path, lambda X, Y: numpy.zeros_like(X),  # flat, zero thickness
            material="Be", density=1.85, f_r_ind=0,
            r_ind=1.0 - 1e-6, r_attenuation=10.0,
            half_width=2e-3, half_height=2e-3,
            **kwargs,
        )
        return oe

    def _beam(self, n=200):
        beam = S4Beam.initialize_as_pencil(n)
        rng = numpy.random.default_rng(0)
        beam.set_column(1, rng.uniform(-1e-3, 1e-3, n))
        beam.set_column(3, rng.uniform(-1e-3, 1e-3, n))
        beam.set_column(11, numpy.full(n, 2 * numpy.pi / 1e-8))
        return beam

    def test_rejects_non_zero_angle_radial(self, tmp_path):
        oe = self._element(tmp_path)
        coordinates = ElementCoordinates(angle_radial=0.01, angle_radial_out=numpy.pi)
        element = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates,
                                                        input_beam=self._beam())
        with pytest.raises(NotImplementedError):
            element.trace_beam()

    def test_rejects_angle_radial_out_not_pi(self, tmp_path):
        oe = self._element(tmp_path)
        coordinates = ElementCoordinates(angle_radial=0.0, angle_radial_out=0.0)
        element = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates,
                                                        input_beam=self._beam())
        with pytest.raises(NotImplementedError):
            element.trace_beam()

    def test_rejects_azimuthal_rotation(self, tmp_path):
        oe = self._element(tmp_path)
        coordinates = ElementCoordinates(angle_radial=0.0, angle_radial_out=numpy.pi, angle_azimuthal=0.1)
        element = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates,
                                                        input_beam=self._beam())
        with pytest.raises(NotImplementedError):
            element.trace_beam()

    def test_rejects_active_movements(self, tmp_path):
        oe = self._element(tmp_path)
        coordinates = ElementCoordinates(angle_radial=0.0, angle_radial_out=numpy.pi)
        movements = S4BeamlineElementMovements(f_move=1, offset_x=1e-6)
        element = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates,
                                                        movements=movements, input_beam=self._beam())
        with pytest.raises(NotImplementedError):
            element.trace_beam()

    def test_rectangle_boundary_crops_rays(self, tmp_path):
        oe = self._element(tmp_path, boundary_shape=Rectangle(x_left=-0.5e-3, x_right=0.5e-3,
                                                                y_bottom=-0.5e-3, y_top=0.5e-3))
        coordinates = ElementCoordinates(p=1.0, q=1.0, angle_radial=0.0, angle_radial_out=numpy.pi)
        beam = self._beam(500)
        element = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates, input_beam=beam)
        out, footprint = element.trace_beam()

        alive = out.get_column(10) > 0
        assert 0 < numpy.sum(alive) < 500  # some cropped, not all
        # `footprint` is in the element-local frame (post rotate/translate): the
        # boundary's two transverse axes are local columns 1 (X, world X) and
        # 2 (Y_local, == world Z) -- NOT column 3, which is the local depth
        # axis and is ~0 for every ray that reached the element.
        lost_x = footprint.get_column(1)[~alive]
        lost_y_local = footprint.get_column(2)[~alive]
        assert numpy.all((numpy.abs(lost_x) > 0.5e-3 - 1e-12) | (numpy.abs(lost_y_local) > 0.5e-3 - 1e-12))
        assert numpy.allclose(footprint.get_column(3), 0.0, atol=1e-9)  # local depth axis
        assert numpy.all(out.get_column(10)[~alive] < 0)
        assert numpy.allclose(out.get_column(1), footprint.get_column(1))  # X unaffected by the swap
        assert numpy.allclose(out.get_column(3)[alive], footprint.get_column(2)[alive], atol=1e-6)

    def test_zero_thickness_map_is_transparent(self, tmp_path):
        oe = self._element(tmp_path)  # flat, zero thickness everywhere
        coordinates = ElementCoordinates(p=0.0, q=0.0, angle_radial=0.0, angle_radial_out=numpy.pi)
        beam = self._beam(50)
        element = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates, input_beam=beam)
        out, footprint = element.trace_beam()

        assert numpy.all(out.get_column(10) > 0)
        assert numpy.allclose(out.get_column(23), beam.get_column(23))  # intensity unchanged
        assert numpy.allclose(out.get_column(13), 0.0)  # no OPD added (p=q=0, t=0)
        assert numpy.allclose(out.get_column(1), beam.get_column(1))  # positions unchanged (p=q=0)

    def test_collimated_beam_through_p_and_q_keeps_transverse_position(self, tmp_path):
        # Regression test: rays must be advanced from local Z=p to the element
        # plane (Z=0) before the thickness kernel runs, and back out over q.
        oe = self._element(tmp_path)  # flat, zero thickness everywhere
        coordinates = ElementCoordinates(p=3.0, q=2.0, angle_radial=0.0, angle_radial_out=numpy.pi)
        n = 30
        beam = S4Beam.initialize_as_pencil(n)
        rng = numpy.random.default_rng(1)
        beam.set_column(1, rng.uniform(-1e-4, 1e-4, n))
        beam.set_column(3, rng.uniform(-1e-4, 1e-4, n))
        beam.set_column(11, numpy.full(n, 2 * numpy.pi / 1e-8))

        element = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates, input_beam=beam)
        out, footprint = element.trace_beam()

        assert numpy.all(out.get_column(10) > 0)
        assert numpy.allclose(out.get_column(1), beam.get_column(1))  # X unchanged (no divergence)
        assert numpy.allclose(out.get_column(3), beam.get_column(3))  # Z unchanged (no divergence)
        assert numpy.allclose(out.get_column(13), 3.0 + 2.0)  # OPD == p + q exactly (vacuum, t=0)

    def test_prism_kick_through_full_pipeline_uses_only_the_q_leg(self, tmp_path):
        # Sharper regression test than the one above: a zero-thickness map
        # cannot tell apart "kick applied at the element, drifted over q"
        # (correct) from "kick applied at the source, drifted over p+q" (a
        # real bug found while building this pipeline: rays were never
        # advanced from local Z=p to Z=0 before the kick was computed/applied).
        a = 3.0  # dt/dx
        oe, x_axis, y_axis, z_nx_ny = _make_phase_deflector(
            tmp_path, lambda X, Y: 1e-2 + a * X, nx=9, ny=9,
            half_width=2e-3, half_height=2e-3,
            material="Be", f_r_ind=0,
            r_ind=1.0 - 5e-6, r_attenuation=0.0,
        )
        p, q = 3.0, 2.0
        coordinates = ElementCoordinates(p=p, q=q, angle_radial=0.0, angle_radial_out=numpy.pi)
        beam = S4Beam.initialize_as_pencil(1)  # collimated, X=Z=0, at the optical axis
        beam.set_column(11, numpy.full(1, 2 * numpy.pi / 1e-8))
        element = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates, input_beam=beam)
        out, footprint = element.trace_beam()

        delta = 5e-6
        expected_kick_x = -delta * a
        expected_x_output = 0.0 + expected_kick_x * q  # drift only over the q leg, not p+q
        assert numpy.isclose(out.get_column(1)[0], expected_x_output, atol=1e-12)
        wrong_x_output = expected_kick_x * (p + q)
        assert not numpy.isclose(out.get_column(1)[0], wrong_x_output, atol=1e-12)

    def test_uniform_slab_opd_through_full_pipeline_matches_manual_formula(self, tmp_path):
        oe, x_axis, y_axis, z_nx_ny = _make_phase_deflector(
            tmp_path, lambda X, Y: numpy.full_like(X, 4e-6), nx=5, ny=5,
            half_width=2e-3, half_height=2e-3,
            material="Be", f_r_ind=0,
            r_ind=1.0 - 6e-6, r_attenuation=0.0,
        )
        coordinates = ElementCoordinates(p=3.0, q=2.0, angle_radial=0.0, angle_radial_out=numpy.pi)
        beam = S4Beam.initialize_as_pencil(4)
        beam.set_column(11, numpy.full(4, 2 * numpy.pi / 1e-8))
        element = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates, input_beam=beam)
        out, footprint = element.trace_beam()

        t = 4e-6
        delta = 6e-6
        expected_opd = 3.0 + 2.0 - delta * t  # p + q from free space, -delta*t from the slab
        assert numpy.allclose(out.get_column(13), expected_opd)

    def test_apply_to_lost_defaults_from_optical_element_setting(self, tmp_path):
        # Reproduces old S4ThinTransmission behavior: apply_to_lost is a
        # persisted optical-element setting, used as the trace_beam default,
        # but still overridable per call.
        oe = self._element(tmp_path, apply_to_lost=True)
        coordinates = ElementCoordinates(p=0.0, q=0.0, angle_radial=0.0, angle_radial_out=numpy.pi)
        beam = self._beam(20)
        beam.rays[0, 9] = -5.0  # mark one ray as already lost before the element
        element = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates, input_beam=beam)

        out_default, _ = element.trace_beam()  # optical_element default apply_to_lost=True
        assert out_default.rays[0, 9] == -5.0  # flag preserved, but amplitude touched
        assert numpy.isclose(out_default.rays[0, 6], numpy.exp(-1.0 * 0.0))  # zero thickness -> amplitude 1

        element2 = S4NumericalMeshPhaseDeflectorElement(optical_element=oe, coordinates=coordinates,
                                                          input_beam=S4Beam(array=beam.rays.copy()))
        out_false, _ = element2.trace_beam(apply_to_lost=False)
        assert out_false.rays[0, 9] == -5.0
        assert out_false.rays[0, 6] == beam.rays[0, 6]  # completely untouched


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
