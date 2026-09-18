# SPDX-License-Identifier: CECILL-2.1
# Copyright (c) 2026 ESRF - the European Synchrotron

"""Shadow4 thin phase/transmission optical element and numerical helpers."""

from __future__ import annotations

from numpy.typing import ArrayLike
from syned.beamline.element_coordinates import ElementCoordinates
from syned.beamline.shape import Ellipse, Rectangle

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.s4_beamline_element import S4BeamlineElement
from shadow4.beamline.s4_beamline_element_movements import S4BeamlineElementMovements
from shadow4.beamline.s4_optical_element_decorators import S4NumericalMeshOpticalElementDecorator
from shadow4.beamline.optical_elements.refractors.s4_interface import S4Interface
import numpy

# S4Interface.f_r_ind supports 10 modes because it allows a different material
# on each side of a (possibly curved) interface. A thin transmission element
# only ever has one material (the slab), surrounded by vacuum on both sides,
# so only the "object side" modes are meaningful here; the image side is
# always pinned to vacuum (r_ind_ima=1.0, r_attenuation_ima=0.0).
#   ri_calculation_mode -> f_r_ind: 0=user -> 0, 1=prerefl -> 1, 2=xraylib -> 4, 3=dabax -> 7
_RI_CALCULATION_MODE_TO_F_R_IND = {0: 0, 1: 1, 2: 4, 3: 7}


class S4ThinTransmission(S4Interface, S4NumericalMeshOpticalElementDecorator):
    """
    Shadow4 thin refractive/transmission element.

    Follows the same base-class pattern as the other Shadow4 refractors
    (e.g. :class:`~shadow4.beamline.optical_elements.refractors.s4_numerical_mesh_interface.S4NumericalMeshInterface`):
    it inherits ``S4Interface`` for material/optical-constant storage and
    ``S4NumericalMeshOpticalElementDecorator`` for the mesh. ``thickness_mesh_file``
    is the only supported mesh source: it must be an HDF5 (or 3-column ASCII)
    file with the standard OASYS surface convention (a 2D projected thickness
    map in metres). Passing raw arrays is intentionally not supported, since
    they cannot be recovered by ``to_python_code()``.

    Unlike a real refractive interface, this element does not compute a
    Snell's-law ray/surface intersection: it stays flat at normal incidence
    and only uses the mesh as a projected thickness map to apply attenuation,
    OPD/phase and phase-gradient angular kicks. Accordingly only the
    "object side" of ``S4Interface`` is used (the slab material); the image
    side is always pinned to vacuum since a thin slab has no second material.
    """

    def __init__(
        self,
        name="Undefined",
        boundary_shape=None,
        thickness_mesh_file=None,
        material="",
        density=1.0,
        ri_calculation_mode=0,
        prerefl_file=None,
        refraction_index=1.0,
        attenuation_coefficient=0.0,
        dabax=None,
        apply_to_lost=True,
        shift_thickness_to_zero=False,
        thickness_scaling=1.0,
        coordinate_scaling=1.0,
        invert_surface=False,
    ):
        if ri_calculation_mode not in _RI_CALCULATION_MODE_TO_F_R_IND:
            raise ValueError("ri_calculation_mode must be 0, 1, 2 or 3.")

        S4NumericalMeshOpticalElementDecorator.__init__(self, xx=None, yy=None, zz=None, surface_data_file=thickness_mesh_file)
        S4Interface.__init__(
            self,
            name=name,
            boundary_shape=boundary_shape,
            surface_shape=self.get_surface_shape_instance(),
            material_object=material,
            material_image=None,
            density_object=density,
            density_image=1.0,
            f_r_ind=_RI_CALCULATION_MODE_TO_F_R_IND[ri_calculation_mode],
            r_ind_obj=refraction_index,
            r_ind_ima=1.0,
            r_attenuation_obj=attenuation_coefficient,
            r_attenuation_ima=0.0,
            file_r_ind_obj=prerefl_file if prerefl_file is not None else "",
            file_r_ind_ima="",
            dabax=dabax,
        )

        if not numpy.isfinite(thickness_scaling):
            raise ValueError("thickness_scaling must be finite.")
        if not numpy.isfinite(coordinate_scaling) or coordinate_scaling <= 0.0:
            raise ValueError("coordinate_scaling must be finite and strictly positive.")

        self._apply_to_lost = apply_to_lost
        self._shift_thickness_to_zero = shift_thickness_to_zero
        self._thickness_scaling = float(thickness_scaling)
        self._coordinate_scaling = float(coordinate_scaling)
        self._invert_surface = bool(invert_surface)

        self.__inputs = {
            "name": name,
            "boundary_shape": boundary_shape,
            "thickness_mesh_file": thickness_mesh_file,
            "material": material,
            "density": density,
            "ri_calculation_mode": ri_calculation_mode,
            "prerefl_file": prerefl_file,
            "refraction_index": refraction_index,
            "attenuation_coefficient": attenuation_coefficient,
            "dabax": self._get_dabax_txt(),
            "apply_to_lost": apply_to_lost,
            "shift_thickness_to_zero": shift_thickness_to_zero,
            "thickness_scaling": thickness_scaling,
            "coordinate_scaling": coordinate_scaling,
            "invert_surface": invert_surface,
        }

    def get_info(self):
        txt = "\n\n"
        txt += "THIN TRANSMISSION ELEMENT\n"
        txt += "  Thin projected-thickness map with refraction and attenuation\n"
        txt += "  Normal-incidence element: angle_radial=0, angle_radial_out=pi\n"
        txt += "  Material: %s\n" % self.get_material_object()
        txt += "  Density: %g g/cm^3\n" % self._density_object
        txt += "  ri_calculation_mode (f_r_ind): %d\n" % self._f_r_ind
        txt += "  apply_to_lost: %s\n" % self._apply_to_lost
        txt += "  shift_thickness_to_zero: %s\n" % self._shift_thickness_to_zero
        txt += "  Thickness scaling: %g\n" % self._thickness_scaling
        txt += "  Coordinate scaling: %g\n" % self._coordinate_scaling
        txt += "  Invert surface: %s\n" % self._invert_surface

        thickness_mesh = self.get_surface_shape_instance()
        if thickness_mesh is not None and thickness_mesh.has_surface_data_file():
            txt += "  Thickness data file: %s\n" % thickness_mesh._surface_data_file
        else:
            txt += "  Thickness data: undefined\n"

        boundary = self.get_boundary_shape()
        if not isinstance(boundary, (Rectangle, Ellipse)):
            txt += "  Boundaries not considered (infinite)\n"
        else:
            txt += "  Boundary shape: %s\n" % boundary.__class__.__name__
            txt += "    Limits: " + repr(boundary.get_boundaries()) + "\n"
        return txt

    def get_thickness_map(self, validate=True):
        """
        Return the internal thickness dictionary used by the numerical kernel.

        Set ``validate=False`` to inspect the transformed map without rejecting
        negative values. This is intended for previews; tracing always uses the
        validated default.
        """
        thickness_mesh = self.get_surface_shape_instance()
        if thickness_mesh is None or not thickness_mesh.has_surface_data_file():
            raise ValueError("No thickness_mesh_file supplied.")

        mesh = self.get_optical_surface_instance()
        x_axis, y_axis = mesh.get_mesh_x_y()
        # S4Mesh stores its z grid as (len(x_axis), len(y_axis)); the thin-element
        # kernel below (and numpy.gradient(profile, y_axis, x_axis)) expects the
        # transposed (len(y_axis), len(x_axis)) convention.
        profile = mesh.get_mesh_z().T

        x_axis = numpy.asarray(x_axis, dtype=float) * self._coordinate_scaling
        y_axis = numpy.asarray(y_axis, dtype=float) * self._coordinate_scaling
        profile = numpy.asarray(profile, dtype=float) * self._thickness_scaling
        if self._invert_surface:
            profile = -profile

        if self._shift_thickness_to_zero:
            finite = numpy.isfinite(profile)
            if numpy.any(finite):
                profile = profile - numpy.nanmin(profile)
        elif validate and numpy.any(profile[numpy.isfinite(profile)] < 0.0):
            raise ValueError(
                "Thickness profile contains negative values. Use "
                "shift_thickness_to_zero=True to subtract the finite minimum."
            )

        return {
            "profile": profile,
            "x_axis": x_axis,
            "y_axis": y_axis,
        }

    def get_optical_constants(self, photon_energy_eV, k):
        """
        Return per-ray ``delta`` and ``beta`` for the material.

        Delegates to the inherited ``S4Interface.get_refraction_indices()``/
        ``get_attenuation_coefficients()`` (constant/prerefl/xraylib/dabax,
        all already handled there) and keeps only the object-side result,
        since the image side is always pinned to vacuum for this element.
        """
        n_real, _ = self.get_refraction_indices(photon_energy_eV)
        mu, _ = self.get_attenuation_coefficients(photon_energy_eV)  # already in m^-1

        n_real = numpy.asarray(n_real, dtype=float)
        mu = numpy.asarray(mu, dtype=float)
        delta = 1.0 - n_real
        beta = mu / (2.0 * numpy.asarray(k, dtype=float))
        return delta, beta

    def to_python_code(self, **kwargs):
        txt = self.to_python_code_boundary_shape()
        txt += "\nfrom shadow4.beamline.optical_elements.refractors.s4_thin_element import S4ThinTransmission"
        txt += "\noptical_element = S4ThinTransmission(name='%s', boundary_shape=boundary_shape," % self.__inputs["name"]
        txt += "\n    thickness_mesh_file=%r," % self.__inputs["thickness_mesh_file"]
        txt += "\n    material='%s', density=%g," % (self.__inputs["material"], self.__inputs["density"])
        txt += "\n    ri_calculation_mode=%d, # 0=user, 1=prerefl, 2=xraylib, 3=dabax" % self.__inputs["ri_calculation_mode"]
        txt += "\n    prerefl_file=%r," % self.__inputs["prerefl_file"]
        txt += "\n    refraction_index=%.10g," % self.__inputs["refraction_index"]
        txt += "\n    attenuation_coefficient=%g," % self.__inputs["attenuation_coefficient"]
        txt += "\n    dabax=%s," % self.__inputs["dabax"]
        txt += "\n    apply_to_lost=%s," % repr(self.__inputs["apply_to_lost"])
        txt += "\n    shift_thickness_to_zero=%s," % repr(self.__inputs["shift_thickness_to_zero"])
        txt += "\n    thickness_scaling=%g," % self.__inputs["thickness_scaling"]
        txt += "\n    coordinate_scaling=%g," % self.__inputs["coordinate_scaling"]
        txt += "\n    invert_surface=%s," % repr(self.__inputs["invert_surface"])
        txt += "\n    )"
        return txt


class S4ThinTransmissionElement(S4BeamlineElement):
    """
    Shadow4 beamline element wrapper for :class:`S4ThinTransmission`.
    """

    def __init__(
        self,
        optical_element: S4ThinTransmission = None,
        coordinates: ElementCoordinates = None,
        movements: S4BeamlineElementMovements = None,
        input_beam: S4Beam = None,
    ):
        super().__init__(
            optical_element=optical_element if optical_element is not None else S4ThinTransmission(),
            coordinates=coordinates if coordinates is not None else ElementCoordinates(),
            movements=movements,
            input_beam=input_beam,
        )
        self.__stored_optical_constants = None

    def get_stored_optical_constants(self):
        return self.__stored_optical_constants

    def _apply_thin_transmission(self, beam, delta, beta, *, flag_lost_value=-1.0, apply_to_lost=True):
        """
        Apply the thin-element projected-thickness kick to a beam already
        expressed in the element's local reference frame.

        This plays the same role as ``S4Interface._apply_interface_refraction``
        (called the same way, from :meth:`trace_beam`, right after the beam
        has been put into the element reference system), but instead of a
        Snell's-law ray/surface intersection with a (possibly curved) mesh,
        it advances each ray to the flat element plane (Z=0) and then applies
        the local thickness map's attenuation, OPD/phase and phase-gradient
        kicks (see :func:`apply_thin_transmission_element`). Kept here (on
        the beamline element, not on ``S4ThinTransmission``) since everything
        it needs -- the optical element and the beam already in the local
        frame -- is already local to :meth:`trace_beam`.
        """
        soe = self.get_optical_element()
        k = beam.get_column(11) * 100.0
        _advance_rays_to_local_surface(beam)
        return apply_thin_transmission_element(
            beam,
            soe.get_thickness_map(),
            k,
            delta,
            beta,
            flag_lost_value=flag_lost_value,
            apply_to_lost=apply_to_lost,
            duplicate=False,
        )

    def trace_beam(self, **params):
        """
        Runs (ray tracing) the input beam through the element.

        Follows the same structure as :meth:`S4InterfaceElement.trace_beam`
        (put the beam in the element reference system, apply the element's
        physics, apply boundaries, transform back to the image plane); the
        Snell's-law refraction step (``_apply_interface_refraction``) is
        replaced by :meth:`_apply_thin_transmission`.

        Parameters
        ----------
        **params

        Returns
        -------
        tuple
            (output_beam, footprint) instances of S4Beam.
        """
        flag_lost_value = params.get("flag_lost_value", -1.0)
        apply_to_lost = params.get("apply_to_lost", None)
        reused_stored_optical_constants = params.get("reused_stored_optical_constants", None)
        angle_tolerance = params.get("angle_tolerance", 1e-9)

        p = self.get_coordinates().p()
        q = self.get_coordinates().q()
        angle_radial = self.get_coordinates().angle_radial()
        angle_radial_out = self.get_coordinates().angle_radial_out()
        alpha1 = self.get_coordinates().angle_azimuthal()

        if not numpy.isclose(angle_radial, 0.0, rtol=0.0, atol=angle_tolerance):
            raise NotImplementedError("S4ThinTransmission supports only angle_radial=0.")
        if not numpy.isclose(angle_radial_out, numpy.pi, rtol=0.0, atol=angle_tolerance):
            raise NotImplementedError("S4ThinTransmission supports only angle_radial_out=pi.")
        if not numpy.isclose(alpha1, 0.0, rtol=0.0, atol=angle_tolerance):
            raise NotImplementedError("S4ThinTransmission angle_azimuthal rotations are not implemented yet.")

        theta_grazing1 = numpy.pi / 2 - angle_radial
        theta_grazing2 = numpy.pi / 2 - angle_radial_out

        #
        input_beam = self.get_input_beam().duplicate()
        soe = self.get_optical_element()

        # retrieve and store optical constants
        if reused_stored_optical_constants is not None:
            delta, beta = reused_stored_optical_constants
        else:
            energy1 = input_beam.get_photon_energy_eV()
            k1 = input_beam.get_column(11) * 100.0
            delta, beta = soe.get_optical_constants(energy1, k1)
            self.__stored_optical_constants = (delta, beta)

        #
        # put beam in element reference system
        #
        input_beam.rotate(alpha1, axis=2)
        input_beam.rotate(theta_grazing1, axis=1)
        input_beam.translation([0.0, -p * numpy.cos(theta_grazing1), p * numpy.sin(theta_grazing1)])

        # element movement:
        movements = self.get_movements()
        if movements is not None:
            if movements.f_move:
                raise NotImplementedError("S4ThinTransmission movements are not implemented yet.")

        #
        # apply the thin-element projected-thickness kick
        #
        footprint = self._apply_thin_transmission(
            input_beam, delta, beta,
            flag_lost_value=flag_lost_value,
            apply_to_lost=soe._apply_to_lost if apply_to_lost is None else apply_to_lost,
        )

        #
        # apply element boundaries
        #
        footprint.apply_boundaries_syned(soe.get_boundary_shape(), flag_lost_value=flag_lost_value)

        #
        # from element reference system to image plane
        #
        output_beam = footprint.duplicate()
        output_beam.change_to_image_reference_system(theta_grazing2, q, refraction_index=1.0, apply_attenuation=0)

        return output_beam, footprint

    def to_python_code(self, **kwargs):
        txt = "\n\n# optical element number XX"
        txt += self.get_optical_element().to_python_code()
        txt += self.to_python_code_coordinates()
        txt += self.to_python_code_movements()
        txt += "\nfrom shadow4.beamline.optical_elements.refractors.s4_thin_element import S4ThinTransmissionElement"
        txt += "\nbeamline_element = S4ThinTransmissionElement(optical_element=optical_element, coordinates=coordinates, movements=movements, input_beam=beam)"
        txt += "\n\nbeam, footprint = beamline_element.trace_beam()"
        return txt


def _advance_rays_to_local_surface(beam: S4Beam) -> None:
    """
    Propagate every ray, in place and using its own direction, from its
    current position in the element's local reference frame to the flat
    element plane at Z=0 (column 3), accumulating the traversed distance
    into the optical path length (column 13).

    This is the flat-plane analogue of the ray/surface intersection that
    ``S4Mesh.apply_refraction_on_beam`` performs for curved interfaces inside
    ``S4Interface._apply_interface_refraction``: after
    ``S4InterfaceElement``-style rotate/translate, every ray starts at local
    Z = p (the source-to-element distance) and travels in -Z; this brings it
    to Z = 0 so :func:`apply_thin_transmission_element` can interpolate the
    thickness map at the ray's actual (possibly divergence-shifted) position.
    """
    rays = beam.rays
    tof = -rays[:, 2] / rays[:, 5]
    rays[:, 0] += tof * rays[:, 3]
    rays[:, 1] += tof * rays[:, 4]
    rays[:, 2] += tof * rays[:, 5]
    rays[:, 12] += tof


def apply_thin_transmission_element(
    beam: S4Beam,
    thickness: dict,
    k: ArrayLike,
    delta: ArrayLike,
    beta: ArrayLike,
    *,
    flag_lost_value: float = -1.0,
    apply_to_lost: bool = True,
    duplicate: bool = True,
) -> S4Beam:
    """
    Apply a thin refractive/transmitting phase element to a Shadow4 beam
    already expressed in the element's local reference frame, i.e. after the
    same "put beam in element reference system" transform that
    ``S4InterfaceElement.trace_beam`` applies before refraction: the
    transverse plane is (X, Y) [columns 1, 2] and the beam travels along
    -Z [column 3], reaching the (flat) element at Z=0.

    Parameters
    ----------
    beam : shadow4.beam.s4_beam.S4Beam
        Beam in the element's local reference frame. Columns 1 and 2 (X, Y)
        are not changed; column 3 (Z) is expected to be (numerically) zero
        for rays landing on the element.
    thickness : dict
        Dictionary with keys ``"profile"``, ``"x_axis"`` and ``"y_axis"``.
        ``profile`` is a 2D array ``profile[y, x]`` containing the projected
        thickness in meters, sampled on the element's local (X, Y) plane
        (``x_axis`` samples local column 1, ``y_axis`` samples local column 2).
    k : scalar or array_like
        Vacuum wavenumber, ``2*pi/wavelength``, in m^-1. If array-like, it
        must have one value per ray.
    delta, beta : scalar or array_like
        Refractive index decrement and absorption index per ray, with
        ``n = 1 - delta + 1j * beta``. If array-like, they must have one value
        per ray.
    flag_lost_value : float, optional
        Negative Shadow4 flag assigned to rays newly lost by this element.
    apply_to_lost : bool, optional
        If True, apply kicks, OPD, phase and attenuation to all stored rays
        with valid thickness support, including rays already marked lost
        with a negative Shadow4 flag. Their existing negative flags are
        preserved. If False, already-lost rays are left unchanged.
    duplicate : bool, optional
        If True, operate on a duplicate of ``beam``. If False, modify ``beam``
        in place.

    Returns
    -------
    shadow4.beam.s4_beam.S4Beam
        Beam with columns 4, 5, 6, 7, 8, 9, 10, 13, 14, 15, 16, 17 and 18
        updated as needed. Columns 1, 2 and 3 are left untouched.

    Notes
    -----
    The complex field factor is

    ``exp(-k * beta * t) * exp(-1j * k * delta * t)``.

    Therefore the relative OPD added to column 13 is ``-delta * t`` and the
    phase added to columns 14 and 15 is ``-k * delta * t``.

    The angular kick is the thin-phase gradient, applied to the local
    transverse direction components:

    ``dVx = -delta * dt/dx`` and ``dVy = -delta * dt/dy``.

    The local "depth" direction component (Vz, column 6) is recomputed to
    preserve a unit direction vector, keeping the sign of the incoming ray
    (rays reaching the element travel in -Z).

    Shadow4 lost-ray flags are used throughout: column 10 > 0 means alive,
    column 10 < 0 means lost. This function does not use the barc4beams
    ``0 = alive, 1 = lost`` convention. Rays that are already lost are still
    propagated through the thin-element calculation when they land on valid
    thickness support if ``apply_to_lost`` is True; their negative flag is
    preserved.
    """
    if flag_lost_value >= 0.0:
        raise ValueError("flag_lost_value must be negative for Shadow4 beams.")

    profile, x_axis, y_axis = _validate_thickness_map(thickness)

    out = beam.duplicate() if duplicate else beam
    rays = out.rays
    n_rays = rays.shape[0]

    k_ray = _as_ray_array("k", k, n_rays)
    delta_ray = _as_ray_array("delta", delta, n_rays)
    beta_ray = _as_ray_array("beta", beta, n_rays)

    valid_support = numpy.isfinite(profile) & (profile >= 0.0)
    if not numpy.any(valid_support):
        raise ValueError("No valid thickness support remains after masking.")

    thickness_map = numpy.where(valid_support, profile, numpy.nan)
    dt_dy, dt_dx = numpy.gradient(thickness_map, y_axis, x_axis)

    x_ray = rays[:, 0]
    y_ray = rays[:, 1]

    thickness_ray = _interp2d_regular(x_axis, y_axis, thickness_map, x_ray, y_ray)
    dt_dx_ray = _interp2d_regular(x_axis, y_axis, dt_dx, x_ray, y_ray)
    dt_dy_ray = _interp2d_regular(x_axis, y_axis, dt_dy, x_ray, y_ray)

    inside = (
        (x_ray >= min(x_axis[0], x_axis[-1]))
        & (x_ray <= max(x_axis[0], x_axis[-1]))
        & (y_ray >= min(y_axis[0], y_axis[-1]))
        & (y_ray <= max(y_axis[0], y_axis[-1]))
    )

    phase = -k_ray * delta_ray * thickness_ray
    opd = -delta_ray * thickness_ray
    amplitude = numpy.exp(-k_ray * beta_ray * thickness_ray)

    kick_x = -delta_ray * dt_dx_ray
    kick_y = -delta_ray * dt_dy_ray

    vx = rays[:, 3].copy()
    vy = rays[:, 4].copy()
    vz = rays[:, 5].copy()

    vx_new = vx + kick_x
    vy_new = vy + kick_y
    vz2_new = 1.0 - vx_new * vx_new - vy_new * vy_new

    valid_through_element = (
        inside
        & numpy.isfinite(thickness_ray)
        & numpy.isfinite(dt_dx_ray)
        & numpy.isfinite(dt_dy_ray)
        & numpy.isfinite(k_ray)
        & (k_ray > 0.0)
        & numpy.isfinite(delta_ray)
        & numpy.isfinite(beta_ray)
        & (beta_ray >= 0.0)
        & numpy.isfinite(phase)
        & numpy.isfinite(opd)
        & numpy.isfinite(amplitude)
        & (amplitude > 0.0)
        & numpy.isfinite(vx_new)
        & numpy.isfinite(vy_new)
        & (vz2_new >= 0.0)
    )

    alive_in = rays[:, 9] > 0.0
    rays_to_apply = numpy.ones(n_rays, dtype=bool) if apply_to_lost else alive_in
    valid = rays_to_apply & valid_through_element
    newly_lost = alive_in & ~valid_through_element
    invalid = rays_to_apply & ~valid_through_element

    sign_vz = numpy.where(vz >= 0.0, 1.0, -1.0)
    vz_new = sign_vz * numpy.sqrt(numpy.maximum(vz2_new, 0.0))

    rays[valid, 3] = vx_new[valid]
    rays[valid, 4] = vy_new[valid]
    rays[valid, 5] = vz_new[valid]

    rays[valid, 12] += opd[valid]
    rays[valid, 13] += phase[valid]
    rays[valid, 14] += phase[valid]

    for column_index in (6, 7, 8, 15, 16, 17):
        rays[valid, column_index] *= amplitude[valid]
        rays[invalid, column_index] = 0.0

    rays[newly_lost, 9] = flag_lost_value

    return out


def _validate_thickness_map(thickness: dict):
    """
    Validate and unpack a regular thin-element thickness map.

    Parameters
    ----------
    thickness : dict
        Dictionary with keys ``"profile"``, ``"x_axis"`` and ``"y_axis"``.
        ``profile`` must be a 2D array with shape ``(ny, nx)`` in metres.
        The axes must be finite, one-dimensional and strictly monotonic.

    Returns
    -------
    profile, x_axis, y_axis
        NumPy arrays containing the validated thickness profile and axes.

    Raises
    ------
    KeyError
        If one of the required keys is missing.
    ValueError
        If the profile shape, axis dimensions, monotonicity or finite-value
        requirements are not satisfied.
    """
    for key in ("profile", "x_axis", "y_axis"):
        if key not in thickness:
            raise KeyError(f"thickness missing key: {key!r}")

    profile = numpy.asarray(thickness["profile"], dtype=float)
    x_axis = numpy.asarray(thickness["x_axis"], dtype=float)
    y_axis = numpy.asarray(thickness["y_axis"], dtype=float)

    if profile.ndim != 2:
        raise ValueError("thickness['profile'] must be 2D with shape (ny, nx).")
    if x_axis.ndim != 1 or y_axis.ndim != 1:
        raise ValueError("thickness['x_axis'] and thickness['y_axis'] must be 1D arrays.")
    if x_axis.size != profile.shape[1] or y_axis.size != profile.shape[0]:
        raise ValueError("Axis lengths must match thickness profile shape (ny, nx).")
    if x_axis.size < 2 or y_axis.size < 2:
        raise ValueError("thickness axes must contain at least two points.")
    if numpy.any(~numpy.isfinite(x_axis)) or numpy.any(~numpy.isfinite(y_axis)):
        raise ValueError("thickness axes must contain finite values.")
    if not (numpy.all(numpy.diff(x_axis) > 0.0) or numpy.all(numpy.diff(x_axis) < 0.0)):
        raise ValueError("thickness['x_axis'] must be strictly monotonic.")
    if not (numpy.all(numpy.diff(y_axis) > 0.0) or numpy.all(numpy.diff(y_axis) < 0.0)):
        raise ValueError("thickness['y_axis'] must be strictly monotonic.")
    if numpy.any(profile[numpy.isfinite(profile)] < 0.0):
        raise ValueError("thickness['profile'] must be non-negative where finite.")

    return profile, x_axis, y_axis


def _as_ray_array(name: str, values, n_rays: int):
    """
    Convert a scalar or per-ray value sequence to a ray-aligned array.

    Parameters
    ----------
    name : str
        Name used in validation error messages.
    values : scalar or array_like
        Scalar value broadcast to all rays, or one-dimensional array with one
        value per ray.
    n_rays : int
        Number of rays expected in the output array.

    Returns
    -------
    numpy.ndarray
        One-dimensional float array with length ``n_rays``.

    Raises
    ------
    ValueError
        If ``values`` is not scalar or one-dimensional, or if its length does
        not match ``n_rays``.
    """
    arr = numpy.asarray(values, dtype=float)
    if arr.ndim == 0:
        return numpy.full(n_rays, float(arr), dtype=float)
    if arr.ndim != 1:
        raise ValueError(f"{name!r} must be scalar or a 1D array.")
    if arr.size != n_rays:
        raise ValueError(f"{name!r} must have length equal to the number of rays ({n_rays}).")
    return arr


def _interp2d_regular(x_axis, y_axis, values, xs, ys):
    """
    Bilinearly interpolate values on a regular rectilinear grid.

    Parameters
    ----------
    x_axis, y_axis : array_like
        One-dimensional grid axes. Either increasing or decreasing order is
        accepted.
    values : array_like
        2D array with shape ``(len(y_axis), len(x_axis))``.
    xs, ys : array_like
        Query coordinates in the same units as ``x_axis`` and ``y_axis``.

    Returns
    -------
    numpy.ndarray
        Interpolated values at ``(xs, ys)``. Points outside the grid, or cells
        touching non-finite corner values, are returned as ``nan``.

    Raises
    ------
    ValueError
        If ``values`` does not have shape ``(len(y_axis), len(x_axis))``.
    """
    x_axis = numpy.asarray(x_axis, dtype=float)
    y_axis = numpy.asarray(y_axis, dtype=float)
    values = numpy.asarray(values, dtype=float)
    xs = numpy.asarray(xs, dtype=float)
    ys = numpy.asarray(ys, dtype=float)

    if values.shape != (y_axis.size, x_axis.size):
        raise ValueError("values shape must be (len(y_axis), len(x_axis)).")

    if x_axis[0] > x_axis[-1]:
        x_axis = x_axis[::-1]
        values = values[:, ::-1]
    if y_axis[0] > y_axis[-1]:
        y_axis = y_axis[::-1]
        values = values[::-1, :]

    ix1 = numpy.searchsorted(x_axis, xs, side="right")
    iy1 = numpy.searchsorted(y_axis, ys, side="right")

    ix1 = numpy.clip(ix1, 1, x_axis.size - 1)
    iy1 = numpy.clip(iy1, 1, y_axis.size - 1)

    ix0 = ix1 - 1
    iy0 = iy1 - 1

    x0 = x_axis[ix0]
    x1 = x_axis[ix1]
    y0 = y_axis[iy0]
    y1 = y_axis[iy1]

    tx = (xs - x0) / (x1 - x0)
    ty = (ys - y0) / (y1 - y0)

    v00 = values[iy0, ix0]
    v10 = values[iy0, ix1]
    v01 = values[iy1, ix0]
    v11 = values[iy1, ix1]

    out = (
        (1.0 - tx) * (1.0 - ty) * v00
        + tx * (1.0 - ty) * v10
        + (1.0 - tx) * ty * v01
        + tx * ty * v11
    )

    bad = (
        (xs < x_axis[0])
        | (xs > x_axis[-1])
        | (ys < y_axis[0])
        | (ys > y_axis[-1])
        | ~numpy.isfinite(v00)
        | ~numpy.isfinite(v10)
        | ~numpy.isfinite(v01)
        | ~numpy.isfinite(v11)
    )
    out[bad] = numpy.nan
    return out


if __name__ == "__main__":
    import numpy as np
    from dabax.dabax_xraylib import DabaxXraylib
    from shadow4.beamline.s4_beamline import S4Beamline

    beamline = S4Beamline()

    #
    #
    #
    from shadow4.sources.source_geometrical.source_geometrical import SourceGeometrical

    light_source = SourceGeometrical(name='Geometrical Source', nrays=5000, seed=5676561)
    light_source.set_spatial_type_point()
    light_source.set_depth_distribution_off()
    light_source.set_angular_distribution_gaussian(sigdix=1e-05, sigdiz=1e-05)
    light_source.set_energy_distribution_uniform(value_min=9950, value_max=10050, unit='eV')
    light_source.set_polarization(polarization_degree=1, phase_diff=0, coherent_beam=0)
    beam = light_source.get_beam()

    beamline.set_light_source(light_source)

    # optical element number XX
    boundary_shape = None
    from shadow4.beamline.optical_elements.refractors.s4_thin_element import S4ThinTransmission

    optical_element = S4ThinTransmission(name='Thin Element SHADOW', boundary_shape=boundary_shape,
                                         thickness_mesh_file='/home/srio/Oasys2/lens_interface_1.h5',
                                         material='Be', density=1.85,
                                         ri_calculation_mode=3,  # 0=user, 1=prerefl, 2=xraylib, 3=dabax
                                         prerefl_file='<none>',
                                         refraction_index=1,
                                         attenuation_coefficient=0,
                                         dabax=DabaxXraylib(file_f1f2="f1f2_Windt.dat",
                                                            file_CrossSec="CrossSec_EPDL97.dat"),
                                         apply_to_lost=True,
                                         shift_thickness_to_zero=False,
                                         thickness_scaling=1,
                                         coordinate_scaling=1,
                                         invert_surface=False,
                                         )
    from syned.beamline.element_coordinates import ElementCoordinates

    coordinates = ElementCoordinates(p=10, q=0, angle_radial=0, angle_azimuthal=0, angle_radial_out=3.141592654)
    movements = None
    from shadow4.beamline.optical_elements.refractors.s4_thin_element import S4ThinTransmissionElement

    beamline_element = S4ThinTransmissionElement(optical_element=optical_element, coordinates=coordinates,
                                                 movements=movements, input_beam=beam)

    beam, footprint = beamline_element.trace_beam()

    beamline.append_beamline_element(beamline_element)

    # test plot
    if 1:
        from srxraylib.plot.gol import plot_scatter

        # plot_scatter(beam.get_photon_energy_eV(nolost=1), beam.get_column(23, nolost=1),
        #              title='(Intensity,Photon Energy)', plot_histograms=0)
        plot_scatter(1e6 * beam.get_column(1, nolost=1), 1e6 * beam.get_column(3, nolost=1), title='(X,Z) in microns')

    print(optical_element.info())
    print(optical_element.get_info())

    print(beamline.to_json())