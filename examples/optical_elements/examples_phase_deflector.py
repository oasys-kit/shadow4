"""Trace a pencil ray through a thin wedge with a known phase gradient."""

from pathlib import Path
from tempfile import TemporaryDirectory

import h5py
import numpy as np
from syned.beamline.element_coordinates import ElementCoordinates

from shadow4.beam.s4_beam import S4Beam
from shadow4.beamline.optical_elements.phase_deflectors.s4_numerical_mesh_phase_deflector import (
    S4NumericalMeshPhaseDeflector,
    S4NumericalMeshPhaseDeflectorElement,
)


def main():
    # At the optical axis, the wedge is 100 µm thick and dt/dx = 0.05 m/m.
    x = np.linspace(-1e-3, 1e-3, 11)
    y = np.linspace(-0.5e-3, 0.5e-3, 7)
    thickness_at_axis = 100e-6
    slope_x = 0.05
    profile = np.broadcast_to(thickness_at_axis + slope_x * x, (y.size, x.size))

    delta = 5e-6
    attenuation = 500.0  # Intensity attenuation coefficient in m^-1.
    p, q = 3.0, 2.0  # Source-to-element and element-to-image distances in m.
    beam = S4Beam.initialize_as_pencil(1)
    beam.set_photon_energy_eV(10_000.0)

    with TemporaryDirectory() as directory:
        mesh = Path(directory) / "wedge.h5"
        with h5py.File(mesh, "w") as handle:
            handle["surface_file/X"] = x
            handle["surface_file/Y"] = y
            handle["surface_file/Z"] = profile   # OASYS HDF5 order is (y, x).

        optic = S4NumericalMeshPhaseDeflector(
            name="Thin wedge", surface_data_file=str(mesh), f_r_ind=0,
            r_ind=1.0 - delta, r_attenuation=attenuation,
        )
        coordinates = ElementCoordinates(p=p, q=q, angle_radial=0.0, angle_radial_out=np.pi)
        element = S4NumericalMeshPhaseDeflectorElement(
            optical_element=optic, coordinates=coordinates, input_beam=beam,
        )
        output, footprint = element.trace_beam()

    kick = -delta * slope_x  # The ray bends toward the thinner side.
    image_x = output.get_column(1)[0]
    transmission = output.get_column(23)[0] / beam.get_column(23)[0]
    assert output.get_column(10)[0] > 0
    np.testing.assert_allclose(footprint.get_column(4)[0], kick, rtol=1e-8)
    np.testing.assert_allclose(image_x, q * kick / np.sqrt(1.0 - kick**2), rtol=1e-8)
    np.testing.assert_allclose(transmission, np.exp(-attenuation * thickness_at_axis), rtol=1e-10)

    print(f"Horizontal direction change: {kick * 1e6:.3f} µrad (small-angle equivalent)")
    print(f"Image displacement: {image_x * 1e6:.3f} µm")
    print(f"Intensity transmission: {transmission:.6f}")


if __name__ == "__main__":
    main()
