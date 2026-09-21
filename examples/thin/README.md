# Refractor vs. phase deflector: when do they agree?

These four scripts compare a real refractive element (`S4Lens`, exact
Snell's-law ray/surface intersection) against the equivalent
`S4NumericalMeshPhaseDeflector` (paraxial thin-element approximation: advance
each ray to a flat reference plane, then apply a phase-gradient kick).

- `check_rafael_0.py` / `check_rafael_1.py` — a Be parabolic lens (`radius=5e-5`),
  compared against a phase deflector using `thin_lens.h5`, a thickness map
  built specifically to represent that lens.
- `check_srio_0.py` / `check_srio_1.py` — a Be lens built from
  `/home/srio/Oasys2/lens_interface_1.h5` (`S4Lens` with `surface_shape=0`
  plus `flag_add_mesh_surface_entrance=1`, i.e. a plane interface with the
  mesh height added on top), compared against a phase deflector using the
  same file directly as its thickness map.

## Result

The rafael pair agrees essentially exactly: same lost-ray count
(12371/500000 both), X/Z RMS matching to 5 significant figures, Xp/Zp RMS
matching to ~0.01 µrad. Deflection angles are small (Xp RMS ≈ 17.7 µrad).

The srio pair does not: lost-ray counts differ slightly (4309 vs 4306 of
5000), and X/Xp RMS differ by 1-2% (X: 80.2 vs 81.5 µm; Xp: 175.2 vs
172.6 µrad). Deflection angles here are ~10x larger (Xp RMS ≈ 175 µrad).

## Why

It's not a units or scaling bug: `S4Lens`'s mesh-loading path
(`S4AdditionalNumericalMeshInterface._apply_interface_refraction`) adds the
ideal surface's height (zero, for `surface_shape=0`) directly onto the raw,
unscaled mesh Z values before doing the real intersection — so `S4Lens` and
`S4NumericalMeshPhaseDeflector` consume the exact same numbers, with no
hidden conversion on either side. A uniform unit error would affect both
equally and can't explain a *relative* mismatch between them.

The real explanation is that `lens_interface_1.h5` isn't a "thin" element in
the sense the phase-deflector model assumes. Srio built it as a **cumulated
profile for several lenses** (a stack of lenslets collapsed into one
thickness map), not a single weak lens — which is why it fits a paraboloid
essentially exactly but with a very strong implied curvature (R ≈ 1.5 µm),
and why the deflection angles it produces are ~10x larger than rafael's
single-lens case.

`S4NumericalMeshPhaseDeflector` is a paraxial, thin-element approximation:
it advances rays to a flat plane and evaluates the phase gradient there,
which is only accurate when the height variation is small relative to the
aperture and deflection angles stay modest. `S4Lens` does an exact
geometric ray/surface intersection, with no such assumption. Rafael's
`thin_lens.h5` was built to respect that assumption (a single, genuinely
weak lens), so the two methods agree almost exactly. A cumulated multi-lens
profile, collapsed into a single thin element and evaluated at one z
position instead of tracing each lens (and the free space between them)
separately, pushes past where the paraxial approximation holds — so the
1-2% difference from the exact calculation is expected, not a bug.

**Takeaway:** `S4NumericalMeshPhaseDeflector` reproduces `S4Lens` results
closely for genuinely thin, weakly-curved elements. For a cumulated profile
standing in for several lenses at once, some deviation from the exact
Snell's-law result is inherent to the thin-element approximation, growing
with the profile's effective curvature and the resulting deflection angles.
