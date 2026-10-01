"""
Monte Carlo model for mosaic crystals (S4MosaicCrystal with calculation_method=1).

Each ray is traced crystallite by crystallite through the mosaic crystal, as in the code
RTbent (L. Alianelli, 2002), with the corrections and the validation described in the
project monte-carlo-mosaics (https://github.com/srio/monte-carlo-mosaics, HOPG/README.md).

The model:

* The crystal is a slab of thickness T under the surface, with crystallites of thickness t0.
  Their normals n follow the mosaic distribution W: a Gaussian tilt, rms s per component,
  about the local lattice normal N (the surface normal, for a symmetric mosaic crystal).
* At each crystallite ("hop") the ray is reflected with probability

      P_acc(n) = min(1, peak exp(-WID (gamma - theta_B)^2) exp(mu tau)),

  with gamma the grazing angle of the ray on the crystallite planes, theta_B the Bragg angle,
  WID = 4 ln2 / width^2, and (peak, width) a Gaussian approximation of the rocking curve of a
  perfect t0 slab. The hop length is tau = tthi/cos(gamma), tthi = t0/tan(theta_B). Otherwise
  the ray crosses the crystallite in a straight line.
* A reflected ray takes the exact Bragg direction (deflection 2 theta_B in the plane of the
  ray and n) or, with flag_direction=0, the specular direction on the crystallite.
* The ray is traced until it crosses the entry face or the back face. It is reflected if it
  has been reflected an odd number of times, and transmitted otherwise. Absorption is the
  weight exp(-mu path).

Implementation: the uneventful hops are skipped exactly, by thinning. Candidate hops occur with
probability M >= E_W[P_acc] per hop, so the number of plain transmissions before the next
candidate is Geometric(M). At a candidate, the crystallite is drawn from a proposal h
centred where it meets the Bragg condition, parametrised by the tilt components (u, w)
along and across the ray. It then reflects with probability W(n) P_acc(n) / (M h(n)). The
per-hop reflection probability and the law of the reflecting crystallite are then exactly
those of the hop-by-hop process. The cost scales with the number of reflections, not of
hops.

Geometry: each ray sees the tangent plane of the surface at its entry point, with the lattice
normal equal to the surface normal there. This is exact for flat crystals, and a good
approximation for curved ones when the lateral travel of the ray inside the crystal
(~ 0.1-1 mm) is much smaller than the radius of curvature.

Approximations (negligible for mosaic widths << 1 rad): the tilt density is Gaussian in the
tilt vector; the skipped hops have the mean length tthi/cos(gamma0).
"""
import numpy

from crystalpy.diffraction.PerfectCrystalDiffraction import PerfectCrystalDiffraction
from crystalpy.util.ComplexAmplitudePhoton import ComplexAmplitudePhoton
from crystalpy.util.Vector import Vector

import scipy.constants as codata

from shadow4.tools.logger import is_verbose

MARGIN = 1.05             # safety factor of the thinning bound M (violations are counted)
N_INT_MAX = 5_000_000     # safety cap on the number of hops per ray (as in RTbent)
# Area / (peak * FWHM) of a Gaussian = 1.0645.
# A Gaussian of peak A and standard deviation sigma has area A sigma sqrt(2 pi) and
# FWHM = 2 sqrt(2 ln2) sigma, so area = A FWHM sqrt(2 pi) / (2 sqrt(2 ln2)) = A FWHM sqrt(pi / (4 ln2)).
# The crystallite reflectivity is the Gaussian peak exp(-4 ln2 Delta^2 / width^2) (RTbent's
# sgflag='G'), with width = FWHM of the crystalpy rocking curve of the t0 slab. It has the same
# integrated reflectivity as that curve if peak = integrated / (GAUSSIAN_AREA * width).
GAUSSIAN_AREA = numpy.sqrt(numpy.pi / (4 * numpy.log(2)))


#
# crystallite parameters (crystalpy)
#
def primary_extinction_depth(diffraction_setup, energy):
    """
    Amplitude primary extinction depth, normal to the surface, sigma polarization, for a
    symmetric Bragg reflection: lambda sin(theta_B) / (pi |psi_H|).

    Parameters
    ----------
    diffraction_setup : instance of crystalpy DiffractionSetupAbstract
    energy : float or numpy array
        Photon energy in eV.

    Returns
    -------
    float or numpy array
        The depth in m.
    """
    wavelength = codata.h * codata.c / (codata.e * numpy.asarray(energy, dtype=float))
    theta_B = diffraction_setup.angleBragg(energy)
    psi_H = numpy.abs(diffraction_setup.psiH(energy))
    return wavelength * numpy.sin(theta_B) / (numpy.pi * psi_H)


def _fwhm(x, y):
    """Full width at half maximum, from the outermost crossings, linearly interpolated."""
    half = y.max() / 2
    i = numpy.flatnonzero(y >= half)
    lo, hi = max(i[0], 1), min(i[-1], y.size - 2)
    x_lo = x[lo - 1] + (half - y[lo - 1]) * (x[lo] - x[lo - 1]) / (y[lo] - y[lo - 1])
    x_hi = x[hi] + (half - y[hi]) * (x[hi + 1] - x[hi]) / (y[hi + 1] - y[hi])
    return x_hi - x_lo


def _integrated(x, y, tail_fraction=0.3):
    """
    Integral of a rocking curve sampled on [-scan, scan], including the tails beyond the scan.

    The tails of a perfect-crystal curve fall as C / x^2 (kinematic slab and dynamical
    curves alike), so a scan of +-n widths misses about 1/(3n) of the integral. On each side,
    C is the mean of y x^2 over the outer tail_fraction of the scan (averaging the
    oscillations of the thin-slab curve), and C / |x_end| is added.
    """
    total = numpy.trapezoid(y, x)
    for side in (x < -(1 - tail_fraction) * abs(x[0]), x > (1 - tail_fraction) * x[-1]):
        if side.sum() > 2:
            C = numpy.mean(y[side] * x[side] ** 2)
            total += C / numpy.abs(x[side]).max()
    return total


def crystallite_rocking_curve(diffraction_setup, energy, t0, n_widths=50, points_per_width=40):
    """
    Rocking curve of a perfect crystallite of thickness t0 (symmetric Bragg, Guigay's method).

    The scan is centred on theta_B and scaled to the expected width of the curve: the
    kinematic width of the t0 slab, 0.886 lambda / (2 t0 cos theta_B), or the (pi) Darwin width
    if larger (thick crystallites), and limited to theta_B / 2. Use _integrated() to include the
    tails beyond the scan: with n_widths=50 the integrated reflectivity is then accurate to
    ~0.1% (HOPG 002, t0 = 0.05-2 extinction depths).

    Parameters
    ----------
    diffraction_setup : instance of crystalpy DiffractionSetupAbstract
    energy : float
        Photon energy in eV.
    t0 : float
        The crystallite thickness in m.
    n_widths : float, optional
        Half range of the angular scan, in units of the expected width of the curve.
    points_per_width : int, optional
        Number of scan points per (expected) width of the curve.

    Returns
    -------
    tuple
        (deviation [rad], reflectivity S, reflectivity P) numpy arrays.
    """
    theta_B = float(numpy.atleast_1d(diffraction_setup.angleBragg(energy))[0])
    wavelength = codata.h * codata.c / (codata.e * energy)
    # expected FWHM: the kinematic slab width, but not below the (pi) Darwin width
    darwin_p = 2 * float(numpy.atleast_1d(diffraction_setup.darwinHalfwidthS(energy))[0]) * \
        numpy.abs(numpy.cos(2 * theta_B))
    width_est = max(0.886 * wavelength / (2 * t0 * numpy.cos(theta_B)), darwin_p)
    # never scan beyond theta_B/2 from the Bragg angle (far from the reflection, at large
    # deviations, the slab formulas are meaningless)
    scan = min(n_widths * width_est, theta_B / 2)
    npoints = int(2 * scan / width_est * points_per_width) | 1
    deviation = numpy.linspace(-scan, scan, npoints)
    energies = numpy.full(npoints, float(energy))
    photons = ComplexAmplitudePhoton(energies,
                                     Vector(numpy.zeros(npoints), numpy.cos(theta_B + deviation),
                                            -numpy.sin(theta_B + deviation)),
                                     Esigma=numpy.ones(npoints, dtype=complex),
                                     Epi=numpy.ones(npoints, dtype=complex))
    perfect_crystal = PerfectCrystalDiffraction.initializeFromDiffractionSetupAndEnergy(
        diffraction_setup, energies, thickness=t0, calculation_strategy_flag=1)
    coefficients = perfect_crystal.calculateDiffraction(photons, calculation_method=1, is_thick=0,
                                                        use_transfer_matrix=0)
    return deviation, numpy.abs(coefficients["S"]) ** 2, numpy.abs(coefficients["P"]) ** 2


def crystallite_parameters(diffraction_setup, energies, t0_flag, t0_factor, t0_value, max_energies=21,
                           return_profiles=False):
    """
    Crystallite inputs of the Monte Carlo for each energy: t0, and the (peak, width) of the
    Gaussian that has the same FWHM and integrated reflectivity as the rocking curve of the
    t0 slab, for S and P polarizations.

    The rocking curves are computed at the distinct energies of the beam, or, if there are more
    than max_energies of them, on a regular grid between their extremes, and interpolated.

    Parameters
    ----------
    diffraction_setup : instance of crystalpy DiffractionSetupAbstract
    energies : numpy array
        The photon energies of the rays in eV.
    t0_flag : int
        0: t0 = t0_factor times the primary extinction depth, 1: t0 = t0_value.
    t0_factor : float
        See t0_flag.
    t0_value : float
        See t0_flag (in m).
    max_energies : int, optional
        Maximum number of energies where the rocking curve is computed.
    return_profiles : bool, optional
        If True, also return the computed crystallite rocking curves.

    Returns
    -------
    dict, or tuple (dict, list)
        A dict of arrays with one value per ray: t0 [m], peak_s, width_s [rad], peak_p,
        width_p [rad]. With return_profiles=True, also a list with one dict per computed
        energy:

        * energy [eV], theta_B [rad], depth (amplitude primary extinction depth, sigma) [m],
          t0 [m], deviation [rad] (angle from theta_B of the scan points);
        * "s" and "p": dicts with reflectivity (the crystalpy rocking curve of the t0 slab
          on the deviation points), gaussian (the Gaussian used by the Monte Carlo, on the same
          points), width (FWHM) [rad], peak (of the Gaussian), maximum (of the curve),
          integrated [rad] (including the tails beyond the scan) and f_primary (integrated /
          kinematic integrated reflectivity Q t0 / sin(theta_B)).
    """
    energies = numpy.asarray(energies, dtype=float)
    distinct = numpy.unique(energies)
    if distinct.size <= max_energies:
        grid = distinct
    elif max_energies == 1:      # a single curve: at the middle of the energy range
        grid = numpy.array([0.5 * (distinct[0] + distinct[-1])])
    else:
        grid = numpy.linspace(distinct[0], distinct[-1], max_energies)
    rows = []
    details = []   # the computed profiles (returned with return_profiles=True, and printed if verbose)
    for e in grid:
        depth = float(numpy.atleast_1d(primary_extinction_depth(diffraction_setup, e))[0])
        t0 = t0_factor * depth if t0_flag == 0 else t0_value
        x, rs, rp = crystallite_rocking_curve(diffraction_setup, e, t0)
        row = [t0]
        detail = dict(energy=e, depth=depth, t0=t0,
                      theta_B=float(numpy.atleast_1d(diffraction_setup.angleBragg(e))[0]),
                      deviation=x)
        # kinematic integrated reflectivity of the t0 slab, Q t0 / sin(theta_B), for f_primary
        wavelength = codata.h * codata.c / (codata.e * e)
        Q_s = numpy.pi ** 2 * float(numpy.abs(numpy.atleast_1d(diffraction_setup.psiH(e))[0] *
                                              numpy.atleast_1d(diffraction_setup.psiH_bar(e))[0])) / \
            (wavelength * numpy.sin(2 * detail["theta_B"]))
        Q = dict(s=Q_s, p=Q_s * numpy.cos(2 * detail["theta_B"]) ** 2)
        for pol, r in (("s", rs), ("p", rp)):
            width = _fwhm(x, r)
            integrated = _integrated(x, r)
            peak = integrated / (GAUSSIAN_AREA * width)
            row += [peak, width]
            detail[pol] = dict(width=width, peak=peak, maximum=r.max(), integrated=integrated,
                               f_primary=integrated / (Q[pol] * t0 / numpy.sin(detail["theta_B"])),
                               reflectivity=r,
                               gaussian=peak * numpy.exp(-4 * numpy.log(2) * (x / width) ** 2))
        rows.append(row)
        details.append(detail)
    rows = numpy.array(rows)
    keys = ["t0", "peak_s", "width_s", "peak_p", "width_p"]
    if grid.size == 1:
        out = {k: numpy.full(energies.size, rows[0, i]) for i, k in enumerate(keys)}
    else:
        out = {k: numpy.interp(energies, grid, rows[:, i]) for i, k in enumerate(keys)}

    if is_verbose():
        print("\nMonte Carlo mosaic crystal: crystallite parameters")
        if t0_flag == 0:
            print(f"   crystallite thickness t0: automatic, {t0_factor:g} x primary extinction depth")
        else:
            print(f"   crystallite thickness t0: user-defined, {1e6 * t0_value:g} um")
        how = "exact" if grid.size == distinct.size else "interpolated"
        print(f"   rays: {energies.size}, distinct energies: {distinct.size}, "
              f"rocking curves computed at {grid.size} energies ({how}), "
              f"from {grid[0]:.3f} to {grid[-1]:.3f} eV")
        print("   depth: amplitude primary extinction depth (sigma)")
        print("   peak = integrated / (1.0645 width)")
        print("   f_primary: integrated / kinematic (Q t0 / sin theta_B), 1 = no primary extinction")
        print(f"   {'energy [eV]':>12s} {'thetaB [deg]':>12s} {'depth [um]':>10s} {'t0 [um]':>9s} "
              f"{'t0/depth':>8s} {'pol':>3s} {'width [urad]':>12s} {'peak':>9s} {'maximum':>9s} "
              f"{'integr. [urad]':>14s} {'f_primary':>9s}")
        for d in details:
            for pol in ("s", "p"):
                c = d[pol]
                head = (f"   {d['energy']:12.3f} {numpy.degrees(d['theta_B']):12.5f} {1e6 * d['depth']:10.5g} "
                        f"{1e6 * d['t0']:9.5g} {d['t0'] / d['depth']:8.4g}") if pol == "s" else " " * 58
                print(f"{head} {pol.upper():>3s} {1e6 * c['width']:12.5g} {c['peak']:9.5g} {c['maximum']:9.5g} "
                      f"{1e6 * c['integrated']:14.5g} {c['f_primary']:9.4g}")
        print("   per-ray values (t0 [m], width [rad]):")
        for k in keys:
            print(f"   {k:>10s}: mean: {out[k].mean():10.5g}, stdev: {round(out[k].std(), 10):10.5g}")
    if return_profiles:
        return out, details
    return out

#
# Monte Carlo
#
def _hop_length(tthi, gamma):
    """Path inside a crystallite for a grazing angle gamma (RTbent's hop geometry)."""
    return numpy.where(gamma <= numpy.pi / 4, tthi / numpy.cos(gamma), tthi / numpy.sin(gamma))


def _reflect(V, n, c, theta_B, flag_direction):
    """Specular reflection on the crystallite (0), or deflection by exactly 2 theta_B (1)."""
    n_facing = -numpy.sign(c)[:, None] * n
    Vout = V + 2 * numpy.abs(c)[:, None] * n_facing
    if flag_direction == 1:
        U = Vout - numpy.sum(Vout * V, axis=1)[:, None] * V
        U /= numpy.linalg.norm(U, axis=1)[:, None]
        Vout = numpy.cos(2 * theta_B)[:, None] * V + numpy.sin(2 * theta_B)[:, None] * U
    return Vout


def _tangent_frame(N):
    """Unit vectors e1 (along x projected on the surface, or y if N ~ x) and e2 = N x e1."""
    ref = numpy.zeros_like(N)
    ref[:, 0] = 1.0
    parallel = numpy.abs(N[:, 0]) > 0.9
    ref[parallel] = [0.0, 1.0, 0.0]
    e1 = ref - numpy.sum(ref * N, axis=1)[:, None] * N
    e1 /= numpy.linalg.norm(e1, axis=1)[:, None]
    return e1, numpy.cross(N, e1)


def trace_mosaic_monte_carlo(X_in, V_in, N_in, theta_B, mu, thickness, s_tilt, tthi, peak, width,
                             flag_direction=1, rng=None):
    """
    Traces the rays through the mosaic crystal crystallite by crystallite (see the module doc).

    Parameters
    ----------
    X_in : numpy array (nrays, 3)
        The entry points on the surface [m].
    V_in : numpy array (nrays, 3)
        The incident unit directions (going into the crystal: V_in . N_in < 0).
    N_in : numpy array (nrays, 3)
        The unit surface normals at the entry points, pointing out of the crystal (toward the
        incident beam). They define the tangent plane and the mean lattice normal of each ray.
    theta_B : numpy array (nrays,)
        The Bragg angle of each ray [rad].
    mu : numpy array (nrays,)
        The linear absorption coefficient [m^-1].
    thickness : float
        The crystal thickness [m].
    s_tilt : float
        The rms of each component of the crystallite tilt [rad] (mosaic FWHM / sqrt(8 ln 2)).
    tthi : numpy array (nrays,)
        The crystallite hop parameter t0 / tan(theta_B) [m].
    peak, width : numpy arrays (nrays,)
        The peak and FWHM [rad] of the Gaussian crystallite rocking curve.
    flag_direction : int, optional
        Reflected direction: 0=specular on the crystallite, 1=exact Bragg deflection 2 theta_B.
    rng : numpy.random.Generator, optional
        The random generator.

    Returns
    -------
    dict
        X_out, V_out (nrays, 3): exit points (on a crystal face) and directions;
        path (nrays,): path inside the crystal [m]; n_reflections, n_hops (nrays,);
        flag (nrays,): 1 reflected, 0 transmitted, -1 stopped at the safety cap;
        stats: dict with the number of iterations, candidate hops, reflections, and bound
        violations (should be 0).
    """
    rng = numpy.random.default_rng() if rng is None else rng
    n = X_in.shape[0]
    X, V = numpy.array(X_in, dtype=float), numpy.array(V_in, dtype=float)
    N = numpy.array(N_in, dtype=float)
    N /= numpy.linalg.norm(N, axis=1)[:, None]
    V /= numpy.linalg.norm(V, axis=1)[:, None]
    e1_all, e2_all = _tangent_frame(N)
    X0 = X.copy()
    theta_B, mu, tthi = numpy.broadcast_to(theta_B, n), numpy.broadcast_to(mu, n), numpy.broadcast_to(tthi, n)
    peak, WID = numpy.broadcast_to(peak, n), 4 * numpy.log(2) / numpy.broadcast_to(width, n) ** 2
    s2 = s_tilt ** 2

    path = numpy.zeros(n)
    n_hops = numpy.zeros(n, dtype=numpy.int64)
    n_ref = numpy.zeros(n, dtype=numpy.int64)
    flag = numpy.full(n, -1)
    stats = dict(iterations=0, candidates=0, reflections=0, violations=0, max_ratio=0.0)

    def face_distance(i, Vi):
        """Path along Vi from X[i] to the entry face (depth 0) or the back face (-thickness)."""
        depth = numpy.sum((X[i] - X0[i]) * N[i], axis=1)
        vz = numpy.sum(Vi * N[i], axis=1)
        with numpy.errstate(divide="ignore", invalid="ignore"):
            s = numpy.where(vz < 0, (depth + thickness) / -vz, numpy.where(vz > 0, -depth / vz, numpy.inf))
        return numpy.maximum(s, 0.0)

    def geometry(i, Vi):
        """Mean grazing angle gamma0 and the slope k = d gamma / d u at zero tilt."""
        c0 = numpy.sum(Vi * N[i], axis=1)
        g0 = numpy.arcsin(numpy.minimum(numpy.abs(c0), 1.0))
        a, b = numpy.sum(Vi * e1_all[i], axis=1), numpy.sum(Vi * e2_all[i], axis=1)
        g = numpy.hypot(a, b)
        return g0, a, b, g, numpy.sign(c0) * g / numpy.cos(g0)

    def crystallite(i, a, b, g, u, w):
        """Crystallite normal for the tilt (u, w) along/across the ray, about N[i]."""
        tx, ty = (a * u - b * w) / g, (b * u + a * w) / g
        pt = numpy.hypot(tx, ty)
        sinc = numpy.where(pt > 0, numpy.sin(pt) / numpy.where(pt > 0, pt, 1.0), 1.0)
        return (sinc * tx)[:, None] * e1_all[i] + (sinc * ty)[:, None] * e2_all[i] + \
            numpy.cos(pt)[:, None] * N[i]

    def grazing(Vi, nvec):
        c = numpy.sum(Vi * nvec, axis=1)
        return c, numpy.arcsin(numpy.minimum(numpy.abs(c), 1.0))

    def finish(i):
        flag[i] = numpy.where(n_ref[i] % 2 == 1, 1, 0)

    alive = numpy.arange(n)
    while alive.size:
        stats["iterations"] += 1
        i = alive
        Vi = V[i]
        g0, a, b, g, k = geometry(i, Vi)
        tau_bar = _hop_length(tthi[i], g0)
        # the lattice is the same along the straight flight (tangent plane): tight bound
        A2 = 2 * WID[i] * k ** 2 * s2
        D = theta_B[i] - g0
        Z = numpy.exp(-WID[i] * D ** 2 / (1 + A2)) / numpy.sqrt(1 + A2)
        M = numpy.clip(peak[i] * numpy.exp(mu[i] * tau_bar * 1.01) * Z * MARGIN, 1e-12, 1.0)
        K = (rng.geometric(M) - 1).astype(float)

        # plain transmitted hops: exit if the flight reaches a face first
        s_face = face_distance(i, Vi)
        out = K * tau_bar >= s_face
        if out.any():
            j = i[out]
            X[j] += s_face[out, None] * V[j]
            path[j] += s_face[out]
            n_hops[j] += numpy.maximum(numpy.ceil(s_face[out] / tau_bar[out]), 1).astype(numpy.int64)
            finish(j)

        # candidate hop
        c_ = ~out
        j = i[c_]
        if j.size:
            stats["candidates"] += j.size
            X[j] += (K[c_] * tau_bar[c_])[:, None] * V[j]
            path[j] += K[c_] * tau_bar[c_]
            n_hops[j] += K[c_].astype(numpy.int64) + 1
            Vj, tBj, WIDj = V[j], theta_B[j], WID[j]
            g0j, aj, bj, gj, kj = g0[c_], a[c_], b[c_], g[c_], k[c_]
            w = numpy.sqrt(s2) * rng.standard_normal(j.size)
            # Bragg tilt u* for this w: two Newton steps from the linear guess
            u_star = (tBj - g0j) / kj
            for _ in range(2):
                _, gam = grazing(Vj, crystallite(j, aj, bj, gj, u_star, w))
                u_star -= (gam - tBj) / kj
            A = 1 / s2 + 2 * WIDj * kj ** 2
            m = 2 * WIDj * kj ** 2 * u_star / A
            u = m + rng.standard_normal(j.size) / numpy.sqrt(A)
            nvec = crystallite(j, aj, bj, gj, u, w)
            c, gam = grazing(Vj, nvec)
            tau_r = _hop_length(tthi[j], gam)
            P_acc = numpy.minimum(1.0, peak[j] * numpy.exp(-WIDj * (gam - tBj) ** 2 + mu[j] * tau_r))
            log_w = -u ** 2 / (2 * s2) - 0.5 * numpy.log(s2)      # W(u), up to 1/sqrt(2 pi)
            log_h = -A * (u - m) ** 2 / 2 + 0.5 * numpy.log(A)     # h(u|w), same constant
            ratio = numpy.exp(log_w - log_h) * P_acc / M[c_]
            stats["violations"] += int(numpy.sum(ratio > 1))
            stats["max_ratio"] = max(stats["max_ratio"], float(ratio.max()))
            acc = rng.random(j.size) < ratio

            # reflection: new direction, then a hop of length tau_r along it (or to the face)
            r = j[acc]
            if r.size:
                stats["reflections"] += r.size
                Vn = _reflect(Vj[acc], nvec[acc], c[acc], tBj[acc], flag_direction)
                n_ref[r] += 1
                V[r] = Vn
                sf = face_distance(r, Vn)
                step = numpy.minimum(sf, tau_r[acc])
                X[r] += step[:, None] * Vn
                path[r] += step
                finish(r[sf <= tau_r[acc]])
            # rejected candidate: a transmitted hop of mean length
            t = j[~acc]
            if t.size:
                tb = _hop_length(tthi[t], g0j[~acc])
                sf = face_distance(t, V[t])
                step = numpy.minimum(sf, tb)
                X[t] += step[:, None] * V[t]
                path[t] += step
                finish(t[sf <= tb])

        still = alive[flag[alive] == -1]
        alive = still[n_hops[still] < N_INT_MAX]

    return dict(X_out=X, V_out=V, path=path, n_reflections=n_ref, n_hops=n_hops, flag=flag, stats=stats)
