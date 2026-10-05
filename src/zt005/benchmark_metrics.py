"""Known-metric benchmark layer (ZT-006.3, PHASE 4/6 groundwork).

This module provides analytic benchmark metrics with known causal structure,
so that the project's geodesic/CTC machinery can be validated against exact
results BEFORE it is ever applied to device-generated fields.

Convention for this layer (differs from the legacy ZT-005 layer, see
docs/physics.md):

* signature (-, +, +, +)
* timelike interval: ds^2 < 0
* null interval:     ds^2 = 0
* spacelike interval: ds^2 > 0
* geometric units c = G = 1 unless stated otherwise

Metrics implemented
-------------------
* Minkowski            — flat reference, no CTC
* Schwarzschild        — static spherical mass, exterior r > 2M, no CTC
* Kerr (Boyer-Lindquist) — rotating black hole, no CTC outside the horizon
* Gödel (cylindrical)  — rotating dust universe, CTC for r > asinh(1)
* Tipler cylinder      — idealized rotating-cylinder metric, CTC for a*r > 1
* Morris-Thorne        — static traversable wormhole, no CTC in static form

References
----------
Gödel, K. (1949). Rev. Mod. Phys. 21, 447.
Tipler, F. J. (1974). Phys. Rev. D 9, 2203.
van Stockum, W. J. (1937). Proc. R. Soc. Edinb. A 57, 135.
Morris, M. S.; Thorne, K. S. (1988). Am. J. Phys. 56, 395.

The Tipler entry is a pedagogical idealization of the interior of a rigidly
rotating cylinder (van Stockum-type line element); it is NOT a claim that a
finite device realizes this metric.
"""

from __future__ import annotations

import numpy as np

SIGNATURE = (-1, +1, +1, +1)


class Metric:
    """Base class: subclasses implement g(x) -> (4, 4) array."""

    name = "base"
    coords = ("t", "x", "y", "z")
    admits_ctc = False

    def g(self, x):
        raise NotImplementedError

    def in_domain(self, x) -> bool:
        return True

    # ---- convenience -------------------------------------------------
    def ds2(self, x, dx):
        """Line element ds^2 = g_mu_nu dx^mu dx^nu at x along dx."""
        gx = np.asarray(self.g(np.asarray(x, float)), float)
        dx = np.asarray(dx, float)
        return float(dx @ gx @ dx)

    def christoffel_numeric(self, x, h=1e-5):
        """Christoffel symbols by central finite differences.

        Gamma^mu_{alpha beta} =
            1/2 g^{mu sigma} (d_a g_{s b} + d_b g_{s a} - d_s g_{a b})

        Accuracy is O(h^2); analytic Christoffels are preferred where
        available (see geodesic_solver.py).
        """
        x = np.asarray(x, float)
        g0 = np.asarray(self.g(x), float)
        ginv = np.linalg.inv(g0)
        dg = np.zeros((4, 4, 4))  # dg[c, a, b] = d_c g_ab
        for c in range(4):
            xp = x.copy(); xp[c] += h
            xm = x.copy(); xm[c] -= h
            dg[c] = (np.asarray(self.g(xp), float)
                     - np.asarray(self.g(xm), float)) / (2.0 * h)
        gamma = np.zeros((4, 4, 4))
        for mu in range(4):
            for a in range(4):
                for b in range(4):
                    s = 0.0
                    for sig in range(4):
                        s += ginv[mu, sig] * (
                            dg[a, sig, b] + dg[b, sig, a] - dg[sig, a, b])
                    gamma[mu, a, b] = 0.5 * s
        return gamma


class MinkowskiMetric(Metric):
    """Flat spacetime in Cartesian coordinates. No CTC exists."""

    name = "minkowski"
    coords = ("t", "x", "y", "z")
    admits_ctc = False

    def g(self, x):
        return np.diag([-1.0, 1.0, 1.0, 1.0])


class SchwarzschildMetric(Metric):
    """Schwarzschild exterior in (t, r, theta, phi), geometric units.

    No closed timelike curves exist in the exterior region r > 2M.
    """

    name = "schwarzschild"
    coords = ("t", "r", "theta", "phi")
    admits_ctc = False

    def __init__(self, M=1.0):
        if M <= 0:
            raise ValueError("M must be positive")
        self.M = float(M)

    def in_domain(self, x):
        return x[1] > 2.0 * self.M and 0.0 < x[2] < np.pi

    def g(self, x):
        M, r, th = self.M, x[1], x[2]
        f = 1.0 - 2.0 * M / r
        return np.diag([-f, 1.0 / f, r * r, r * r * np.sin(th) ** 2])

    def christoffel_analytic(self, x):
        """Exact Schwarzschild connection (exterior)."""
        M, r, th = self.M, x[1], x[2]
        f = 1.0 - 2.0 * M / r
        s, c = np.sin(th), np.cos(th)
        G = np.zeros((4, 4, 4))
        # Gamma^t_{tr} = 1/2 g^{tt} d_r g_{tt} = M / (r^2 f)
        G[0, 0, 1] = G[0, 1, 0] = M / (r * r * f)
        G[1, 0, 0] = M / (r * r) * f
        G[1, 1, 1] = -M / (r * r) / f
        G[1, 2, 2] = -r * f
        G[1, 3, 3] = -r * f * s * s
        G[2, 1, 2] = G[2, 2, 1] = 1.0 / r
        G[2, 3, 3] = -s * c
        G[3, 1, 3] = G[3, 3, 1] = 1.0 / r
        G[3, 2, 3] = G[3, 3, 2] = c / s if s != 0 else 0.0
        return G


class KerrMetric(Metric):
    """Kerr metric in Boyer-Lindquist (t, r, theta, phi), geometric units.

    The ring singularity's CTC region lies behind the horizon and is not part
    of the benchmarked domain; outside the horizon no CTC exists here.
    """

    name = "kerr"
    coords = ("t", "r", "theta", "phi")
    admits_ctc = False  # outside the outer horizon

    def __init__(self, M=1.0, a=0.7):
        if M <= 0 or abs(a) >= M:
            raise ValueError("require M > 0 and |a| < M")
        self.M, self.a = float(M), float(a)

    @property
    def r_plus(self):
        return self.M + np.sqrt(self.M ** 2 - self.a ** 2)

    def in_domain(self, x):
        return x[1] > self.r_plus and 0.0 < x[2] < np.pi

    def g(self, x):
        M, a = self.M, self.a
        r, th = x[1], x[2]
        s2 = np.sin(th) ** 2
        Sigma = r * r + a * a * np.cos(th) ** 2
        Delta = r * r - 2.0 * M * r + a * a
        g = np.zeros((4, 4))
        g[0, 0] = -(1.0 - 2.0 * M * r / Sigma)
        g[0, 3] = g[3, 0] = -2.0 * M * a * r * s2 / Sigma
        g[1, 1] = Sigma / Delta
        g[2, 2] = Sigma
        g[3, 3] = (r * r + a * a
                   + 2.0 * M * a * a * r * s2 / Sigma) * s2
        return g


class GodelMetric(Metric):
    """Gödel universe in cylindrical (t, r, phi, z), geometric units.

    ds^2 = (2/omega^2) [ -dt^2 - 2*sqrt(2)*sinh^2(r) dt dphi
                         + dr^2 + dz^2 + (sinh^2 r - sinh^4 r) dphi^2 ]

    The phi-coordinate circles are closed; they become timelike when
    sinh^2(r) > 1, i.e. r > asinh(1) = ln(1 + sqrt(2)) ~ 0.8814.
    This is the classic analytic CTC benchmark.
    """

    name = "godel"
    coords = ("t", "r", "phi", "z")
    admits_ctc = True
    CTC_RADIUS = float(np.arcsinh(1.0))  # = ln(1 + sqrt(2))

    def __init__(self, omega=1.0):
        if omega <= 0:
            raise ValueError("omega must be positive")
        self.omega = float(omega)

    def in_domain(self, x):
        return x[1] >= 0.0

    def g(self, x):
        r = x[1]
        w2 = self.omega ** 2
        sh = np.sinh(r)
        g = np.zeros((4, 4))
        g[0, 0] = -2.0 / w2
        g[0, 2] = g[2, 0] = -2.0 * np.sqrt(2.0) * sh * sh / w2
        g[1, 1] = 2.0 / w2
        g[3, 3] = 2.0 / w2
        g[2, 2] = 2.0 / w2 * (sh * sh - sh ** 4)
        return g

    def ctc_radius(self):
        return self.CTC_RADIUS


class TiplerCylinderMetric(Metric):
    """Idealized rotating-cylinder (Tipler/van Stockum-type) line element.

    ds^2 = -(dt + a r^2 dphi)^2 + dr^2 + dz^2 + r^2 dphi^2

    The coordinate phi-circles are timelike when a*r > 1.  This is the
    standard pedagogical idealization used in discussions of Tipler's
    rotating-cylinder causality violation; it is cylindrically symmetric
    and infinite in extent, so it is only ever used as a known-answer test
    for the CTC detector — never as a device model.
    """

    name = "tipler_cylinder"
    coords = ("t", "r", "phi", "z")
    admits_ctc = True

    def __init__(self, a=1.0):
        if a <= 0:
            raise ValueError("a must be positive")
        self.a = float(a)

    def in_domain(self, x):
        return x[1] >= 0.0

    def g(self, x):
        r, a = x[1], self.a
        g = np.zeros((4, 4))
        g[0, 0] = -1.0
        g[0, 2] = g[2, 0] = -a * r * r
        g[1, 1] = 1.0
        g[3, 3] = 1.0
        g[2, 2] = r * r * (1.0 - a * a * r * r)
        return g

    def ctc_radius(self):
        return 1.0 / self.a


class MorrisThorneMetric(Metric):
    """Static Morris-Thorne wormhole with zero redshift function.

    ds^2 = -dt^2 + dl^2 + (l^2 + b0^2)(dtheta^2 + sin^2(theta) dphi^2)

    Coordinates (t, l, theta, phi); the throat is at l = 0 with areal
    radius b0.  The static metric has no CTC; CTCs in the MT literature
    require converting one mouth into a time-shifted frame, which is NOT
    modeled here.
    """

    name = "morris_thorne"
    coords = ("t", "l", "theta", "phi")
    admits_ctc = False

    def __init__(self, b0=1.0):
        if b0 <= 0:
            raise ValueError("b0 must be positive")
        self.b0 = float(b0)

    def in_domain(self, x):
        return 0.0 < x[2] < np.pi

    def g(self, x):
        l, th, b0 = x[1], x[2], self.b0
        r2 = l * l + b0 * b0
        return np.diag([-1.0, 1.0, r2, r2 * np.sin(th) ** 2])


def all_benchmarks():
    """Default benchmark instances used by tests and the runner."""
    return [
        MinkowskiMetric(),
        SchwarzschildMetric(M=1.0),
        KerrMetric(M=1.0, a=0.7),
        GodelMetric(omega=1.0),
        TiplerCylinderMetric(a=1.0),
        MorrisThorneMetric(b0=1.0),
    ]
