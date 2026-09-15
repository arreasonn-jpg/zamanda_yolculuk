"""Geodesic integrator (ZT-006.3, PHASE 5 groundwork).

Integrates

    d^2 x^mu / d tau^2 = - Gamma^mu_{alpha beta}
                           (dx^alpha/dtau) (dx^beta/dtau)

with a fixed-step RK4 scheme.  Christoffel symbols come either from an
analytic provider (preferred where known) or from central finite
differences of the metric (Metric.christoffel_numeric).

Validation strategy: the norm g_mu_nu u^mu u^nu is an exact constant of
motion along any geodesic; drift in that norm is the primary integrator
quality diagnostic (see tests/test_geodesic_solver.py).

No physics conclusion about any device is made by this module.
"""

from __future__ import annotations

import numpy as np


def _make_christoffel_provider(metric, analytic=None, numeric_step=1e-5):
    if analytic is not None:
        return analytic
    return lambda x: metric.christoffel_numeric(x, h=numeric_step)


def christoffel_from_metric(metric, analytic=False, numeric_step=1e-5):
    """Return a callable x -> Gamma^mu_{alpha beta} for a benchmark metric."""
    if analytic:
        fn = getattr(metric, "christoffel_analytic", None)
        if fn is None:
            raise ValueError(f"no analytic Christoffels for {metric.name}")
        return fn
    return _make_christoffel_provider(metric, None, numeric_step)


def _accel(x, u, gamma_fn):
    G = gamma_fn(x)
    return -np.einsum("mab,a,b->m", G, u, u)


def integrate_geodesic(metric, x0, u0, tau_span, n_steps,
                       christoffel=None, analytic=False,
                       numeric_step=1e-5):
    """Integrate one geodesic.

    Parameters
    ----------
    metric      : benchmark_metrics.Metric instance
    x0, u0      : initial position and 4-velocity dx^mu/dtau
    tau_span    : total affine-parameter length to integrate
    n_steps     : number of RK4 steps
    christoffel : optional callable x -> Gamma (overrides ``analytic``)
    analytic    : use metric.christoffel_analytic when available

    Returns dict with trajectory arrays and the norm-drift diagnostic.
    """
    gamma_fn = christoffel or christoffel_from_metric(
        metric, analytic=analytic, numeric_step=numeric_step)

    x = np.asarray(x0, float).copy()
    u = np.asarray(u0, float).copy()
    dt = float(tau_span) / int(n_steps)

    xs = np.zeros((n_steps + 1, 4))
    us = np.zeros((n_steps + 1, 4))
    norms = np.zeros(n_steps + 1)
    xs[0], us[0] = x, u
    norms[0] = float(u @ metric.g(x) @ u)

    for i in range(n_steps):
        k1x = u
        k1u = _accel(x, u, gamma_fn)
        k2x = u + 0.5 * dt * k1u
        k2u = _accel(x + 0.5 * dt * k1x, u + 0.5 * dt * k1u, gamma_fn)
        k3x = u + 0.5 * dt * k2u
        k3u = _accel(x + 0.5 * dt * k2x, u + 0.5 * dt * k2u, gamma_fn)
        k4x = u + dt * k3u
        k4u = _accel(x + dt * k3x, u + dt * k3u, gamma_fn)
        x = x + (dt / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        u = u + (dt / 6.0) * (k1u + 2 * k2u + 2 * k3u + k4u)
        xs[i + 1], us[i + 1] = x, u
        norms[i + 1] = float(u @ metric.g(x) @ u)

    ref = max(abs(norms[0]), 1e-300)
    return {
        "x": xs,
        "u": us,
        "tau": np.linspace(0.0, tau_span, n_steps + 1),
        "norm0": float(norms[0]),
        "norm_final": float(norms[-1]),
        "norm_drift_rel": float(np.max(np.abs(norms - norms[0])) / ref),
    }


def circular_orbit_initial_conditions_schwarzschild(M, r0):
    """Timelike circular equatorial orbit at radius r0 (r0 > 3M required).

    Returns (x0, u0) with u normalized so g(u,u) = -1.
    """
    if r0 <= 3.0 * M:
        raise ValueError("stable circular timelike orbits require r0 > 3M")
    omega = np.sqrt(M / r0 ** 3)          # dphi/dt (Kepler)
    f = 1.0 - 2.0 * M / r0
    ut = 1.0 / np.sqrt(f - r0 * r0 * omega * omega)
    x0 = np.array([0.0, r0, np.pi / 2, 0.0])
    u0 = np.array([ut, 0.0, 0.0, omega * ut])
    return x0, u0
