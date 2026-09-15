"""Retarded linearized-Einstein source solver (ZT-006.3).

This is the module that replaces the frozen/instantaneous Green kernel as
the physically correct source-to-metric transition for time-dependent
sources.  It implements the retarded solution of the linearized Einstein
equation in the Lorenz gauge,

    Box hbar_munu = -16 pi G / c^4 * T_munu        (labeled BLOCKER item:
                                                     "retarded spacetime
                                                     source solver")

    hbar_munu(t, x) = (4G / c^4) * INT T_munu(t - |x-x'|/c, x')
                                     / |x-x'| d^3x'

Note the integrand is T_munu itself, NOT its trace reverse; the
trace-reversal lives on the h-side (hbar_munu = h_munu - 1/2 eta_munu h).
Lorenz-gauge consistency then follows from stress-energy conservation
partial_mu T^{mu nu} = 0 alone (see gauge_checks.py).

Conventions
-----------
* Coordinates x^mu = (c t, x, y, z); signature (-, +, +, +).
* Component values of the trace-reversed tensor are identical in the
  (+---) legacy convention used by the ZT-005 EM layer, so existing
  stress-energy arrays can be fed in unchanged (see docs/physics.md).
* Units: SI.

Honesty constraints
-------------------
* This is LINEARIZED gravity on a fixed Minkowski background.  It does not
  produce a full g_munu(t,x) solving the nonlinear Einstein equations.
* No CTC or time-machine conclusion follows from this module.
"""

from __future__ import annotations

import numpy as np

C = 299_792_458.0
G = 6.67430e-11
KAPPA = 4.0 * G / C ** 4          # prefactor of the retarded integral
ETA = np.diag([-1.0, 1.0, 1.0, 1.0])


def singular_cell_kernel(volume):
    """Cell-averaged 1/R kernel for a cell containing the observer.

    The midpoint approximation vol/R diverges when the observer sits on a
    source point, but the true integral of 1/|x-x'| over the cell is
    finite.  We replace such cells by the spherical average with the same
    volume,  <1/R>_cell = 3 / (2 r_cell),  r_cell = (3 vol / 4 pi)^{1/3}.
    This is the standard singular-cell correction for midpoint quadrature
    of 1/R kernels and removes a ~1e12 spurious spike per coincident cell.
    """
    r_cell = (3.0 * np.asarray(volume, float) / (4.0 * np.pi)) ** (1.0 / 3.0)
    return 1.5 / r_cell


def trace_reverse(T):
    """Tbar_munu = T_munu - 1/2 eta_munu T, with T = eta^{ab} T_ab.

    Component values are convention-independent between (-+++) and (+---);
    see module docstring.
    """
    T = np.asarray(T, float)
    single = T.ndim == 2
    Tb = T if single else None
    if single:
        tr = np.einsum("ab,ab->", ETA, T)
        return T - 0.5 * ETA * tr
    tr = np.einsum("ab,...ab->...", ETA, T)
    return T - 0.5 * ETA * tr[..., None, None]


def retarded_hbar(source_fn, obs_points, obs_times, src_points,
                  src_weights, trace_reverse_source=False,
                  softening_m=0.0):
    """Evaluate the retarded trace-reversed metric perturbation.

    Parameters
    ----------
    source_fn : callable(t_array (Ns,), x_array (Ns, 3)) -> (Ns, 4, 4)
        Stress-energy T^{mu nu} evaluated at per-point retarded times and
        positions.  Trace-reversing the source here would be WRONG for the
        retarded integral; keep the default False.
    obs_points : (No, 3) observation positions [m]
    obs_times  : (Nt,) observation times [s]
    src_points : (Ns, 3) source quadrature points [m]
    src_weights: (Ns,) quadrature volume weights [m^3]

    Returns
    -------
    hbar : (Nt, No, 4, 4) trace-reversed perturbation in the (-+++) layer.
    """
    obs = np.asarray(obs_points, float)
    src = np.asarray(src_points, float)
    w = np.asarray(src_weights, float)
    times = np.atleast_1d(np.asarray(obs_times, float))
    if obs.ndim == 1:
        obs = obs[None, :]

    out = np.zeros((len(times), len(obs), 4, 4))
    for i, xo in enumerate(obs):
        R = np.linalg.norm(src - xo, axis=1)
        if softening_m > 0:
            R = np.sqrt(R * R + softening_m ** 2)
        coincident = R < 1e-9
        kernel = np.where(coincident, 1.0, w / np.maximum(R, 1e-9))
        if np.any(coincident):
            kernel[coincident] = (singular_cell_kernel(w[coincident])
                                  * w[coincident])
        for j, t in enumerate(times):
            t_ret = t - R / C
            T = np.asarray(source_fn(t_ret, src), float)
            if trace_reverse_source:
                T = trace_reverse(T)
            out[j, i] = KAPPA * np.einsum("n,nab->ab", kernel, T)
    return out


def static_hbar(source_T, obs_points, src_points, src_weights,
                trace_reverse_source=False, softening_m=0.0):
    """Instantaneous (t-independent) limit of the retarded integral.

    Used only as a convergence/reference diagnostic: a truly static source
    must reproduce this from :func:`retarded_hbar`.
    """
    obs = np.asarray(obs_points, float)
    if obs.ndim == 1:
        obs = obs[None, :]
    src = np.asarray(src_points, float)
    T = np.asarray(source_T, float)
    if trace_reverse_source:
        T = trace_reverse(T)
    w = np.asarray(src_weights, float)
    out = np.zeros((len(obs), 4, 4))
    for i, xo in enumerate(obs):
        R = np.linalg.norm(src - xo, axis=1)
        if softening_m > 0:
            R = np.sqrt(R * R + softening_m ** 2)
        coincident = R < 1e-9
        kernel = np.where(coincident, 1.0, w / np.maximum(R, 1e-9))
        if np.any(coincident):
            kernel[coincident] = (singular_cell_kernel(w[coincident])
                                  * w[coincident])
        out[i] = KAPPA * np.einsum("n,nab->ab", kernel, T)
    return out
