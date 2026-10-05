"""Observer Geodesic Evolution and Proper-Time / Coordinate-Time Tracker.

Addresses Section 9, 10, and 11:
- Integrates the timelike geodesic equation:
    d^2 x^mu / d tau^2 + Gamma^mu_{alpha beta} (dx^alpha/dtau) (dx^beta/dtau) = 0
- Tracks both coordinate time t = x^0/c and proper time tau
- Tracks elapsed coordinate time Delta t and elapsed proper time Delta tau > 0
- Reconstructs arrival coordinate: x_arrival = (x(tau_f), y(tau_f), z(tau_f))
- Continuously evaluates interval classification:
    ds^2 = g_munu dx^mu dx^nu  ->  timelike (< 0), null (= 0), spacelike (> 0)
- Scenario:
    Earth surface -> active region -> metric generated -> observer enters field
    -> geodesic evolves -> worldline reconstructed.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

C = 299_792_458.0


@dataclass
class WorldlinePoint:
    tau: float               # Proper time [s]
    t: float                 # Coordinate time [s]
    x: np.ndarray            # 4-position (c*t, x, y, z) [m]
    u: np.ndarray            # 4-velocity dx^mu/dtau
    ds2: float               # Local line element ds^2 = g(u, u)
    is_timelike: bool        # Whether ds^2 < 0
    dt_dtau: float           # Coordinate time rate dt/dtau = u^0 / c


def evolve_traveler_geodesic(
    metric,
    x0: np.ndarray,
    u0: np.ndarray,
    tau_span: float,
    n_steps: int = 200,
    gamma_fn=None,
) -> dict:
    """Evolve a traveler's geodesic through spacetime and record full kinematic trajectory.

    Parameters
    ----------
    metric : Metric instance providing g(x)
    x0 : initial 4-position (c*t0, x0, y0, z0) [m]
    u0 : initial 4-velocity dx^mu/dtau (must be timelike: u @ g @ u < 0)
    tau_span : total proper time length [s]
    n_steps : number of RK4 integration steps
    gamma_fn : callable x -> Gamma^mu_{ab} (optional, defaults to metric.christoffel_numeric)
    """
    if gamma_fn is None:
        gamma_fn = lambda pos: metric.christoffel_numeric(pos, h=1e-5)

    x = np.asarray(x0, float).copy()
    u = np.asarray(u0, float).copy()
    dtau = float(tau_span) / int(n_steps)

    g0 = metric.g(x)
    norm0 = float(u @ g0 @ u)
    if norm0 >= 0:
        raise ValueError(f"Initial velocity must be timelike (g(u,u) < 0), got norm = {norm0}")

    def accel(pos, vel):
        G = gamma_fn(pos)
        return -np.einsum("mab,a,b->m", G, vel, vel)

    history: list[WorldlinePoint] = []
    tau = 0.0

    # Initial state
    history.append(WorldlinePoint(
        tau=tau,
        t=float(x[0] / C),
        x=x.copy(),
        u=u.copy(),
        ds2=norm0,
        is_timelike=bool(norm0 < -1e-12),
        dt_dtau=float(u[0] / C),
    ))

    norms = [norm0]

    for step in range(n_steps):
        k1x = u
        k1u = accel(x, u)

        k2x = u + 0.5 * dtau * k1u
        k2u = accel(x + 0.5 * dtau * k1x, u + 0.5 * dtau * k1u)

        k3x = u + 0.5 * dtau * k2u
        k3u = accel(x + 0.5 * dtau * k2x, u + 0.5 * dtau * k2u)

        k4x = u + dtau * k3u
        k4u = accel(x + dtau * k3x, u + dtau * k3u)

        x = x + (dtau / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        u = u + (dtau / 6.0) * (k1u + 2 * k2u + 2 * k3u + k4u)
        tau += dtau

        g_curr = metric.g(x)
        curr_norm = float(u @ g_curr @ u)
        norms.append(curr_norm)

        history.append(WorldlinePoint(
            tau=tau,
            t=float(x[0] / C),
            x=x.copy(),
            u=u.copy(),
            ds2=curr_norm,
            is_timelike=bool(curr_norm < -1e-12),
            dt_dtau=float(u[0] / C),
        ))

    ref_norm = max(abs(norm0), 1e-300)
    norm_drift_rel = float(np.max(np.abs(np.array(norms) - norm0)) / ref_norm)

    t_initial = history[0].t
    t_final = history[-1].t
    delta_t = t_final - t_initial
    delta_tau = tau

    all_timelike = all(pt.is_timelike for pt in history)
    # Check if backward coordinate time evolution occurred (time reversal)
    has_negative_dt = any(pt.dt_dtau < 0 for pt in history)

    return {
        "initial_position": history[0].x.tolist(),
        "arrival_position": history[-1].x.tolist(),
        "initial_coordinate_time_s": t_initial,
        "arrival_coordinate_time_s": t_final,
        "elapsed_coordinate_time_delta_t_s": delta_t,
        "elapsed_proper_time_delta_tau_s": delta_tau,
        "ratio_dt_dtau_mean": float(delta_t / max(delta_tau, 1e-30)),
        "all_timelike": all_timelike,
        "causality_inversion_detected": has_negative_dt,
        "norm_drift_rel": norm_drift_rel,
        "n_steps": n_steps,
        "trajectory_points": len(history),
    }
