"""LQC-modified scalar field evolution in FRW (v8.1-alpha2).

Background:
    H^2 = (8 pi G / 3) rho (1 - rho / rho_c)
    rho = (1/2) phidot^2 + V(phi)
    H = adot / a

Klein-Gordon:
    phiddot + 3 H phidot + dV/dphi = 0

We integrate in cosmic time t and track a(t), phi(t), rho(t), H(t).
"""
from __future__ import annotations
import numpy as np

G = 6.67430e-11
C = 2.99792458e8
HBAR = 1.054571817e-34
RHO_PLANCK = C**5 / (HBAR * G**2)
RHO_C = 0.41 * RHO_PLANCK


def H_squared_lqc(rho, rho_c=RHO_C):
    return (8.0 * np.pi * G / 3.0) * rho * (1.0 - rho / rho_c)


def rho_scalar(phi, phi_dot, V_func):
    return 0.5 * phi_dot**2 + V_func(phi)


def evolve_lqc_scalar(phi0, phi_dot0, V_func, dV_func,
                       rho_c=RHO_C, dt=1e-40, n_steps=1000):
    """Explicit Euler integration of LQC FRW + scalar."""
    phi = phi0
    phi_dot = phi_dot0
    a = 1.0
    t_hist = [0.0]
    phi_hist = [phi]
    rho_hist = [rho_scalar(phi, phi_dot, V_func)]
    H_hist = [np.sqrt(max(H_squared_lqc(rho_hist[-1], rho_c), 0.0))]
    a_hist = [a]

    for step in range(1, n_steps + 1):
        rho = rho_scalar(phi, phi_dot, V_func)
        H2 = H_squared_lqc(rho, rho_c)
        if H2 >= 0:
            H = np.sqrt(H2)
        else:
            # In quantum regime, effective bounce: H -> 0
            H = 0.0
        # Klein-Gordon
        phi_ddot = -3.0 * H * phi_dot - dV_func(phi)
        phi_dot = phi_dot + dt * phi_ddot
        phi = phi + dt * phi_dot
        a = a + dt * H * a
        t_hist.append(step * dt)
        phi_hist.append(phi)
        rho_hist.append(rho_scalar(phi, phi_dot, V_func))
        H_hist.append(H)
        a_hist.append(a)

    return {
        "t": np.array(t_hist), "phi": np.array(phi_hist),
        "rho": np.array(rho_hist), "H": np.array(H_hist),
        "a": np.array(a_hist),
    }


__all__ = ["H_squared_lqc", "rho_scalar", "evolve_lqc_scalar",
           "RHO_C", "RHO_PLANCK"]
