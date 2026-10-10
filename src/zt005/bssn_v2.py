"""v2.0 unified BSSN driver: state + Gamma + gauge + Z4c damping, RK4."""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState, minkowski_state
from .bssn_rhs_grid import bssn_rhs_grid
from .bssn_gamma import dGamma_RHS
from .bssn_gauge import dt_alpha_1pluslog, dt_beta_gamma_driver
from .bssn_z4c import damping_dK, damping_dGamma
from .bssn import BSSN, hamiltonian_constraint, momentum_constraint
from .bssn_grid import d_central


def _rhs(state: BSSNState, h):
    """Full v2.0 RHS: returns dict of derivatives."""
    r = bssn_rhs_grid(state.phi, state.gtilde, state.K, state.Atilde, h)
    dG = dGamma_RHS(state.Gamma, state.gtilde, state.Atilde,
                    state.K, state.phi, state.alpha, h)
    # Gauge
    dalpha = dt_alpha_1pluslog(state.alpha, state.K)
    dbeta = dt_beta_gamma_driver(dG, state.beta, eta=1.0)
    # Z4c damping (pointwise average H, M)
    Nx, Ny, Nz = state.phi.shape
    dK_z4c = np.zeros_like(state.K)
    dG_z4c = np.zeros_like(state.Gamma)
    for ix in range(Nx):
        for iy in range(Ny):
            for iz in range(Nz):
                b = BSSN(float(state.phi[ix,iy,iz]), state.gtilde[ix,iy,iz],
                         float(state.K[ix,iy,iz]), state.Atilde[ix,iy,iz],
                         np.zeros(3))
                H = hamiltonian_constraint(b, 0.0)
                M = momentum_constraint(b, np.zeros(3))
                dK_z4c[ix,iy,iz] = damping_dK(state.alpha[ix,iy,iz], H)
                dG_z4c[ix,iy,iz] = damping_dGamma(M)
    return {
        "dphi": r["dphi"], "dgt": r["dgt"],
        "dK": r["dK"] + dK_z4c, "dAt": r["dAt"],
        "dGamma": dG + dG_z4c, "dalpha": dalpha, "dbeta": dbeta,
    }


def rk4_step_v2(state: BSSNState, h, dt):
    """One RK4 step on the extended state."""
    def add(s, k, scale):
        return BSSNState(
            phi=s.phi + scale*k["dphi"],
            gtilde=s.gtilde + scale*k["dgt"],
            K=s.K + scale*k["dK"],
            Atilde=s.Atilde + scale*k["dAt"],
            Gamma=s.Gamma + scale*k["dGamma"],
            alpha=s.alpha + scale*k["dalpha"],
            beta=s.beta + scale*k["dbeta"],
        )
    k1 = _rhs(state, h)
    k2 = _rhs(add(state, k1, 0.5*dt), h)
    k3 = _rhs(add(state, k2, 0.5*dt), h)
    k4 = _rhs(add(state, k3, dt), h)
    return BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"] + 2*k2["dphi"] + 2*k3["dphi"] + k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"] + 2*k2["dgt"] + 2*k3["dgt"] + k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"] + 2*k2["dK"] + 2*k3["dK"] + k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"] + 2*k2["dAt"] + 2*k3["dAt"] + k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"] + 2*k2["dGamma"] + 2*k3["dGamma"] + k4["dGamma"]),
        alpha=state.alpha + (dt/6)*(k1["dalpha"] + 2*k2["dalpha"] + 2*k3["dalpha"] + k4["dalpha"]),
        beta=state.beta + (dt/6)*(k1["dbeta"] + 2*k2["dbeta"] + 2*k3["dbeta"] + k4["dbeta"]),
    )


__all__ = ["rk4_step_v2"]
