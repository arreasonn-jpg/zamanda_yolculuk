"""RK4 driver using full BSSN RHS (beta1)."""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState, minkowski_state
from .bssn_full import bssn_rhs_beta1
from .bssn_z4c import damping_dK, damping_dGamma
from .bssn_ricci_full import hamiltonian_constraint_full


def rk4_step_beta1(state: BSSNState, h, dt):
    def rhs(s):
        r = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
        # Z4c damping (uses real H grid)
        H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
        r["dK"] = r["dK"] + damping_dK(s.alpha, H, kappa1=1.0)
        return r
    def add(s, k, scale):
        return BSSNState(
            phi=s.phi + scale*k["dphi"],
            gtilde=s.gtilde + scale*k["dgt"],
            K=s.K + scale*k["dK"],
            Atilde=s.Atilde + scale*k["dAt"],
            Gamma=s.Gamma + scale*k["dGamma"],
            alpha=s.alpha.copy(),
            beta=s.beta.copy(),
        )
    k1 = rhs(state)
    k2 = rhs(add(state, k1, 0.5*dt))
    k3 = rhs(add(state, k2, 0.5*dt))
    k4 = rhs(add(state, k3, dt))
    return BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha.copy(),
        beta=state.beta.copy(),
    )


__all__ = ["rk4_step_beta1"]
