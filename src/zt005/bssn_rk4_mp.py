"""RK4 driver with moving-puncture gauge (v3.0-alpha2)."""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_state_mp import BSSNStateMP
from .bssn_rhs_grid import bssn_rhs_grid
from .bssn_shift import shift_terms
from .bssn_gamma import dGamma_RHS
from .bssn_gauge_mp import dt_alpha_mp, dt_beta_mp, dt_B_mp
from .bssn_z4c import damping_dK
from .bssn_ricci_full import hamiltonian_constraint_full


def _rhs_mp(s: BSSNStateMP, h):
    r = bssn_rhs_grid(s.phi, s.gtilde, s.K, s.Atilde, h)
    L = shift_terms(s.phi, s.gtilde, s.K, s.Atilde, s.Gamma, s.beta, h)
    r["dphi"] = r["dphi"] + L["Lphi"]
    r["dgt"]  = r["dgt"]  + L["Lgt"]
    r["dK"]   = r["dK"]   + L["LK"]
    r["dAt"]  = r["dAt"]  + L["LAt"]
    r["dGamma"] = dGamma_RHS(s.Gamma, s.gtilde, s.Atilde, s.K,
                              s.phi, s.alpha, h) + L["LGamma"]

    H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
    r["dK"] = r["dK"] + damping_dK(s.alpha, H, kappa1=1.0)

    r["dalpha"] = dt_alpha_mp(s.alpha, s.K)
    r["dbeta"] = dt_beta_mp(s.B)
    r["dB"] = dt_B_mp(r["dGamma"], s.B, eta=1.0)
    return r


def _add(s: BSSNStateMP, k, scale):
    return BSSNStateMP(
        phi=s.phi + scale*k["dphi"],
        gtilde=s.gtilde + scale*k["dgt"],
        K=s.K + scale*k["dK"],
        Atilde=s.Atilde + scale*k["dAt"],
        Gamma=s.Gamma + scale*k["dGamma"],
        alpha=s.alpha + scale*k["dalpha"],
        beta=s.beta + scale*k["dbeta"],
        B=s.B + scale*k["dB"],
    )


def rk4_step_mp(state: BSSNStateMP, h, dt):
    k1 = _rhs_mp(state, h)
    k2 = _rhs_mp(_add(state, k1, 0.5*dt), h)
    k3 = _rhs_mp(_add(state, k2, 0.5*dt), h)
    k4 = _rhs_mp(_add(state, k3, dt), h)
    return BSSNStateMP(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha + (dt/6)*(k1["dalpha"]+2*k2["dalpha"]+2*k3["dalpha"]+k4["dalpha"]),
        beta=state.beta + (dt/6)*(k1["dbeta"]+2*k2["dbeta"]+2*k3["dbeta"]+k4["dbeta"]),
        B=state.B + (dt/6)*(k1["dB"]+2*k2["dB"]+2*k3["dB"]+k4["dB"]),
    )


__all__ = ["rk4_step_mp"]
