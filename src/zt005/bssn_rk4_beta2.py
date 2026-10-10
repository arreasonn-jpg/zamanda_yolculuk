"""RK4 driver with shift evolution (v2.0-beta2)."""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_full import bssn_rhs_beta1
from .bssn_shift import shift_terms
from .bssn_gauge import dt_alpha_1pluslog, dt_beta_gamma_driver
from .bssn_z4c import damping_dK
from .bssn_ricci_full import hamiltonian_constraint_full


def _rhs(s: BSSNState, h):
    r = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
    L = shift_terms(s.phi, s.gtilde, s.K, s.Atilde, s.Gamma, s.beta, h)
    r["dphi"] = r["dphi"] + L["Lphi"]
    r["dgt"]  = r["dgt"]  + L["Lgt"]
    r["dK"]   = r["dK"]   + L["LK"]
    r["dAt"]  = r["dAt"]  + L["LAt"]
    r["dGamma"] = r["dGamma"] + L["LGamma"]

    H = hamiltonian_constraint_full(s.phi, s.gtilde, s.K, s.Atilde, h)
    r["dK"] = r["dK"] + damping_dK(s.alpha, H, kappa1=1.0)

    r["dalpha"] = dt_alpha_1pluslog(s.alpha, s.K)
    r["dbeta"] = dt_beta_gamma_driver(r["dGamma"], s.beta, eta=1.0)
    return r


def rk4_step_beta2(state: BSSNState, h, dt):
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
    def avg(kk):
        return {key: (kk[0][key] + 2*kk[1][key] + 2*kk[2][key] + kk[3][key])/6.0
                for key in kk[0]}
    a = avg([k1, k2, k3, k4])
    return BSSNState(
        phi=state.phi + dt*a["dphi"],
        gtilde=state.gtilde + dt*a["dgt"],
        K=state.K + dt*a["dK"],
        Atilde=state.Atilde + dt*a["dAt"],
        Gamma=state.Gamma + dt*a["dGamma"],
        alpha=state.alpha + dt*a["dalpha"],
        beta=state.beta + dt*a["dbeta"],
    )


__all__ = ["rk4_step_beta2"]
