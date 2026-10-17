"""v12.7-fix: Z4c-damped vectorized moving puncture RK4."""
from __future__ import annotations
import numpy as np
from .bssn_state_mp import BSSNStateMP
from .bssn_rhs_vec import bssn_rhs_vec
from .bssn_bc_improved import sponge_3d

def _rhs(s, h, sponge, kappa, eta=1.0):
    r = bssn_rhs_vec(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h,
                     sponge=sponge, kappa=kappa, beta=s.beta, Gamma=s.Gamma)
    r["dalpha"] = -2.0 * s.alpha * s.K
    r["dbeta"]  = 0.75 * s.B
    r["dB"]     = r["dGamma"] - eta * s.B
    return r

def _add(s, k, sc):
    return BSSNStateMP(
        phi=s.phi + sc*k["dphi"], gtilde=s.gtilde + sc*k["dgt"],
        K=s.K + sc*k["dK"], Atilde=s.Atilde + sc*k["dAt"],
        Gamma=s.Gamma + sc*k["dGamma"],
        alpha=s.alpha + sc*k["dalpha"],
        beta=s.beta + sc*k["dbeta"],
        B=s.B + sc*k["dB"])

def rk4_step_mp_vec_z4c(state, h, dt, width_cells=3, kappa=0.02, eta=1.0):
    N = state.phi.shape[0]
    sponge = sponge_3d(N, width_cells)
    k1 = _rhs(state, h, sponge, kappa, eta)
    k2 = _rhs(_add(state, k1, 0.5*dt), h, sponge, kappa, eta)
    k3 = _rhs(_add(state, k2, 0.5*dt), h, sponge, kappa, eta)
    k4 = _rhs(_add(state, k3, dt), h, sponge, kappa, eta)
    return BSSNStateMP(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha + (dt/6)*(k1["dalpha"]+2*k2["dalpha"]+2*k3["dalpha"]+k4["dalpha"]),
        beta=state.beta + (dt/6)*(k1["dbeta"]+2*k2["dbeta"]+2*k3["dbeta"]+k4["dbeta"]),
        B=state.B + (dt/6)*(k1["dB"]+2*k2["dB"]+2*k3["dB"]+k4["dB"]))

__all__ = ["rk4_step_mp_vec_z4c"]
