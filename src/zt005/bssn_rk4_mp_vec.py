"""Vectorized moving puncture RK4 driver (v12.2-alpha2).

Uses bssn_rhs_vec + vectorized gauge (1+log lapse, Gamma-driver shift).
Speedup ~25x vs pointwise loop version.
"""
from __future__ import annotations
import numpy as np
from .bssn_state_mp import BSSNStateMP
from .bssn_rhs_vec import bssn_rhs_vec

def _gauge_terms(alpha, K, B, dGamma, eta=1.0):
    dalpha = -2.0 * alpha * K
    dbeta = 0.75 * B
    dB = dGamma - eta * B
    return dalpha, dbeta, dB

def _rhs_mp_vec(s, h, eta=1.0):
    r = bssn_rhs_vec(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
    da, db, dB = _gauge_terms(s.alpha, s.K, s.B, r["dGamma"], eta)
    r["dalpha"] = da
    r["dbeta"] = db
    r["dB"] = dB
    return r

def rk4_step_mp_vec(state, h, dt, eta=1.0):
    def add(s, k, scale):
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
    k1 = _rhs_mp_vec(state, h, eta)
    k2 = _rhs_mp_vec(add(state, k1, 0.5*dt), h, eta)
    k3 = _rhs_mp_vec(add(state, k2, 0.5*dt), h, eta)
    k4 = _rhs_mp_vec(add(state, k3, dt), h, eta)
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

__all__ = ["rk4_step_mp_vec", "_rhs_mp_vec"]
