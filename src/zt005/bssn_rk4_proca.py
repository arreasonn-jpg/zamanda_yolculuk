"""RK4 driver: BSSN metric + Proca vector field (v10.2-alpha2).
State: A_1, A_2, A_3 (vector potential) + E_1, E_2, E_3 (momenta)
Gauge: A_0 = 0 (temporal gauge)
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_full import bssn_rhs_beta1
from .matter.matter_proca import proca_sources
from .bssn_grid import d_central

def _lap(f, h):
    grad = np.stack(d_central(f, h), axis=-1)
    L = np.zeros_like(f)
    for ax in range(3):
        L += d_central(grad[..., ax], h)[ax]
    return L

def _proca_rhs(A1, A2, A3, E1, E2, E3, alpha, K, h, m2=0.0):
    """Maxwell/Proca in temporal gauge:
    d_t A_i = alpha E_i
    d_t E_i = alpha (K E_i + Lap(A_i) - m2 A_i)
    """
    dA1 = alpha * E1; dA2 = alpha * E2; dA3 = alpha * E3
    dE1 = alpha * (K * E1 + _lap(A1, h) - m2 * A1)
    dE2 = alpha * (K * E2 + _lap(A2, h) - m2 * A2)
    dE3 = alpha * (K * E3 + _lap(A3, h) - m2 * A3)
    return dA1, dA2, dA3, dE1, dE2, dE3

def rk4_step_proca(state, A1, A2, A3, E1, E2, E3, h, dt, m2=0.0):
    def rhs_all(s, a1, a2, a3, e1, e2, e3):
        rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
        src = proca_sources(np.zeros_like(a1), a1, a2, a3,
                            e1, e2, e3, s.gtilde, s.phi, h, m2=m2)
        e4phi = np.exp(4.0 * s.phi)
        gamma = e4phi[..., None, None] * s.gtilde
        ginv = np.linalg.inv(gamma)
        rhs["dK"] = rhs["dK"] + 4.0 * np.pi * s.alpha * (src["rho"] + src["S"])
        STF = src["S_ij"] - (1.0/3.0) * src["S"][..., None, None] * gamma
        rhs["dAt"] = rhs["dAt"] + 8.0 * np.pi * s.alpha[..., None, None] * \
            np.exp(-4.0 * s.phi)[..., None, None] * STF
        rhs["dGamma"] = rhs["dGamma"] - 16.0 * np.pi * s.alpha[..., None] * \
            np.einsum('...ij,...j->...i', ginv, src["S_i"])
        m = _proca_rhs(a1, a2, a3, e1, e2, e3, s.alpha, s.K, h, m2=m2)
        return rhs, m
    def add_m(s, k, scale):
        return BSSNState(
            phi=s.phi + scale*k["dphi"],
            gtilde=s.gtilde + scale*k["dgt"],
            K=s.K + scale*k["dK"],
            Atilde=s.Atilde + scale*k["dAt"],
            Gamma=s.Gamma + scale*k["dGamma"],
            alpha=s.alpha.copy(), beta=s.beta.copy())
    k1, m1 = rhs_all(state, A1, A2, A3, E1, E2, E3)
    s2 = add_m(state, k1, 0.5*dt)
    a1_2 = A1+0.5*dt*m1[0]; a2_2 = A2+0.5*dt*m1[1]; a3_2 = A3+0.5*dt*m1[2]
    e1_2 = E1+0.5*dt*m1[3]; e2_2 = E2+0.5*dt*m1[4]; e3_2 = E3+0.5*dt*m1[5]
    k2, m2_ = rhs_all(s2, a1_2, a2_2, a3_2, e1_2, e2_2, e3_2)
    s3 = add_m(state, k2, 0.5*dt)
    a1_3 = A1+0.5*dt*m2_[0]; a2_3 = A2+0.5*dt*m2_[1]; a3_3 = A3+0.5*dt*m2_[2]
    e1_3 = E1+0.5*dt*m2_[3]; e2_3 = E2+0.5*dt*m2_[4]; e3_3 = E3+0.5*dt*m2_[5]
    k3, m3 = rhs_all(s3, a1_3, a2_3, a3_3, e1_3, e2_3, e3_3)
    s4 = add_m(state, k3, dt)
    a1_4 = A1+dt*m3[0]; a2_4 = A2+dt*m3[1]; a3_4 = A3+dt*m3[2]
    e1_4 = E1+dt*m3[3]; e2_4 = E2+dt*m3[4]; e3_4 = E3+dt*m3[5]
    k4, m4 = rhs_all(s4, a1_4, a2_4, a3_4, e1_4, e2_4, e3_4)
    new_state = BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha.copy(), beta=state.beta.copy())
    na1 = A1 + (dt/6)*(m1[0]+2*m2_[0]+2*m3[0]+m4[0])
    na2 = A2 + (dt/6)*(m1[1]+2*m2_[1]+2*m3[1]+m4[1])
    na3 = A3 + (dt/6)*(m1[2]+2*m2_[2]+2*m3[2]+m4[2])
    ne1 = E1 + (dt/6)*(m1[3]+2*m2_[3]+2*m3[3]+m4[3])
    ne2 = E2 + (dt/6)*(m1[4]+2*m2_[4]+2*m3[4]+m4[4])
    ne3 = E3 + (dt/6)*(m1[5]+2*m2_[5]+2*m3[5]+m4[5])
    return new_state, na1, na2, na3, ne1, ne2, ne3

__all__ = ["rk4_step_proca"]
