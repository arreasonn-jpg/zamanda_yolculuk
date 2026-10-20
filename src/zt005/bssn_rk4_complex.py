"""RK4 driver coupling complex scalar to BSSN (v6.0-alpha2).

Two independent Klein-Gordon equations (U(1) symmetric):
  d_t Psi_R = alpha Pi_R
  d_t Pi_R  = alpha K Pi_R + alpha Lap(Psi_R) - alpha dV/dPsi_R
  (same for I)
with V = V(|Psi|^2), so dV/dPsi_R = 2 Psi_R V'(|Psi|^2).
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_full import bssn_rhs_beta1
from .matter.matter_complex import complex_sources
from .bssn_grid import d_central


def _lap(f, h):
    grad = np.stack(d_central(f, h), axis=-1)
    L = np.zeros_like(f)
    for ax in range(3):
        dd = d_central(grad[..., ax], h)
        L += dd[ax]
    return L


def _complex_field_rhs(Psi_R, Psi_I, Pi_R, Pi_I, alpha, K, h,
                       dV_func=None):
    dPsiR_t = alpha * Pi_R
    dPsiI_t = alpha * Pi_I
    lapR = _lap(Psi_R, h)
    lapI = _lap(Psi_I, h)
    if dV_func is None:
        dVR = np.zeros_like(Psi_R)
        dVI = np.zeros_like(Psi_I)
    else:
        dVR, dVI = dV_func(Psi_R, Psi_I)
    dPiR_t = alpha * K * Pi_R + alpha * lapR - alpha * dVR
    dPiI_t = alpha * K * Pi_I + alpha * lapI - alpha * dVI
    return dPsiR_t, dPsiI_t, dPiR_t, dPiI_t


def rk4_step_complex(state, Psi_R, Psi_I, Pi_R, Pi_I, h, dt,
                     dV_func=None):
    def metric_rhs(s, PR, PI, PiR, PiI):
        rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
        src = complex_sources(PR, PI, PiR, PiI, s.gtilde, s.phi, h)
        rho = src["rho"]; S_i = src["S_i"]; S_ij = src["S_ij"]; S = src["S"]
        e4phi = np.exp(4.0 * s.phi)
        gamma = e4phi[..., None, None] * s.gtilde
        ginv = np.linalg.inv(gamma)
        rhs["dK"] = rhs["dK"] + 4.0 * np.pi * s.alpha * (rho + S)
        STF = S_ij - (1.0/3.0) * S[..., None, None] * gamma
        rhs["dAt"] = rhs["dAt"] + 8.0 * np.pi * s.alpha[..., None, None] * \
                     np.exp(-4.0*s.phi)[..., None, None] * STF
        rhs["dGamma"] = rhs["dGamma"] - 16.0 * np.pi * s.alpha[..., None] * \
                        np.einsum('...ij,...j->...i', ginv, S_i)
        return rhs

    def add_m(s, k, scale):
        return BSSNState(
            phi=s.phi + scale*k["dphi"],
            gtilde=s.gtilde + scale*k["dgt"],
            K=s.K + scale*k["dK"],
            Atilde=s.Atilde + scale*k["dAt"],
            Gamma=s.Gamma + scale*k["dGamma"],
            alpha=s.alpha.copy(), beta=s.beta.copy(),
        )

    def field_rhs(s, PR, PI, PiR, PiI):
        return _complex_field_rhs(PR, PI, PiR, PiI, s.alpha, s.K, h, dV_func)

    k1 = metric_rhs(state, Psi_R, Psi_I, Pi_R, Pi_I)
    m1 = field_rhs(state, Psi_R, Psi_I, Pi_R, Pi_I)
    s2 = add_m(state, k1, 0.5*dt)
    PR2 = Psi_R + 0.5*dt*m1[0]; PI2 = Psi_I + 0.5*dt*m1[1]
    PiR2 = Pi_R + 0.5*dt*m1[2]; PiI2 = Pi_I + 0.5*dt*m1[3]

    k2 = metric_rhs(s2, PR2, PI2, PiR2, PiI2)
    m2 = field_rhs(s2, PR2, PI2, PiR2, PiI2)
    s3 = add_m(state, k2, 0.5*dt)
    PR3 = Psi_R + 0.5*dt*m2[0]; PI3 = Psi_I + 0.5*dt*m2[1]
    PiR3 = Pi_R + 0.5*dt*m2[2]; PiI3 = Pi_I + 0.5*dt*m2[3]

    k3 = metric_rhs(s3, PR3, PI3, PiR3, PiI3)
    m3 = field_rhs(s3, PR3, PI3, PiR3, PiI3)
    s4 = add_m(state, k3, dt)
    PR4 = Psi_R + dt*m3[0]; PI4 = Psi_I + dt*m3[1]
    PiR4 = Pi_R + dt*m3[2]; PiI4 = Pi_I + dt*m3[3]

    k4 = metric_rhs(s4, PR4, PI4, PiR4, PiI4)
    m4 = field_rhs(s4, PR4, PI4, PiR4, PiI4)

    new_state = BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha.copy(), beta=state.beta.copy(),
    )
    new_PR = Psi_R + (dt/6)*(m1[0]+2*m2[0]+2*m3[0]+m4[0])
    new_PI = Psi_I + (dt/6)*(m1[1]+2*m2[1]+2*m3[1]+m4[1])
    new_PiR = Pi_R + (dt/6)*(m1[2]+2*m2[2]+2*m3[2]+m4[2])
    new_PiI = Pi_I + (dt/6)*(m1[3]+2*m2[3]+2*m3[3]+m4[3])
    return new_state, new_PR, new_PI, new_PiR, new_PiI


__all__ = ["rk4_step_complex"]
