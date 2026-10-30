"""RK4 driver coupling complex scalar field to BSSN (v6.0-alpha2).

Complex scalar with U(1) symmetric potential:
  V = m^2 |psi|^2 + lambda |psi|^4

Field equations:
  d_t psi_R = alpha Pi_R
  d_t psi_I = alpha Pi_I
  d_t Pi_R = alpha K Pi_R + alpha Lap(psi_R) - alpha dV/dpsi_R
  d_t Pi_I = alpha K Pi_I + alpha Lap(psi_I) - alpha dV/dpsi_I
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_full import bssn_rhs_beta1
from .matter.complex_scalar import (
    complex_scalar_sources,
    complex_scalar_field_rhs,
)


def rk4_step_complex(state: BSSNState, psi_R, psi_I, Pi_R, Pi_I, h, dt,
                     m2=1.0, lam=1.0):
    """One RK4 step on (metric + complex scalar field)."""
    def metric_rhs(s, psi_R, psi_I, Pi_R, Pi_I):
        rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
        src = complex_scalar_sources(psi_R, psi_I, Pi_R, Pi_I,
                                      s.gtilde, s.phi, h, m2=m2, lam=lam)
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

    def field_rhs(s, psi_R, psi_I, Pi_R, Pi_I):
        return complex_scalar_field_rhs(psi_R, psi_I, Pi_R, Pi_I,
                                         s.alpha, s.K, h, m2=m2, lam=lam)

    def add_metric(s, k, scale):
        return BSSNState(
            phi=s.phi + scale*k["dphi"],
            gtilde=s.gtilde + scale*k["dgt"],
            K=s.K + scale*k["dK"],
            Atilde=s.Atilde + scale*k["dAt"],
            Gamma=s.Gamma + scale*k["dGamma"],
            alpha=s.alpha.copy(), beta=s.beta.copy(),
        )

    # k1
    k1 = metric_rhs(state, psi_R, psi_I, Pi_R, Pi_I)
    m1 = field_rhs(state, psi_R, psi_I, Pi_R, Pi_I)
    s2 = add_metric(state, k1, 0.5*dt)
    psiR2 = psi_R + 0.5*dt*m1[0]; psiI2 = psi_I + 0.5*dt*m1[1]
    PiR2 = Pi_R + 0.5*dt*m1[2]; PiI2 = Pi_I + 0.5*dt*m1[3]

    # k2
    k2 = metric_rhs(s2, psiR2, psiI2, PiR2, PiI2)
    m2_ = field_rhs(s2, psiR2, psiI2, PiR2, PiI2)
    s3 = add_metric(state, k2, 0.5*dt)
    psiR3 = psi_R + 0.5*dt*m2_[0]; psiI3 = psi_I + 0.5*dt*m2_[1]
    PiR3 = Pi_R + 0.5*dt*m2_[2]; PiI3 = Pi_I + 0.5*dt*m2_[3]

    # k3
    k3 = metric_rhs(s3, psiR3, psiI3, PiR3, PiI3)
    m3 = field_rhs(s3, psiR3, psiI3, PiR3, PiI3)
    s4 = add_metric(state, k3, dt)
    psiR4 = psi_R + dt*m3[0]; psiI4 = psi_I + dt*m3[1]
    PiR4 = Pi_R + dt*m3[2]; PiI4 = Pi_I + dt*m3[3]

    # k4
    k4 = metric_rhs(s4, psiR4, psiI4, PiR4, PiI4)
    m4 = field_rhs(s4, psiR4, psiI4, PiR4, PiI4)

    new_state = BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha.copy(), beta=state.beta.copy(),
    )
    new_psiR = psi_R + (dt/6)*(m1[0]+2*m2_[0]+2*m3[0]+m4[0])
    new_psiI = psi_I + (dt/6)*(m1[1]+2*m2_[1]+2*m3[1]+m4[1])
    new_PiR = Pi_R + (dt/6)*(m1[2]+2*m2_[2]+2*m3[2]+m4[2])
    new_PiI = Pi_I + (dt/6)*(m1[3]+2*m2_[3]+2*m3[3]+m4[3])
    return new_state, new_psiR, new_psiI, new_PiR, new_PiI


__all__ = ["rk4_step_complex"]
