"""RK4 driver for Brans-Dicke gravity (v8.0-alpha2).

Evolves BSSN metric + BD scalar phi with conjugate momentum Pi.
G_eff = G / phi modifies the effective gravitational constant.
"""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_full import bssn_rhs_beta1
from .matter_brans_dicke import brans_dicke_sources
from .bssn_grid import d_central


def _lap(f, h):
    grad = np.stack(d_central(f, h), axis=-1)
    L = np.zeros_like(f)
    for ax in range(3):
        dd = d_central(grad[..., ax], h)
        L += dd[ax]
    return L


def _bd_field_rhs(phi_BD, Pi_BD, alpha, K, h, omega=40000.0):
    """BD scalar field equation:
    Box phi = (1/(3 + 2 omega)) [8 pi T_matter - dV/dphi]
    In weak field / slow motion, reduce to wave equation:
    d_t Pi = alpha K Pi + alpha Lap(phi)
    """
    return alpha * Pi_BD, alpha * K * Pi_BD + alpha * _lap(phi_BD, h)


def rk4_step_brans_dicke(state, phi_BD, Pi_BD, h, dt, omega=40000.0):
    """One RK4 step on (metric + BD scalar)."""
    def rhs_all(s, ph, Pph):
        rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
        src = brans_dicke_sources(ph, Pph, s.gtilde, s.phi, h, omega=omega)
        e4phi = np.exp(4.0 * s.phi)
        gamma = e4phi[..., None, None] * s.gtilde
        ginv = np.linalg.inv(gamma)
        rhs["dK"] = rhs["dK"] + 4.0 * np.pi * s.alpha * (src["rho"] + src["S"])
        STF = src["S_ij"] - (1.0/3.0) * src["S"][..., None, None] * gamma
        rhs["dAt"] = rhs["dAt"] + 8.0 * np.pi * s.alpha[..., None, None] * \
                     np.exp(-4.0*s.phi)[..., None, None] * STF
        rhs["dGamma"] = rhs["dGamma"] - 16.0 * np.pi * s.alpha[..., None] * \
                        np.einsum('...ij,...j->...i', ginv, src["S_i"])
        m = _bd_field_rhs(ph, Pph, s.alpha, s.K, h, omega)
        return rhs, m

    def add_metric(s, k, scale):
        return BSSNState(
            phi=s.phi + scale*k["dphi"],
            gtilde=s.gtilde + scale*k["dgt"],
            K=s.K + scale*k["dK"],
            Atilde=s.Atilde + scale*k["dAt"],
            Gamma=s.Gamma + scale*k["dGamma"],
            alpha=s.alpha.copy(), beta=s.beta.copy(),
        )

    k1, m1 = rhs_all(state, phi_BD, Pi_BD)
    s2 = add_metric(state, k1, 0.5*dt)
    ph2 = phi_BD + 0.5*dt*m1[0]; Pph2 = Pi_BD + 0.5*dt*m1[1]

    k2, m2 = rhs_all(s2, ph2, Pph2)
    s3 = add_metric(state, k2, 0.5*dt)
    ph3 = phi_BD + 0.5*dt*m2[0]; Pph3 = Pi_BD + 0.5*dt*m2[1]

    k3, m3 = rhs_all(s3, ph3, Pph3)
    s4 = add_metric(state, k3, dt)
    ph4 = phi_BD + dt*m3[0]; Pph4 = Pi_BD + dt*m3[1]

    k4, m4 = rhs_all(s4, ph4, Pph4)

    new_state = BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha.copy(), beta=state.beta.copy(),
    )
    new_phi_BD = phi_BD + (dt/6)*(m1[0]+2*m2[0]+2*m3[0]+m4[0])
    new_Pi_BD = Pi_BD + (dt/6)*(m1[1]+2*m2[1]+2*m3[1]+m4[1])
    return new_state, new_phi_BD, new_Pi_BD


__all__ = ["rk4_step_brans_dicke"]
