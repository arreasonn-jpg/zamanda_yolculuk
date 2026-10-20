"""RK4 driver coupling PTUM fields (Phi, psi) to BSSN (v5.0-alpha2)."""
from __future__ import annotations
import numpy as np
from .bssn_state import BSSNState
from .bssn_full import bssn_rhs_beta1
from .matter.matter_ptum import ptum_sources
from .matter.ptum_field import V_fuel, V_proxy, m_eff_squared
from .bssn_grid import d_central


def _ptum_field_rhs(Phi, Pi_Phi, psi, Pi_psi, alpha, K, h,
                    V0=1.0, lam=0.05, m2=1.0, alpha_p=1.0, Eth=1.0,
                    lam_psi=1.0):
    """Field equations:
      d_t Phi = alpha Pi_Phi
      d_t Pi_Phi = alpha K Pi_Phi + alpha Lap(Phi) - alpha dV_fuel/dPhi
                    - alpha dV_proxy/dPhi
      d_t psi = alpha Pi_psi
      d_t Pi_psi = alpha K Pi_psi + alpha Lap(psi)
                    - alpha dV_proxy/dpsi
    (flat laplacian; shift terms omitted for v5.0-alpha2)
    """
    def lap(f):
        grad = np.stack(d_central(f, h), axis=-1)
        L = np.zeros_like(f)
        for ax in range(3):
            dd = d_central(grad[..., ax], h)
            L += dd[ax]
        return L

    # potentials and derivatives
    V_f = V_fuel(Phi, V0=V0, lam=lam)
    eps = 1e-5
    dV_f = (V_fuel(Phi + eps, V0=V0, lam=lam) - V_fuel(Phi - eps, V0=V0, lam=lam)) / (2*eps)
    # dV_proxy/dPhi = 1/2 psi^2 dm_eff^2/dPhi
    #   dm_eff^2/dPhi = -2 m^2 alpha Phi / Eth^2
    meff2 = m_eff_squared(Phi, m2=m2, alpha=alpha_p, Eth=Eth)
    dmeff2_dPhi = -2.0 * m2 * alpha_p * Phi / (Eth**2)
    dVp_dPhi = 0.5 * psi**2 * dmeff2_dPhi
    # dV_proxy/dpsi = meff2 psi + lam_psi psi^3
    dVp_dpsi = meff2 * psi + lam_psi * psi**3

    dPhi_t = alpha * Pi_Phi
    dPiPhi_t = alpha * K * Pi_Phi + alpha * lap(Phi) - alpha * (dV_f + dVp_dPhi)
    dpsi_t = alpha * Pi_psi
    dPipsi_t = alpha * K * Pi_psi + alpha * lap(psi) - alpha * dVp_dpsi
    return dPhi_t, dPiPhi_t, dpsi_t, dPipsi_t


def rk4_step_ptum(state: BSSNState, Phi, Pi_Phi, psi, Pi_psi, h, dt,
                  V0=1.0, lam=0.05, m2=1.0, alpha_p=1.0, Eth=1.0,
                  lam_psi=1.0):
    """One RK4 step on (metric + 2 PTUM fields)."""
    def metric_rhs(s, Phi, Pi_Phi, psi, Pi_psi):
        rhs = bssn_rhs_beta1(s.phi, s.gtilde, s.K, s.Atilde, s.alpha, h)
        src = ptum_sources(Phi, Pi_Phi, psi, Pi_psi, s.gtilde, s.phi, h,
                            V0=V0, lam=lam, m2=m2, alpha=alpha_p, Eth=Eth,
                            lam_psi=lam_psi)
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

    def field_rhs(s, Phi, Pi_Phi, psi, Pi_psi):
        return _ptum_field_rhs(Phi, Pi_Phi, psi, Pi_psi, s.alpha, s.K, h,
                                V0=V0, lam=lam, m2=m2, alpha_p=alpha_p,
                                Eth=Eth, lam_psi=lam_psi)

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
    k1 = metric_rhs(state, Phi, Pi_Phi, psi, Pi_psi)
    m1 = field_rhs(state, Phi, Pi_Phi, psi, Pi_psi)
    s2 = add_metric(state, k1, 0.5*dt)
    Phi2 = Phi + 0.5*dt*m1[0]; PiPhi2 = Pi_Phi + 0.5*dt*m1[1]
    psi2 = psi + 0.5*dt*m1[2]; Pipsi2 = Pi_psi + 0.5*dt*m1[3]

    # k2
    k2 = metric_rhs(s2, Phi2, PiPhi2, psi2, Pipsi2)
    m2_ = field_rhs(s2, Phi2, PiPhi2, psi2, Pipsi2)
    s3 = add_metric(state, k2, 0.5*dt)
    Phi3 = Phi + 0.5*dt*m2_[0]; PiPhi3 = Pi_Phi + 0.5*dt*m2_[1]
    psi3 = psi + 0.5*dt*m2_[2]; Pipsi3 = Pi_psi + 0.5*dt*m2_[3]

    # k3
    k3 = metric_rhs(s3, Phi3, PiPhi3, psi3, Pipsi3)
    m3 = field_rhs(s3, Phi3, PiPhi3, psi3, Pipsi3)
    s4 = add_metric(state, k3, dt)
    Phi4 = Phi + dt*m3[0]; PiPhi4 = Pi_Phi + dt*m3[1]
    psi4 = psi + dt*m3[2]; Pipsi4 = Pi_psi + dt*m3[3]

    # k4
    k4 = metric_rhs(s4, Phi4, PiPhi4, psi4, Pipsi4)
    m4 = field_rhs(s4, Phi4, PiPhi4, psi4, Pipsi4)

    new_state = BSSNState(
        phi=state.phi + (dt/6)*(k1["dphi"]+2*k2["dphi"]+2*k3["dphi"]+k4["dphi"]),
        gtilde=state.gtilde + (dt/6)*(k1["dgt"]+2*k2["dgt"]+2*k3["dgt"]+k4["dgt"]),
        K=state.K + (dt/6)*(k1["dK"]+2*k2["dK"]+2*k3["dK"]+k4["dK"]),
        Atilde=state.Atilde + (dt/6)*(k1["dAt"]+2*k2["dAt"]+2*k3["dAt"]+k4["dAt"]),
        Gamma=state.Gamma + (dt/6)*(k1["dGamma"]+2*k2["dGamma"]+2*k3["dGamma"]+k4["dGamma"]),
        alpha=state.alpha.copy(), beta=state.beta.copy(),
    )
    new_Phi = Phi + (dt/6)*(m1[0]+2*m2_[0]+2*m3[0]+m4[0])
    new_PiPhi = Pi_Phi + (dt/6)*(m1[1]+2*m2_[1]+2*m3[1]+m4[1])
    new_psi = psi + (dt/6)*(m1[2]+2*m2_[2]+2*m3[2]+m4[2])
    new_Pipsi = Pi_psi + (dt/6)*(m1[3]+2*m2_[3]+2*m3[3]+m4[3])
    return new_state, new_Phi, new_PiPhi, new_psi, new_Pipsi


__all__ = ["rk4_step_ptum"]
