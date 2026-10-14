"""Direct GW initial data injection (v10.1-alpha1).
Inject a TT-gauge linearized GW packet propagating along z.
Tests the BSSN -> NP extraction pipeline directly, independent
of source physics.
"""
from __future__ import annotations
import numpy as np

def gw_packet(N, L, A_amp, sigma, z0, polarization="plus"):
    """Return BSSN-consistent metric perturbation.
    gamma_ij = delta_ij + h_ij (linearized)
    """
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    envelope = A_amp * np.exp(-((Z - z0)**2) / (2.0 * sigma**2))
    if polarization == "plus":
        e_ij = np.zeros((3, 3))
        e_ij[0, 0] = 1.0
        e_ij[1, 1] = -1.0
    elif polarization == "cross":
        e_ij = np.zeros((3, 3))
        e_ij[0, 1] = 1.0
        e_ij[1, 0] = 1.0
    else:
        raise ValueError(f"Unknown polarization: {polarization}")
    h = np.einsum('...ij,...->...ij',
                  np.broadcast_to(e_ij, (N, N, N, 3, 3)),
                  envelope)
    gamma = np.broadcast_to(np.eye(3), (N, N, N, 3, 3)).copy() + h
    detg = np.linalg.det(gamma)
    phi = (1.0/12.0) * np.log(np.maximum(detg, 1e-30))
    gtilde = np.exp(-4.0 * phi[..., None, None]) * gamma
    # Time derivative: for +z propagation, h_t = -h_z (speed c=1)
    dh_dz = np.gradient(h, x[1] - x[0], axis=2)
    K_ij = -0.5 * dh_dz
    K_tr = np.einsum('...ij,...ij->...', np.linalg.inv(gamma), K_ij)
    Atilde = np.exp(-4.0 * phi[..., None, None]) * (
        K_ij - (1.0/3.0) * gamma * K_tr[..., None, None])
    return {"gtilde": gtilde, "phi": phi, "Atilde": Atilde,
            "K": K_tr, "h_ij": h, "envelope": envelope}

__all__ = ["gw_packet"]
