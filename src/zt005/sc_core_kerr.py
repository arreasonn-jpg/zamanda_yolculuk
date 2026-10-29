"""Rotating SC regular core (de Sitter + Israel) IC (v3.0-alpha6).

Combines the SC regular core (sc_core) with frame-dragging shift:

  Interior (r < R0): de Sitter f = 1 - (r/ell)^2
  Exterior (r > R0): Schwarzschild f = 1 - Rs/r
  Matching at R0 = ell

  beta^phi ~ 2 a M_eff / r^3   (leading-order slow rotation)
  where M_eff(r) = M * (1 - f(r)) * r / (2M) for r > R0
  and M_eff(r) = M * (r/ell)^2 for r < R0 (de Sitter interior)

Only valid for |a| <= 0.5 M (slow-rotation regime).
"""
from __future__ import annotations
import numpy as np
from .bssn_state_mp import BSSNStateMP
from .sc_core import f_SC, core_parameters, schwarzschild_radius


def sc_core_kerr_state(N, L, M=1.0, rho0=1e-3, a=0.1, eps=1e-6):
    """SC regular core (de Sitter + Israel) with frame-dragging shift."""
    if abs(a) > 0.5 * M:
        raise ValueError("a must satisfy |a| <= M/2 for slow rotation")
    
    x = np.linspace(-L/2, L/2, N)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    R = np.sqrt(X**2 + Y**2 + Z**2) + eps
    
    params = core_parameters(M, rho0)
    Rs = params["Rs"]
    R0 = params["R0"]
    ell = params["ell"]
    
    # Metric function f(r)
    f = f_SC(R, M, R0, ell)
    
    # Conformal factor and lapse
    psi = np.where(f > 0, f**(-0.25), 1.0)
    phi = np.log(psi)
    alpha = np.sqrt(np.maximum(f, 1e-10))
    
    # Effective mass for frame-dragging
    # For r > R0: M_eff = M (Schwarzschild)
    # For r < R0: M_eff = M * (r/ell)^2 (de Sitter interior mass)
    M_eff = np.where(R < R0, M * (R / ell)**2, M)
    
    # Frame-dragging shift: beta^phi ~ 2 a M_eff / r^3
    mag = 2.0 * a * M_eff / (R**3)
    rho_xy = np.sqrt(X**2 + Y**2) + eps
    bx = -Y / rho_xy * mag
    by =  X / rho_xy * mag
    bz = np.zeros_like(bx)
    
    return BSSNStateMP(
        phi=phi,
        gtilde=np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy(),
        K=np.zeros((N,N,N)),
        Atilde=np.zeros((N,N,N,3,3)),
        Gamma=np.zeros((N,N,N,3)),
        alpha=alpha,
        beta=np.stack([bx, by, bz], axis=-1),
        B=np.zeros((N,N,N,3)),
    ), params


__all__ = ["sc_core_kerr_state"]