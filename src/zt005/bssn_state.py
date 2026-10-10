"""Extended BSSN state for v2.0: includes Gamma, alpha, beta.

v1.x state = {phi, gtilde, K, Atilde}
v2.0 state = above + {Gamma, alpha, beta}
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class BSSNState:
    phi: np.ndarray       # (N,N,N)
    gtilde: np.ndarray    # (N,N,N,3,3), det = 1
    K: np.ndarray         # (N,N,N)
    Atilde: np.ndarray    # (N,N,N,3,3), traceless
    Gamma: np.ndarray     # (N,N,N,3)  -- conformal connection
    alpha: np.ndarray     # (N,N,N)    -- lapse
    beta: np.ndarray      # (N,N,N,3)  -- shift

    def copy(self) -> "BSSNState":
        return BSSNState(
            phi=self.phi.copy(), gtilde=self.gtilde.copy(),
            K=self.K.copy(), Atilde=self.Atilde.copy(),
            Gamma=self.Gamma.copy(), alpha=self.alpha.copy(),
            beta=self.beta.copy(),
        )

    def shape(self) -> tuple:
        return self.phi.shape


def minkowski_state(N: int) -> BSSNState:
    """Flat space: alpha=1, beta=0, phi=0, gtilde=I, K=0, Atilde=0, Gamma=0."""
    return BSSNState(
        phi=np.zeros((N,N,N)),
        gtilde=np.broadcast_to(np.eye(3)[None,None,None,:,:], (N,N,N,3,3)).copy(),
        K=np.zeros((N,N,N)),
        Atilde=np.zeros((N,N,N,3,3)),
        Gamma=np.zeros((N,N,N,3)),
        alpha=np.ones((N,N,N)),
        beta=np.zeros((N,N,N,3)),
    )


__all__ = ["BSSNState", "minkowski_state"]
