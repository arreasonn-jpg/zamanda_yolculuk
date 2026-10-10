"""BSSNState with auxiliary B^i field for moving-puncture gauge."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .bssn_state import BSSNState


@dataclass
class BSSNStateMP:
    phi: np.ndarray
    gtilde: np.ndarray
    K: np.ndarray
    Atilde: np.ndarray
    Gamma: np.ndarray
    alpha: np.ndarray
    beta: np.ndarray
    B: np.ndarray          # auxiliary field (N,N,N,3)

    def to_base(self) -> BSSNState:
        return BSSNState(phi=self.phi, gtilde=self.gtilde, K=self.K,
                         Atilde=self.Atilde, Gamma=self.Gamma,
                         alpha=self.alpha, beta=self.beta)

    @classmethod
    def from_base(cls, s: BSSNState, B=None):
        if B is None:
            B = np.zeros(s.phi.shape + (3,))
        return cls(phi=s.phi, gtilde=s.gtilde, K=s.K, Atilde=s.Atilde,
                   Gamma=s.Gamma, alpha=s.alpha, beta=s.beta, B=B)


def minkowski_mp_state(N):
    from .bssn_state import minkowski_state
    return BSSNStateMP.from_base(minkowski_state(N))


__all__ = ["BSSNStateMP", "minkowski_mp_state"]
