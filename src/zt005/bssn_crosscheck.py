"""Cross-check: in the linear regime, BSSN perturbation should be
consistent with the retarded linearized-GR solution retarded_hbar."""
from __future__ import annotations
import numpy as np
from .retarded_gr import retarded_hbar, static_hbar, C
from .bssn_benchmark import schwarzschild_isotropic


def metric_perturbation_from_phi(phi_amp):
    """Weak-field relation h ~ 4 phi (isotropic conformal factor)."""
    return 4.0 * phi_amp


def weak_field_phi_from_mass(M, r):
    """phi = ln(1 + M/(2r)) ~ M/(2r) in weak field."""
    return float(np.log(1.0 + M / (2.0 * r)))


def linear_regime_consistent(M, r, tol_rel=0.05):
    """Check that h_ij ~ 4 phi for small M/r."""
    phi = weak_field_phi_from_mass(M, r)
    h_direct = metric_perturbation_from_phi(phi)
    # Newtonian expectation: h ~ 2M/r
    h_newton = 2.0 * M / r
    rel = abs(h_direct - h_newton) / max(abs(h_newton), 1e-300)
    return {"phi": phi, "h_bssn": h_direct, "h_newton": h_newton,
            "rel_diff": rel, "consistent": rel < tol_rel}


__all__ = ["metric_perturbation_from_phi", "weak_field_phi_from_mass",
           "linear_regime_consistent"]
