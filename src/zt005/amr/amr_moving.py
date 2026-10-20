"""Moving fine-patch AMR window logic (v13).

Pure helpers for COM-tracked recentering of the fine window:
 - binary center of mass from |Atilde|^2 peaks (amr_bh_track),
 - even-offset window start so the factor-2 prolongation grid
   stays aligned,
 - fine-grid diagnosis mask (window-local coordinates).
"""
from __future__ import annotations
import numpy as np


def com_drift_cells(com, h_c):
    """Coarse-cell drift (float) of the binary COM from x=0."""
    return com / h_c


def target_lo(com, h_c, n_c, n_sub):
    """Even-offset window start snapped to the COM drift.

    Even offsets keep the factor-2 prolongation grid aligned with the
    coarse nodes; the start is clipped to the valid range.
    """
    center_idx = (n_c - n_sub) // 2
    shift_cells = int(round(com / h_c))
    lo = center_idx + 2 * (shift_cells // 2)
    return int(np.clip(lo, 0, n_c - n_sub))


def fine_grid_mask(n_f, L_f, r_min=1.5):
    """Window-local diagnosis mask: r > r_min on the fine grid."""
    x = np.linspace(-L_f / 2, L_f / 2, n_f)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    return np.sqrt(X ** 2 + Y ** 2 + Z ** 2) > r_min


__all__ = ["com_drift_cells", "target_lo", "fine_grid_mask"]
