"""Unit tests for the moving fine-patch window logic (zt005.amr.amr_moving)."""
import numpy as np

from zt005.amr.amr_moving import com_drift_cells, fine_grid_mask, target_lo

N_C, N_SUB = 17, 9
H_C = 20.0 / (N_C - 1)
CENTER_IDX = (N_C - N_SUB) // 2


def test_com_drift_cells():
    assert com_drift_cells(0.0, H_C) == 0.0
    assert np.isclose(com_drift_cells(1.25, H_C), 1.0)


def test_target_lo_centered_zero():
    assert target_lo(0.0, H_C, N_C, N_SUB) == CENTER_IDX


def test_target_lo_even_alignment():
    # even coarse-cell shifts keep the factor-2 prolongation aligned;
    # odd 1-cell shifts snap to the even multiple (floor by 2)
    for com_cells in (0.5, 1.0, 1.5, 2.0, 3.0):
        lo = target_lo(com_cells * H_C, H_C, N_C, N_SUB)
        shift = lo - CENTER_IDX
        assert shift % 2 == 0
        assert shift in (0, 2, 4, -2, -4)


def test_target_lo_clipped():
    lo = target_lo(50.0 * H_C, H_C, N_C, N_SUB)
    assert lo == N_C - N_SUB
    lo = target_lo(-50.0 * H_C, H_C, N_C, N_SUB)
    assert lo == 0


def test_fine_grid_mask_shape_and_threshold():
    n_f, L_f = 17, 10.0
    mask = fine_grid_mask(n_f, L_f, 1.5)
    assert mask.shape == (n_f, n_f, n_f)
    assert mask[0, 0, 0]          # corner is far
    assert mask[n_f // 2, n_f // 2, n_f // 2] == 0   # center excluded
    # window-local r threshold: corner distance of the window is
    # sqrt(3) * (L_f/2) = 8.66 > 1.5
    assert mask[-1, -1, -1]
