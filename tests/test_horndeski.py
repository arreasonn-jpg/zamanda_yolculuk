import numpy as np
import pytest
from zt005.horndeski import (
    L2_kessence, L4_nonminimal,
    horndeski_GR_limit, horndeski_BD_limit,
    ptum_qft_window_in_horndeski,
)


def test_L2_canonical():
    """Canonical kinetic: L2 = X."""
    assert L2_kessence(1.0, 0.5) == 0.5


def test_L4_GR_limit():
    """GR: G4 = 1/(16 pi G)."""
    R = 1.0
    L = L4_nonminimal(R, 1.0)
    assert L == pytest.approx(1.0 / (16.0 * np.pi * 6.6743e-11), rel=1e-4)


def test_GR_limit_structure():
    lim = horndeski_GR_limit()
    assert lim["G3"] == 0
    assert lim["G5"] == 0


def test_BD_limit_structure():
    lim = horndeski_BD_limit(40000.0)
    assert "phi" in lim["G4"]
    assert "40000" in lim["G2"]


def test_qft_window_invariant():
    """PTUM QFT window is G4-independent."""
    w1 = ptum_qft_window_in_horndeski(1.0)
    w2 = ptum_qft_window_in_horndeski(5.0)
    assert w1["A_min"] == w2["A_min"] == 1.35
    assert w1["A_max"] == w2["A_max"] == 1.87
    assert w1["invariant"] is True


def test_multiple_G4_same_window():
    """Across G4 in [0.1, 10]: window identical."""
    windows = [ptum_qft_window_in_horndeski(g) for g in [0.1, 1.0, 10.0]]
    assert all(w["A_min"] == 1.35 for w in windows)
    assert all(w["A_max"] == 1.87 for w in windows)
