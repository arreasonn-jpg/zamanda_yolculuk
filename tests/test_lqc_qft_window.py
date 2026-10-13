import numpy as np
import pytest
from zt005.lqc_qft_window import (
    lqc_correction_factor, ptum_rho_physical,
    qft_window_shift_due_to_lqc, RHO_PTUM_UNIT,
)
from zt005.lqc_scalar_evolution import RHO_C


def test_ptum_rho_physical_order():
    """PTUM density ~ 10^18-10^19 kg/m^3."""
    rho = ptum_rho_physical(1.8)
    assert 1e17 < rho < 1e21


def test_lqc_correction_negligible():
    """LQC correction for PTUM is essentially 1."""
    rho = ptum_rho_physical(1.8)
    corr = lqc_correction_factor(rho)
    assert corr == pytest.approx(1.0, abs=1e-10)


def test_lqc_shift_zero_machine_precision():
    """QFT window shift due to LQC is below machine precision."""
    result = qft_window_shift_due_to_lqc(1.8)
    assert result["shift"] < 1e-70


def test_ratio_scale():
    """rho_ptum / rho_c ~ 1e-77."""
    rho = ptum_rho_physical(1.8)
    ratio = rho / RHO_C
    assert 1e-80 < ratio < 1e-70


def test_a_to_rho_monotonic():
    """Higher A -> larger |rho_min|."""
    rho_14 = ptum_rho_physical(1.4)
    rho_20 = ptum_rho_physical(2.0)
    assert rho_20 > rho_14


def test_a_not_in_table_raises():
    with pytest.raises(ValueError):
        ptum_rho_physical(5.0)
