
"""ZT-009.4 tests."""
import numpy as np
from pathlib import Path
from zt005.run_zt009_4 import (
    fig1_unified_mechanisms, fig2_em_scaling, fig3_best_case_tiers,
    fig4_monte_carlo, fig5_sensitivity_tornado, fig6_engineering_envelope,
    FIGDIR,
)


def _ensure_dir():
    FIGDIR.mkdir(exist_ok=True)


def test_fig1_generates():
    _ensure_dir()
    p = fig1_unified_mechanisms(FIGDIR)
    assert p is not None
    assert p.exists()
    assert p.stat().st_size > 1000


def test_fig2_generates():
    _ensure_dir()
    p = fig2_em_scaling(FIGDIR)
    assert p is not None
    assert p.exists()


def test_fig3_generates():
    _ensure_dir()
    p = fig3_best_case_tiers(FIGDIR)
    assert p is not None
    assert p.exists()


def test_fig4_generates():
    _ensure_dir()
    p = fig4_monte_carlo(FIGDIR)
    assert p is not None
    assert p.exists()


def test_fig5_generates():
    _ensure_dir()
    p = fig5_sensitivity_tornado(FIGDIR)
    assert p is not None
    assert p.exists()


def test_fig6_generates():
    _ensure_dir()
    p = fig6_engineering_envelope(FIGDIR)
    assert p is not None
    assert p.exists()
