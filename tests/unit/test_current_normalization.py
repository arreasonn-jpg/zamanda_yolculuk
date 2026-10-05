"""Unit tests for current normalization and classification."""

import pytest
from zt005.current_normalization import (
    CurrentSpec,
    make_current_spec,
    get_geometry_turns,
    multi_resolution_current_audit,
)
from zt005.model import Backpack


def test_current_spec_definitions():
    spec = CurrentSpec(i_peak=100.0, n_turns=6)
    assert spec.i_parent == 100.0
    assert spec.i_parent_peak == 100.0
    assert pytest.approx(spec.i_parent_rms, rel=1e-12) == 100.0 / (2.0 ** 0.5)
    assert spec.i_total_peak == 600.0
    assert pytest.approx(spec.i_total_rms, rel=1e-12) == 600.0 / (2.0 ** 0.5)
    assert spec.i_peak == 100.0
    assert pytest.approx(spec.i_rms, rel=1e-12) == 100.0 / (2.0 ** 0.5)

    d = spec.as_dict()
    assert d["I_parent"] == 100.0
    assert d["I_total"] == 600.0
    assert d["I_peak"] == 100.0
    assert pytest.approx(d["I_rms"], rel=1e-12) == 100.0 / (2.0 ** 0.5)


def test_geometry_turn_counts():
    assert get_geometry_turns("G1") == 6
    assert get_geometry_turns("G2") == 4
    assert get_geometry_turns("G3") == 2
    assert get_geometry_turns("G4") == 3


def test_multi_resolution_current_conservation():
    bp = Backpack()
    for geom in ["G1", "G2", "G3", "G4"]:
        audit = multi_resolution_current_audit(
            geom, bp, radius_m=0.0025, rules=[(1, 6), (2, 8), (3, 10)]
        )
        assert audit["all_pass"] is True
        assert audit["max_resolution_discrepancy"] < 1e-12
        for ev in audit["evaluations"]:
            assert ev["pass"] is True
            assert ev["current_residual"] < 1e-12
