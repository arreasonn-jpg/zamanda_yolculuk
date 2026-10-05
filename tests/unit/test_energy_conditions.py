"""Unit tests for stress-energy physical consistency and energy conditions."""

import numpy as np
import pytest
from zt005.energy_conditions import (
    verify_maxwell_energy_density,
    evaluate_energy_conditions,
    classify_conservation_gate,
)
from zt005.validated_stress_energy import em_stress_energy


def test_verify_maxwell_energy_density():
    E = np.array([[100.0, 0.0, 0.0], [0.0, 200.0, 0.0]])
    B = np.array([[0.0, 0.01, 0.0], [0.005, 0.0, 0.0]])
    se = em_stress_energy(E, B)
    report = verify_maxwell_energy_density(E, B, se["T"][:, 0, 0])
    assert report["pass"] is True
    assert report["all_positive_T00"] is True
    assert report["max_rel_err"] < 1e-12


def test_evaluate_energy_conditions():
    E = np.array([[100.0, 50.0, 0.0]])
    B = np.array([[0.0, 0.01, 0.005]])
    se = em_stress_energy(E, B)
    cond = evaluate_energy_conditions(se["T"])
    assert cond["overall_energy_conditions_pass"] is True
    assert cond["t00_positivity"]["pass"] is True
    assert cond["WEC"]["pass"] is True
    assert cond["NEC"]["pass"] is True
    assert cond["DEC"]["pass"] is True


def test_classify_conservation_gate():
    # Residuals below 0.05 -> PASS
    res_pass = classify_conservation_gate([0.01, 0.02, 0.01, 0.03], scale=1.0)
    assert res_pass["overall_gate"] == "PASS"

    # Residual between 0.05 and 0.20 -> WARN
    res_warn = classify_conservation_gate([0.01, 0.12, 0.01, 0.03], scale=1.0)
    assert res_warn["overall_gate"] == "WARN"

    # Residual above 0.20 -> FAIL
    res_fail = classify_conservation_gate([0.01, 0.25, 0.01, 0.03], scale=1.0)
    assert res_fail["overall_gate"] == "FAIL"
