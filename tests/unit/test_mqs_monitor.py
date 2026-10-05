"""Unit tests for MQS regime monitor and uncertainty budget engine."""

import pytest
from zt005.mqs_monitor import evaluate_mqs_validity
from zt005.uncertainty_budget import compute_uncertainty_budget, ErrorBreakdown


def test_mqs_regime_classification():
    # 1 MHz with 1m size -> ka ~ 0.021 -> MQS VALID
    v1 = evaluate_mqs_validity(1e6, 1.0)
    assert v1["status"] == "MQS VALID"
    assert v1["is_mqs_valid"] is True
    assert v1["ka"] < 0.1

    # 10 MHz with 2.5m size -> ka ~ 0.52 -> FULL-WAVE REQUIRED
    v2 = evaluate_mqs_validity(10e6, 2.5)
    assert v2["status"] == "FULL-WAVE REQUIRED"
    assert v2["is_mqs_valid"] is False


def test_error_breakdown_rss():
    eb = ErrorBreakdown(
        grid_error=3.0,
        quadrature_error=4.0,
        source_discretization_error=0.0,
        finite_conductor_error=0.0,
        floating_point_error=0.0,
        model_approximation_error=0.0,
    )
    # sqrt(3^2 + 4^2) = 5.0
    assert pytest.approx(eb.total_uncertainty) == 5.0


def test_compute_uncertainty_budget():
    budget = compute_uncertainty_budget(1.98e-46, h_prev_grid=1.975e-46, h_solver_b=1.985e-46)
    assert budget["total_uncertainty"] > 0
    assert budget["relative_uncertainty_percent"] < 1.0
    assert "h = (" in budget["formatted_result"]
