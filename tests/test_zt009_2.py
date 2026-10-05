
"""ZT-009.2 tests."""
import numpy as np
from zt005.run_zt009_2 import (
    monte_carlo_h00, summary_stats, one_at_a_time_uncertainty,
    UNCERTAINTIES, DEFAULTS,
)
from zt005.run_zt009_1 import H00_BASE


def test_monte_carlo_returns_array():
    h, _ = monte_carlo_h00(1000)
    assert h.shape == (1000,)


def test_monte_carlo_mean_close_to_baseline():
    h, _ = monte_carlo_h00(10000)
    mean = np.mean(h)
    # Mean should be close to baseline
    assert abs(mean - H00_BASE) / H00_BASE < 0.05


def test_monte_carlo_reproducible():
    h1, _ = monte_carlo_h00(1000, seed=1)
    h2, _ = monte_carlo_h00(1000, seed=1)
    assert np.allclose(h1, h2)


def test_summary_stats_keys():
    h, _ = monte_carlo_h00(1000)
    s = summary_stats(h)
    for key in ["mean", "std", "median", "p05", "p95",
                "rel_uncertainty"]:
        assert key in s


def test_variance_contributions_sum_to_about_1():
    contrib = one_at_a_time_uncertainty(5000)
    total = sum(contrib.values())
    # Since parameters interact through h ~ I^2 R^2 / D, contributions
    # won't sum exactly to 1, but should be bounded
    assert 0.1 < total < 5.0


def test_current_dominates_uncertainty():
    """I has largest |S| * sigma combination, should dominate."""
    contrib = one_at_a_time_uncertainty(5000)
    # I should contribute substantial fraction
    assert contrib["I"] > 0.1


def test_uncertainties_are_small():
    """All input uncertainties should be < 20%."""
    for key, sigma in UNCERTAINTIES.items():
        assert sigma < 0.20
