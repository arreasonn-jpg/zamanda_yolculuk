
"""ZT-007B.4 tests - unified comparison."""
import numpy as np
from zt005.run_zt007_b4 import (
    collect_mechanism_data, find_closest_mechanism,
)


def test_all_4_mechanisms_present():
    d = collect_mechanism_data()
    assert set(d.keys()) == {"EM", "Kerr", "Tipler", "Wormhole"}


def test_each_mechanism_has_gap_field():
    d = collect_mechanism_data()
    for key, m in d.items():
        assert "gap" in m
        assert m["gap"] > 0


def test_all_gaps_are_astronomical():
    """All 4 mechanisms should be at least 10^6 above device scale."""
    d = collect_mechanism_data()
    for key, m in d.items():
        assert m["gap"] > 1e6, f"{key} gap = {m['gap']}"


def test_only_wormhole_uses_exotic():
    d = collect_mechanism_data()
    exotic = [k for k, m in d.items() if m["uses_exotic"]]
    assert exotic == ["Wormhole"]


def test_tipler_not_physically_allowed():
    d = collect_mechanism_data()
    assert d["Tipler"]["physically_allowed"] is False


def test_closest_mechanism_is_allowed():
    d = collect_mechanism_data()
    best_allowed, best_overall = find_closest_mechanism(d)
    assert d[best_allowed[0]]["physically_allowed"] is True


def test_em_gap_exceeds_expected():
    """EM current gap is ~10^20; energy-density gap is ~10^40."""
    d = collect_mechanism_data()
    assert d["EM"]["gap"] > 1e19
    assert d["EM"]["rho_gap"] > 1e40


def test_wormhole_gap_exceeds_1e35():
    d = collect_mechanism_data()
    assert d["Wormhole"]["gap"] > 1e35
