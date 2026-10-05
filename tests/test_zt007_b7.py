
"""ZT-007B.7 tests - final unified bound."""
import numpy as np
from zt005.run_zt007_b7 import (
    mechanism_table, closest_mechanism,
)


def test_six_mechanisms():
    mechs = mechanism_table()
    assert len(mechs) == 6


def test_all_have_required_fields():
    mechs = mechanism_table()
    for m in mechs:
        for key in ["name", "key", "param_name", "param_device",
                    "param_required", "gap", "rho_device", "rho_required",
                    "rho_gap", "exotic", "physically_allowed", "reference"]:
            assert key in m, f"{m['key']} missing {key}"


def test_all_gaps_positive():
    mechs = mechanism_table()
    for m in mechs:
        assert m["gap"] > 0
        assert m["rho_gap"] > 0


def test_gap_range_wide():
    mechs = mechanism_table()
    log_gaps = [np.log10(m["gap"]) for m in mechs]
    assert max(log_gaps) - min(log_gaps) > 5  # at least 5 decades spread


def test_tipler_not_physically_allowed():
    mechs = mechanism_table()
    tipler = [m for m in mechs if m["key"] == "Tipler"][0]
    assert tipler["physically_allowed"] is False


def test_closest_allowed_excludes_tipler():
    mechs = mechanism_table()
    allowed, _ = closest_mechanism(mechs)
    assert allowed[0] != "Tipler"


def test_wormhole_and_alcubierre_are_exotic():
    mechs = mechanism_table()
    for key in ["Wormhole", "Alcubierre"]:
        m = [x for x in mechs if x["key"] == key][0]
        assert m["exotic"] is True
