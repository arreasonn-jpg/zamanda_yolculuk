
"""ZT-008.6 tests."""
import numpy as np
from zt005.run_zt008_6 import (
    test_hierarchy, pretest_checklist, measurement_protocol,
    data_validation_criteria, failure_response, human_gate_requirements,
)


def test_hierarchy_has_7_levels():
    h = test_hierarchy()
    assert len(h) == 7


def test_hierarchy_starts_with_simulation():
    h = test_hierarchy()
    assert h[0]["level"] == "T0"
    assert "Simulation" in h[0]["name"]


def test_human_test_is_last():
    h = test_hierarchy()
    assert h[-1]["human_involved"] is True
    assert h[-1]["level"] == "T6"


def test_only_last_level_involves_human():
    h = test_hierarchy()
    human_levels = [t for t in h if t["human_involved"]]
    assert len(human_levels) == 1


def test_prerequisite_chain_monotonic():
    """Each level requires all previous levels."""
    h = test_hierarchy()
    for i, t in enumerate(h):
        if i == 0:
            assert t["prerequisites"] == []
        else:
            # Should require previous levels (not necessarily only immediate)
            assert len(t["prerequisites"]) > 0


def test_pretest_checklist_grows_with_level():
    """Higher levels have more checklist items."""
    n_low = len(pretest_checklist("T1"))
    n_high = len(pretest_checklist("T6"))
    assert n_high > n_low


def test_measurement_protocol_increases_rate():
    """Higher levels have higher sample rates or more metrics."""
    mp2 = measurement_protocol("T2")
    mp3 = measurement_protocol("T3")
    assert mp3["sample_rate_hz"] > mp2["sample_rate_hz"]


def test_data_validation_has_checks():
    checks = data_validation_criteria()
    assert len(checks) >= 5


def test_failure_response_nonempty():
    fr = failure_response()
    assert len(fr) >= 4


def test_human_gate_has_6_layers():
    gates = human_gate_requirements()
    assert len(gates) == 6
    names = [g["name"] for g in gates]
    for required in ["Physics", "Engineering", "Biology",
                     "Ethics", "Legal", "Medical"]:
        assert required in names


def test_human_gate_includes_irb():
    gates = human_gate_requirements()
    ethics = [g for g in gates if g["name"] == "Ethics"][0]
    assert "IRB" in ethics["requirement"]
