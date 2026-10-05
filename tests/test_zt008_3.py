
"""ZT-008.3 tests."""
import numpy as np
from zt005.run_zt008_3 import (
    fault_modes, hardware_cutoff_chain, state_machine,
    check_safety_layers, worst_case_response_time, redundancy_analysis,
)


def test_fault_modes_nonempty():
    fm = fault_modes()
    assert len(fm) >= 6


def test_all_faults_have_id_and_severity():
    for fm in fault_modes():
        assert "id" in fm
        assert "severity" in fm
        assert fm["severity"] in ["low", "medium", "high"]


def test_all_faults_have_mitigation():
    for fm in fault_modes():
        assert "mitigation" in fm
        assert len(fm["mitigation"]) > 0


def test_cutoff_chain_has_multiple():
    chain = hardware_cutoff_chain()
    assert len(chain) >= 4


def test_cutoff_response_times_positive():
    for c in hardware_cutoff_chain():
        assert c["response_ms"] > 0


def test_redundancy_sufficient():
    red = redundancy_analysis()
    assert red["redundant"] >= 3


def test_state_machine_has_off_state():
    sm = state_machine()
    assert "OFF" in sm


def test_state_machine_transitions_valid():
    sm = state_machine()
    all_states = set(sm.keys())
    for state, trans in sm.items():
        for action, target in trans.items():
            assert target in all_states, f"{state}->{action}->{target} invalid"


def test_safety_layers_no_issues():
    issues = check_safety_layers()
    assert len(issues) == 0


def test_worst_response_bounded():
    wc = worst_case_response_time()
    assert 1 < wc < 5000
