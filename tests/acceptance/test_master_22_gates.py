"""Official Acceptance Test: All 22 Scientific & Causal Pipeline Gates."""

import json
from pathlib import Path
import pytest
from zt005.run_all_gates import main as run_all_gates_main


def test_master_22_acceptance_gates():
    """Assert that every single gate in the official 22-gate table evaluates to PASS."""
    from zt005.run_all_gates import OUT
    import sys

    # Run quick mode to verify all gates pass
    sys.argv = ["run_all_gates", "--mode", "quick"]
    run_all_gates_main()

    assert OUT.exists()
    data = json.loads(OUT.read_text(encoding="utf-8"))

    assert data["all_pass"] is True
    assert data["n_gates"] == 22
    assert data["n_passed"] == 22

    # Check that each gate individually passed
    for gate_name, gate_info in data["gates"].items():
        assert gate_info["pass"] is True, f"Gate '{gate_name}' failed: {gate_info}"
