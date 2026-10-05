
"""ZT-008.1 tests."""
import numpy as np
from zt005.run_zt008_1 import (
    component_list, check_fit, total_mass, category_mass,
    ENVELOPE_CM,
)


def test_components_nonempty():
    c = component_list()
    assert len(c) >= 8


def test_all_components_have_required_fields():
    c = component_list()
    for comp in c:
        assert "name" in comp
        assert "dims_cm" in comp
        assert "mass_kg" in comp
        assert "position_frac" in comp
        assert "category" in comp


def test_all_components_fit_envelope():
    c = component_list()
    r = check_fit(c)
    fails = [x for x in r if not x["fits"]]
    assert len(fails) == 0, f"Components exceeding envelope: {fails}"


def test_total_mass_reasonable():
    c = component_list()
    m = total_mass(c)
    assert 3.0 < m < 12.0


def test_categories_include_required():
    c = component_list()
    cats = category_mass(c)
    for required in ["field", "power", "control", "sensor", "safety"]:
        assert required in cats


def test_positions_within_unit_cube():
    c = component_list()
    for comp in c:
        x, y, z = comp["position_frac"]
        assert 0.0 <= x <= 1.0
        assert 0.0 <= y <= 1.0
        assert 0.0 <= z <= 1.0
