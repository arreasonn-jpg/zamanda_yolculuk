import numpy as np
import pytest
from zt005.brans_dicke_observables import (
    gr_perihelion_precession, bd_perihelion_precession,
    bd_deviation, mercury_test, gw_speed_ratio,
)


def test_gr_mercury():
    """GR Mercury precession ~ 43 arcsec/century."""
    out = mercury_test(omega=1e10)  # effectively GR
    assert out["gr_arcsec_per_century"] == pytest.approx(43.0, abs=0.5)


def test_bd_deviation_scaling():
    """BD deviation ~ 3/(4 omega) for large omega."""
    eps_40000 = bd_deviation(40000.0)
    eps_80000 = bd_deviation(80000.0)
    assert abs(eps_40000 / eps_80000 - 2.0) < 0.01


def test_cassini_bound():
    """omega = 40000 -> deviation < 2e-5."""
    eps = bd_deviation(40000.0)
    assert abs(eps) < 2e-5


def test_gw_speed_gr_limit():
    """omega large -> c_GW/c -> 1."""
    r = gw_speed_ratio(1e10)
    assert r == pytest.approx(1.0, abs=1e-9)


def test_gw_speed_deviation():
    """omega=40000 -> c_GW/c ~ 1 - 2.5e-5."""
    r = gw_speed_ratio(40000.0)
    assert 0.9999 < r < 1.0


def test_mercury_gr_bd_close():
    """At Cassini bound, BD and GR Mercury precession differ by < 1e-3 arcsec."""
    out = mercury_test(omega=40000.0)
    assert abs(out["deviation_arcsec"]) < 1e-3
