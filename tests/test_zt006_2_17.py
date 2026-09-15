import numpy as np
from zt005.metric_validation import metric_from_hbar, metric_invariants, timelike_check

def test_minkowski():
    g = metric_from_hbar(np.zeros((4,4)))
    assert np.allclose(g, np.diag([-1.,1.,1.,1.]))
    assert np.isclose(metric_invariants(g)["det"], -1.0)
    assert metric_invariants(g)["lorentzian_signature"]
    assert timelike_check(g)["timelike"]

def test_tiny_perturbation():
    h=np.zeros((4,4)); h[0,0]=1e-12; h[1,1]=2e-12
    g=metric_from_hbar(h)
    assert metric_invariants(g)["lorentzian_signature"]
    assert timelike_check(g)["timelike"]
