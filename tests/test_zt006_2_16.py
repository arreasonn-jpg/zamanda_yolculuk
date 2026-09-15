import numpy as np
from zt005.metric_validation import metric_invariants,timelike_check

def test_metric_identity():
    g=np.diag([-1.,1.,1.,1.]); inv=metric_invariants(g)
    assert np.isclose(inv['det'],-1.)

def test_rest_timelike():
    g=np.diag([-1.,1.,1.,1.]); r=timelike_check(g,[1,0,0,0])
    assert r['timelike'] and r['interval']<0
