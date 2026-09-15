
import numpy as np
from zt005.robust_convergence import robust_compare

def test_robust_compare_zero():
    x=np.ones((5,4,4))
    r=robust_compare(x,x)
    assert r["mean"]==0 and r["max"]==0

def test_robust_compare_positive():
    x=np.ones((5,3)); y=2*np.ones((5,3))
    r=robust_compare(x,y)
    assert r["mean"]>0 and r["p95"]>0
