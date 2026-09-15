
import numpy as np
from zt005.frozen_source_gr import relative_metric_error

def test_metric_error_zero():
    x=np.zeros((3,4,4)); r=relative_metric_error(x,x)
    assert r["mean"]==0 and r["max"]==0

def test_metric_error_positive():
    ref=np.ones((3,4,4)); x=2*ref
    r=relative_metric_error(x,ref)
    assert r["mean"]>0
