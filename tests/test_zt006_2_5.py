
import numpy as np
from zt005.source_conditioning import cancellation_metrics

def test_cancellation_no_cancel():
    c=np.array([[[1.,0,0]],[[2.,0,0]]])
    m=cancellation_metrics(c,c.sum(axis=0))
    assert np.isclose(m["kappa_mean"],1.)

def test_cancellation_detects_opposite_sources():
    c=np.array([[[1.,0,0]],[[-0.5,0,0]]])
    m=cancellation_metrics(c,c.sum(axis=0))
    assert m["kappa_mean"]>1.5
