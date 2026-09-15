
import numpy as np
from zt005.wire_regularization import regularization_report

def test_regularization_report():
    pts=np.array([[0.,0.,0.],[1.,0.,0.],[2.,0.,0.]])
    src=[(np.array([[1.,0.,0.],[1.,0.,1.]]),np.zeros((2,3)))]
    r=regularization_report(pts,src,(.5,1.))
    assert len(r)==2
    assert all(0<=x["excluded_fraction"]<=1 for x in r)

def test_cutoff_ordering():
    pts=np.array([[0.,0.,0.],[1.,0.,0.]])
    src=[(np.array([[1.,0.,0.],[1.,0.,1.]]),np.zeros((2,3)))]
    r=regularization_report(pts,src,(.01,.1,.5))
    assert r[0]["excluded_fraction"]<=r[1]["excluded_fraction"]<=r[2]["excluded_fraction"]
