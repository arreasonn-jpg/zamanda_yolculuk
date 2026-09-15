import numpy as np
from zt005.convergence_exclusion import nearest_wire_distance, exclusion_report

def test_nearest_wire_distance():
    line=np.array([[1.,0.,0.],[1.,0.,1.]])
    p=np.array([[0.,0.,0.],[1.,0.,.5]])
    d=nearest_wire_distance(p,[line])
    assert np.all(d >= 0)
    assert d[1] < d[0]

def test_exclusion_report():
    pts=np.array([[0.,0.,0.],[1.,0.,0.],[2.,0.,0.]])
    line=np.array([[1.,-1.,0.],[1.,1.,0.]])
    rows=exclusion_report(pts,[line],.5,(1.,2.))
    assert len(rows)==2
    assert rows[0]["threshold_m"] < rows[1]["threshold_m"]
    assert all(0<=r["excluded_fraction"]<=1 for r in rows)
