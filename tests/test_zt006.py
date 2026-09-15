
import numpy as np
from zt005.linear_gr import trace_reversed_source,proper_time_fractional_shift_from_h00

def test_trace_reverse_shape():
    T=np.zeros((2,4,4)); T[:,0,0]=1.
    R=trace_reversed_source(T)
    assert R.shape==(2,4,4) and np.isfinite(R).all()

def test_proxy():
    assert np.isclose(proper_time_fractional_shift_from_h00(2e-6),1e-6)
