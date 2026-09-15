
import numpy as np
from zt005.finite_wire_convergence import relative_l2

def test_relative_zero():
    x=np.ones((3,4,4))
    assert relative_l2(x,x)==0.0

def test_relative_positive():
    x=np.ones((3,4,4)); y=2*x
    assert relative_l2(x,y)>0
