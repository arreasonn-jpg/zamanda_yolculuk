
import numpy as np
from zt005.local_conditioning import gradient_norm

def test_linear_gradient():
    p=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[0,0,1.],[1,1,1.],[2,0,0],[0,2,0]],float)
    v=np.column_stack([2*p[:,0],-3*p[:,1],p[:,2]])
    g=gradient_norm(v,p)
    assert np.nanmedian(g)>3.0
