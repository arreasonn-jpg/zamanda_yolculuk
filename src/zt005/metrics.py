import numpy as np

def rel_l2(a,b,mask=None):
    if mask is not None:
        a=a[mask]; b=b[mask]
    den=np.linalg.norm(b.ravel())
    return float(np.linalg.norm((a-b).ravel())/max(den,1e-30))

def rms(a,mask=None):
    if mask is not None: a=a[mask]
    return float(np.sqrt(np.mean(a*a)))

def max_abs(a,mask=None):
    if mask is not None: a=a[mask]
    return float(np.max(np.abs(a)))
