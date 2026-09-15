
import numpy as np

C=299792458.0
G=6.67430e-11
ETA=np.diag([1.,-1.,-1.,-1.])

def inverse_trace_reverse(hbar):
    tr=np.einsum("ab,ab",ETA,hbar)
    return hbar-0.5*ETA*tr

def spatial_hbar_from_source(obs_points, source_points, T_source, volume_m3, softening_m=0.0):
    obs=np.asarray(obs_points,float)
    src=np.asarray(source_points,float)
    T=np.asarray(T_source,float)
    w=volume_m3/len(src)
    d=np.linalg.norm(obs[:,None,:]-src[None,:,:],axis=2)
    if softening_m>0:
        d=np.sqrt(d*d+softening_m**2)
    d=np.maximum(d,1e-12)
    out=np.zeros((len(obs),4,4),float)
    # trace reverse in lower-index (+---) convention
    tr=np.einsum("ab,nab->n",ETA,T)
    trrev=T-0.5*ETA[None,:,:]*tr[:,None,None]
    for i in range(len(obs)):
        weights=w/d[i]
        out[i]=(4*G/C**4)*np.einsum("n,nab->ab",weights,trrev)
    return out

def proper_time_proxy(h00):
    return 0.5*np.asarray(h00,float)
