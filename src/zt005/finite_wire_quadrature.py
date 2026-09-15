
from __future__ import annotations
import numpy as np

def min_distance_to_filament(points, sources):
    p=np.asarray(points,float)
    d=np.full(len(p),np.inf)
    for src in sources:
        wire=np.asarray(src[0],float)
        d=np.minimum(d,np.min(np.linalg.norm(p[:,None,:]-wire[None,:,:],axis=2),axis=1))
    return d

def masked_volume_integral(obs, points, weights, T, keep, green_factor, softening_m):
    p=np.asarray(points,float)
    w=np.asarray(weights,float)
    vals=np.asarray(T,float)
    keep=np.asarray(keep,bool)
    p=p[keep]; w=w[keep]; vals=vals[keep]
    d=np.linalg.norm(p-np.asarray(obs,float)[None,:],axis=1)
    d=np.sqrt(d*d+softening_m**2)
    d=np.maximum(d,1e-15)
    raw=np.sum((w/d)[:,None,None]*vals,axis=0)
    return green_factor*raw
