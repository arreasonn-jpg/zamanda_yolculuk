
from __future__ import annotations
import numpy as np

def regularized_green_integral(obs, points, weights, T, green_factor, rho_m=0.02):
    obs=np.asarray(obs,float); p=np.asarray(points,float)
    w=np.asarray(weights,float); vals=np.asarray(T,float)
    d=np.linalg.norm(p-obs[None,:],axis=1)
    near=d<rho_m
    out=np.zeros((4,4),float)

    if np.any(near):
        wn=w[near]
        tloc=np.average(vals[near],axis=0,weights=wn)
        # Exact spherical-cell identity: ∫_{r<rho} dV/r = 2π rho².
        out += green_factor*(2*np.pi*rho_m**2)*tloc

    far=~near
    if np.any(far):
        df=np.maximum(d[far],1e-15)
        out += green_factor*np.sum((w[far]/df)[:,None,None]*vals[far],axis=0)
    return out

def exclusion_fraction(obs, points, rho_m):
    d=np.linalg.norm(np.asarray(points)-np.asarray(obs)[None,:],axis=1)
    return float(np.mean(d<rho_m))
