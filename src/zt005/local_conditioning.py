
from __future__ import annotations
import numpy as np
from .geometry import build_sources
from .maxwell import vector_potential
from .pointwise_curl import curl_point_4

def field_values(name, points, backpack, drive, h):
    sources=build_sources(name,backpack)
    omega=2*np.pi*drive.frequency_hz
    def A_fn(p):
        return vector_potential(np.atleast_2d(p),sources,drive.peak_current_a)[0]
    p=np.asarray(points,float)
    A=np.vstack([A_fn(x) for x in p])
    B=np.vstack([curl_point_4(A_fn,x,h) for x in p])
    E=omega*A
    return E,B

def gradient_norm(samples, points):
    p=np.asarray(points,float); v=np.asarray(samples,float)
    out=[]
    for i in range(len(p)):
        d=p-p[i]
        r=np.linalg.norm(d,axis=1)
        idx=np.argsort(r)[1:7]
        if len(idx)<3:
            out.append(np.nan); continue
        M=d[idx]
        J=np.zeros((3,3))
        for c in range(3):
            g,*_=np.linalg.lstsq(M,v[idx,c]-v[i,c],rcond=None)
            J[c,:]=g
        out.append(np.linalg.norm(J,'fro'))
    return np.asarray(out)

def local_conditioning(E,B,points):
    ge=gradient_norm(E,points); gb=gradient_norm(B,points)
    em=np.linalg.norm(E,axis=1); bm=np.linalg.norm(B,axis=1)
    re=ge/np.maximum(em,1e-30); rb=gb/np.maximum(bm,1e-30)
    return {
        "gradE_mean":float(np.nanmean(ge)),
        "gradB_mean":float(np.nanmean(gb)),
        "gradE_over_E_mean":float(np.nanmean(re)),
        "gradB_over_B_mean":float(np.nanmean(rb)),
        "gradE_over_E_p95":float(np.nanquantile(re,0.95)),
        "gradB_over_B_p95":float(np.nanquantile(rb,0.95)),
    }
