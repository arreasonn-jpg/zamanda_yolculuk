
from __future__ import annotations
import numpy as np

def source_sample_points_on_quadrature(cyl, nr,nphi,nz):
    from .cyl_quadrature import cylindrical_midpoints
    return cylindrical_midpoints(cyl.radius_m,cyl.height_m,nr,nphi,nz)

def min_distance_to_sources(points, sources):
    p=np.asarray(points,float)
    d=np.full(len(p),np.inf)
    for src in sources:
        wire=np.asarray(src[0],float)
        # Dense point approximation; segment-aware distance is expensive.
        dd=np.min(np.linalg.norm(p[:,None,:]-wire[None,:,:],axis=2),axis=1)
        d=np.minimum(d,dd)
    return d

def regularization_report(points,sources,cutoffs=(0.002,0.005,0.01,0.02)):
    d=min_distance_to_sources(points,sources)
    out=[]
    for c in cutoffs:
        keep=d>=c
        out.append({
            "cutoff_m":float(c),
            "kept_fraction":float(keep.mean()),
            "excluded_fraction":float((~keep).mean()),
            "min_distance_m":float(d.min()),
            "p01_distance_m":float(np.quantile(d,0.01)),
            "median_distance_m":float(np.median(d))
        })
    return out

def masked_mean(values,keep):
    v=np.asarray(values)
    k=np.asarray(keep)
    if not np.any(k):
        return float("nan")
    return float(np.mean(v[k]))
