
from __future__ import annotations
import numpy as np
from .cyl_quadrature import cylindrical_midpoints, quadrature_integral_many
from .spatial_gr import inverse_trace_reverse

def frozen_green_metric(obs_points, source_points, weights, T, softening_m=0.0025):
    hb = quadrature_integral_many(obs_points, source_points, T, weights, softening_m)
    return np.stack([inverse_trace_reverse(x) for x in hb])

def relative_metric_error(metric, reference):
    a=np.asarray(metric); b=np.asarray(reference)
    d=np.linalg.norm((a-b).reshape(len(a),-1),axis=1)
    r=np.linalg.norm(b.reshape(len(b),-1),axis=1)
    rel=d/np.maximum(r,1e-300)
    return {
        "mean":float(np.mean(rel)),
        "median":float(np.median(rel)),
        "p95":float(np.quantile(rel,0.95)),
        "max":float(np.max(rel)),
    }
