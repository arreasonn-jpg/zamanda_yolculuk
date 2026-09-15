
from __future__ import annotations
import numpy as np
from .geometry import build_sources
from .maxwell import vector_potential, biot_savart
from .validated_stress_energy import em_stress_energy
from .pointwise_curl import curl_point_4

def source_components(name, points, backpack, drive):
    sources=build_sources(name,backpack)
    A=[]; B=[]
    for src in sources:
        A.append(vector_potential(points,[src],drive.peak_current_a))
        B.append(biot_savart(points,[src],drive.peak_current_a))
    return sources,np.asarray(A),np.asarray(B)

def cancellation_metrics(components,total,eps=1e-300):
    cn=np.linalg.norm(components,axis=-1)
    tn=np.linalg.norm(total,axis=-1)
    ratio=cn.sum(axis=0)/np.maximum(tn,eps)
    finite=ratio[np.isfinite(ratio)]
    return {
        "kappa_mean":float(np.mean(finite)) if len(finite) else float("inf"),
        "kappa_p95":float(np.quantile(finite,0.95)) if len(finite) else float("inf"),
        "kappa_max":float(np.max(ratio)) if len(ratio) else float("inf"),
        "component_sum_norm_mean":float(cn.sum(axis=0).mean()),
        "total_norm_mean":float(tn.mean()),
        "near_zero_total_fraction":float(np.mean(tn < 1e-12*np.maximum(cn.sum(axis=0),eps))),
    }

def component_energy_metrics(name,points,backpack,drive,h=0.005):
    sources=build_sources(name,backpack)
    omega=2*np.pi*drive.frequency_hz
    rows=[]
    for i,src in enumerate(sources):
        def A_fn(p):
            return vector_potential(np.atleast_2d(p),[src],drive.peak_current_a)[0]
        A=np.vstack([A_fn(p) for p in points])
        B=np.vstack([curl_point_4(A_fn,p,h) for p in points])
        E=omega*A
        T=em_stress_energy(E,B)["T"]
        rows.append({
            "source_index":i,
            "B_mean":float(np.linalg.norm(B,axis=1).mean()),
            "E_mean":float(np.linalg.norm(E,axis=1).mean()),
            "T00_mean":float(T[:,0,0].mean()),
            "T0i_mean_norm":float(np.linalg.norm(T[:,0,1:],axis=1).mean()),
        })
    return rows
