
from __future__ import annotations
import numpy as np
from .maxwell import vector_potential
from .finite_wire_sources import build_finite_sources_grouped
from .pointwise_curl import curl_point_4
from .validated_stress_energy import em_stress_energy

def fields_with_groups(name,points,backpack,drive,radius_m,cross_samples=6,h=0.005):
    groups=build_finite_sources_grouped(name,backpack,radius_m,cross_samples)
    sources=[s for g in groups for s in g]
    currents=[drive.peak_current_a/len(g) for g in groups for _ in g]
    omega=2*np.pi*drive.frequency_hz

    def A_fn(p):
        p=np.atleast_2d(p)
        A=np.zeros((len(p),3))
        for src,I in zip(sources,currents):
            A += vector_potential(p,[src],I)
        return A[0]

    p=np.asarray(points,float)
    A=np.vstack([A_fn(x) for x in p])
    B=np.vstack([curl_point_4(A_fn,x,h) for x in p])
    E=omega*A
    T=em_stress_energy(E,B)["T"]
    return E,B,T,sources,currents,groups

def fields_at_points(name,points,backpack,drive,radius_m,cross_samples=6,h=0.005):
    E,B,T,sources,currents,groups=fields_with_groups(
        name,points,backpack,drive,radius_m,cross_samples,h
    )
    # Backward-compatible return: one representative current per subfilament.
    return E,B,T,sources,(drive.peak_current_a/max(len(groups[0]),1))

def field_stats(E,B):
    en=np.linalg.norm(E,axis=1); bn=np.linalg.norm(B,axis=1)
    return {"E_mean":float(en.mean()),"E_p95":float(np.quantile(en,0.95)),
            "B_mean":float(bn.mean()),"B_p95":float(np.quantile(bn,0.95))}
