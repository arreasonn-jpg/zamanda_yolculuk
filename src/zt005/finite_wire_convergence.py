
from __future__ import annotations
import numpy as np
from .finite_wire_study import fields_with_groups

def relative_l2(a,b):
    a=np.asarray(a); b=np.asarray(b)
    return float(np.linalg.norm((a-b).ravel())/max(np.linalg.norm(b.ravel()),1e-300))

def convergence_rows(name,points,backpack,drive,radius_m,sample_counts,h=0.005):
    ref_n=max(sample_counts)
    Eref,Bref,Tref,*_=fields_with_groups(name,points,backpack,drive,radius_m,ref_n,h)
    rows=[]
    for n in sample_counts:
        E,B,T,*_=fields_with_groups(name,points,backpack,drive,radius_m,n,h)
        rows.append({
            "cross_samples":int(n),
            "E_rel_vs_ref":relative_l2(E,Eref),
            "B_rel_vs_ref":relative_l2(B,Bref),
            "T_rel_vs_ref":relative_l2(T,Tref),
            "T00_mean":float(T[:,0,0].mean()),
            "T0i_mean_norm":float(np.linalg.norm(T[:,0,1:],axis=1).mean())
        })
    return rows
