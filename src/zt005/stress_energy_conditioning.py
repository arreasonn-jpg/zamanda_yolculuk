
from __future__ import annotations
import numpy as np

def relative_l2(a,b):
    a=np.asarray(a); b=np.asarray(b)
    return float(np.linalg.norm((a-b).ravel())/max(np.linalg.norm(b.ravel()),1e-300))

def componentwise_errors(T,Tref):
    out={}
    for name,(sl1,sl2) in {
        "T00":((slice(None),0,0),None),
        "T0i":((slice(None),0,slice(1,4)),None),
        "Tij":((slice(None),slice(1,4),slice(1,4)),None)
    }.items():
        x=T[sl1]; y=Tref[sl1]
        out[name]={
            "absolute_rms":float(np.sqrt(np.mean((x-y)**2))),
            "reference_rms":float(np.sqrt(np.mean(y**2))),
            "relative_rms":float(np.sqrt(np.mean((x-y)**2))/max(np.sqrt(np.mean(y**2)),1e-300)),
            "reference_mean_abs":float(np.mean(np.abs(y))),
        }
    return out

def field_error_components(E,B,Eref,Bref):
    out={}
    for name,X,Y in [("E",E,Eref),("B",B,Bref)]:
        out[name]={
            "relative_l2":relative_l2(X,Y),
            "mean_norm":float(np.mean(np.linalg.norm(Y,axis=1))),
            "rms_abs_error":float(np.sqrt(np.mean((X-Y)**2))),
        }
    return out
