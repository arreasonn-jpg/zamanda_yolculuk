from __future__ import annotations
import numpy as np
from .gauss_cyl_quadrature import gauss_cylindrical_nodes
from .batch_em_source import compute_T_batch
from .cyl_quadrature import EINSTEIN_GREEN_FACTOR
from .spatial_gr import inverse_trace_reverse

def external_observers():
    return np.array([[0.,0.,2.],[0.,0.,3.],[0.,0.,5.],[2.,0.,0.],[3.,0.,0.]])

def compute_external_metric(name, observers, cyl, backpack, drive, grid):
    pts,w=gauss_cylindrical_nodes(cyl.radius_m,cyl.height_m,*grid)
    T=compute_T_batch(name,pts,backpack,drive,h=0.005)
    out=[]
    for obs in observers:
        d=np.linalg.norm(pts-np.asarray(obs)[None,:],axis=1)
        hbar=EINSTEIN_GREEN_FACTOR*np.sum((w/d)[:,None,None]*T,axis=0)
        out.append(inverse_trace_reverse(hbar))
    return np.asarray(out),len(pts)

def component_errors(a,b):
    return {
        "h00_abs":float(abs(a[0,0]-b[0,0])),
        "h00_rel":float(abs(a[0,0]-b[0,0])/max(abs(b[0,0]),1e-300)),
        "h0i_abs":float(np.linalg.norm(a[0,1:]-b[0,1:])),
        "tensor_rel":float(np.linalg.norm((a-b).ravel())/max(np.linalg.norm(b.ravel()),1e-300))
    }
