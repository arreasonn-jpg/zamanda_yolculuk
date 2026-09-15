from __future__ import annotations
import numpy as np
from .geometry import build_sources
from .maxwell import vector_potential
from .pointwise_curl import curl_point_4
from .validated_stress_energy import em_stress_energy
from .spatial_gr import spatial_hbar_from_source, inverse_trace_reverse

def make_source_points(cyl, n, seed=42):
    rng=np.random.default_rng(seed); pts=[]
    while len(pts)<n:
        q=rng.uniform([-0.95*cyl.radius_m,-0.95*cyl.height_m/2,-0.95*cyl.radius_m],[0.95*cyl.radius_m,0.95*cyl.height_m/2,0.95*cyl.radius_m])
        if q[0]**2+q[2]**2 <= (0.95*cyl.radius_m)**2: pts.append(q)
    return np.asarray(pts,float)

def compute_T(name,cyl,backpack,drive,points,h=0.005):
    sources=build_sources(name,backpack); omega=2*np.pi*drive.frequency_hz
    def A_fn(p): return vector_potential(np.atleast_2d(p),sources,drive.peak_current_a)[0]
    A=np.vstack([A_fn(p) for p in points]); B=np.vstack([curl_point_4(A_fn,p,h) for p in points]); E=omega*A
    return em_stress_energy(E,B)['T']

def metric_at_center(T,points,cyl,softening=0.0025):
    V=np.pi*cyl.radius_m**2*cyl.height_m
    hb=spatial_hbar_from_source(np.array([[0.,0.,0.]]),points,T,V,softening_m=softening)[0]
    return inverse_trace_reverse(hb)

def source_convergence(name,cyl,backpack,drive,counts=(200,400,800,1600),h=0.005):
    allp=make_source_points(cyl,max(counts)); prev=None; rows=[]
    for n in counts:
        p=allp[:n]; T=compute_T(name,cyl,backpack,drive,p,h); m=metric_at_center(T,p,cyl)
        h00=float(m[0,0]); h0=float(np.linalg.norm(m[0,1:])); ch=None if prev is None else abs(h00-prev)/max(abs(h00),1e-300)
        rows.append({'source_points':n,'h00':h00,'h0i_norm':h0,'relative_change_h00':ch}); prev=h00
    return rows

def current_scaling(name,cyl,backpack,drive,currents=(25,50,100,200,400),source_points=400,h=0.005):
    p=make_source_points(cyl,source_points); rows=[]
    for I in currents:
        d=type(drive)(frequency_hz=drive.frequency_hz,peak_current_a=float(I)); T=compute_T(name,cyl,backpack,d,p,h); m=metric_at_center(T,p,cyl)
        h00=float(m[0,0]); h0=float(np.linalg.norm(m[0,1:])); rows.append({'current_A':float(I),'h00':h00,'h0i_norm':h0,'h00_over_I2':h00/I**2,'h0i_over_I2':h0/I**2})
    I=np.array([r['current_A'] for r in rows]); a=np.abs(np.array([r['h00'] for r in rows])); b=np.abs(np.array([r['h0i_norm'] for r in rows]))
    return rows,float(np.polyfit(np.log(I),np.log(a),1)[0]),float(np.polyfit(np.log(I),np.log(b),1)[0])
