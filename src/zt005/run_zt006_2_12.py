
import json,argparse
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .finite_wire_study import fields_at_points,field_stats
from .source_conditioning import cancellation_metrics
from .maxwell import biot_savart,vector_potential
from .geometry import build_sources
from .pointwise_curl import curl_point_4
from .validated_stress_energy import em_stress_energy

OUT=Path("zt006_2_12_results.json")

def sample_points(cyl,n=160,seed=909):
    rng=np.random.default_rng(seed); pts=[]
    while len(pts)<n:
        q=rng.uniform([-0.8*cyl.radius_m,-0.8*cyl.height_m/2,-0.8*cyl.radius_m],
                      [0.8*cyl.radius_m,0.8*cyl.height_m/2,0.8*cyl.radius_m])
        if q[0]**2+q[2]**2 <= (0.8*cyl.radius_m)**2: pts.append(q)
    return np.asarray(pts)

def baseline(name,pts,b,d,h=0.005):
    sources=build_sources(name,b); omega=2*np.pi*d.frequency_hz
    def A_fn(p): return vector_potential(np.atleast_2d(p),sources,d.peak_current_a)[0]
    A=np.vstack([A_fn(x) for x in pts])
    B=np.vstack([curl_point_4(A_fn,x,h) for x in pts])
    E=omega*A
    T=em_stress_energy(E,B)["T"]
    return E,B,T

def rel_l2(a,b):
    return float(np.linalg.norm((a-b).ravel())/max(np.linalg.norm(b.ravel()),1e-300))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    a=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    pts=sample_points(cyl,120 if a.mode=="quick" else 250)
    radii=[0.,.001,.0025,.005] if a.mode=="quick" else [0.,.001,.0025,.005,.0075]
    cross=4 if a.mode=="quick" else 8
    print("="*104)
    print("ZAMANDA YOLCULUK — ZT-006.2.12")
    print("CURRENT-CONSERVING FINITE-CROSS-SECTION CONDUCTOR MODEL")
    print("="*104)
    print(f"Points={len(pts)} radii={radii} cross={cross}")
    results={}
    for name in ["G1","G2","G3","G4"]:
        E0,B0,T0=baseline(name,pts,b,d)
        rows=[]; print(f"\n[{name}]")
        for r in radii:
            E,B,T,sources,Ieach=fields_at_points(name,pts,b,d,r,cross,.005)
            st=field_stats(E,B)
            comps=np.asarray([biot_savart(pts,[s],Ieach) for s in sources])
            cm=cancellation_metrics(comps,comps.sum(axis=0))
            row={"radius_m":r,"source_count":len(sources),
                 "current_per_filament_A":Ieach,**st,
                 "T00_mean":float(T[:,0,0].mean()),
                 "T0i_mean_norm":float(np.linalg.norm(T[:,0,1:],axis=1).mean()),
                 "E_vs_r0_relL2":rel_l2(E,E0),
                 "B_vs_r0_relL2":rel_l2(B,B0),
                 "T_vs_r0_relL2":rel_l2(T,T0),
                 "kappa_mean":cm["kappa_mean"],"kappa_p95":cm["kappa_p95"],"kappa_max":cm["kappa_max"]}
            rows.append(row)
            print(f"  r={r:.4f} N={len(sources):3d} Ieach={Ieach:.4f} "
                  f"B={row['B_mean']:.6e} E={row['E_mean']:.6e} T00={row['T00_mean']:.6e} "
                  f"BΔ={row['B_vs_r0_relL2']:.3e} TΔ={row['T_vs_r0_relL2']:.3e}")
        results[name]=rows
    OUT.write_text(json.dumps({"stage":"ZT-006.2.12","mode":a.mode,"results":results},indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: current-conserving finite-wire study; no GR/CTC conclusion.")
if __name__=="__main__": main()
