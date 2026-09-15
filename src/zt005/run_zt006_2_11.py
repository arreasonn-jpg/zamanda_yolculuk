
import argparse,json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .finite_wire_study import fields_at_points,field_stats
from .validated_stress_energy import em_stress_energy
from .source_conditioning import cancellation_metrics

OUT=Path("zt006_2_11_results.json")

def sample_points(cyl,n=160,seed=909):
    rng=np.random.default_rng(seed); pts=[]
    while len(pts)<n:
        q=rng.uniform([-0.8*cyl.radius_m,-0.8*cyl.height_m/2,-0.8*cyl.radius_m],
                      [0.8*cyl.radius_m,0.8*cyl.height_m/2,0.8*cyl.radius_m])
        if q[0]**2+q[2]**2 <= (0.8*cyl.radius_m)**2: pts.append(q)
    return np.asarray(pts)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    args=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    pts=sample_points(cyl,120 if args.mode=="quick" else 250)
    radii=[0.0,0.001,0.0025,0.005] if args.mode=="quick" else [0.0,0.001,0.0025,0.005,0.0075]
    cross=4 if args.mode=="quick" else 8

    print("="*102)
    print("ZAMANDA YOLCULUK — ZT-006.2.11")
    print("FINITE-CROSS-SECTION CONDUCTOR MODEL")
    print("="*102)
    print(f"Validation points: {len(pts)}")
    print(f"Conductor radii: {radii}")
    print(f"Cross-section filaments per source: {cross}")
    results={}
    for name in ["G1","G2","G3","G4"]:
        rows=[]; print(f"\n[{name}]")
        for r in radii:
            E,B,sources=fields_at_points(name,pts,b,d,r,cross_samples=cross,h=.005)
            T=em_stress_energy(E,B)["T"]; st=field_stats(E,B)
            from .maxwell import biot_savart
            current_per= d.peak_current_a/max(len(sources),1)
            comps=np.asarray([biot_savart(pts,[src],current_per) for src in sources])
            total=comps.sum(axis=0)
            cm=cancellation_metrics(comps,total)
            row={"radius_m":float(r),"source_count":len(sources),**st,
                 "T00_mean":float(T[:,0,0].mean()),
                 "T0i_mean_norm":float(np.linalg.norm(T[:,0,1:],axis=1).mean()),
                 "kappa_mean":cm["kappa_mean"],"kappa_p95":cm["kappa_p95"],"kappa_max":cm["kappa_max"]}
            rows.append(row)
            print(f"  r={r:.4f} m sources={len(sources):4d} B={row['B_mean']:.6e} "
                  f"E={row['E_mean']:.6e} T00={row['T00_mean']:.6e} kappa={row['kappa_mean']:.3f}")
        results[name]=rows
    OUT.write_text(json.dumps({"stage":"ZT-006.2.11","mode":args.mode,"results":results},indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: finite-cross-section conductor study; no GR/CTC conclusion.")

if __name__=="__main__": main()
