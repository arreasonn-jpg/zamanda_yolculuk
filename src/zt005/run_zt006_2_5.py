
import argparse,json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder, Backpack, Drive
from .source_conditioning import source_components,cancellation_metrics,component_energy_metrics

OUT=Path("zt006_2_5_results.json")

def sample_points(cyl,n=500,seed=123):
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
    pts=sample_points(cyl,200 if args.mode=="quick" else 500)
    print("="*96)
    print("ZAMANDA YOLCULUK — ZT-006.2.5")
    print("SOURCE CANCELLATION / COMPONENT CONDITIONING STUDY")
    print("="*96)
    print(f"Validation points: {len(pts)}")
    results={}
    for name in ["G1","G2","G3","G4"]:
        sources,Acomp,Bcomp=source_components(name,pts,b,d)
        Atotal=Acomp.sum(axis=0); Btotal=Bcomp.sum(axis=0)
        m=cancellation_metrics(Bcomp,Btotal)
        comp=component_energy_metrics(name,pts,b,d)
        super_err=float(np.linalg.norm(Btotal-Bcomp.sum(axis=0))/max(np.linalg.norm(Btotal),1e-300))
        results[name]={"source_count":len(sources),"B_cancellation":m,
                       "component_energy":comp,"superposition_relative_error":super_err}
        print(f"\n{name} sources={len(sources)} kappa_mean={m['kappa_mean']:.6e} "
              f"kappa_p95={m['kappa_p95']:.6e} kappa_max={m['kappa_max']:.6e} "
              f"near_zero={m['near_zero_total_fraction']:.3%} "
              f"superposition={super_err:.6e}")
        for r in comp:
            print(f"  src={r['source_index']} B={r['B_mean']:.6e} E={r['E_mean']:.6e} "
                  f"T00={r['T00_mean']:.6e} T0i={r['T0i_mean_norm']:.6e}")
    OUT.write_text(json.dumps({"stage":"ZT-006.2.5","mode":args.mode,"results":results},indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: source conditioning only; no GR/CTC conclusion.")
if __name__=="__main__": main()
