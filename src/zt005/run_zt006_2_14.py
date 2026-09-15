
import argparse,json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .finite_wire_convergence import convergence_rows

OUT=Path("zt006_2_14_results.json")

def sample_points(cyl,n=120,seed=111):
    rng=np.random.default_rng(seed); pts=[]
    while len(pts)<n:
        q=rng.uniform([-0.8*cyl.radius_m,-0.8*cyl.height_m/2,-0.8*cyl.radius_m],
                      [0.8*cyl.radius_m,0.8*cyl.height_m/2,0.8*cyl.radius_m])
        if q[0]**2+q[2]**2 <= (0.8*cyl.radius_m)**2:
            pts.append(q)
    return np.asarray(pts)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    a=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    pts=sample_points(cyl,80 if a.mode=="quick" else 180)
    radii=[.001,.0025,.005] if a.mode=="quick" else [.001,.0025,.005,.0075]
    counts=[2,4,6,8] if a.mode=="quick" else [2,4,6,8,12]

    print("="*104)
    print("ZAMANDA YOLCULUK — ZT-006.2.14")
    print("FINITE-WIRE CROSS-SECTION CONVERGENCE")
    print("="*104)
    print(f"Validation points: {len(pts)}")
    print(f"Radii: {radii}")
    print(f"Cross-section sample counts: {counts}")

    results={}
    for name in ["G1","G2","G3","G4"]:
        print(f"\\n[{name}]")
        results[name]={}
        for r in radii:
            rows=convergence_rows(name,pts,b,d,r,counts,.005)
            results[name][str(r)]=rows
            print(f"  radius={r:.4f} m")
            for row in rows:
                print(f"    n={row['cross_samples']:2d} Eerr={row['E_rel_vs_ref']:.3e} "
                      f"Berr={row['B_rel_vs_ref']:.3e} Terr={row['T_rel_vs_ref']:.3e}")

    OUT.write_text(json.dumps({"stage":"ZT-006.2.14","mode":a.mode,"results":results},indent=2),encoding="utf-8")
    print(f"\\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: finite-wire source-distribution convergence; no GR/CTC conclusion.")

if __name__=="__main__":
    main()
