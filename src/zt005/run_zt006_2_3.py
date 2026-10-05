
import argparse, json
from pathlib import Path
import numpy as np

from .model import ActiveCylinder, Backpack, Drive
from .error_budget import validation_points, T_at_fixed_points, relative_tensor_l2, kernel_resolution_with_frozen_T, softening_sensitivity
from .cyl_quadrature import cylindrical_midpoints

OUT=Path("results/exploratory/zt006_2_3_results.json")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    args=ap.parse_args()

    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    pts=validation_points()
    obs=np.array([[0.,0.,0.],[.5,0,0],[0,.5,0],[0,0,.5],[.9,0,0]],float)

    hvals=[0.02,0.01,0.005]
    grids=[(4,8,4),(6,12,6),(8,16,8)] if args.mode=="quick" else [(8,16,8),(12,24,12),(16,32,16)]
    soft=[0.0,0.001,0.0025,0.005,0.01]

    print("="*94)
    print("ZAMANDA YOLCULUK — ZT-006.2.3")
    print("Tμν / GREEN-KERNEL / SOFTENING ERROR BUDGET")
    print("="*94)
    print(f"Fixed validation points: {len(pts)}")
    print(f"Field h values: {hvals}")
    print(f"Green grids: {grids}")

    results={}

    for name in ["G1","G2","G3","G4"]:
        print(f"\n[{name}] T_mu_nu field sensitivity")
        Tref=T_at_fixed_points(name,pts,b,d,hvals[-1])
        field_rows=[]
        for h in hvals:
            T=T_at_fixed_points(name,pts,b,d,h)
            field_rows.append({
                "h_m":h,
                "T00_mean":float(T[:,0,0].mean()),
                "T0i_mean_norm":float(np.linalg.norm(T[:,0,1:],axis=1).mean()),
                "relative_tensor_error_vs_h005":relative_tensor_l2(T,Tref),
            })
        for r in field_rows:
            print(f"  h={r['h_m']:.4f} relT={r['relative_tensor_error_vs_h005']:.6e}")

        print(f"[{name}] Green-kernel resolution sensitivity")
        kernel_rows=kernel_resolution_with_frozen_T(name,cyl,b,d,0.005,grids,obs)
        for r in kernel_rows:
            print(f"  grid={r['resolution']} N={r['points']} h00={r['center_h00']:.6e}")

        print(f"[{name}] softening sensitivity")
        src_pts,w=cylindrical_midpoints(cyl.radius_m,cyl.height_m,8,16,8)
        T=T_at_fixed_points(name,src_pts,b,d,0.005)
        soft_rows=softening_sensitivity(obs,src_pts,T,w,soft)

        results[name]={
            "field_sensitivity":field_rows,
            "green_resolution":kernel_rows,
            "softening":soft_rows
        }

    OUT.write_text(json.dumps({
        "stage":"ZT-006.2.3",
        "mode":args.mode,
        "results":results
    },indent=2),encoding="utf-8")

    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: numerical error-budget study; no Earth/CTC conclusion.")

if __name__=="__main__":
    main()
