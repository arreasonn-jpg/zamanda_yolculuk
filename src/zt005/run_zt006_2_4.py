
import argparse,json
from pathlib import Path
import numpy as np

from .model import ActiveCylinder, Backpack, Drive
from .geometry import build_sources
from .batch_em_source import compute_T_batch
from .cyl_quadrature import cylindrical_midpoints,quadrature_integral_many,EINSTEIN_GREEN_FACTOR
from .spatial_gr import inverse_trace_reverse
from .wire_regularization import min_distance_to_sources,regularization_report

OUT=Path("results/exploratory/zt006_2_4_results.json")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    args=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    obs=np.array([[0.,0.,0.],[.5,0,0],[0,.5,0],[0,0,.5],[.9,0,0]],float)
    grids=[(4,8,4),(6,12,6),(8,16,8)] if args.mode=="quick" else [(8,16,8),(12,24,12),(16,32,16)]
    cutoffs=[0.002,0.005,0.01,0.02]

    print("="*96)
    print("ZAMANDA YOLCULUK — ZT-006.2.4")
    print("FINITE-WIRE REGULARIZATION / SOURCE-SINGULARITY DIAGNOSTICS")
    print("="*96)
    print(f"Mode: {args.mode}")
    print(f"Einstein Green factor: {EINSTEIN_GREEN_FACTOR:.6e}")
    print("IMPORTANT: zero-thickness filaments are not trusted arbitrarily close to the source wire.")

    results={}
    for name in ["G1","G2","G3","G4"]:
        sources=build_sources(name,b)
        print(f"\n[{name}]")
        geo={"wire_report":None,"resolutions":[]}

        # Use the finest grid for a geometry-aware wire proximity report.
        fpts,_=cylindrical_midpoints(cyl.radius_m,cyl.height_m,8,16,8)
        geo["wire_report"]=regularization_report(fpts,sources,cutoffs)
        for row in geo["wire_report"]:
            print(f"  cutoff={row['cutoff_m']:.3f} m excluded={100*row['excluded_fraction']:.3f}% "
                  f"p01={row['p01_distance_m']:.4e} m")

        for res in grids:
            nr,nphi,nz=res
            pts,w=cylindrical_midpoints(cyl.radius_m,cyl.height_m,nr,nphi,nz)
            T=compute_T_batch(name,pts,b,d,h=.005)
            dmin=min_distance_to_sources(pts,sources)
            for cutoff in cutoffs:
                keep=dmin>=cutoff
                if not np.any(keep):
                    continue
                Tm=T[keep]
                wm=w[keep]
                pm=pts[keep]
                hb=quadrature_integral_many(obs,pm,Tm,wm,softening_m=0.0025)
                hp=np.stack([inverse_trace_reverse(x) for x in hb])
                geo["resolutions"].append({
                    "resolution":list(res),
                    "cutoff_m":float(cutoff),
                    "kept_points":int(keep.sum()),
                    "kept_fraction":float(keep.mean()),
                    "center_h00":float(hp[0,0,0]),
                    "center_h0i_norm":float(np.linalg.norm(hp[0,0,1:]))
                })
            print(f"  grid={res} N={len(pts)} computed")

        results[name]=geo

    OUT.write_text(json.dumps({"stage":"ZT-006.2.4","mode":args.mode,"results":results},indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: source singularity / finite-wire regularization diagnostics; no GR/CTC conclusion.")

if __name__=="__main__":
    main()
