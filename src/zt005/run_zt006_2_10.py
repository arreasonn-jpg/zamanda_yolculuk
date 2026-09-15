
import argparse,json
from pathlib import Path
import numpy as np

from .model import ActiveCylinder,Backpack,Drive
from .geometry import build_sources
from .batch_em_source import compute_T_batch
from .gauss_cyl_quadrature import gauss_cylindrical_nodes
from .finite_wire_quadrature import min_distance_to_filament,masked_volume_integral
from .spatial_gr import inverse_trace_reverse
from .cyl_quadrature import EINSTEIN_GREEN_FACTOR

OUT=Path("zt006_2_10_results.json")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    a=ap.parse_args()

    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    obs=np.array([[0,0,0],[.5,0,0],[0,.5,0],[0,0,.5],[.9,0,0]],float)

    grids=[(6,12,6),(8,16,8),(12,24,12)] if a.mode=="quick" else [(8,16,8),(12,24,12),(16,32,16),(20,40,20)]
    cutoffs=[0.005,0.01,0.02,0.03]

    print("="*102)
    print("ZAMANDA YOLCULUK — ZT-006.2.10")
    print("FINITE-WIRE DISTANCE-CUTOFF CONVERGENCE")
    print("="*102)
    print(f"Einstein Green factor: {EINSTEIN_GREEN_FACTOR:.6e}")
    print(f"Grids: {grids}")
    print(f"Cutoffs: {cutoffs}")
    print("Purpose: determine whether Green-integral instability is dominated by ideal-filament near-source samples.")

    results={}
    for name in ["G1","G2","G3","G4"]:
        sources=build_sources(name,b)
        rows=[]
        print(f"\n[{name}]")
        for res in grids:
            nr,nphi,nz=res
            pts,w=gauss_cylindrical_nodes(cyl.radius_m,cyl.height_m,nr,nphi,nz)
            T=compute_T_batch(name,pts,b,d,h=.005)
            dist=min_distance_to_filament(pts,sources)
            for cutoff in cutoffs:
                keep=dist>=cutoff
                vals=[]
                for ob in obs:
                    hb=masked_volume_integral(ob,pts,w,T,keep,EINSTEIN_GREEN_FACTOR,softening_m=.0025)
                    hp=inverse_trace_reverse(hb)
                    vals.append((float(hp[0,0]),float(np.linalg.norm(hp[0,1:]))))
                center=vals[0]
                rows.append({
                    "resolution":list(res),
                    "cutoff_m":float(cutoff),
                    "kept_fraction":float(keep.mean()),
                    "center_h00":center[0],
                    "center_h0i_norm":center[1],
                    "all_observations":vals
                })
            print(f"  grid={res} N={len(pts)} min_dist={dist.min():.4e} m")

        # Group by cutoff and report fine-grid pairwise stabilization.
        summary={}
        for c in cutoffs:
            rr=[r for r in rows if r["cutoff_m"]==c]
            summary[str(c)]={
                "center_h00_by_grid":[r["center_h00"] for r in rr],
                "center_h0i_by_grid":[r["center_h0i_norm"] for r in rr],
            }
        results[name]={"rows":rows,"summary":summary}

    OUT.write_text(json.dumps({
        "stage":"ZT-006.2.10","mode":a.mode,
        "results":results
    },indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: finite-wire/cutoff sensitivity; no GR/CTC conclusion.")

if __name__=="__main__":
    main()
