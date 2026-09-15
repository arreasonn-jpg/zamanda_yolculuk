
import argparse,json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .batch_em_source import compute_T_batch
from .gauss_cyl_quadrature import gauss_cylindrical_nodes
from .cyl_quadrature import EINSTEIN_GREEN_FACTOR
from .spatial_gr import inverse_trace_reverse
from .singular_green import regularized_green_integral,exclusion_fraction

OUT=Path("zt006_2_18_results.json")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    a=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    obs=np.array([[0.,0.,0.],[.25,0,0],[.5,0,0],[0,.25,0],[0,0,.5],[.75,0,0]])
    grids=[(6,12,6),(8,16,8),(12,24,12)] if a.mode=="quick" else [(12,24,12),(16,32,16),(20,40,20)]
    rhos=[.005,.01,.02,.04] if a.mode=="quick" else [.0025,.005,.01,.02,.04]

    print("="*108)
    print("ZAMANDA YOLCULUK — ZT-006.2.18")
    print("SINGULAR-CELL REGULARIZED GREEN-KERNEL VALIDATION")
    print("="*108)
    print(f"Einstein Green factor: {EINSTEIN_GREEN_FACTOR:.6e}")
    print(f"rho values: {rhos}")

    results={}
    for name in ["G1","G2","G3","G4"]:
        rows=[]; print(f"\n[{name}]")
        for grid in grids:
            pts,w=gauss_cylindrical_nodes(cyl.radius_m,cyl.height_m,*grid)
            T=compute_T_batch(name,pts,b,d,h=.005)
            for rho in rhos:
                hs=[]; ex=[]
                for ob in obs:
                    hb=regularized_green_integral(ob,pts,w,T,EINSTEIN_GREEN_FACTOR,rho)
                    hs.append(inverse_trace_reverse(hb))
                    ex.append(exclusion_fraction(ob,pts,rho))
                hs=np.asarray(hs)
                rows.append({
                    "grid":list(grid),"rho_m":rho,
                    "mean_exclusion_fraction":float(np.mean(ex)),
                    "center_h00":float(hs[0,0,0]),
                    "center_h0i_norm":float(np.linalg.norm(hs[0,0,1:])),
                    "center_h_fro":float(np.linalg.norm(hs[0])
                )})
            print(f"  grid={grid} N={len(pts)} done")
        results[name]=rows

    OUT.write_text(json.dumps({"stage":"ZT-006.2.18","mode":a.mode,"results":results},indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: singular-cell Green-kernel validation; no CTC conclusion.")

if __name__=="__main__": main()
