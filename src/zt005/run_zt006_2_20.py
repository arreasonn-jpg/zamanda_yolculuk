import argparse,json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .external_green_validation import external_observers,compute_external_metric,component_errors

OUT=Path("zt006_2_20_results.json")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    a=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    obs=external_observers()
    grids=[(6,12,6),(8,16,8),(12,24,12)] if a.mode=="quick" else [(8,16,8),(12,24,12),(16,32,16),(20,40,20)]
    print("="*108)
    print("ZAMANDA YOLCULUK — ZT-006.2.20")
    print("EXTERNAL-OBSERVER GREEN-KERNEL CONVERGENCE")
    print("="*108)
    print("No interior 1/r singularity: observers are outside the active cylinder.")
    results={}
    for name in ["G1","G2","G3","G4"]:
        print(f"\n[{name}]")
        fields={}
        for grid in grids:
            h,N=compute_external_metric(name,obs,cyl,b,d,grid)
            fields[str(grid)]=h
            print(f"  grid={grid} N={N} |h|F:", " ".join(f"{np.linalg.norm(x):.3e}" for x in h))
        fine=fields[str(grids[-1])]
        rows=[]
        for grid in grids:
            arr=fields[str(grid)]
            errs=[component_errors(arr[i],fine[i]) for i in range(len(obs))]
            row={"grid":list(grid),
                 "max_tensor_rel":float(max(e["tensor_rel"] for e in errs)),
                 "max_h00_rel":float(max(e["h00_rel"] for e in errs)),
                 "max_h00_abs":float(max(e["h00_abs"] for e in errs)),
                 "max_h0i_abs":float(max(e["h0i_abs"] for e in errs)),
                 "observer_errors":errs}
            rows.append(row)
            print(f"    grid={grid} tensorRel={row['max_tensor_rel']:.3e} "
                  f"h00Rel={row['max_h00_rel']:.3e} |Δh00|={row['max_h00_abs']:.3e} "
                  f"|Δh0i|={row['max_h0i_abs']:.3e}")
        results[name]=rows
    OUT.write_text(json.dumps({"stage":"ZT-006.2.20","mode":a.mode,
                               "observers":obs.tolist(),"results":results},indent=2),encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: external-observer Green convergence; no CTC conclusion.")

if __name__=="__main__": main()
