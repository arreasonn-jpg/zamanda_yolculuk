
import json,argparse
from pathlib import Path
import numpy as np
from .model import ActiveCylinder,Backpack,Drive
from .robust_convergence import sample_points,field_T,robust_compare,wire_distance

OUT=Path("results/exploratory/zt006_2_8_results.json")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--mode",choices=["quick","full"],default="quick")
    args=ap.parse_args()
    cyl,b,d=ActiveCylinder(),Backpack(),Drive()
    seeds=[11,22] if args.mode=="quick" else [11,22,33,44,55]
    n=150 if args.mode=="quick" else 400
    hs=[0.02,0.01,0.005]
    cutoffs=[0.0,0.01,0.02]
    print("="*100)
    print("ZAMANDA YOLCULUK — ZT-006.2.8")
    print("ROBUST CROSS-SEED / DISTANCE-AWARE FIELD CONVERGENCE")
    print("="*100)
    print(f"seeds={seeds}, points/seed={n}, h={hs}, cutoffs={cutoffs}")
    results={}
    for name in ["G1","G2","G3","G4"]:
        seed_rows=[]
        for seed in seeds:
            pts=sample_points(cyl,n,seed)
            Eref,Bref,Tref,sources=field_T(name,pts,b,d,0.005)
            dist=wire_distance(pts,sources)
            rows=[]
            for h in hs[:-1]:
                E,B,T,_=field_T(name,pts,b,d,h)
                for c in cutoffs:
                    keep=dist>=c
                    if keep.sum()<10: continue
                    rows.append({
                        "h_m":h,
                        "cutoff_m":c,
                        "kept_fraction":float(keep.mean()),
                        "B_error":robust_compare(B[keep],Bref[keep]),
                        "T_error":robust_compare(T[keep],Tref[keep]),
                    })
            seed_rows.append({"seed":seed,"rows":rows})
        results[name]=seed_rows
        # Compact console summary: worst median T error at h=.01, cutoff 0.
        vals=[]
        for sr in seed_rows:
            for r in sr["rows"]:
                if r["h_m"]==0.01 and r["cutoff_m"]==0.0: vals.append(r["T_error"]["median"])
        print(f"{name} h=.01 T-median-error across seeds: {vals}")
    OUT.write_text(json.dumps({"stage":"ZT-006.2.8","mode":args.mode,"results":results},indent=2),encoding="utf-8")
    print(f"Results written to {OUT}")
    print("SCIENTIFIC STATUS: robustness/reproducibility study; no GR/CTC conclusion.")

if __name__=="__main__": main()
