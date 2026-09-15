
import argparse, json, time
from pathlib import Path
from .model import ActiveCylinder, Backpack, Drive
from .production import evaluate_geometry

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--grids",default="25,41,61",
                    help="grid side counts; y uses floor(0.56*n) to retain aspect ratio")
    ap.add_argument("--circle-segments",default="128,256",
                    help="comma-separated filament segment counts")
    ap.add_argument("--helix-segments",default="240,480")
    ap.add_argument("--chunk",type=int,default=4096)
    args=ap.parse_args()

    cyl,b,drive=ActiveCylinder(),Backpack(),Drive()
    grids=[int(x) for x in args.grids.split(",")]
    cs=[int(x) for x in args.circle_segments.split(",")]
    hs=[int(x) for x in args.helix_segments.split(",")]
    if len(cs)!=len(hs):
        raise SystemExit("--circle-segments and --helix-segments must have equal lengths")

    print("="*76)
    print("ZAMANDA YOLCULUK — ZT-005.2")
    print("PRODUCTION-RESOLUTION FIELD CONVERGENCE STUDY")
    print("="*76)
    print("Hard constraints: active D=2.00 m, H=2.50 m; backpack 0.32 x 0.15 x 0.45 m")
    print(f"Grids: {grids}")
    print(f"Filament pairs: {list(zip(cs,hs))}")
    print()

    results=[]
    for n,c,h in zip(grids,cs,hs):
        ny=max(17,int(round(n*2.5/2.0)))
        shape=(n,ny,n)
        for name in ["G1","G2","G3","G4"]:
            t=time.perf_counter()
            r=evaluate_geometry(name,cyl,b,drive.peak_current_a,shape,c,h,args.chunk)
            r["runtime_s"]=time.perf_counter()-t
            results.append(r)
            print(f"{name} grid={n}x{ny}x{n} seg={c}/{h} "
                  f"curlA-vs-BS={r['B_curlA_vs_BS_relL2']:.6e} "
                  f"divB={r['divB_relative']:.6e} "
                  f"time={r['runtime_s']:.1f}s")

    with open("zt005_2_results.json","w",encoding="utf-8") as f:
        json.dump(results,f,indent=2)

    print("\nResults written to zt005_2_results.json")
    print("SCIENTIFIC STATUS: production convergence only; no GR/CTC conclusion.")

if __name__=="__main__":
    main()
