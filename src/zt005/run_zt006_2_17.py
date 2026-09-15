import argparse, json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder, Backpack, Drive
from .batch_em_source import compute_T_batch
from .gauss_cyl_quadrature import gauss_cylindrical_nodes
from .cyl_quadrature import EINSTEIN_GREEN_FACTOR, quadrature_integral_many
from .spatial_gr import inverse_trace_reverse
from .metric_validation import diagnostics

OUT = Path("zt006_2_17_results.json")

def observers():
    return np.array([[0,0,0],[0.25,0,0],[0.5,0,0],[0,0.25,0],[0,0,0.5],[0.75,0,0]], float)

def calc(name, obs, cyl, b, d, grid):
    pts, w = gauss_cylindrical_nodes(cyl.radius_m, cyl.height_m, *grid)
    T = compute_T_batch(name, pts, b, d, h=0.005)
    raw = quadrature_integral_many(obs, pts, T, w, softening_m=0.0025)
    trace = np.stack([inverse_trace_reverse(x) for x in raw])
    return raw, trace, len(pts)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["quick","full"], default="quick")
    args = ap.parse_args()
    cyl, b, d = ActiveCylinder(), Backpack(), Drive()
    obs = observers()
    grids = [(8,16,8),(12,24,12)] if args.mode=="quick" else [(12,24,12),(16,32,16),(20,40,20)]
    results = {}
    print("="*108)
    print("ZAMANDA YOLCULUK — ZT-006.2.17")
    print("METRIC-CONVENTION / WEAK-FIELD SANITY GATE")
    print("="*108)
    print(f"Einstein Green factor: {EINSTEIN_GREEN_FACTOR:.6e}")
    for name in ["G1","G2","G3","G4"]:
        print(f"\n[{name}]")
        rows=[]
        for grid in grids:
            raw, trace, n = calc(name, obs, cyl, b, d, grid)
            dr = diagnostics(raw[0])
            dt = diagnostics(trace[0])
            row = {
                "grid": list(grid),
                "source_points": n,
                "raw_hbar": raw[0].tolist(),
                "trace_reversed_hbar": trace[0].tolist(),
                "raw": dr,
                "trace_reversed": dt,
            }
            rows.append(row)
            print(f"  grid={grid} N={n} raw||h||={dr['h_frobenius']:.6e} "
                  f"trace||h||={dt['h_frobenius']:.6e} "
                  f"raw det={dr['metric_invariants']['det']:.6e} "
                  f"trace det={dt['metric_invariants']['det']:.6e} "
                  f"raw timelike={dr['rest_frame_timelike']['timelike']} "
                  f"trace timelike={dt['rest_frame_timelike']['timelike']}")
        results[name]=rows
    OUT.write_text(json.dumps({"stage":"ZT-006.2.17","mode":args.mode,"results":results}, indent=2), encoding="utf-8")
    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: metric-convention/weak-field sanity gate; no CTC conclusion.")

if __name__ == "__main__":
    main()
