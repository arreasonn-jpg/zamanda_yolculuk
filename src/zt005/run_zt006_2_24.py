"""ZT-006.2.24 - external-observer convergence with effective a_wire."""
import argparse, json
from pathlib import Path
import numpy as np
from .model import ActiveCylinder, Backpack, Drive
from .green_convergence import focused_source_nodes, fast_regularized_em_fields
from .green_solvers import solver_a_tensor_product
from .geometry import build_sources

OUT = Path("results/validated/zt006_2_24_results.json")
A_WIRE_EFF = 0.08  # numerik etkin yumusatma (m); fiziksel tel 5 mm

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["quick","full"], default="quick")
    a = ap.parse_args()
    cyl, bk, dv = ActiveCylinder(), Backpack(), Drive()
    obs = np.array([[0.,0.,2.],[0.,0.,3.],[0.,0.,5.],[2.,0.,0.],[3.,0.,0.]])
    grids = [(6,12,6),(8,16,8),(12,24,12)] if a.mode=="quick" else \
            [(8,16,8),(12,24,12),(16,32,16),(20,40,20)]
    out = {"stage":"ZT-006.2.24","mode":a.mode,
           "a_wire_effective_m":A_WIRE_EFF,"physical_wire_radius_m":0.005,
           "results":{}}
    all_ok = True
    for name in ["G1","G2","G3","G4"]:
        srcs = build_sources(name, bk)
        fields = {}
        for g in grids:
            pts, w = focused_source_nodes(cyl.radius_m, cyl.height_m, *g)
            T = fast_regularized_em_fields(pts, srcs, dv.peak_current_a,
                                           dv.frequency_hz, a_wire=A_WIRE_EFF)
            fields[str(g)] = solver_a_tensor_product(obs, pts, w, T)
        rows = []
        for i in range(1, len(grids)):
            hp, hc = fields[str(grids[i-1])], fields[str(grids[i])]
            rel = float(np.linalg.norm(hc-hp)/max(np.linalg.norm(hc),1e-300))
            rows.append({"from":list(grids[i-1]),"to":list(grids[i]),"rel":rel})
            if rel >= 0.05: all_ok = False
        out["results"][name] = rows
    out["convergence_passed"] = all_ok
    out["criterion"] = "rel < 0.05 between consecutive grids"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print("\nCONVERGENCE:", "PASS" if all_ok else "FAIL")

if __name__=="__main__": main()
