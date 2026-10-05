
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from .convergence import benchmark_single_loop

OUT=Path("results/exploratory/zt005_1_results")
OUT.mkdir(parents=True, exist_ok=True)

def main():
    data=benchmark_single_loop()

    print("="*76)
    print("ZAMANDA YOLCULUK — ZT-005.1")
    print("ANALYTIC LOOP BENCHMARK + DISCRETIZATION CONVERGENCE")
    print("="*76)
    print("Benchmark: circular filament loop, on-axis Biot-Savart analytic solution")
    print("B_axis(z) = mu0 I R^2 / [2(R^2+z^2)^(3/2)]")
    print()
    print("FILAMENT CONVERGENCE")
    for row in data["segment_convergence"]:
        print(f"Nseg={row['segments']:4d}  relL2={row['relative_L2_error']:.6e}")
    print()
    print("GRID CONVERGENCE OF curl(A)")
    for row in data["grid_convergence"]:
        print(f"N={row['grid_n']:3d}^3  relL2={row['relative_L2_error_curlA_axis']:.6e}")

    with open("zt005_1_results.json","w",encoding="utf-8") as f:
        json.dump(data,f,indent=2)

    s=data["segment_convergence"]
    g=data["grid_convergence"]

    plt.figure(figsize=(6,4))
    plt.loglog([r["segments"] for r in s],[r["relative_L2_error"] for r in s],"o-")
    plt.xlabel("Filament segments")
    plt.ylabel("Relative L2 error")
    plt.title("Biot-Savart filament convergence")
    plt.grid(True,which="both",alpha=.25)
    plt.tight_layout()
    plt.savefig(OUT/"filament_convergence.png",dpi=160)
    plt.close()

    plt.figure(figsize=(6,4))
    plt.loglog([r["grid_n"] for r in g],[r["relative_L2_error_curlA_axis"] for r in g],"o-")
    plt.xlabel("Grid points per axis")
    plt.ylabel("Relative L2 error")
    plt.title("curl(A) grid convergence")
    plt.grid(True,which="both",alpha=.25)
    plt.tight_layout()
    plt.savefig(OUT/"grid_convergence.png",dpi=160)
    plt.close()

    print()
    print("Results written to zt005_1_results.json and zt005_1_results/*.png")
    print("SCIENTIFIC STATUS: benchmark/convergence stage; no GR/CTC conclusion.")

if __name__=="__main__":
    main()
