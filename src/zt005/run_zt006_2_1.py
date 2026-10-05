
import json
from pathlib import Path
import numpy as np

from .model import ActiveCylinder, Backpack, Drive
from .geometry import build_sources
from .maxwell import vector_potential
from .pointwise_curl import curl_point_4
from .validated_stress_energy import em_stress_energy
from .cyl_quadrature import cylindrical_midpoints, quadrature_integral_many
from .spatial_gr import inverse_trace_reverse

OUT=Path("results/exploratory/zt006_2_1_results.json")
G = 6.67430e-11
C = 299792458.0

def compute_T(name, points, cyl, backpack, drive, h=0.005):
    sources = build_sources(name, backpack)
    omega = 2*np.pi*drive.frequency_hz

    def A_fn(p):
        return vector_potential(np.atleast_2d(p), sources, drive.peak_current_a)[0]

    A = np.vstack([A_fn(p) for p in points])
    B = np.vstack([curl_point_4(A_fn, p, h) for p in points])
    E = omega*A
    return em_stress_energy(E, B)["T"]

def main():
    cyl = ActiveCylinder()
    b = Backpack()
    drive = Drive()

    observers = np.array([
        [0.,0.,0.],
        [0.5,0.,0.],
        [0.,0.5,0.],
        [0.,0.,0.5],
        [0.9,0.,0.],
    ])

    resolutions = [
        (8, 16, 8),
        (12, 24, 12),
        (16, 32, 16),
    ]

    print("="*88)
    print("ZAMANDA YOLCULUK — ZT-006.2.1")
    print("DETERMINISTIC CYLINDRICAL QUADRATURE")
    print("="*88)
    print(f"Active volume: {np.pi*cyl.radius_m**2*cyl.height_m:.6f} m^3")
    print("Resolutions: (Nr,Nphi,Nz) =", resolutions)
    print("Source: EM only; Earth/device background excluded.")

    results = {}

    for name in ["G1","G2","G3","G4"]:
        print(f"\n[{name}]")
        rows = []

        # Build source T on the finest source grid once, then reuse prefixes
        # only if dimensions differ through exact nested construction.
        for nr, nphi, nz in resolutions:
            src_pts, weights = cylindrical_midpoints(
                cyl.radius_m, cyl.height_m, nr, nphi, nz
            )

            T = compute_T(name, src_pts, cyl, b, drive, h=0.005)
            # Retain a softened kernel only as a numerical regularizer.
            hbar = quadrature_integral_many(
                observers, src_pts, T, weights, softening_m=0.0025
            )

            hphys = np.stack([inverse_trace_reverse(x) for x in hbar])
            h00 = hphys[:,0,0]
            h0 = np.linalg.norm(hphys[:,0,1:], axis=1)

            row = {
                "resolution": [nr,nphi,nz],
                "source_points": int(len(src_pts)),
                "h00": h00.tolist(),
                "h0i_norm": h0.tolist(),
                "center_h00": float(h00[0]),
                "center_h0i_norm": float(h0[0]),
            }
            rows.append(row)

            print(
                f"  ({nr:2d},{nphi:2d},{nz:2d}) "
                f"N={len(src_pts):5d} "
                f"center h00={h00[0]:.6e} "
                f"|h0i|={h0[0]:.6e}"
            )

        # Relative change of center h00 across resolutions
        for i in range(1, len(rows)):
            prev = rows[i-1]["center_h00"]
            cur = rows[i]["center_h00"]
            rows[i]["relative_change_center_h00"] = abs(cur-prev)/max(abs(cur),1e-300)
        rows[0]["relative_change_center_h00"] = None

        results[name] = rows

    OUT.write_text(json.dumps({
        "stage": "ZT-006.2.1",
        "observers": observers.tolist(),
        "results": results,
    }, indent=2), encoding="utf-8")

    print(f"\nResults written to {OUT}")
    print("SCIENTIFIC STATUS: deterministic source-integration convergence; no GR/CTC conclusion.")

if __name__ == "__main__":
    main()
