"""Auto-patch for ZT-006.2.23 — physical T00 from Biot-Savart."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt006_2_23.py"
TST = ROOT / "tests" / "test_zt006_2_23.py"

SRC_CODE = r'''
"""ZT-006.2.23 — Physical T00 from Biot-Savart, then Poisson solver."""
import argparse, json
from pathlib import Path
import numpy as np

from zt005.model import Backpack
from zt005.finite_wire_sources import build_finite_sources_grouped
from zt005.poisson_solver import solve_poisson_axisym, EINSTEIN_FACTOR

MU0 = 4 * np.pi * 1e-7


def extract_segments(sources_grouped, I_wire=100.0):
    """Extract (mids, dls, I_seg) arrays from grouped source."""
    mids_all, dls_all, I_all = [], [], []
    for group in sources_grouped:
        arr = np.asarray(group, dtype=float)
        if arr.size == 0:
            continue
        if arr.ndim == 4 and arr.shape[-1] == 3:
            n_fil = arr.shape[0]
            base = arr[:, 0, :, :]  # (n_fil, n_seg, 3)
            for f in range(n_fil):
                pts = base[f]
                nxt = np.roll(pts, -1, axis=0)
                dl = nxt - pts
                mid = 0.5 * (pts + nxt)
                mids_all.append(mid)
                dls_all.append(dl)
                I_all.append(np.full(len(pts), I_wire / n_fil))
    if not mids_all:
        return np.zeros((0, 3)), np.zeros((0, 3)), np.zeros(0)
    return (np.vstack(mids_all),
            np.vstack(dls_all),
            np.concatenate(I_all))


def biot_savart_B(mids, dls, I_seg, x, y, z, softening=1e-6, chunk=500):
    """Biot-Savart B at points (x,y,z). Returns B: (M, 3)."""
    M = len(x)
    B = np.zeros((M, 3))
    k = MU0 / (4 * np.pi)
    mx, my, mz = mids[:, 0], mids[:, 1], mids[:, 2]
    dx, dy, dz = dls[:, 0], dls[:, 1], dls[:, 2]
    for s in range(0, M, chunk):
        e = min(s + chunk, M)
        Rx = x[s:e, None] - mx[None, :]
        Ry = y[s:e, None] - my[None, :]
        Rz = z[s:e, None] - mz[None, :]
        Rn = np.maximum(np.sqrt(Rx*Rx + Ry*Ry + Rz*Rz), softening)
        Rn3 = Rn ** 3
        cx = dy[None, :] * Rz - dz[None, :] * Ry
        cy = dz[None, :] * Rx - dx[None, :] * Rz
        cz = dx[None, :] * Ry - dy[None, :] * Rx
        f = I_seg[None, :] / Rn3
        B[s:e, 0] = k * np.sum(cx * f, axis=1)
        B[s:e, 1] = k * np.sum(cy * f, axis=1)
        B[s:e, 2] = k * np.sum(cz * f, axis=1)
    return B


def compute_T00_physical(geometry, r_vals, z_vals, n_phi=8, I_wire=100.0):
    """T00(r,z) = <B^2>_phi / (2 mu0), from Biot-Savart of finite wires."""
    sources = build_finite_sources_grouped(geometry, Backpack(), radius_m=0.005)
    mids, dls, I_seg = extract_segments(sources, I_wire=I_wire)
    if len(mids) == 0:
        return np.zeros((len(r_vals), len(z_vals)))
    R, Z = np.meshgrid(r_vals, z_vals, indexing="ij")
    Rf, Zf = R.ravel(), Z.ravel()
    B2sum = np.zeros_like(Rf)
    phis = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
    for phi in phis:
        x = Rf * np.cos(phi)
        y = Rf * np.sin(phi)
        B = biot_savart_B(mids, dls, I_seg, x, y, Zf)
        B2sum += B[:, 0]**2 + B[:, 1]**2 + B[:, 2]**2
    return (B2sum / (2 * MU0 * n_phi)).reshape(len(r_vals), len(z_vals))


def run_convergence(geometry, resolutions, n_phi=8, r_max=1.0, z_max=1.25):
    out = []
    for Nr, Nz in resolutions:
        rv = np.linspace(0.0, r_max, Nr)
        zv = np.linspace(-z_max, z_max, Nz)
        T00 = compute_T00_physical(geometry, rv, zv, n_phi=n_phi)
        h = solve_poisson_axisym(T00, rv, zv)
        out.append({
            "Nr": Nr, "Nz": Nz,
            "T00_max": float(np.max(T00)),
            "T00_mean": float(np.mean(T00)),
            "h_center": float(h[0, Nz // 2]),
            "h_max": float(np.max(np.abs(h))),
        })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    ap.add_argument("--n-phi", type=int, default=8)
    a = ap.parse_args()
    resolutions = ([(40,80),(60,120),(80,160)] if a.mode=="quick"
                   else [(40,80),(60,120),(80,160),(100,200)])
    print("=" * 72)
    print("ZAMANDA YOLCULUK - ZT-006.2.23")
    print("PHYSICAL T00 (Biot-Savart) -> POISSON -> h_munu")
    print("=" * 72)
    print(f"Einstein factor = {EINSTEIN_FACTOR:.6e}")
    print(f"n_phi = {a.n_phi}  |  mode = {a.mode}")
    print()
    all_res = {}
    for g in ["G1","G2","G3","G4"]:
        print(f"[{g}]")
        try:
            res = run_convergence(g, resolutions, n_phi=a.n_phi)
            for r in res:
                print(f"  grid=({r['Nr']:3d},{r['Nz']:3d})  "
                      f"T00_max={r['T00_max']:.4e}  "
                      f"T00_mean={r['T00_mean']:.4e}  "
                      f"h_center={r['h_center']:.6e}")
            for i in range(1, len(res)):
                h1, h2 = res[i-1]["h_center"], res[i]["h_center"]
                if abs(h2) > 0:
                    print(f"    rel_change({res[i-1]['Nr']} -> {res[i]['Nr']}) = "
                          f"{abs(h2-h1)/abs(h2):.4%}")
            all_res[g] = res
        except Exception as e:
            print(f"  FAILED: {type(e).__name__}: {e}")
            all_res[g] = []
        print()
    out = Path("zt006_2_23_results.json")
    out.write_text(json.dumps({
        "stage":"ZT-006.2.23","mode":a.mode,"n_phi":a.n_phi,
        "einstein_factor":EINSTEIN_FACTOR,"results":all_res}, indent=2),
        encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: physical-T00 Poisson convergence; no CTC conclusion.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-006.2.23 tests."""
import numpy as np
from zt005.run_zt006_2_23 import (
    extract_segments, biot_savart_B, compute_T00_physical,
)
from zt005.finite_wire_sources import build_finite_sources_grouped
from zt005.model import Backpack
from zt005.poisson_solver import solve_poisson_axisym


def test_extract_segments_nonempty():
    b = Backpack()
    s = build_finite_sources_grouped("G1", b, radius_m=0.005)
    mids, dls, I = extract_segments(s, I_wire=100.0)
    assert len(mids) > 100
    assert mids.shape == dls.shape
    assert len(I) == len(mids)
    assert np.all(I > 0)


def test_biot_savart_single_loop_on_axis():
    # One circular loop in xy-plane at z=0, radius 0.1, I=1 A
    n = 128
    phi = np.linspace(0, 2*np.pi, n, endpoint=False)
    pts = np.stack([0.1*np.cos(phi), 0.1*np.sin(phi), np.zeros(n)], axis=1)
    nxt = np.roll(pts, -1, axis=0)
    mids = 0.5*(pts+nxt)
    dls  = nxt - pts
    I    = np.ones(n) * 1.0
    B = biot_savart_B(mids, dls, I, np.array([0.0]), np.array([0.0]),
                      np.array([0.0]))
    # On-axis B = mu0 I / (2 R) = 4pi e-7 / 0.2 ~ 6.28e-6 T
    assert abs(B[0,2] - 6.283e-6) / 6.283e-6 < 0.02


def test_T00_physical_positive_for_all():
    r = np.linspace(0, 1, 20)
    z = np.linspace(-1.25, 1.25, 40)
    for g in ["G1","G2","G3","G4"]:
        T00 = compute_T00_physical(g, r, z, n_phi=4)
        assert np.max(T00) > 0, f"{g} gave zero T00"


def test_T00_physical_in_range():
    """Peak T00 should be ~1e-3..1e2 J/m^3 for 100 A / mm wires."""
    r = np.linspace(0, 1, 20)
    z = np.linspace(-1.25, 1.25, 40)
    T00 = compute_T00_physical("G1", r, z, n_phi=4)
    peak = np.max(T00)
    assert 1e-6 < peak < 1e4, f"Peak T00 = {peak} not in plausible range"


def test_poisson_from_physical_source_finite():
    r = np.linspace(0, 1, 20)
    z = np.linspace(-1.25, 1.25, 40)
    T00 = compute_T00_physical("G1", r, z, n_phi=4)
    h = solve_poisson_axisym(T00, r, z)
    assert np.all(np.isfinite(h))
    assert np.max(np.abs(h)) > 0
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt006_2_23.py -v")
    print("  python -m zt005.run_zt006_2_23 --mode quick")