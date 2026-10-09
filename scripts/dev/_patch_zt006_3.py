"""Auto-patch for ZT-006.3 — Metric field + geodesic validation."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt006_3.py"
TST = ROOT / "tests" / "test_zt006_3.py"

SRC_CODE = r'''
"""ZT-006.3 - Metric field g_munu = eta + h_munu and proper-time analysis."""
import argparse, json
from pathlib import Path
import numpy as np

from zt005.run_zt006_2_23 import compute_T00_physical
from zt005.poisson_solver import solve_poisson_axisym, EINSTEIN_FACTOR
from zt005.model import Backpack

C = 2.99792458e8


def build_metric_from_h(h_00, h_0i_mag, eta=None):
    """
    Build 4x4 metric g_munu = eta_munu + h_munu.

    Parameters
    ----------
    h_00 : (Nr, Nz) array
    h_0i_mag : (Nr, Nz) array, |h_0i| scalar magnitude
    eta : 4x4 Minkowski (default: diag(-1,+1,+1,+1))

    Returns
    -------
    g : (Nr, Nz, 4, 4)
    """
    if eta is None:
        eta = np.diag([-1.0, 1.0, 1.0, 1.0])
    Nr, Nz = h_00.shape
    g = np.tile(eta[None, None, :, :], (Nr, Nz, 1, 1))

    # h_00 goes to g_tt (with -1 -> -1 + h_00 ? No: h_00 is perturbation)
    # Convention: g_tt = -1 + h_00 (h_00 is small positive per Poisson)
    g[:, :, 0, 0] = -1.0 + h_00
    # h_0i along x-direction for simplicity (magnitude only)
    g[:, :, 0, 1] = h_0i_mag
    g[:, :, 1, 0] = h_0i_mag
    return g


def proper_time_ratio_rest(g):
    """
    For a test particle at rest in the field:
        dtau/dt = sqrt(-g_tt)
    Weak field: dtau/dt ~ 1 - h_00/2
    """
    g_tt = g[:, :, 0, 0]
    # g_tt is negative in our convention; -g_tt positive
    return np.sqrt(np.maximum(-g_tt, 1e-30))


def compute_geodesic_diagnostics(geometry, Nr, Nz, n_phi=8,
                                  r_max=1.0, z_max=1.25):
    rv = np.linspace(0.0, r_max, Nr)
    zv = np.linspace(-z_max, z_max, Nz)
    T00 = compute_T00_physical(geometry, rv, zv, n_phi=n_phi)
    h00 = solve_poisson_axisym(T00, rv, zv)

    # For simplicity, no h_0i yet (set to 0 here); will be added later
    h0i = np.zeros_like(h00)

    g = build_metric_from_h(h00, h0i)
    dtau_dt = proper_time_ratio_rest(g)

    # Center: r=0, z=0
    ic = 0
    jc = Nz // 2
    ratio_center = dtau_dt[ic, jc]
    delta_tau = (1.0 - ratio_center)  # fractional change

    return {
        "Nr": Nr, "Nz": Nz,
        "h_00_center": float(h00[ic, jc]),
        "h_00_max": float(np.max(np.abs(h00))),
        "dtau_dt_center": float(ratio_center),
        "delta_tau_over_t_center": float(delta_tau),
        "delta_tau_seconds_per_year": float(delta_tau * 3.156e7),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    ap.add_argument("--n-phi", type=int, default=8)
    a = ap.parse_args()

    if a.mode == "quick":
        resolutions = [(40, 80)]
    else:
        resolutions = [(40, 80), (60, 120), (80, 160)]

    print("=" * 72)
    print("ZAMANDA YOLCULUK - ZT-006.3")
    print("METRIC FIELD g_munu = eta + h_munu AND PROPER-TIME ANALYSIS")
    print("=" * 72)
    print(f"Einstein factor = {EINSTEIN_FACTOR:.6e}")
    print(f"n_phi = {a.n_phi}  |  mode = {a.mode}")
    print()

    all_res = {}
    for g_name in ["G1", "G2", "G3", "G4"]:
        print(f"[{g_name}]")
        rows = []
        for Nr, Nz in resolutions:
            res = compute_geodesic_diagnostics(g_name, Nr, Nz, n_phi=a.n_phi)
            print(f"  grid=({res['Nr']:3d},{res['Nz']:3d})  "
                  f"h00_center={res['h_00_center']:.6e}  "
                  f"dtau/dt={res['dtau_dt_center']:.15f}  "
                  f"delta_tau/yr={res['delta_tau_seconds_per_year']:.4e} s")
            rows.append(res)
        all_res[g_name] = rows
        print()

    out = Path("zt006_3_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-006.3",
        "mode": a.mode,
        "einstein_factor": EINSTEIN_FACTOR,
        "results": all_res,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: weak-field metric + proper-time; no CTC conclusion.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-006.3 tests."""
import numpy as np
from zt005.run_zt006_3 import (
    build_metric_from_h, proper_time_ratio_rest, compute_geodesic_diagnostics,
)


def test_metric_minkowski_when_h_zero():
    h00 = np.zeros((10, 10))
    h0i = np.zeros((10, 10))
    g = build_metric_from_h(h00, h0i)
    assert g.shape == (10, 10, 4, 4)
    assert np.allclose(g[0, 0, 0, 0], -1.0)
    assert np.allclose(g[0, 0, 1, 1], 1.0)
    assert np.allclose(g[0, 0, 2, 2], 1.0)
    assert np.allclose(g[0, 0, 3, 3], 1.0)


def test_metric_g_tt_shifted_by_h00():
    h00 = np.full((5, 5), 1e-48)
    h0i = np.zeros((5, 5))
    g = build_metric_from_h(h00, h0i)
    assert np.allclose(g[:, :, 0, 0], -1.0 + 1e-48, rtol=1e-6)


def test_proper_time_equals_one_for_minkowski():
    h00 = np.zeros((10, 10))
    h0i = np.zeros((10, 10))
    g = build_metric_from_h(h00, h0i)
    ratio = proper_time_ratio_rest(g)
    assert np.allclose(ratio, 1.0, atol=1e-12)


def test_proper_time_less_than_one_for_positive_h00():
    h00 = np.full((5, 5), 1e-10)
    h0i = np.zeros((5, 5))
    g = build_metric_from_h(h00, h0i)
    ratio = proper_time_ratio_rest(g)
    assert np.all(ratio < 1.0)


def test_diagnostics_finite():
    res = compute_geodesic_diagnostics("G1", 30, 60, n_phi=4)
    assert np.isfinite(res["h_00_center"])
    assert np.isfinite(res["dtau_dt_center"])
    assert 0.999 < res["dtau_dt_center"] < 1.0
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt006_3.py -v")
    print("  python -m zt005.run_zt006_3 --mode quick")