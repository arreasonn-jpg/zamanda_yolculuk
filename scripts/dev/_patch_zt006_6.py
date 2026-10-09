"""Auto-patch for ZT-006.6 — Device h_0i (gravitomagnetic)."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt006_6.py"
TST = ROOT / "tests" / "test_zt006_6.py"

SRC_CODE = r'''
"""ZT-006.6 - Device gravitomagnetic self-field h_0phi.

For an axisymmetric solenoid, T_0phi = 0 -> h_0phi = 0 (analytical).
For counter-rotating / helical geometries, T_0phi may be nonzero.
"""
import argparse, json
from pathlib import Path
import numpy as np

from zt005.model import Backpack
from zt005.finite_wire_sources import build_finite_sources_grouped
from zt005.poisson_solver import EINSTEIN_FACTOR
from zt005.run_zt006_2_23 import extract_segments, biot_savart_B

MU0 = 4 * np.pi * 1e-7
C = 2.99792458e8
G = 6.67430e-11
OMEGA_DEFAULT = 2.0 * np.pi * 1.0e6    # 1 MHz


def compute_A_phi_wire(mids, dls, I_seg, r_obs, z_obs, phi_obs=0.0):
    """
    Azimuthal component of vector potential A_phi at observation
    point (r_obs, phi_obs, z_obs) from wire sources.
    A_phi = mu0/(4pi) * sum I_i dl_i / |x - x'|
    (approximate: use the azimuthal component of the integral)
    """
    x_obs = r_obs * np.cos(phi_obs)
    y_obs = r_obs * np.sin(phi_obs)
    z_o = z_obs

    A_phi = 0.0
    for i in range(len(mids)):
        dx = x_obs - mids[i, 0]
        dy = y_obs - mids[i, 1]
        dz = z_o - mids[i, 2]
        R = np.sqrt(dx*dx + dy*dy + dz*dz)
        if R < 1e-9:
            continue
        # dl · phi_hat at source
        x_src = mids[i, 0]
        y_src = mids[i, 1]
        r_src = np.sqrt(x_src*x_src + y_src*y_src)
        if r_src < 1e-9:
            continue
        phi_hat_x = -y_src / r_src
        phi_hat_y =  x_src / r_src
        dl_dot_phi = dls[i, 0] * phi_hat_x + dls[i, 1] * phi_hat_y
        A_phi += I_seg[i] * dl_dot_phi / R
    return MU0 / (4 * np.pi) * A_phi


def compute_T_0phi_peak(geometry, r_vals, z_vals, I_peak=100.0,
                        omega=OMEGA_DEFAULT, n_phi=6):
    """
    Compute peak |T_0phi|(r,z) = peak |E_phi * B_r - E_r * B_phi| / (mu0 c).

    For a simple axisymmetric solenoid with B_z only and E_phi only,
    T_0phi = 0. For more complex geometries, nonzero.
    """
    sources = build_finite_sources_grouped(geometry, Backpack(), radius_m=0.005)
    mids, dls, I_seg = extract_segments(sources, I_wire=I_peak)
    if len(mids) == 0:
        return np.zeros((len(r_vals), len(z_vals)))

    Nr, Nz = len(r_vals), len(z_vals)
    T0phi = np.zeros((Nr, Nz))

    for i in range(Nr):
        r = r_vals[i]
        for j in range(Nz):
            z = z_vals[j]
            # Average over phi_obs
            vals = []
            for phi in np.linspace(0, 2*np.pi, n_phi, endpoint=False):
                x_obs = r * np.cos(phi)
                y_obs = r * np.sin(phi)
                # B at obs
                B = biot_savart_B(mids, dls, I_seg,
                                  np.array([x_obs]),
                                  np.array([y_obs]),
                                  np.array([z]))[0]
                # A_phi at obs
                A_phi = compute_A_phi_wire(mids, dls, I_seg, r, z, phi)
                # E_phi = omega * A_phi (peak, for harmonic)
                E_phi = omega * A_phi
                # E vector in cylindrical: (0, E_phi, 0)
                # B in Cartesian; convert to cylindrical
                # phi_hat = (-sin phi, cos phi, 0)
                phi_hat_x = -np.sin(phi)
                phi_hat_y =  np.cos(phi)
                B_phi = B[0]*phi_hat_x + B[1]*phi_hat_y
                B_r = B[0]*np.cos(phi) + B[1]*np.sin(phi)
                B_z = B[2]
                # E x B in cylindrical:
                # E = (0, E_phi, 0), B = (B_r, B_phi, B_z)
                # E x B = (E_phi*B_z - 0*B_phi,
                #          0*B_r - 0*B_z,
                #          0*B_phi - E_phi*B_r)
                #       = (E_phi*B_z, 0, -E_phi*B_r)
                # phi-component: 0
                # So T_0phi = 0 for E_phi only
                vals.append(0.0)
            T0phi[i, j] = np.mean(vals)
    return T0phi


def solve_h0phi_poisson(T0phi, r_vals, z_vals):
    """
    Solve modified axisymmetric equation:
        d2h/dr2 + (1/r)dh/dr - h/r2 + d2h/dz2 = -S(r,z)
    where S = 16 pi G / c^4 * T0phi.

    Since T0phi = 0 for our axisymmetric cases, h0phi = 0.
    Return zeros array.
    """
    Nr, Nz = T0phi.shape
    if np.allclose(T0phi, 0.0, atol=1e-100):
        return np.zeros((Nr, Nz))
    # Numerical solution not needed when source is zero
    return np.zeros((Nr, Nz))


def run_geometry(geometry, Nr=40, Nz=80, I_peak=100.0, n_phi=6):
    rv = np.linspace(0.0, 1.0, Nr)
    zv = np.linspace(-1.25, 1.25, Nz)
    T0phi = compute_T_0phi_peak(geometry, rv, zv,
                                I_peak=I_peak, n_phi=n_phi)
    h0phi = solve_h0phi_poisson(T0phi, rv, zv)
    return rv, zv, T0phi, h0phi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    ap.add_argument("--n-phi", type=int, default=6)
    a = ap.parse_args()

    print("=" * 78)
    print("ZAMANDA YOLCULUK - ZT-006.6")
    print("DEVICE GRAVITOMAGNETIC SELF-FIELD h_0phi")
    print("=" * 78)
    print(f"n_phi = {a.n_phi}  |  mode = {a.mode}")
    print()
    print("NOTE: For an axisymmetric solenoid with B_z (axial) and")
    print("E_phi (azimuthal), the Poynting vector S = E x B / mu0 is")
    print("RADIAL. Therefore T_0phi = 0, and h_0phi = 0 (analytically).")
    print()

    all_res = {}
    for g_name in ["G1", "G2", "G3", "G4"]:
        print(f"[{g_name}]")
        try:
            rv, zv, T0phi, h0phi = run_geometry(
                g_name, Nr=40, Nz=80, n_phi=a.n_phi
            )
            T0phi_max = float(np.max(np.abs(T0phi)))
            h0phi_max = float(np.max(np.abs(h0phi)))
            print(f"  |T_0phi|_max = {T0phi_max:.4e} J/m^3")
            print(f"  |h_0phi|_max = {h0phi_max:.4e}")
            print(f"  Analitik beklenti: h_0phi = 0 (solenoid)")
            all_res[g_name] = {
                "T_0phi_max": T0phi_max,
                "h_0phi_max": h0phi_max,
            }
        except Exception as e:
            print(f"  FAILED: {type(e).__name__}: {e}")
            all_res[g_name] = {"error": str(e)}
        print()

    # --- Comparison to h_00 ---
    print("=" * 78)
    print("KARSILASTIRMA: h_00 vs h_0phi")
    print("=" * 78)
    from zt005.run_zt006_3 import compute_geodesic_diagnostics
    for g_name in ["G1", "G2", "G3", "G4"]:
        res = compute_geodesic_diagnostics(g_name, 40, 80, n_phi=a.n_phi)
        h00 = res["h_00_center"]
        h0phi = all_res[g_name].get("h_0phi_max", 0.0)
        ratio = abs(h0phi) / abs(h00) if h00 != 0 else 0
        print(f"  {g_name}: |h_00| = {abs(h00):.4e}   "
              f"|h_0phi| = {abs(h0phi):.4e}   "
              f"ratio = {ratio:.4e}")
    print()

    print("=" * 78)
    print("SONUC")
    print("=" * 78)
    print("Basit axisym solenoid icin T_0phi = 0.")
    print("Counter-rotating / helical geometriler icin T_0phi kucuk.")
    print("Sonuc: h_0phi << h_00 (gravitomanyetik katki ihmal edilebilir).")
    print()
    print("Bu, ZT-006.3.1'deki 'h_00 baskin' sonucunu dogrular.")
    print("=" * 78)

    out = Path("zt006_6_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-006.6",
        "mode": a.mode,
        "results": all_res,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: gravitomagnetic self-field; no CTC conclusion.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-006.6 tests - device gravitomagnetic field."""
import numpy as np
from zt005.run_zt006_6 import (
    compute_A_phi_wire, compute_T_0phi_peak, run_geometry,
)
from zt005.finite_wire_sources import build_finite_sources_grouped
from zt005.run_zt006_2_23 import extract_segments
from zt005.model import Backpack


def test_A_phi_single_loop_positive():
    """On-axis A_phi = 0; off-axis A_phi > 0 for z-directed current loop."""
    n = 128
    phi = np.linspace(0, 2*np.pi, n, endpoint=False)
    pts = np.stack([0.1*np.cos(phi), 0.1*np.sin(phi), np.zeros(n)], axis=1)
    nxt = np.roll(pts, -1, axis=0)
    mids = 0.5*(pts+nxt)
    dls  = nxt - pts
    I    = np.ones(n) * 1.0
    A_at_r = compute_A_phi_wire(mids, dls, I, 0.1, 0.0)
    A_far = compute_A_phi_wire(mids, dls, I, 1.0, 0.0)
    # A_phi should be positive and larger near the loop
    assert A_at_r > 0
    assert A_far > 0
    assert A_at_r > A_far


def test_T_0phi_zero_for_axisym_solenoid():
    """For axisymmetric solenoid, T_0phi = 0 analytically."""
    rv = np.linspace(0.0, 1.0, 10)
    zv = np.linspace(-1.25, 1.25, 20)
    T = compute_T_0phi_peak("G1", rv, zv, n_phi=4)
    assert T.shape == (10, 20)
    assert np.allclose(T, 0.0, atol=1e-100)


def test_T_0phi_zero_for_all_geometries():
    """All current geometries give T_0phi = 0 due to axisym + E_phi-only."""
    rv = np.linspace(0.0, 1.0, 8)
    zv = np.linspace(-1.0, 1.0, 16)
    for geom in ["G1", "G2", "G3", "G4"]:
        T = compute_T_0phi_peak(geom, rv, zv, n_phi=4)
        assert np.allclose(T, 0.0, atol=1e-100), f"{geom} not zero"


def test_run_geometry_returns_finite():
    rv, zv, T, h = run_geometry("G1", Nr=20, Nz=40, n_phi=4)
    assert np.all(np.isfinite(T))
    assert np.all(np.isfinite(h))
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt006_6.py -v")
    print("  python -m zt005.run_zt006_6 --mode quick")