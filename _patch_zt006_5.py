"""Auto-patch for ZT-006.5 — Geodesic integration."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt006_5.py"
TST = ROOT / "tests" / "test_zt006_5.py"

SRC_CODE = r'''
"""ZT-006.5 - Geodesic integration in the weak-field metric."""
import argparse, json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from zt005.run_zt006_2_23 import compute_T00_physical
from zt005.poisson_solver import solve_poisson_axisym
from zt005.run_zt006_3 import build_metric_from_h

C = 2.99792458e8


def build_metric_grid(geometry, Nr=40, Nz=80, n_phi=8,
                      r_max=1.0, z_max=1.25):
    """Return (g_grid, r_vals, z_vals) for a geometry."""
    rv = np.linspace(0.0, r_max, Nr)
    zv = np.linspace(-z_max, z_max, Nz)
    T00 = compute_T00_physical(geometry, rv, zv, n_phi=n_phi)
    h00 = solve_poisson_axisym(T00, rv, zv)
    h0i = np.zeros_like(h00)
    g = build_metric_from_h(h00, h0i)
    return g, rv, zv


def interpolate_metric(g_grid, r_vals, z_vals, r, z):
    """
    Bilinear interpolation of g_munu at (r, z).
    Returns 4x4 array.
    """
    Nr, Nz = len(r_vals), len(z_vals)
    ri = np.clip(np.searchsorted(r_vals, r) - 1, 0, Nr - 2)
    zi = np.clip(np.searchsorted(z_vals, z) - 1, 0, Nz - 2)
    dr = (r_vals[ri+1] - r_vals[ri])
    dz = (z_vals[zi+1] - z_vals[zi])
    fr = (r - r_vals[ri]) / dr
    fz = (z - z_vals[zi]) / dz
    g = (g_grid[ri,   zi  ] * (1-fr) * (1-fz)
       + g_grid[ri+1, zi  ] * fr     * (1-fz)
       + g_grid[ri,   zi+1] * (1-fr) * fz
       + g_grid[ri+1, zi+1] * fr     * fz)
    return g


def christoffel_numerical(g_grid, r_vals, z_vals, coords, eps=None):
    """
    Compute Gamma^mu_{alpha beta} at (r, z) via finite differences.

    Coordinates used: x^mu = (t, x, y, z) but we work in cylindrical
    (r, z) with phi=0, so x = r, y = 0. g_munu depends only on (r, z).

    We use central differences on g_ab in (r, z).
    """
    Nr, Nz = len(r_vals), len(z_vals)
    r, z = coords
    if eps is None:
        eps = min((r_vals[1] - r_vals[0]) * 4,
                  (z_vals[1] - z_vals[0]) * 4)

    # Clamp to interior
    r = np.clip(r, r_vals[0] + eps, r_vals[-1] - eps)
    z = np.clip(z, z_vals[0] + eps, z_vals[-1] - eps)

    def g_at(rr, zz):
        return interpolate_metric(g_grid, r_vals, z_vals, rr, zz)

    # Derivatives: dg/dr, dg/dz
    dg_dr = (g_at(r + eps, z) - g_at(r - eps, z)) / (2 * eps)
    dg_dz = (g_at(r, z + eps) - g_at(r, z - eps)) / (2 * eps)

    # For our metric, depend only on r (axisym). So d/dx = d/dr,
    # d/dy = 0, d/dz = dg_dz.
    # In Cartesian x, y, z: r = sqrt(x^2 + y^2), dr/dx = x/r, dr/dy = y/r.
    # For simplicity, we'll work in cylindrical (r, phi, z) then convert.
    # BUT: geodesics in axisym metric are usually simpler in Cartesian
    # (x, y, z) with the metric having only (r, z) dependence.

    # For an axisym, static metric with g_ab = g_ab(r, z):
    # Gamma^t_{t a} = (1/2) g^{t t} d_a g_tt
    # Gamma^i_{t t} = -(1/2) g^{i i} d_i g_tt
    # Others by symmetry.

    g = g_at(r, z)
    g_inv = np.linalg.inv(g)

    # Build 4x4x4 Christoffel in Cartesian-like coordinates.
    # Coordinate labels: 0=t, 1=x, 2=y, 3=z.
    # g depends on x via r = sqrt(x^2 + y^2), so dg/dx = dg/dr * x/r.
    # We'll treat x = r (equivalent to phi=0 slice; y = 0).
    # In that slice, d/dx = d/dr exactly, d/dy = 0.

    # Build the derivative array: partial_a g_{bc} for a, b, c in 0..3
    dg = np.zeros((4, 4, 4))  # dg[a, b, c] = d_a g_{bc}

    # d/dx = d/dr, d/dy = 0, d/dz = d/dz
    dg[1, :, :] = dg_dr   # d/dx
    dg[2, :, :] = 0.0     # d/dy (axisym)
    dg[3, :, :] = dg_dz   # d/dz

    # Gamma^a_{bc} = (1/2) g^{ad} (d_b g_{dc} + d_c g_{db} - d_d g_{bc})
    Gamma = np.zeros((4, 4, 4))
    for a in range(4):
        for b in range(4):
            for c in range(4):
                s = 0.0
                for d in range(4):
                    s += g_inv[a, d] * (
                        dg[b, d, c] + dg[c, d, b] - dg[d, b, c]
                    )
                Gamma[a, b, c] = 0.5 * s
    return Gamma


def geodesic_rhs(tau, y, g_grid, r_vals, z_vals):
    """
    y = [t, x, y_coord, z, dt/dtau, dx/dtau, dy/dtau, dz/dtau]
    Actually 'x' here is x-coordinate (Cartesian), not time.
    Let's use y = [X^0, X^1, X^2, X^3, U^0, U^1, U^2, U^3]
    where X^0 = t, X^1 = x, X^2 = y, X^3 = z.
    """
    X = y[:4]
    U = y[4:]
    t, x, yc, z = X

    r = np.sqrt(x*x + yc*yc)
    Gamma = christoffel_numerical(g_grid, r_vals, z_vals, (r, z))

    # dU^a/dtau = -Gamma^a_{bc} U^b U^c
    dU = np.zeros(4)
    for a in range(4):
        s = 0.0
        for b in range(4):
            for c in range(4):
                s += Gamma[a, b, c] * U[b] * U[c]
        dU[a] = -s
    return np.concatenate([U, dU])


def integrate_geodesic(g_grid, r_vals, z_vals,
                       x0=(0.0, 0.0, 0.0),
                       U0=(1.0, 0.0, 0.0, 0.0),
                       tau_max=1.0e-3, n_steps=200):
    """
    Integrate geodesic from x0 with 4-velocity U0 over tau in [0, tau_max].

    Default: rest particle at origin (U0 = (1, 0, 0, 0)).
    """
    y0 = np.concatenate([np.asarray(x0, dtype=float),
                         np.asarray(U0, dtype=float)])
    taus = np.linspace(0.0, tau_max, n_steps)
    sol = solve_ivp(
        lambda tau, y: geodesic_rhs(tau, y, g_grid, r_vals, z_vals),
        (0.0, tau_max), y0, t_eval=taus,
        method="RK45", rtol=1e-10, atol=1e-12
    )
    return sol


def geodesic_norm(g_grid, r_vals, z_vals, X, U):
    """Compute g_munu U^mu U^nu along a trajectory."""
    N = X.shape[1]
    norms = np.zeros(N)
    for i in range(N):
        t, x, yc, z = X[:, i]
        r = np.sqrt(x*x + yc*yc)
        g = interpolate_metric(g_grid, r_vals, z_vals, r, z)
        norms[i] = U[:, i] @ g @ U[:, i]
    return norms


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    ap.add_argument("--n-phi", type=int, default=8)
    a = ap.parse_args()

    print("=" * 72)
    print("ZAMANDA YOLCULUK - ZT-006.5")
    print("GEODESIC INTEGRATION IN WEAK-FIELD METRIC")
    print("=" * 72)
    print(f"n_phi = {a.n_phi}  |  mode = {a.mode}")
    print()

    all_res = {}
    for g_name in ["G1", "G2", "G3", "G4"]:
        print(f"[{g_name}]")
        try:
            g_grid, rv, zv = build_metric_grid(g_name, Nr=40, Nz=80,
                                               n_phi=a.n_phi)
            # Rest particle at origin
            sol_rest = integrate_geodesic(g_grid, rv, zv,
                                          x0=(0, 0, 0),
                                          U0=(1, 0, 0, 0),
                                          tau_max=1e-3)
            X_rest = sol_rest.y[:4]
            U_rest = sol_rest.y[4:]
            norms = geodesic_norm(g_grid, rv, zv, X_rest, U_rest)

            # Check: g_tt at origin
            g_origin = interpolate_metric(g_grid, rv, zv, 0.0, 0.0)
            g_tt = g_origin[0, 0]
            U_t_final = U_rest[0, -1]
            print(f"  grid=({len(rv)},{len(zv)})  "
                  f"g_tt(0,0)={g_tt:.6e}  "
                  f"dt/dtau_final={U_t_final:.15f}  "
                  f"||U||^2_mean={norms.mean():.4e}  "
                  f"||U||^2_var={norms.var():.4e}")
            all_res[g_name] = {
                "g_tt_origin": float(g_tt),
                "U_t_final": float(U_t_final),
                "norm_mean": float(norms.mean()),
                "norm_var": float(norms.var()),
            }
        except Exception as e:
            print(f"  FAILED: {type(e).__name__}: {e}")
            all_res[g_name] = {"error": str(e)}
        print()

    print("=" * 72)
    print("YORUM")
    print("=" * 72)
    print("Norm (g_munu U^mu U^nu) should be ~ -1 for a timelike geodesic.")
    print("Variation across trajectory should be << 1 if integration is good.")
    print("U_t_final ~ 1.000... if proper time ~ coordinate time.")
    print("=" * 72)

    out = Path("zt006_5_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-006.5",
        "mode": a.mode,
        "n_phi": a.n_phi,
        "results": all_res,
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: geodesic integration; no CTC conclusion.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-006.5 tests - geodesic integration."""
import numpy as np
from zt005.run_zt006_5 import (
    build_metric_grid, interpolate_metric,
    christoffel_numerical, integrate_geodesic, geodesic_norm,
)


def test_metric_interpolation_returns_4x4():
    g, rv, zv = build_metric_grid("G1", Nr=20, Nz=40, n_phi=4)
    g_interp = interpolate_metric(g, rv, zv, 0.5, 0.0)
    assert g_interp.shape == (4, 4)


def test_christoffel_shape():
    g, rv, zv = build_metric_grid("G1", Nr=20, Nz=40, n_phi=4)
    Gamma = christoffel_numerical(g, rv, zv, (0.3, 0.0))
    assert Gamma.shape == (4, 4, 4)


def test_rest_geodesic_norm_is_minus_one():
    """A particle at rest should have g_uu U^u U^v ~ -1."""
    g, rv, zv = build_metric_grid("G1", Nr=20, Nz=40, n_phi=4)
    sol = integrate_geodesic(g, rv, zv,
                             x0=(0, 0, 0),
                             U0=(1.0, 0, 0, 0),
                             tau_max=1e-4)
    X = sol.y[:4]; U = sol.y[4:]
    norms = geodesic_norm(g, rv, zv, X, U)
    assert np.abs(norms.mean() + 1.0) < 1e-3


def test_rest_geodesic_stays_at_origin():
    """Without velocity, spatial coordinates should stay at 0."""
    g, rv, zv = build_metric_grid("G1", Nr=20, Nz=40, n_phi=4)
    sol = integrate_geodesic(g, rv, zv,
                             x0=(0, 0, 0),
                             U0=(1.0, 0, 0, 0),
                             tau_max=1e-4)
    X = sol.y[:4]
    # Spatial coords should not have moved appreciably
    assert np.max(np.abs(X[1:])) < 1e-3


def test_proper_time_close_to_coordinate_time():
    """For weak field, dt/dtau ~ 1."""
    g, rv, zv = build_metric_grid("G1", Nr=20, Nz=40, n_phi=4)
    sol = integrate_geodesic(g, rv, zv,
                             x0=(0, 0, 0),
                             U0=(1.0, 0, 0, 0),
                             tau_max=1e-4)
    U_t = sol.y[4]
    assert np.allclose(U_t, 1.0, atol=1e-6)
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt006_5.py -v")
    print("  python -m zt005.run_zt006_5 --mode quick")