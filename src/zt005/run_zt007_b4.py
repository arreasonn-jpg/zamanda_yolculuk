
"""ZT-007B.4 - Unified comparative feasibility analysis."""
import argparse, json
from pathlib import Path
import numpy as np

from zt005.poisson_solver import EINSTEIN_FACTOR
from zt005.run_zt006_3 import compute_geodesic_diagnostics
from zt005.run_zt007_b1 import (
    lab_device_angular_momentum, required_J_for_earth_drag,
)
from zt005.run_zt007_b2 import ctc_threshold, critical_omega
from zt005.run_zt007_b3 import (
    exotic_mass_required, casimir_energy_density,
)

G = 6.67430e-11
C = 2.99792458e8
HBAR = 1.054571817e-34
SECONDS_PER_YEAR = 3.15576e7


def collect_mechanism_data():
    """
    Return a unified dict of {mechanism_name: {...}} with all the
    input/output parameters needed for the comparative table.
    """
    # --- EM (Path A) ---
    res_g1 = compute_geodesic_diagnostics("G1", 40, 80, n_phi=8)
    h00_G1 = res_g1["h_00_center"]
    h00_target = 2.0 / SECONDS_PER_YEAR
    I_device = 100.0
    I_required = I_device * np.sqrt(h00_target / h00_G1)
    # Effective energy density scaling: rho ~ I^2
    rho_device = 5e-4    # J/m^3 (G1 typical)
    rho_required = rho_device * (I_required / I_device)**2

    em = {
        "name": "EM field (Path A)",
        "parameter_name": "Current I (A)",
        "parameter_device": I_device,
        "parameter_required": float(I_required),
        "gap": float(I_required / I_device),
        "rho_device": rho_device,
        "rho_required": float(rho_required),
        "rho_gap": float(rho_required / rho_device),
        "uses_exotic": False,
        "physically_allowed": True,
        "note": "Energy-density gap ~10^41",
    }

    # --- Kerr (Path B.1) ---
    J_device = lab_device_angular_momentum(100.0, np.pi*0.07**2, 0.07, 6)
    J_required = required_J_for_earth_drag(1e-6)  # human-perceptible drag
    kerr = {
        "name": "Kerr frame dragging",
        "parameter_name": "Angular momentum J (kg m^2/s)",
        "parameter_device": float(J_device),
        "parameter_required": float(J_required),
        "gap": float(J_required / J_device),
        "uses_exotic": False,
        "physically_allowed": True,
        "note": "Requires astrophysical J scale",
    }

    # --- Tipler (Path B.2) ---
    R_target = 1.0
    omega_device = 1e3   # rad/s (typical lab spin)
    rho_device_tip = 1e4  # kg/m^3 (steel)
    # Critical density for CTC at R=1 m, omega=1e3:
    rho_crit_tip = (C**2 / G) / (2.0 * omega_device * R_target**3)  # kg/m^3
    tipler = {
        "name": "Tipler cylinder",
        "parameter_name": "Density rho (kg/m^3)",
        "parameter_device": rho_device_tip,
        "parameter_required": float(rho_crit_tip),
        "gap": float(rho_crit_tip / rho_device_tip),
        "uses_exotic": False,
        "physically_allowed": False,  # infinite length
        "note": "Requires infinite cylinder (Hawking 1992)",
    }

    # --- Morris-Thorne (Path B.3) ---
    b0 = 1e-3   # 1 mm throat
    M_exotic = exotic_mass_required(b0)
    V_throat = 4.0 / 3.0 * np.pi * b0**3
    rho_required_worm = M_exotic * C**2 / V_throat
    rho_casimir_max = abs(casimir_energy_density(1e-9))  # 1 nm
    worm = {
        "name": "Morris-Thorne wormhole",
        "parameter_name": "Exotic density |rho| (J/m^3)",
        "parameter_device": float(rho_casimir_max),
        "parameter_required": float(rho_required_worm),
        "gap": float(rho_required_worm / rho_casimir_max),
        "uses_exotic": True,
        "physically_allowed": True,  # in classical GR
        "note": "Only classical-GR-consistent CTC mechanism",
    }

    return {
        "EM": em,
        "Kerr": kerr,
        "Tipler": tipler,
        "Wormhole": worm,
    }


def print_unified_table(data):
    print("=" * 92)
    print("UNIFIED FEASIBILITY TABLE (device-scale assessment)")
    print("=" * 92)
    header = (f"{'Mechanism':<26s} {'Param':<14s} "
              f"{'Device':>12s} {'Required':>14s} {'Gap':>12s} "
              f"{'Exotic?':>8s} {'Allowed?':>10s}")
    print(header)
    print("-" * 92)
    for key, d in data.items():
        pname = d["parameter_name"][:13]
        dev = f"{d['parameter_device']:.2e}"
        req = f"{d['parameter_required']:.2e}"
        gap = f"{d['gap']:.2e}"
        exo = "yes" if d["uses_exotic"] else "no"
        allow = "yes" if d["physically_allowed"] else "NO"
        print(f"{d['name']:<26s} {pname:<14s} "
              f"{dev:>12s} {req:>14s} {gap:>12s} "
              f"{exo:>8s} {allow:>10s}")
    print()


def print_energy_density_summary(data):
    """Convert all mechanisms to a comparable 'energy density gap'."""
    print("=" * 92)
    print("ENERGY-DENSITY-EQUIVALENT GAP (log10 scale)")
    print("=" * 92)
    print(f"{'Mechanism':<26s} {'log10(rho_req/rho_dev)':>26s} {'Short note':<35s}")
    print("-" * 92)

    log_gaps = {}
    for key, d in data.items():
        if "rho_gap" in d:
            logg = np.log10(d["rho_gap"])
            log_gaps[key] = logg
            print(f"{d['name']:<26s} {logg:>26.2f} {d['note']:<35s}")
        else:
            # For Kerr and Tipler, use proxy: gap in fundamental units
            logg = np.log10(d["gap"])
            log_gaps[key] = logg
            print(f"{d['name']:<26s} {logg:>26.2f} {d['note']:<35s}")
    print()
    return log_gaps


def find_closest_mechanism(data):
    """Which mechanism is closest to feasibility?"""
    gaps = {k: d["gap"] for k, d in data.items()}
    # Only consider physically allowed mechanisms
    allowed = {k: d["gap"] for k, d in data.items() if d["physically_allowed"]}
    best_allowed = min(allowed.items(), key=lambda x: x[1])
    best_overall = min(gaps.items(), key=lambda x: x[1])
    return best_allowed, best_overall


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    a = ap.parse_args()

    print("=" * 92)
    print("ZAMANDA YOLCULUK - ZT-007B.4")
    print("UNIFIED COMPARATIVE FEASIBILITY ANALYSIS")
    print("=" * 92)
    print()

    data = collect_mechanism_data()

    print_unified_table(data)
    log_gaps = print_energy_density_summary(data)

    best_allowed, best_overall = find_closest_mechanism(data)
    print("=" * 92)
    print("CLOSEST MECHANISM")
    print("=" * 92)
    print(f"Best physically-allowed:  {best_allowed[0]:<12s}  "
          f"gap = {best_allowed[1]:.2e}")
    print(f"Best overall (incl. non-physical):  {best_overall[0]:<12s}  "
          f"gap = {best_overall[1]:.2e}")
    print()

    # --- Common gap range ---
    gaps = [d["gap"] for d in data.values()]
    log_gaps_list = [np.log10(g) for g in gaps]
    print(f"Common gap range (all 4 mechanisms): "
          f"10^{min(log_gaps_list):.1f} to 10^{max(log_gaps_list):.1f}")
    print()

    # --- Paper summary ---
    print("=" * 92)
    print("PAPER SUMMARY")
    print("=" * 92)
    print("All four physically-motivated CTC mechanisms require resources")
    print("at least 10^35 (Kerr) to 10^41 (EM/Wormhole) times beyond")
    print("device-scale feasibility. No mechanism is remotely accessible.")
    print()
    print("Key insight for the paper:")
    print("  - EM (Path A): 10^41 energy density gap")
    print("  - Kerr:        10^35 angular momentum gap")
    print("  - Tipler:      requires infinite cylinder (Hawking 1992)")
    print("  - Wormhole:    10^40 exotic energy gap (Casimir-based)")
    print()
    print("This is a UNIFIED negative result, not mechanism-specific.")
    print("=" * 92)

    out = Path("zt007_b4_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-007B.4",
        "mode": a.mode,
        "data": {k: {kk: (float(vv) if isinstance(vv, (int,float)) else vv)
                     for kk, vv in v.items()} for k, v in data.items()},
        "log_gaps": {k: float(v) for k, v in log_gaps.items()},
        "best_allowed": {"mechanism": best_allowed[0], "gap": float(best_allowed[1])},
        "best_overall":  {"mechanism": best_overall[0],  "gap": float(best_overall[1])},
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: unified feasibility comparison; no CTC construction.")


if __name__ == "__main__":
    main()
