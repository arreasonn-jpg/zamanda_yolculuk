"""Auto-patch for ZT-007B.1 — Kerr frame-dragging feasibility."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src" / "zt005" / "run_zt007_b1.py"
TST = ROOT / "tests" / "test_zt007_b1.py"

SRC_CODE = r'''
"""ZT-007B.1 - Kerr frame dragging: how much J is needed?"""
import argparse, json
from pathlib import Path
import numpy as np

G = 6.67430e-11
C = 2.99792458e8
MSUN = 1.98892e30


def kerr_frame_dragging(M, a_star, r, theta=np.pi/2):
    """
    Frame-dragging angular velocity of Kerr metric (Boyer-Lindquist).

    Omega = 2 M a r / (r^3 + a^2 r + 2 M a^2)

    with a = a_star * M (geometric units, G = c = 1).

    Parameters
    ----------
    M : float  mass in kg
    a_star : float  dimensionless spin (0 to 1)
    r : array  radius in meters
    theta : float  polar angle (rad)

    Returns
    -------
    Omega : array  angular velocity in rad/s
    """
    # Convert to geometric units (G = c = 1)
    M_geo = G * M / C**2       # meters
    a_geo = a_star * M_geo

    # Formula (SI): Omega_phys = c^3/(G M) * (M_geo a_geo r_geo) / (...)
    # Actually easier: use dimensionless formula and multiply
    r_geo = r  # in meters, since we're working with geometric coordinates
    num = 2.0 * M_geo * a_geo * r_geo
    den = r_geo**3 + a_geo**2 * r_geo + 2.0 * M_geo * a_geo**2
    Omega_geo = num / den   # units: 1/meter (since geometrized)
    Omega_phys = C * Omega_geo  # convert 1/m -> rad/s
    return Omega_phys


def ergosphere_radius(M, a_star, theta):
    """r_ergo = M + sqrt(M^2 - a^2 cos^2 theta) in geometric units (m)."""
    M_geo = G * M / C**2
    a_geo = a_star * M_geo
    cos2 = np.cos(theta)**2
    disc = M_geo**2 - a_geo**2 * cos2
    return M_geo + np.sqrt(np.maximum(disc, 0.0))


def lab_device_angular_momentum(I, A_coil, R_coil, N_turns=1):
    """
    Approximate J of a current-carrying coil.
    J = N * I * A (magnetic moment) * mu0 ... actually J = I * A * N for
    angular momentum when the current carries mass. In our case, angular
    momentum of EM field: L ~ (1/mu0) * I * N * pi R^2 ... order-of-magnitude.
    """
    # Use magnetic dipole moment m = N I A; J_field ~ m/c^2 (SI, since
    # field energy = mc^2 relation for angular momentum is subtle).
    # Simpler: for EM field, J_EM = (1/c^2) integral r x (E x B/mu0) dV
    # For a uniform solenoid, J_EM ~ I^2 R^3 mu0 / c^2 order-of-magnitude.
    mu0 = 4 * np.pi * 1e-7
    J_em = mu0 * (I * N_turns)**2 * R_coil**3 / C**2
    return J_em


def compare_sources():
    """Compare frame dragging of known astrophysical sources."""
    # Earth
    M_E = 5.972e24
    J_E = 5.86e33
    a_E = C * J_E / (G * M_E**2)
    Omega_E = kerr_frame_dragging(M_E, a_E, 6.371e6)

    # Neutron star (typical)
    M_NS = 1.4 * MSUN
    J_NS = 1e40
    a_NS = C * J_NS / (G * M_NS**2)
    a_NS = min(a_NS, 0.99)
    Omega_NS = kerr_frame_dragging(M_NS, a_NS, 1e4)

    # Kerr BH (a* = 0.998)
    M_BH = 10 * MSUN
    a_BH = 0.998
    M_geo_BH = G * M_BH / C**2
    r_horizon = M_geo_BH * (1 + np.sqrt(1 - a_BH**2))
    Omega_BH = kerr_frame_dragging(M_BH, a_BH, r_horizon * 1.01)

    return {
        "Earth": {
            "M_kg": M_E, "a_star": float(a_E),
            "Omega_drag_surface_rad_s": float(Omega_E),
        },
        "NeutronStar": {
            "M_kg": M_NS, "a_star": float(a_NS),
            "Omega_drag_surface_rad_s": float(Omega_NS),
        },
        "KerrBH_10Msun": {
            "M_kg": M_BH, "a_star": a_BH,
            "Omega_drag_horizon_rad_s": float(Omega_BH),
        },
    }


def required_J_for_earth_drag(target_Omega):
    """
    Given the Kerr formula Omega ~ (G J) / (c^2 r^3),
    invert to find required J for target Omega at r = 1 m.
    This is an approximation for a_star << 1.
    """
    r = 1.0
    J_required = target_Omega * C**2 * r**3 / G
    return J_required


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick","full"])
    a = ap.parse_args()

    print("=" * 72)
    print("ZAMANDA YOLCULUK - ZT-007B.1")
    print("KERR FRAME DRAGGING: MINIMUM REQUIRED ANGULAR MOMENTUM")
    print("=" * 72)
    print()

    # ---- Known sources ----
    print("Known astrophysical frame dragging (at surface):")
    print(f"{'Source':<22s} {'M (kg)':>12s} {'a*':>8s} {'Omega (rad/s)':>18s}")
    print("-" * 72)
    for name, data in compare_sources().items():
        M = data["M_kg"]
        a_s = data["a_star"]
        Om = data.get("Omega_drag_surface_rad_s",
                       data.get("Omega_drag_horizon_rad_s", 0))
        print(f"{name:<22s} {M:>12.4e} {a_s:>8.3f} {Om:>18.6e}")
    print()

    # ---- Our device ----
    print("Our device (G1, 100 A, 6 loops, R=7 cm):")
    I = 100.0
    R = 0.07
    N = 6
    J_device = lab_device_angular_momentum(I, np.pi*R**2, R, N)
    print(f"  Current I          = {I} A")
    print(f"  Coil radius R      = {R} m")
    print(f"  Number of loops N  = {N}")
    print(f"  Device J (EM field) = {J_device:.4e} kg*m^2/s")
    print()

    # ---- Required J for interesting drag ----
    print("=" * 72)
    print("GEREKLI J: 1 m yaricapta anlamli drag")
    print("=" * 72)
    targets = [
        ("GPS-level drag (1e-12 rad/s)", 1e-12),
        ("Human-perceptible (1e-6 rad/s)", 1e-6),
        ("1 second/hour drag (5e-4 rad/s)", 5e-4),
        ("Extreme (1 rad/s)", 1.0),
        ("Horizon-like (1e5 rad/s)", 1e5),
    ]
    print(f"{'Target':<40s} {'Omega (rad/s)':>15s} {'J_required (kg m^2/s)':>22s}")
    print("-" * 80)
    for name, Omega in targets:
        J_req = required_J_for_earth_drag(Omega)
        print(f"{name:<40s} {Omega:>15.4e} {J_req:>22.4e}")
    print()

    # ---- Gap ----
    print("=" * 72)
    print("SONUC")
    print("=" * 72)
    J_earth = 5.86e33
    ratio_device = J_device / J_earth
    print(f"Earth J = {J_earth:.3e} kg m^2/s")
    print(f"Device J = {J_device:.3e} kg m^2/s")
    print(f"Device / Earth orani = {ratio_device:.3e}")
    print()

    J_required_gps = required_J_for_earth_drag(1e-12)
    print(f"GPS-level drag requires J ~ {J_required_gps:.3e} kg m^2/s")
    print(f"Bu, Dunya'nin acisal momentumunun "
          f"{J_required_gps/J_earth:.3f} kati")
    print(f"Veya Guney'in acisal momentumunun "
          f"{J_required_gps/1e41:.3e} kati")
    print()

    # ---- Required device scale ----
    # Our J_device scales as I^2 R^3 N
    # To reach J_required_gps:  I_new = I * sqrt(J_req / J_device)
    I_req = I * np.sqrt(J_required_gps / J_device)
    print(f"Gerekli akim (ayni geometri) = {I_req:.4e} A")
    print(f"  Yildirim akiminin {I_req/1e5:.2e} kati")
    print(f"  Magnetar akiminin {I_req/1e24:.2e} kati")
    print()
    print("=" * 72)

    out = Path("zt007_b1_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-007B.1",
        "mode": a.mode,
        "device_J": float(J_device),
        "J_earth": J_earth,
        "required_J_gps_level": float(J_required_gps),
        "required_current_A": float(I_req),
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: Kerr feasibility; no CTC conclusion.")


if __name__ == "__main__":
    main()
'''

TST_CODE = r'''
"""ZT-007B.1 tests."""
import numpy as np
from zt005.run_zt007_b1 import (
    kerr_frame_dragging, compare_sources,
    lab_device_angular_momentum, required_J_for_earth_drag,
)
from zt005.run_zt007_b1 import G, C, MSUN


def test_earth_frame_dragging_reasonable():
    """Earth's Lense-Thirring at surface ~ 10^-14 rad/s."""
    M_E = 5.972e24
    J_E = 5.86e33
    a_E = C * J_E / (G * M_E**2)
    Om = kerr_frame_dragging(M_E, a_E, 6.371e6)
    assert 1e-16 < Om < 1e-12


def test_kerr_bh_dragging_is_large():
    """Kerr BH horizon angular velocity should be ~ 10^3-10^5 rad/s."""
    M = 10 * MSUN
    a_star = 0.998
    M_geo = G * M / C**2
    r_H = M_geo * (1 + np.sqrt(1 - a_star**2))
    Om = kerr_frame_dragging(M, a_star, r_H * 1.01)
    assert 1e2 < Om < 1e6


def test_sources_comparison_complete():
    s = compare_sources()
    assert "Earth" in s
    assert "NeutronStar" in s
    assert "KerrBH_10Msun" in s


def test_device_J_is_tiny():
    """Our lab device J ~ 10^-3 kg m^2/s."""
    J = lab_device_angular_momentum(100.0, np.pi * 0.07**2, 0.07, 6)
    assert 1e-6 < J < 1.0


def test_required_J_earth_drag_huge():
    """For 1 rad/s drag at r=1 m, J_required ~ 1e40."""
    J = required_J_for_earth_drag(1.0)
    assert 1e38 < J < 1e42


def test_device_to_target_ratio():
    """Device J is >20 orders below Earth J."""
    J_dev = lab_device_angular_momentum(100.0, np.pi*0.07**2, 0.07, 6)
    J_earth = 5.86e33
    assert J_dev / J_earth < 1e-30
'''

if __name__ == "__main__":
    SRC.write_text(SRC_CODE, encoding="utf-8")
    print(f"[OK] Wrote {SRC}")
    TST.write_text(TST_CODE, encoding="utf-8")
    print(f"[OK] Wrote {TST}")
    print()
    print("Now run:")
    print("  python -m pytest tests\\test_zt007_b1.py -v")
    print("  python -m zt005.run_zt007_b1 --mode quick")