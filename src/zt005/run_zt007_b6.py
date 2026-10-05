
"""ZT-007B.6 - Godel metric: CTC radius vs rotation rate."""
import argparse, json
from pathlib import Path
import numpy as np

G = 6.67430e-11
C = 2.99792458e8
HBAR = 1.054571817e-34
MSUN = 1.98892e30
MEARTH = 5.972e24
PC = 3.0857e16


def godel_density(omega):
    """
    Godel universe energy density:
        rho = omega^2 / (2 pi G)
    For rho > 0 (dust), omega real.
    """
    return omega**2 / (2.0 * np.pi * G)


def godel_ctc_radius(omega):
    """
    Critical radius for closed timelike curves in Godel metric:
        r_CTC = (c / (omega * sqrt(2))) * ln(1 + sqrt(2))
              ~ 0.6232 * c / omega
    """
    return C / (omega * np.sqrt(2)) * np.log(1.0 + np.sqrt(2))


def godel_omega_for_ctc_at_R(R):
    """Inverse: given target radius R, what omega is needed?"""
    return C * np.log(1.0 + np.sqrt(2)) / (R * np.sqrt(2))


def godel_comparison_table():
    """Known objects vs Godel-like rotation."""
    rows = []
    # Godel universe itself
    omega_godel = 1.0 / (1.0 * PC / C)   # ~ 1/(pc/c) rotation
    rows.append(("Godel universe (rho ~ 1e-28)", omega_godel,
                 godel_ctc_radius(omega_godel), godel_density(omega_godel)))
    # Slow rotation: omega = 1e-15 rad/s (galactic-scale)
    for omega, name in [
        (1e-20, "Cosmic rotation limit"),
        (1e-15, "Galactic rotation"),
        (1e-9,  "Solar orbital"),
        (1e-3,  "Slow lab rotation"),
        (1.0,   "1 rad/s"),
        (1e3,   "Millisecond pulsar"),
        (1e6,   "Neutron star surface"),
        (1e12,  "Exotic (hypothetical)"),
        (1e21,  "Planck rotation"),
    ]:
        rows.append((name, omega, godel_ctc_radius(omega),
                     godel_density(omega)))
    return rows


def mass_of_uniform_sphere(rho, R):
    """Mass of uniform sphere."""
    return rho * (4.0/3.0) * np.pi * R**3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="quick", choices=["quick", "full"])
    a = ap.parse_args()

    print("=" * 78)
    print("ZAMANDA YOLCULUK - ZT-007B.6")
    print("GODEL METRIC: CTC RADIUS vs ROTATION RATE")
    print("=" * 78)
    print()
    print("Godel universe CTC radius:")
    print("    r_CTC = c * ln(1 + sqrt(2)) / (omega * sqrt(2))")
    print("          ~ 0.6232 * c / omega")
    print("Energy density:")
    print("    rho = omega^2 / (2 pi G)")
    print()

    # --- Table of r_CTC vs omega ---
    print("=" * 78)
    print("CTC RADIUS vs ROTATION RATE")
    print("=" * 78)
    print(f"{'Scenario':<30s} {'omega (rad/s)':>16s} {'r_CTC (m)':>16s} "
          f"{'rho (kg/m^3)':>18s}")
    print("-" * 78)
    for name, omega, r_ctc, rho in godel_comparison_table():
        print(f"{name:<30s} {omega:>16.4e} {r_ctc:>16.4e} {rho:>18.4e}")
    print()

    # --- Lab-scale CTC: what omega? ---
    print("=" * 78)
    print("INVERSE PROBLEM: omega for CTC at R = 1 m")
    print("=" * 78)
    R = 1.0
    omega_required = godel_omega_for_ctc_at_R(R)
    rho_required = godel_density(omega_required)
    mass_required = mass_of_uniform_sphere(rho_required, R)
    print(f"  R          = {R} m")
    print(f"  omega_req  = {omega_required:.4e} rad/s")
    print(f"  rho_req    = {rho_required:.4e} kg/m^3")
    print(f"  Mass (V)   = {mass_required:.4e} kg  "
          f"= {mass_required/MEARTH:.4e} Earth masses")
    print(f"             = {mass_required/MSUN:.4e} Solar masses")
    print()

    # --- Compare to known matter densities ---
    print("=" * 78)
    print("KARSILASTIRMA: rho_required vs known densities")
    print("=" * 78)
    refs = {
        "Water": 1e3,
        "Steel": 7.9e3,
        "Earth core": 1.1e4,
        "Sun core": 1.5e5,
        "White dwarf": 1e9,
        "Nuclear": 2.3e17,
        "Neutron star core": 5e17,
        "QGP": 1e18,
        "Planck density": 5.2e96,
    }
    for name, rho_ref in refs.items():
        ratio = rho_required / rho_ref
        print(f"  {name:<20s} rho = {rho_ref:.4e}  "
              f"ratio(required/ref) = {ratio:.4e}")
    print()

    # --- Quantum inequality (Tipler bound analog) ---
    print("=" * 78)
    print("FORD-ROMAN QUANTUM INEQUALITY")
    print("=" * 78)
    print("Godel's classical solution assumes perfect fluid + Lambda.")
    print("With quantum fields, negative energy must be bounded.")
    print()
    print(f"  rho_required (Godel CTC at 1 m) = {rho_required:.4e} kg/m^3")
    print(f"  Equivalent energy density (c^2) = {rho_required*C**2:.4e} J/m^3")
    print(f"  Casimir max (1 nm)              = {np.pi**2*HBAR*C/(720*1e-9**4):.4e} J/m^3")
    gap = rho_required * C**2 / (np.pi**2 * HBAR * C / (720 * 1e-9**4))
    print(f"  Gap                             = {gap:.4e}")
    print()

    # --- Conclusion ---
    print("=" * 78)
    print("SONUC")
    print("=" * 78)
    print("Godel metricinde CTC icin:")
    print(f"  1 m yaricap: omega ~ {omega_required:.2e} rad/s, "
          f"rho ~ {rho_required:.2e} kg/m^3")
    print()
    print(f"Bu, 1 m yaricapli bir kure icinde:")
    print(f"  {mass_required/MSUN:.2e} Gunes kutlesi yogunlugu gerektirir")
    print()
    print("Karsilastirma:")
    print(f"  Proton Compton frekansi       ~ 1.4e23 Hz")
    print(f"  Planck frekansi               ~ 1.9e43 Hz")
    print(f"  Gerekli omega                 ~ {omega_required:.2e} rad/s")
    print()
    if omega_required > 1e43:
        print("Gerekli omega, Planck frekansinin UZERINDE.")
        print("Kuantum gravite rejimi disinda tanimli degil.")
    else:
        print("Gerekli omega, pre-Planck fakat fiziksel olarak ulasilamaz.")
    print()
    print("Godel metrici, klasik GR'de tam cozum olmasina ragmen,")
    print("device-scale'de UYGULANAMAZ.")
    print("=" * 78)

    out = Path("zt007_b6_results.json")
    out.write_text(json.dumps({
        "stage": "ZT-007B.6",
        "mode": a.mode,
        "ctc_coefficient": float(np.log(1+np.sqrt(2))/np.sqrt(2)),
        "omega_required_1m": float(omega_required),
        "rho_required_1m": float(rho_required),
        "mass_required_1m_solar": float(mass_required / MSUN),
    }, indent=2), encoding="utf-8")
    print(f"Results written to {out}")
    print("SCIENTIFIC STATUS: Godel feasibility; no CTC construction.")


if __name__ == "__main__":
    main()
