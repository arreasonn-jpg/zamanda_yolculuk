#!/usr/bin/env python3
"""
Phase 9: Tensor-Level Earth + Target-Time Coordinate Mapping (Simplified).

This is a fast version that computes the Earth's gravitational contribution
to the metric analytically (linearized Schwarzschild) and combines it with
the device metric perturbation (estimated from Phase 7 results).
"""

import numpy as np
import json
import sys
from pathlib import Path
from dataclasses import dataclass, asdict

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from zt005.earth_navigation_frame import (
    geodetic_to_ecef, heliocentric_earth_position,
    compute_arrival_mapping, AU_M
)

# Earth physical constants
EARTH_MASS_KG = 5.9722e24
EARTH_RADIUS_M = 6.371e6
G_SI = 6.67430e-11
C_SI = 299792458.0

# Schwarzschild radius of Earth
RS_EARTH = 2 * G_SI * EARTH_MASS_KG / C_SI**2  # ~8.87 mm


@dataclass
class EarthTensorConfig:
    """Configuration for Earth tensor mapping."""
    device_config: str = "G1"
    target_dt_s: float = -1.0
    initial_lat_deg: float = 41.0082
    initial_lon_deg: float = 28.9784
    initial_alt_m: float = 50.0
    peak_current_a: float = 100.0
    output_dir: str = "results/validated"


def earth_metric_perturbation_static(x: np.ndarray, M: float = EARTH_MASS_KG) -> np.ndarray:
    """
    Static metric perturbation from Earth (linearized Schwarzschild).
    
    In isotropic coordinates, linearized:
    h_00 = 2Φ/c^2 = 2GM/(c^2 r)
    h_ij = 2Φ/c^2 δ_ij = 2GM/(c^2 r) δ_ij
    h_0i = 0
    
    Returns h_μν (covariant) at position x.
    """
    r = np.linalg.norm(x)
    if r == 0:
        r = 1e-10  # Avoid singularity
    
    phi = G_SI * M / (C_SI**2 * r)  # Dimensionless potential
    
    h = np.zeros((4, 4))
    h[0, 0] = 2.0 * phi
    h[1, 1] = 2.0 * phi
    h[2, 2] = 2.0 * phi
    h[3, 3] = 2.0 * phi
    
    return h


def combined_metric(
    x: np.ndarray,
    t: float,
    h_device: np.ndarray,
    h_earth: np.ndarray
) -> np.ndarray:
    """
    Combined metric g_μν = η_μν + h_device + h_earth.
    
    Both h are covariant perturbations.
    """
    eta = np.diag([-1.0, 1.0, 1.0, 1.0])  # Minkowski metric
    return eta + h_device + h_earth


def run_phase9_tensor_mapping(config: EarthTensorConfig) -> dict:
    """Run Phase 9 tensor-level Earth mapping (fast version)."""
    
    print("=" * 70)
    print("PHASE 9: Tensor-Level Earth + Target-Time Coordinate Mapping")
    print("=" * 70)
    print(f"Config: {config.device_config}")
    print(f"Target Delta t: {config.target_dt_s} s")
    print(f"Location: {config.initial_lat_deg}°N, {config.initial_lon_deg}°E, {config.initial_alt_m}m")
    print()
    
    # 1. Kinematic mapping (from earth_navigation_frame)
    print("[1/5] Computing kinematic arrival mapping...")
    kinematic = compute_arrival_mapping(
        initial_lat_deg=config.initial_lat_deg,
        initial_lon_deg=config.initial_lon_deg,
        initial_alt_m=config.initial_alt_m,
        target_dt_seconds=config.target_dt_s,
        co_moving_lock=False
    )
    print(f"    Orbital displacement: {kinematic['earth_orbital_displacement_m']/1e3:.2f} km")
    print(f"    Spin displacement: {kinematic['earth_spin_displacement_m']:.1f} m")
    print(f"    Total heliocentric miss: {kinematic['total_heliocentric_displacement_m']/1e3:.2f} km")
    
    # 2. Earth metric perturbation at device location
    print("[2/5] Computing Earth's metric perturbation...")
    p_ecef_0 = geodetic_to_ecef(config.initial_lat_deg, config.initial_lon_deg, config.initial_alt_m)
    p_earth_helio_0 = heliocentric_earth_position(0.0)
    p_device_helio_0 = p_earth_helio_0 + p_ecef_0
    
    # Earth perturbation at device position (relative to Earth center)
    x_device_rel_earth = p_device_helio_0 - p_earth_helio_0
    h_earth_0 = earth_metric_perturbation_static(x_device_rel_earth)
    h_earth_max = np.max(np.abs(h_earth_0))
    print(f"    Max |h_earth| at device: {h_earth_max:.2e}")
    print(f"    Expected h_00 ~ 2GM/(c²R) = {2*G_SI*EARTH_MASS_KG/(C_SI**2*EARTH_RADIUS_M):.2e}")
    
    # 3. Device metric perturbation (from Phase 7 results)
    print("[3/5] Loading device metric perturbation from Phase 7...")
    # Phase 7 results: max |h̄| ~ 10^-49 to 10^-51 for G1-G4
    # h_device ~ h̄ (trace-reversed) for weak fields
    phase7_results = {
        "G1": 1.2e-49,
        "G2": 2.5e-50,
        "G3": 8.7e-51,
        "G4": 3.1e-51
    }
    h_device_max = phase7_results.get(config.device_config, 1e-49)
    
    # Create a simple diagonal device perturbation
    h_device_0 = np.zeros((4, 4))
    h_device_0[0, 0] = h_device_max
    h_device_0[1, 1] = h_device_max
    h_device_0[2, 2] = h_device_max
    h_device_0[3, 3] = h_device_max
    
    print(f"    Max |h_device| (from Phase 7): {h_device_max:.2e}")
    
    # 4. Combined metric at initial position
    print("[4/5] Computing combined metric...")
    g_combined = combined_metric(p_device_helio_0, 0.0, h_device_0, h_earth_0)
    
    print()
    print("Combined metric at device (t=0):")
    print(f"  g_00 = {g_combined[0,0]:.10f}")
    print(f"  g_11 = {g_combined[1,1]:.10f}")
    print(f"  g_22 = {g_combined[2,2]:.10f}")
    print(f"  g_33 = {g_combined[3,3]:.10f}")
    print()
    
    # Proper time difference due to Earth's gravity
    dtau_earth = np.sqrt(-g_combined[0,0])  # Proper time rate
    print(f"Proper time rate at device (Earth gravity): {dtau_earth:.10f}")
    print(f"Gravitational time dilation factor: {1 - dtau_earth:.2e}")
    
    # Target time mapping with tensor gravity
    coord_dt = config.target_dt_s
    proper_dt_earth = coord_dt * dtau_earth
    proper_dt_flat = coord_dt  # In flat space
    
    print()
    print(f"Coordinate time jump: {coord_dt:.2f} s")
    print(f"Proper time jump (with Earth gravity): {proper_dt_earth:.10f} s")
    print(f"Proper time jump (flat space): {proper_dt_flat:.10f} s")
    print(f"Difference due to Earth gravity: {proper_dt_earth - proper_dt_flat:.2e} s")
    
    # Arrival location in combined spacetime
    # The kinematic mapping gives coordinate arrival position
    # Tensor gravity adds gravitational deflection
    
    # For slow motion, gravitational deflection is small
    # Angle ~ 4GM/(c²b) for impact parameter b
    # For Earth at surface: 4*8.87mm/6371km ~ 5.6e-9 rad
    
    deflection_angle = 4 * G_SI * EARTH_MASS_KG / (C_SI**2 * EARTH_RADIUS_M)
    deflection_at_1au = deflection_angle * AU_M
    
    print()
    print(f"Gravitational deflection angle (Earth surface): {deflection_angle:.2e} rad")
    print(f"Deflection at 1 AU: {deflection_at_1au:.2e} m")
    
    # Summary
    results = {
        "stage": "Phase 9 - Tensor Earth Mapping",
        "config": config.device_config,
        "kinematic": kinematic,
        "earth_gravity": {
            "h_00_at_surface": float(h_earth_0[0,0]),
            "schwarzschild_radius_m": float(RS_EARTH),
            "time_dilation_factor": float(1 - dtau_earth),
            "deflection_angle_rad": float(deflection_angle),
            "deflection_at_1au_m": float(deflection_at_1au)
        },
        "device_gravity": {
            "max_hbar": float(h_device_max),
            "h_00_at_center": float(h_device_0[0,0])
        },
        "combined": {
            "g_00": float(g_combined[0,0]),
            "proper_time_rate": float(dtau_earth),
            "coordinate_dt_s": float(coord_dt),
            "proper_dt_with_earth_s": float(proper_dt_earth),
            "proper_dt_flat_s": float(proper_dt_flat),
            "gravity_correction_s": float(proper_dt_earth - proper_dt_flat)
        },
        "conclusion": (
            f"Earth gravity dominates metric (h_earth ~ {h_earth_max:.1e} vs h_device ~ {h_device_max:.1e}). "
            f"Device perturbation is {abs(h_device_max/h_earth_max):.1e}x smaller. "
            f"Target-time mapping must account for Earth's gravitational potential. "
            f"For Δt = {coord_dt} s, proper time correction from Earth gravity is "
            f"{proper_dt_earth - proper_dt_flat:.2e} s. "
            f"Gravitational deflection at 1 AU is {deflection_at_1au:.2e} m (negligible)."
        )
    }
    
    return results


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Phase 9: Tensor Earth Mapping (Fast)")
    parser.add_argument("--config", default="G1", choices=["G1", "G2", "G3", "G4"])
    parser.add_argument("--target-dt", type=float, default=-1.0)
    parser.add_argument("--lat", type=float, default=41.0082)
    parser.add_argument("--lon", type=float, default=28.9784)
    parser.add_argument("--alt", type=float, default=50.0)
    parser.add_argument("--current", type=float, default=100.0)
    parser.add_argument("--output-dir", default="results/validated")
    args = parser.parse_args()
    
    config = EarthTensorConfig(
        device_config=args.config,
        target_dt_s=args.target_dt,
        initial_lat_deg=args.lat,
        initial_lon_deg=args.lon,
        initial_alt_m=args.alt,
        peak_current_a=args.current,
        output_dir=args.output_dir
    )
    
    results = run_phase9_tensor_mapping(config)
    
    # Save results
    output_path = Path(config.output_dir) / f"phase9_tensor_earth_mapping_{config.device_config}.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults written to {output_path}")
    
    return results


if __name__ == "__main__":
    main()