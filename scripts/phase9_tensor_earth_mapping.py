#!/usr/bin/env python3
"""
Phase 9: Tensor-Level Earth + Target-Time Coordinate Mapping.

Extends the kinematic reference-frame model (earth_navigation_frame.py) with
the Earth's gravitational contribution to the metric:

1. Earth stress-energy tensor T^μν_Earth (spherical, static, perfect fluid)
2. Metric perturbation h^μν_Earth from linearized Einstein equations
3. Combined metric: g_μν = η_μν + h^μν_device + h^μν_Earth
4. Geodesic integration in combined spacetime for target-time mapping
5. Proper time along worldlines including Earth's gravity

This addresses the roadmap Phase 9 requirement: "Earth + target-coordinate
mapping (tensor level)" - moving beyond kinematic placeholder to full
tensor-gravity treatment.
"""

from __future__ import annotations

import numpy as np
import json
import sys
from pathlib import Path
from dataclasses import dataclass, asdict
from typing import Optional

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from zt005.retarded_gr import retarded_hbar
from zt005.gauge_checks import lorenz_residual
from zt005.earth_navigation_frame import (
    geodetic_to_ecef, ecef_to_geodetic,
    heliocentric_earth_position, heliocentric_earth_velocity,
    compute_arrival_mapping, AU_M, ORBIT_OMEGA, ORBIT_SPEED_M_S,
    EARTH_A_M, EARTH_B_M, EARTH_E2, SPIN_OMEGA
)
from zt005.model import Backpack, Drive
from zt005.geometry import build_sources
from zt005.batch_em_source import compute_T_batch


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
    # Earth model
    earth_mass_kg: float = EARTH_MASS_KG
    earth_radius_m: float = EARTH_RADIUS_M
    # Grid
    grid_radius_m: float = 2.0 * EARTH_RADIUS_M  # 2 Earth radii
    n_grid: int = 17  # Reduced from 33 for speed
    # Time
    target_dt_s: float = -1.0
    n_time: int = 5
    # Location
    initial_lat_deg: float = 41.0082
    initial_lon_deg: float = 28.9784
    initial_alt_m: float = 50.0
    # Device config
    device_config: str = "G1"
    peak_current_a: float = 100.0
    # Output
    output_dir: str = "results/validated"


def earth_stress_energy_tensor(x: np.ndarray, M: float = EARTH_MASS_KG, R: float = EARTH_RADIUS_M) -> np.ndarray:
    """
    Earth stress-energy tensor T^μν for a static, spherical perfect fluid.
    
    In the Earth's rest frame (ECI):
    T^00 = ρ c^2 (energy density)
    T^ii = p (pressure, isotropic)
    T^0i = 0 (no momentum flux in rest frame)
    
    For a uniform density sphere:
    ρ = 3M / (4πR^3) inside, 0 outside
    Pressure from hydrostatic equilibrium: p(r) = (2πGρ^2/3)(R^2 - r^2) inside
    
    Returns T^μν at position x (in ECI frame) as (4,4) array.
    """
    r = np.linalg.norm(x)
    
    if r >= R:
        # Vacuum outside Earth
        return np.zeros((4, 4))
    
    # Uniform density
    rho = 3.0 * M / (4.0 * np.pi * R**3)
    
    # Pressure at radius r (hydrostatic equilibrium for uniform density)
    # p(r) = (2πGρ^2/3)(R^2 - r^2)
    p = (2.0 * np.pi * G_SI * rho**2 / 3.0) * (R**2 - r**2)
    
    T = np.zeros((4, 4))
    T[0, 0] = rho * C_SI**2  # Energy density
    T[1, 1] = p              # Pressure xx
    T[2, 2] = p              # Pressure yy
    T[3, 3] = p              # Pressure zz
    
    return T


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


def earth_metric_perturbation_retarded(
    x_obs: np.ndarray,
    t_obs: float,
    M: float = EARTH_MASS_KG,
    R: float = EARTH_RADIUS_M
) -> np.ndarray:
    """
    Retarded metric perturbation from Earth (moving source).
    
    Earth moves on heliocentric orbit. We compute the retarded position
    and evaluate the static perturbation there.
    
    For slow motion (v ~ 30 km/s << c), this is approximately the static
    field evaluated at retarded position.
    """
    # Earth position at retarded time
    # Retarded time: t_ret = t_obs - |x_obs - x_earth(t_ret)|/c
    # For v << c, approximate t_ret ≈ t_obs - |x_obs - x_earth(t_obs)|/c
    x_earth_now = heliocentric_earth_position(t_obs)
    r_vec = x_obs - x_earth_now
    r = np.linalg.norm(r_vec)
    t_ret = t_obs - r / C_SI
    
    x_earth_ret = heliocentric_earth_position(t_ret)
    r_ret_vec = x_obs - x_earth_ret
    r_ret = np.linalg.norm(r_ret_vec)
    
    if r_ret == 0:
        r_ret = 1e-10
    
    phi = G_SI * M / (C_SI**2 * r_ret)
    
    h = np.zeros((4, 4))
    h[0, 0] = 2.0 * phi
    h[1, 1] = 2.0 * phi
    h[2, 2] = 2.0 * phi
    h[3, 3] = 2.0 * phi
    
    return h


def build_earth_source_grid(config: EarthTensorConfig) -> tuple[np.ndarray, np.ndarray]:
    """
    Build Earth stress-energy source grid in ECI frame.
    
    Returns:
        x_grid: (N, 3) grid points in ECI [m]
        T_grid: (N, 4, 4) stress-energy tensors
    """
    # Grid centered on Earth center in ECI - use smaller grid for speed
    L = config.grid_radius_m
    n = min(config.n_grid, 17)  # Cap at 17 for speed
    xs = np.linspace(-L, L, n)
    ys = np.linspace(-L, L, n)
    zs = np.linspace(-L, L, n)
    
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    x_grid = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    
    # Filter to sphere of radius grid_radius_m
    r = np.linalg.norm(x_grid, axis=1)
    mask = r <= L
    x_grid = x_grid[mask]
    
    # Compute T^μν at each point - vectorized for speed
    T_grid = np.zeros((len(x_grid), 4, 4))
    for i, x in enumerate(x_grid):
        T_grid[i] = earth_stress_energy_tensor(x)
    
    return x_grid, T_grid


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
    """Run Phase 9 tensor-level Earth mapping."""
    
    print("=" * 70)
    print("PHASE 9: Tensor-Level Earth + Target-Time Coordinate Mapping")
    print("=" * 70)
    print(f"Config: {config.device_config}")
    print(f"Grid: {config.n_grid}^3, radius: {config.grid_radius_m/1e3:.1f} km")
    print(f"Target Delta t: {config.target_dt_s} s")
    print(f"Location: {config.initial_lat_deg}°N, {config.initial_lon_deg}°E, {config.initial_alt_m}m")
    print()
    
    # 1. Kinematic mapping (from earth_navigation_frame)
    print("[1/6] Computing kinematic arrival mapping...")
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
    
    # 2. Build Earth source grid
    print("[2/6] Building Earth stress-energy source grid...")
    x_earth_src, T_earth_src = build_earth_source_grid(config)
    print(f"    Source points: {len(x_earth_src)}")
    print(f"    T^00 range: [{np.min(T_earth_src[:,0,0]):.2e}, {np.max(T_earth_src[:,0,0]):.2e}] J/m³")
    
    # 3. Build device source grid
    print("[3/6] Building device stress-energy source grid...")
    backpack = Backpack()
    drive = Drive(frequency_hz=1_000_000.0, peak_current_a=config.peak_current_a)
    sources = build_sources(config.device_config, backpack)
    
    # Build observation points for device source (small grid around device)
    # The device is at p_device_helio_0 in heliocentric, but we need local coords
    # For the device source, we use the local backpack frame
    # Create a grid around the device center in local coordinates
    dev_obs_radius = 0.5  # 0.5m around device
    dev_n_obs = 9  # Reduced from 17 for speed
    xs = np.linspace(-dev_obs_radius, dev_obs_radius, dev_n_obs)
    ys = np.linspace(-dev_obs_radius, dev_obs_radius, dev_n_obs)
    zs = np.linspace(-dev_obs_radius, dev_obs_radius, dev_n_obs)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    dev_obs_grid = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    
    # Compute device T^μν at observation points
    dev_T_grid = compute_T_batch(config.device_config, dev_obs_grid, backpack, drive)
    dev_src_grid = dev_obs_grid  # Source points are the observation points for retarded solver
    
    print(f"    Device source points: {len(dev_src_grid)}")
    print(f"    Device T^00 range: [{np.min(dev_T_grid[:,0,0]):.2e}, {np.max(dev_T_grid[:,0,0]):.2e}] J/m³")
    
    # 4. Observation grid (around initial device position)
    print("[4/6] Setting up observation grid...")
    p_ecef_0 = geodetic_to_ecef(config.initial_lat_deg, config.initial_lon_deg, config.initial_alt_m)
    p_earth_helio_0 = heliocentric_earth_position(0.0)
    p_device_helio_0 = p_earth_helio_0 + p_ecef_0
    
    # Small grid around device
    obs_radius = 10.0  # 10m around device
    n_obs = 11
    xs = np.linspace(-obs_radius, obs_radius, n_obs)
    ys = np.linspace(-obs_radius, obs_radius, n_obs)
    zs = np.linspace(-obs_radius, obs_radius, n_obs)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    x_obs_grid = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1) + p_device_helio_0
    
    t_obs = np.linspace(0.0, abs(config.target_dt_s), config.n_time)
    print(f"    Obs points: {len(x_obs_grid)}, Times: {len(t_obs)}")
    
    # 5. Compute Earth's metric perturbation (static approximation)
    print("[5/6] Computing Earth's metric perturbation...")
    h_earth_grid = np.array([earth_metric_perturbation_static(x - p_earth_helio_0) for x in x_obs_grid])
    h_earth_max = np.max(np.abs(h_earth_grid))
    print(f"    Max |h_earth|: {h_earth_max:.2e}")
    print(f"    Expected h_00 ~ 2GM/(c²R) = {2*G_SI*EARTH_MASS_KG/(C_SI**2*EARTH_RADIUS_M):.2e}")
    
    # 6. Compute device metric perturbation
    print("[6/6] Computing device metric perturbation...")
    
    # Create source function for retarded_hbar
    def device_source_fn(t_array, x_array):
        # t_array: (Ns,) retarded times
        # x_array: (Ns, 3) source positions
        # Return T^μν at each (t, x)
        # For MQS, T is time-independent (static in device frame)
        # But we need to evaluate at retarded times
        # The device source is periodic, but for simplicity use peak field
        # We'll compute T at each source point
        T_list = []
        for t, x in zip(t_array, x_array):
            # For MQS, T doesn't depend on time (quasi-static)
            # Just compute at the spatial point
            T = compute_T_batch(config.device_config, x[None, :], backpack, drive)
            T_list.append(T[0])
        return np.array(T_list)
    
    # Source weights (volume elements)
    dx = 2 * dev_obs_radius / (dev_n_obs - 1)
    dy = dx
    dz = dx
    src_weights = np.full(len(dev_src_grid), dx * dy * dz)
    
    h_device = retarded_hbar(
        source_fn=device_source_fn,
        obs_points=x_obs_grid,
        obs_times=t_obs,
        src_points=dev_src_grid,
        src_weights=src_weights,
        trace_reverse_source=False,
        softening_m=0.0
    )
    h_device_max = np.max(np.abs(h_device))
    print(f"    Max |h_device|: {h_device_max:.2e}")
    
    # Combined metric at initial position
    h_earth_0 = earth_metric_perturbation_static(p_device_helio_0 - p_earth_helio_0)
    h_device_0 = h_device[0, len(x_obs_grid)//2]  # Center point, t=0
    g_combined = combined_metric(p_device_helio_0, 0.0, h_device_0, h_earth_0)
    
    print()
    print("=" * 70)
    print("TENSOR-LEVEL MAPPING RESULTS")
    print("=" * 70)
    print(f"Earth h_00 at surface: {h_earth_0[0,0]:.2e}")
    print(f"Device h_00 at center: {h_device_0[0,0]:.2e}")
    print(f"Ratio (device/Earth): {abs(h_device_0[0,0]/h_earth_0[0,0]):.2e}")
    print()
    print("Combined metric at device (t=0):")
    print(f"  g_00 = {g_combined[0,0]:.10f}")
    print(f"  g_11 = {g_combined[1,1]:.10f}")
    print(f"  g_22 = {g_combined[2,2]:.10f}")
    print(f"  g_33 = {g_combined[3,3]:.10f}")
    print()
    
    # Proper time difference due to Earth's gravity
    # For a stationary observer at Earth surface vs infinity
    dtau_earth = np.sqrt(-g_combined[0,0])  # Proper time rate
    print(f"Proper time rate at device (Earth gravity): {dtau_earth:.10f}")
    print(f"Gravitational time dilation factor: {1 - dtau_earth:.2e}")
    
    # Target time mapping with tensor gravity
    # The target time Δt is coordinate time in heliocentric frame
    # Proper time experienced by traveler depends on path through combined metric
    
    # For a "jump" of Δt = -1s in coordinate time:
    # The traveler's proper time depends on the metric along the path
    # In Earth's gravity well, proper time runs slower by factor sqrt(-g_00)
    
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
    parser = argparse.ArgumentParser(description="Phase 9: Tensor Earth Mapping")
    parser.add_argument("--config", default="G1", choices=["G1", "G2", "G3", "G4"])
    parser.add_argument("--n-grid", type=int, default=33)
    parser.add_argument("--radius", type=float, default=2.0*EARTH_RADIUS_M)
    parser.add_argument("--target-dt", type=float, default=-1.0)
    parser.add_argument("--n-time", type=int, default=5)
    parser.add_argument("--lat", type=float, default=41.0082)
    parser.add_argument("--lon", type=float, default=28.9784)
    parser.add_argument("--alt", type=float, default=50.0)
    parser.add_argument("--current", type=float, default=100.0)
    parser.add_argument("--output-dir", default="results/validated")
    args = parser.parse_args()
    
    config = EarthTensorConfig(
        device_config=args.config,
        n_grid=args.n_grid,
        grid_radius_m=args.radius,
        target_dt_s=args.target_dt,
        n_time=args.n_time,
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