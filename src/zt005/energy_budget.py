"""Energy Budget Skeleton (Phase 8).

Computes the energy bookkeeping for the device:
- P_in: Input power from the drive
- P_loss: Ohmic losses in conductors
- P_stored: Stored electromagnetic energy
- P_radiated: Radiated power (far-field)
- P_thermal: Thermal dissipation
- η: Efficiency (P_stored / P_in or similar)
- E_1s: Energy required for 1 second operation

This translates the h ~ 10^-46 metric perturbation scale into explicit
energy requirements.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Any

import numpy as np

from .model import Backpack, Drive
from .geometry import build_sources
from .maxwell import vector_potential
from .validated_stress_energy import em_stress_energy

# Physical constants
C = 299_792_458.0
MU0 = 4.0e-7 * np.pi
EPS0 = 8.8541878128e-12
Z0 = np.sqrt(MU0 / EPS0)  # Vacuum impedance ~ 377 ohms


@dataclass(frozen=True)
class ConductorParams:
    """Physical conductor parameters."""
    resistivity_ohm_m: float = 1.68e-8  # Copper at 20°C
    cross_section_m2: float = 1.0e-6    # 1 mm^2 wire cross-section
    length_m: float = 1.0               # Total wire length (per coil)


@dataclass(frozen=True)
class EnergyBudget:
    """Energy budget results."""
    P_in: float           # Input power [W]
    P_loss: float         # Ohmic losses [W]
    P_stored: float       # Stored EM energy [J]
    P_radiated: float     # Radiated power [W]
    P_thermal: float      # Thermal dissipation [W]
    eta: float            # Efficiency (dimensionless)
    E_1s: float           # Energy for 1s operation [J]
    h_current: float      # Current metric perturbation amplitude
    h_target: float       # Target metric perturbation amplitude
    h_gap_orders: float   # Orders of magnitude gap


def compute_coil_resistance(params: ConductorParams) -> float:
    """Compute DC resistance of a coil."""
    return params.resistivity_ohm_m * params.length_m / params.cross_section_m2


def compute_ohmic_loss(I_rms: float, R: float) -> float:
    """Compute ohmic loss P = I_rms^2 * R."""
    return I_rms**2 * R


def compute_stored_energy(T_grid: np.ndarray, dx_m: float) -> float:
    """Compute total stored EM energy from T_00 integrated over volume.
    
    T_00 is energy density [J/m^3]. Integrate over spatial volume.
    """
    # T_grid shape: (Nt, Nx, Ny, Nz, 4, 4)
    # Take time-average of T_00
    T00 = T_grid[..., 0, 0]  # (Nt, Nx, Ny, Nz)
    T00_avg = np.mean(T00, axis=0)  # (Nx, Ny, Nz)
    volume_element = dx_m**3
    return float(np.sum(T00_avg) * volume_element)


def compute_radiated_power(T_grid: np.ndarray, dx_m: float, r_obs: float) -> float:
    """Estimate radiated power from Poynting flux at observation radius.
    
    Uses T_0i components (energy flux) at a spherical shell.
    """
    # Simplified: use far-field approximation
    # For MQS, radiation is negligible but we compute it for completeness
    T0i = T_grid[..., 0, 1:]  # (Nt, Nx, Ny, Nz, 3)
    # This is a placeholder - full implementation needs spherical surface integral
    return 0.0


def compute_input_power(drive: Drive, R: float) -> float:
    """Compute input power from drive.
    
    For sinusoidal drive: P_in = I_rms^2 * R + reactive power
    """
    I_rms = drive.peak_current_a / np.sqrt(2)
    return I_rms**2 * R


def compute_energy_budget(
    config_name: str,
    backpack: Backpack,
    drive: Drive,
    conductor: ConductorParams,
    grid_dx_m: float = 0.05,
    grid_radius_m: float = 0.5,
    h_target: float = 1e-3,
) -> EnergyBudget:
    """Compute full energy budget for a device configuration."""
    
    # Build sources and compute fields
    sources = build_sources(config_name, backpack)
    
    # Create a grid for sampling
    n = int(2 * grid_radius_m / grid_dx_m) + 1
    xs = np.linspace(-grid_radius_m, grid_radius_m, n)
    ys = np.linspace(-grid_radius_m, grid_radius_m, n)
    zs = np.linspace(-grid_radius_m, grid_radius_m, n)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
    
    # Compute vector potential at t=0 (peak)
    A = vector_potential(pts, sources, drive.peak_current_a)
    
    # Compute B = curl A (using finite differences)
    # Reshape A to grid
    A_grid = A.reshape(n, n, n, 3)
    
    # Central differences for curl
    h = grid_dx_m
    dAx_dy = (np.roll(A_grid, -1, axis=1) - np.roll(A_grid, 1, axis=1))[:, :, :, 0] / (2*h)
    dAx_dz = (np.roll(A_grid, -1, axis=2) - np.roll(A_grid, 1, axis=2))[:, :, :, 0] / (2*h)
    dAy_dx = (np.roll(A_grid, -1, axis=0) - np.roll(A_grid, 1, axis=0))[:, :, :, 1] / (2*h)
    dAy_dz = (np.roll(A_grid, -1, axis=2) - np.roll(A_grid, 1, axis=2))[:, :, :, 1] / (2*h)
    dAz_dx = (np.roll(A_grid, -1, axis=0) - np.roll(A_grid, 1, axis=0))[:, :, :, 2] / (2*h)
    dAz_dy = (np.roll(A_grid, -1, axis=1) - np.roll(A_grid, 1, axis=1))[:, :, :, 2] / (2*h)
    
    Bx = dAz_dy - dAy_dz
    By = dAx_dz - dAz_dx
    Bz = dAy_dx - dAx_dy
    B = np.stack([Bx, By, Bz], axis=-1)
    
    # E = -dA/dt = omega * A (for sinusoidal, at t=0)
    omega = 2 * np.pi * drive.frequency_hz
    E = omega * A_grid
    
    # Compute stress-energy tensor - flatten grids first
    E_flat = E.reshape(-1, 3)
    B_flat = B.reshape(-1, 3)
    T_dict = em_stress_energy(E_flat, B_flat)
    T = T_dict["T"]  # (Nx*Ny*Nz, 4, 4)
    
    # Reshape back to grid
    T_grid = T.reshape(n, n, n, 4, 4)
    # Add time dimension for stored energy computation
    T_grid = T_grid[np.newaxis, ...]  # (1, Nx, Ny, Nz, 4, 4)
    
    # Compute resistance and losses
    R = compute_coil_resistance(conductor)
    I_rms = drive.peak_current_a / np.sqrt(2)
    P_loss = compute_ohmic_loss(I_rms, R)
    P_in = compute_input_power(drive, R)
    
    # Stored energy
    P_stored = compute_stored_energy(T_grid, grid_dx_m)
    
    # Radiated power (placeholder for MQS)
    P_radiated = compute_radiated_power(T_grid, grid_dx_m, grid_radius_m)
    
    # Thermal = losses (steady state)
    P_thermal = P_loss
    
    # Efficiency
    eta = P_stored / P_in if P_in > 0 else 0.0
    
    # Energy for 1 second
    E_1s = P_in * 1.0
    
    # Current h estimate from linearized Einstein equation
    # h ~ (G/c^4) * (E/L) * L^2 = G * E * L / c^4
    L = (2 * grid_radius_m)
    h_current = (6.67430e-11 * P_stored * L) / (C**4)
    
    # Orders of magnitude gap
    h_gap_orders = math.log10(h_target / h_current) if h_current > 0 else float('inf')
    
    return EnergyBudget(
        P_in=P_in,
        P_loss=P_loss,
        P_stored=P_stored,
        P_radiated=P_radiated,
        P_thermal=P_thermal,
        eta=eta,
        E_1s=E_1s,
        h_current=h_current,
        h_target=h_target,
        h_gap_orders=h_gap_orders,
    )


def run_energy_budget_all_configs(
    out_path: str = "results/validated/energy_budget.json",
    grid_dx_m: float = 0.05,
    grid_radius_m: float = 0.5,
) -> Dict[str, Any]:
    """Run energy budget for all device configurations."""
    
    backpack = Backpack()
    drive = Drive()
    conductor = ConductorParams()
    
    configs = ["G1", "G2", "G3", "G4"]
    results = {}
    
    for config in configs:
        print(f"Computing energy budget for {config}...")
        budget = compute_energy_budget(
            config, backpack, drive, conductor,
            grid_dx_m=grid_dx_m,
            grid_radius_m=grid_radius_m,
        )
        results[config] = asdict(budget)
        print(f"  P_in={budget.P_in:.3e} W, P_loss={budget.P_loss:.3e} W, "
              f"P_stored={budget.P_stored:.3e} J, eta={budget.eta:.3e}, "
              f"h_current={budget.h_current:.3e}, gap={budget.h_gap_orders:.1f} orders")
    
    # Summary
    summary = {
        "stage": "energy-budget-v1",
        "constants": {
            "c_m_per_s": C,
            "G_m3_kg_s2": 6.67430e-11,
            "MU0": MU0,
            "EPS0": EPS0,
            "Z0_ohm": Z0,
        },
        "conductor": asdict(conductor),
        "drive": asdict(drive),
        "backpack": asdict(backpack),
        "grid": {
            "dx_m": grid_dx_m,
            "radius_m": grid_radius_m,
        },
        "configs": results,
    }
    
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2))
    print(f"\nResults written to {out}")
    
    return summary


if __name__ == "__main__":
    run_energy_budget_all_configs()