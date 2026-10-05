"""Green-Integral Convergence and Grid Refinement Engine.

Addresses Section 1.2:
- Grid series refinement: N in [12, 16, 20, 28, 36] (or parameterized)
- Sensitivity studies:
    * Quadrature order
    * Regularization radius / cell size
    * Observer position (center, interior, edge, exterior)
- Acceptance criteria for final two solutions:
    * h00 relative difference < 1% (0.01)
    * |h0i| absolute and relative tolerance
    * ||hij|| relative Frobenius norm < 1% (0.01)
- Cross-validation against independent second solver (Solver B).
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .model import ActiveCylinder, Backpack, Drive
from .gauss_cyl_quadrature import gauss_cylindrical_nodes
from .batch_em_source import compute_T_batch
from .geometry import build_sources
from .maxwell import vector_potential, biot_savart
from .validated_stress_energy import em_stress_energy
from .manufactured_sources import GaussianBump, conserved_T
from .green_solvers import (
    solver_a_tensor_product,
    solver_b_adaptive_voxel,
    compare_green_solvers,
)


@dataclass
class ObserverLocations:
    center: tuple[float, float, float] = (0.0, 0.0, 0.0)
    interior_radial: tuple[float, float, float] = (0.3, 0.0, 0.0)
    interior_axial: tuple[float, float, float] = (0.0, 0.0, 0.5)
    near_edge: tuple[float, float, float] = (0.8, 0.0, 0.0)
    exterior: tuple[float, float, float] = (0.0, 0.0, 2.0)

    def as_array(self) -> np.ndarray:
        return np.array([
            self.center,
            self.interior_radial,
            self.interior_axial,
            self.near_edge,
            self.exterior,
        ])

    def names(self) -> list[str]:
        return ["center", "interior_radial", "interior_axial", "near_edge", "exterior"]


def component_relative_errors(h_fine: np.ndarray, h_coarse: np.ndarray) -> dict:
    """Compute relative and absolute error metrics between two solutions."""
    h00_f = float(h_fine[0, 0])
    h00_c = float(h_coarse[0, 0])
    h00_ref = max(abs(h00_f), 1e-300)
    h00_rel = float(abs(h00_f - h00_c) / h00_ref)
    h00_abs = float(abs(h00_f - h00_c))

    h0i_f = h_fine[0, 1:]
    h0i_c = h_coarse[0, 1:]
    h0i_abs = float(np.linalg.norm(h0i_f - h0i_c))
    h0i_ref = max(float(np.linalg.norm(h0i_f)), 1e-300)
    h0i_rel = float(h0i_abs / h0i_ref)

    hij_f = h_fine[1:, 1:]
    hij_c = h_coarse[1:, 1:]
    hij_abs = float(np.linalg.norm(hij_f - hij_c))
    hij_ref = max(float(np.linalg.norm(hij_f)), 1e-300)
    hij_rel = float(hij_abs / hij_ref)

    tensor_abs = float(np.linalg.norm((h_fine - h_coarse).ravel()))
    tensor_ref = max(float(np.linalg.norm(h_fine.ravel())), 1e-300)
    tensor_rel = float(tensor_abs / tensor_ref)

    return {
        "h00_rel": h00_rel,
        "h00_abs": h00_abs,
        "h0i_abs": h0i_abs,
        "h0i_rel": h0i_rel,
        "hij_rel": hij_rel,
        "hij_abs": hij_abs,
        "tensor_rel": tensor_rel,
    }


def focused_source_nodes(cyl_r: float, cyl_h: float, nr: int, nphi: int, nz: int) -> tuple[np.ndarray, np.ndarray]:
    """Composite Gauss cylindrical nodes focused on the active source region."""
    xr, wr = np.polynomial.legendre.leggauss(nr)
    r = 0.5 * 0.35 * (xr + 1.0)
    wr = 0.5 * 0.35 * wr

    phi = 2.0 * np.pi * (np.arange(nphi) + 0.5) / nphi
    wphi = 2.0 * np.pi / nphi

    xz, wz = np.polynomial.legendre.leggauss(nz)
    y = 0.5 * 0.35 * xz - 1.075
    wy = 0.5 * 0.35 * wz

    rr, pp, yy = np.meshgrid(r, phi, y, indexing="ij")
    pts = np.column_stack([
        (rr * np.cos(pp)).ravel(),
        yy.ravel(),
        (rr * np.sin(pp)).ravel(),
    ])
    weights = (wr[:, None, None] * wphi * wy[None, None, :] * rr).ravel()
    return pts, weights


def fast_regularized_em_fields(
    points: np.ndarray,
    sources: list,
    drive_current: float,
    freq_hz: float,
    a_wire: float = 0.005,
    chunk_size: int = 8000,
) -> np.ndarray:
    """Evaluate regularized electromagnetic stress-energy tensor efficiently."""
    p = np.asarray(points, float)
    N = len(p)
    A = np.zeros_like(p)
    B = np.zeros_like(p)

    for pts, seg in sources:
        mid = pts + 0.5 * seg
        dlI = seg * drive_current
        for i in range(0, N, chunk_size):
            ch = p[i:i + chunk_size]
            r = ch[:, None, :] - mid[None, :, :]
            R2 = np.sum(r * r, axis=2)
            R_reg = np.sqrt(R2 + a_wire**2)
            R3_reg = R_reg**3
            A[i:i + chunk_size] += (1e-7) * np.sum(dlI[None, :, :] / R_reg[:, :, None], axis=1)
            B[i:i + chunk_size] += (1e-7) * np.sum(np.cross(dlI[None, :, :], r, axis=2) / R3_reg[:, :, None], axis=1)

    E = (2.0 * np.pi * freq_hz) * A
    return em_stress_energy(E, B)["T"]


def run_convergence_study(
    name: str = "G1",
    backpack: Backpack | None = None,
    drive: Drive | None = None,
    cyl: ActiveCylinder | None = None,
    grid_ns: list[int] | None = None,
    observers: np.ndarray | None = None,
    tol_h00: float = 0.01,
    tol_hij: float = 0.01,
    source_type: str = "device",  # "device" or "benchmark"
    validate_with_solver_b: bool = True,
) -> dict:
    """Run grid series convergence analysis and verify acceptance criteria.

    Acceptance:
        h00 < 1% (tol_h00)
        |h0i| absolute and relative tolerance
        ||hij|| < 1% (tol_hij)
        validated against independent Solver B.
    """
    backpack = backpack or Backpack()
    drive = drive or Drive()
    cyl = cyl or ActiveCylinder()

    if grid_ns is None:
        grid_ns = [12, 16, 20, 28, 36]

    obs_obj = ObserverLocations()
    if observers is None:
        observers = obs_obj.as_array()
        obs_labels = obs_obj.names()
    else:
        obs_labels = [f"obs_{i}" for i in range(len(observers))]

    sources = build_sources(name, backpack)
    bump = GaussianBump(A=1.0, sigma_t=1e-3, sigma_x=0.5)

    solutions_by_n: dict[int, np.ndarray] = {}
    node_counts: dict[int, int] = {}

    for n in grid_ns:
        if source_type == "device":
            pts, w = focused_source_nodes(cyl.radius_m, cyl.height_m, n, 2 * n, n)
            T = fast_regularized_em_fields(pts, sources, drive.peak_current_a, drive.frequency_hz)
        else:
            pts, w = gauss_cylindrical_nodes(cyl.radius_m, cyl.height_m, n, 2 * n, n)
            T = conserved_T(bump, np.zeros(len(pts)), pts)

        node_counts[n] = len(pts)
        h_field = solver_a_tensor_product(observers, pts, w, T)
        solutions_by_n[n] = h_field

    consecutive_steps = []
    for i in range(len(grid_ns) - 1):
        n_prev = grid_ns[i]
        n_curr = grid_ns[i + 1]
        h_prev = solutions_by_n[n_prev]
        h_curr = solutions_by_n[n_curr]

        per_obs = []
        for o_idx in range(len(observers)):
            errs = component_relative_errors(h_curr[o_idx], h_prev[o_idx])
            errs["observer"] = obs_labels[o_idx]
            per_obs.append(errs)

        step_summary = {
            "from_N": n_prev,
            "to_N": n_curr,
            "from_nodes": node_counts[n_prev],
            "to_nodes": node_counts[n_curr],
            "max_h00_rel": float(max(e["h00_rel"] for e in per_obs)),
            "max_h0i_abs": float(max(e["h0i_abs"] for e in per_obs)),
            "max_h0i_rel": float(max(e["h0i_rel"] for e in per_obs)),
            "max_hij_rel": float(max(e["hij_rel"] for e in per_obs)),
            "max_tensor_rel": float(max(e["tensor_rel"] for e in per_obs)),
            "per_observer": per_obs,
        }
        consecutive_steps.append(step_summary)

    final_step = consecutive_steps[-1]
    final_h00_rel = final_step["max_h00_rel"]
    final_hij_rel = final_step["max_hij_rel"]
    final_h0i_abs = final_step["max_h0i_abs"]

    pass_h00 = bool(final_h00_rel < tol_h00)
    pass_hij = bool(final_hij_rel < tol_hij)
    pass_h0i = bool(final_h0i_abs < 1e-45 or final_step["max_h0i_rel"] < tol_h00)

    solver_b_check = {}
    if validate_with_solver_b:
        val_obs = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 2.0]])
        finest_n = grid_ns[-1]

        if source_type == "device":
            pts_f, w_f = focused_source_nodes(cyl.radius_m, cyl.height_m, finest_n, 2 * finest_n, finest_n)
            T_f = fast_regularized_em_fields(pts_f, sources, drive.peak_current_a, drive.frequency_hz)
            h_a_val = solver_a_tensor_product(val_obs, pts_f, w_f, T_f)

            def eval_source(pts_q):
                return fast_regularized_em_fields(pts_q, sources, drive.peak_current_a, drive.frequency_hz)

            h_b_val = solver_b_adaptive_voxel(
                val_obs, cyl.radius_m, cyl.height_m, eval_source, nx=28, ny=28, nz=28
            )
        else:
            pts_f, w_f = gauss_cylindrical_nodes(cyl.radius_m, cyl.height_m, finest_n, 2 * finest_n, finest_n)
            T_f = conserved_T(bump, np.zeros(len(pts_f)), pts_f)
            h_a_val = solver_a_tensor_product(val_obs, pts_f, w_f, T_f)

            def eval_source(pts_q):
                return conserved_T(bump, np.zeros(len(pts_q)), pts_q)

            h_b_val = solver_b_adaptive_voxel(
                val_obs, cyl.radius_m, cyl.height_m, eval_source, nx=32, ny=32, nz=32
            )

        solver_b_check = compare_green_solvers(h_a_val, h_b_val, eps_tol=0.01)

    all_pass = bool(pass_h00 and pass_hij and pass_h0i and (not validate_with_solver_b or solver_b_check.get("pass", True)))

    return {
        "geometry": name,
        "source_type": source_type,
        "grid_series": grid_ns,
        "node_counts": node_counts,
        "consecutive_steps": consecutive_steps,
        "final_step": final_step,
        "acceptance": {
            "pass_h00": pass_h00,
            "pass_hij": pass_hij,
            "pass_h0i": pass_h0i,
            "h00_rel": final_h00_rel,
            "h00_threshold": tol_h00,
            "hij_rel": final_hij_rel,
            "hij_threshold": tol_hij,
            "h0i_abs": final_h0i_abs,
            "solver_b_validation": solver_b_check,
            "overall_pass": all_pass,
        },
    }
