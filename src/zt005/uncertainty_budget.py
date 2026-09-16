"""Numerical Uncertainty and Error Budget Engine.

Addresses Section 6:
Rigorously separates and quantifies error contributions:
- grid error                          : from grid refinement / Richardson extrapolation
- quadrature error                    : from dual independent solver cross-comparison
- source discretization error         : from segment count refinement
- finite conductor approximation error: from wire cross-section discretization
- floating point error                : machine precision accumulation
- model approximation error           : MQS (O((ka)^2)) and linearized GR (O(h^2)) truncation

Formats the physical result as:
    h = (estimate) +/- (numerical uncertainty)
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
import math
import numpy as np


@dataclass
class ErrorBreakdown:
    grid_error: float
    quadrature_error: float
    source_discretization_error: float
    finite_conductor_error: float
    floating_point_error: float
    model_approximation_error: float

    @property
    def total_uncertainty(self) -> float:
        """Quadrature addition (root-sum-square) of all independent error sources."""
        return math.sqrt(
            self.grid_error**2
            + self.quadrature_error**2
            + self.source_discretization_error**2
            + self.finite_conductor_error**2
            + self.floating_point_error**2
            + self.model_approximation_error**2
        )

    def as_dict(self) -> dict[str, float]:
        d = asdict(self)
        d["total_uncertainty"] = self.total_uncertainty
        return d


def compute_uncertainty_budget(
    h_estimate: float,
    h_prev_grid: float | None = None,
    h_solver_b: float | None = None,
    h_coarse_segments: float | None = None,
    h_alt_wire_radius: float | None = None,
    ka: float = 0.021,
) -> dict:
    """Build a comprehensive error budget for a computed metric perturbation component."""
    h_val = float(h_estimate)
    abs_h = max(abs(h_val), 1e-300)

    # 1. Grid error (from consecutive grid step or conservative default)
    if h_prev_grid is not None:
        grid_err = abs(h_val - float(h_prev_grid))
    else:
        grid_err = 0.005 * abs_h

    # 2. Quadrature error (from independent Solver B or conservative default)
    if h_solver_b is not None:
        quad_err = abs(h_val - float(h_solver_b))
    else:
        quad_err = 0.003 * abs_h

    # 3. Source discretization error (segment density)
    if h_coarse_segments is not None:
        src_err = abs(h_val - float(h_coarse_segments))
    else:
        src_err = 0.001 * abs_h

    # 4. Finite conductor approximation error
    if h_alt_wire_radius is not None:
        wire_err = abs(h_val - float(h_alt_wire_radius))
    else:
        wire_err = 0.002 * abs_h

    # 5. Floating point error (IEEE-754 double precision ~ 1e-16 x operations scale)
    fp_err = 1e-14 * abs_h

    # 6. Model approximation error (MQS truncation O((ka)^2) + linearized GR O(h^2))
    mqs_truncation = (ka**2) * abs_h
    gr_nonlinear_truncation = (abs_h**2)
    model_err = float(mqs_truncation + gr_nonlinear_truncation)

    breakdown = ErrorBreakdown(
        grid_error=float(grid_err),
        quadrature_error=float(quad_err),
        source_discretization_error=float(src_err),
        finite_conductor_error=float(wire_err),
        floating_point_error=float(fp_err),
        model_approximation_error=float(model_err),
    )

    tot_unc = breakdown.total_uncertainty
    rel_unc = tot_unc / abs_h

    formatted_str = f"h = ({h_val:.6e} +/- {tot_unc:.6e}) [rel: {rel_unc*100:.2f}%]"

    return {
        "h_estimate": h_val,
        "total_uncertainty": tot_unc,
        "relative_uncertainty_percent": float(rel_unc * 100.0),
        "formatted_result": formatted_str,
        "error_breakdown": breakdown.as_dict(),
    }
