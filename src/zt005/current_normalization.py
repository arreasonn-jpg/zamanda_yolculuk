"""Current Normalization and Classification Module.

Addresses Section 1.1 of the scientific maturity specification:
- Disentangles and strictly defines:
    * I_parent : current per parent conductor turn [A]
    * I_total  : total current across all turns/loops in the assembly [A]
    * I_peak   : peak amplitude of the AC waveform [A]
    * I_rms    : root-mean-square amplitude I_peak / sqrt(2) [A]
- Guarantees current conservation across cross-sectional quadrature refinement:
    Sum_k w_k == 1.0  =>  Sum_k I_k == I_parent
    Sum_all I_k == I_total
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from .model import Backpack
from .area_weighted_finite_wire import build_area_weighted_groups
from .geometry import build_sources


@dataclass(frozen=True)
class CurrentSpec:
    """Strictly distinguishes parent, total, peak, and rms current representations."""
    i_peak: float = 100.0  # Peak current per parent conductor turn [A]
    n_turns: int = 1       # Number of physical loops / turns

    def __post_init__(self):
        if self.i_peak < 0:
            raise ValueError(f"i_peak must be non-negative, got {self.i_peak}")
        if self.n_turns < 1:
            raise ValueError(f"n_turns must be at least 1, got {self.n_turns}")

    @property
    def i_parent(self) -> float:
        """Alias for i_parent_peak: peak current in one physical conductor."""
        return self.i_peak

    @property
    def i_rms(self) -> float:
        """RMS current in one parent conductor turn [A]."""
        return self.i_peak / math.sqrt(2.0)

    @property
    def i_parent_peak(self) -> float:
        """Peak current in one parent conductor turn [A]."""
        return self.i_peak

    @property
    def i_parent_rms(self) -> float:
        """RMS current in one parent conductor turn [A] for sinusoidal AC."""
        return self.i_peak / math.sqrt(2.0)

    @property
    def i_total_peak(self) -> float:
        """Total peak current across all turns in the geometry [A]."""
        return self.n_turns * self.i_peak

    @property
    def i_total_rms(self) -> float:
        """Total RMS current across all turns in the geometry [A]."""
        return self.n_turns * self.i_parent_rms

    def as_dict(self) -> dict[str, float]:
        return {
            "I_parent": self.i_parent,
            "I_parent_peak": self.i_parent_peak,
            "I_parent_rms": self.i_parent_rms,
            "I_total": self.i_total_peak,
            "I_total_peak": self.i_total_peak,
            "I_total_rms": self.i_total_rms,
            "I_peak": self.i_peak,
            "I_rms": self.i_parent_rms,
            "n_turns": self.n_turns,
        }


def get_geometry_turns(name: str) -> int:
    """Return the number of parent loops for a given geometry."""
    counts = {
        "G1": 6,  # 6 circular loops along y
        "G2": 4,  # 4 circular loops along x
        "G3": 2,  # 2 counter-wound helical coils
        "G4": 3,  # 3 circular loops (2 along y, 1 along x)
    }
    if name not in counts:
        raise ValueError(f"Unknown geometry {name}")
    return counts[name]


def make_current_spec(name: str, i_peak: float = 100.0) -> CurrentSpec:
    """Create a CurrentSpec correctly calibrated to the named geometry."""
    n_turns = get_geometry_turns(name)
    return CurrentSpec(i_peak=i_peak, n_turns=n_turns)


def verify_current_conservation(
    groups: list[list[tuple[tuple[np.ndarray, np.ndarray], float]]],
    current_spec: CurrentSpec,
    tol: float = 1e-12,
) -> dict:
    """Verify that conductor subdivision preserves parent and total current.

    Parameters
    ----------
    groups : list of parent turns, each containing list of ((pts, seg), weight)
    current_spec : CurrentSpec defining I_parent and I_total
    tol : numerical tolerance for weight sum comparison
    """
    n_groups = len(groups)
    group_weight_sums = []
    group_current_sums = []

    for grp in groups:
        w_sum = sum(w for _, w in grp)
        group_weight_sums.append(float(w_sum))
        i_sum = sum(w * current_spec.i_parent_peak for _, w in grp)
        group_current_sums.append(float(i_sum))

    total_weight_sum = float(sum(group_weight_sums))
    total_current_calc = float(sum(group_current_sums))
    expected_total_current = current_spec.i_total_peak

    weights_normalized = all(abs(ws - 1.0) < tol for ws in group_weight_sums)
    total_current_matches = abs(total_current_calc - expected_total_current) < tol * max(expected_total_current, 1.0)

    max_weight_residual = max(abs(ws - 1.0) for ws in group_weight_sums)
    total_current_residual = abs(total_current_calc - expected_total_current)

    is_passed = weights_normalized and total_current_matches and (n_groups == current_spec.n_turns)

    return {
        "pass": bool(is_passed),
        "n_groups": n_groups,
        "n_subfilaments_total": sum(len(grp) for grp in groups),
        "expected_turns": current_spec.n_turns,
        "I_parent_peak": current_spec.i_parent_peak,
        "I_parent_rms": current_spec.i_parent_rms,
        "I_total_peak": current_spec.i_total_peak,
        "I_total_rms": current_spec.i_total_rms,
        "calculated_total_current": total_current_calc,
        "max_group_weight_residual": float(max_weight_residual),
        "total_current_residual": float(total_current_residual),
        "group_weight_sums": group_weight_sums,
    }


def multi_resolution_current_audit(
    name: str,
    backpack: Backpack,
    radius_m: float = 0.0025,
    rules: list[tuple[int, int]] | None = None,
    i_peak: float = 100.0,
) -> dict:
    """Check that different cross-section quadrature resolutions yield identical total currents.

    Acceptance criterion: Same physical configuration produces identical total current
    across all quadrature resolutions.
    """
    if rules is None:
        rules = [(1, 6), (2, 8), (3, 10), (4, 12), (5, 16)]

    spec = make_current_spec(name, i_peak=i_peak)
    evaluations = []

    for nr, nphi in rules:
        groups = build_area_weighted_groups(
            name, backpack, radius_m=radius_m, nr=nr, nphi=nphi
        )
        report = verify_current_conservation(groups, spec)
        evaluations.append({
            "nr": nr,
            "nphi": nphi,
            "subfilaments_per_turn": nr * nphi,
            "total_subfilaments": report["n_subfilaments_total"],
            "calculated_total_current": report["calculated_total_current"],
            "current_residual": report["total_current_residual"],
            "pass": report["pass"],
        })

    all_pass = all(e["pass"] for e in evaluations)
    # Check pairwise consistency across resolutions
    currents = [e["calculated_total_current"] for e in evaluations]
    max_resolution_discrepancy = max(abs(c - currents[0]) for c in currents)

    return {
        "geometry": name,
        "current_spec": spec.as_dict(),
        "n_resolutions_tested": len(rules),
        "max_resolution_discrepancy": float(max_resolution_discrepancy),
        "all_pass": bool(all_pass and max_resolution_discrepancy < 1e-12),
        "evaluations": evaluations,
    }
