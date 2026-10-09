"""Quantitative gap between the current device and CTC-relevant scales.

Reads the device baseline, the EM source magnitude, and the exotic-matter
models, then produces a JSON report with all the "orders of magnitude away"
numbers. This is the honest deliverable: a bound, not a device.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Dict

import numpy as np

from .exotic_matter import (
    casimir_parallel_plates,
    morris_thorne_shell,
    tipler_cylinder,
    compare_to_em_source,
)

C = 2.99792458e8
G = 6.67430e-11
HBAR = 1.054571817e-34
PLANCK_ENERGY_J = 1.956e9
SUN_MASS_KG = 1.989e30
SUN_ENERGY_J = SUN_MASS_KG * C ** 2


def planck_scale_energy_J() -> float:
    return PLANCK_ENERGY_J


def sun_total_energy_J() -> float:
    return SUN_ENERGY_J


def energy_for_metric_perturbation(h_target: float,
                                   volume_m3: float = 1.0) -> float:
    """Rough energy scale needed for a metric perturbation of amplitude h
    over a volume V, from the linearized Einstein equation G ~ 8 pi G T / c^4.

    h ~ (G / c^4) * (E / V) * L^2  with L ~ V^{1/3}
    => E ~ h * c^4 * V / (G * V^{2/3}) = h * c^4 * V^{1/3} / G
    """
    L = volume_m3 ** (1.0 / 3.0)
    return h_target * (C ** 4) * L ** 3 / (G * L)  # = h c^4 L^2 / G


def build_report(em_rho_max: float = 1e-30,
                 h_current: float = 1e-46,
                 h_target: float = 1e-3) -> Dict:
    """Produce a JSON-serialisable report."""

    casimir = casimir_parallel_plates(separation_m=1e-6)
    mt = morris_thorne_shell(b0_m=1.0)
    tipler = tipler_cylinder(radius_m=1.0, omega_rad_s=1e6)

    report = {
        "stage": "energy-gap-report-v1",
        "constants": {
            "c_m_per_s": C,
            "G_m3_kg_s2": G,
            "hbar_J_s": HBAR,
        },
        "current_device": {
            "h_current": h_current,
            "em_rho_max_J_per_m3": em_rho_max,
        },
        "targets": {
            "h_target": h_target,
            "planck_energy_J": planck_scale_energy_J(),
            "sun_total_energy_J": sun_total_energy_J(),
        },
        "exotic_vs_em": {
            "casimir_d1um": compare_to_em_source(em_rho_max, casimir, 1e-6),
            "morris_thorne_b0_1m": compare_to_em_source(em_rho_max, mt, 1.0),
            "tipler_1m_1e6rad_s": compare_to_em_source(em_rho_max, tipler, 1.0),
        },
        "energy_to_perturbation": {},
    }

    for V in (1e-3, 1.0, 1e3):
        E = energy_for_metric_perturbation(h_target, V)
        report["energy_to_perturbation"][f"V={V}_m3"] = {
            "energy_J": E,
            "log10_E": math.log10(E) if E > 0 else None,
            "ratio_to_planck": E / planck_scale_energy_J(),
            "ratio_to_sun_total": E / sun_total_energy_J(),
        }

    # Orders of magnitude to close the h gap
    report["h_gap_orders"] = {
        "current_log10": math.log10(h_current),
        "target_log10": math.log10(h_target),
        "gap_orders": math.log10(h_target / h_current),
    }

    return report


def main(out_path: str = "results/validated/energy_gap_report.json") -> Dict:
    rep = build_report()
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rep, indent=2))
    print(json.dumps(rep, indent=2))
    return rep


if __name__ == "__main__":
    main()
