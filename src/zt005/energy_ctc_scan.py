"""Energy-to-CTC Threshold Scan Engine.

Addresses Section 18:
Establishes the quantitative physical chain:
    I -> (E, B) -> T_munu -> g_munu -> CTC Search

Performs the multi-decade current and energy threshold search:
    100 A   -> h ~ 1.98e-46 -> no CTC
    1 kA    -> h ~ 1.98e-44 -> no CTC
    10 kA   -> h ~ 1.98e-42 -> no CTC
    100 kA  -> h ~ 1.98e-40 -> no CTC
    ...
    I_CTC   -> h ~ O(1)     -> CTC threshold

Calculates:
- I_CTC [A] : Current required for metric inversion / closed timelike loops
- E_CTC [J] : Total stored electromagnetic energy required
- Orders-of-magnitude distance between laboratory apparatus and physical CTC formation.
"""

from __future__ import annotations

import math
import numpy as np

C = 299_792_458.0
G_CONST = 6.67430e-11


def run_energy_to_ctc_scan(
    baseline_current_a: float = 100.0,
    baseline_h00: float = 1.98e-46,
    baseline_energy_j: float = 0.0136,
    current_scan_levels: list[float] | None = None,
) -> dict:
    """Execute the multi-decade energy-to-CTC threshold scan.

    Because Maxwell stress-energy scales quadratically with current:
        T_munu ~ I^2
        h_munu ~ (4G / c^4) * integral(T/R) ~ I^2
    The metric perturbation reaches the causal-inversion regime (light-cone tipping, g_phiphi < 0)
    when h_scale ~ 1.0.
    """
    if current_scan_levels is None:
        current_scan_levels = [
            100.0,          # 100 A (Laboratory baseline)
            1.0e3,          # 1 kA
            1.0e4,          # 10 kA (Superconducting magnet scale)
            1.0e5,          # 100 kA
            1.0e6,          # 1 MA (Z-machine pulse scale)
            1.0e7,          # 10 MA
            1.0e9,          # 1 GA
            1.0e15,         # 1 PA
            1.0e20,         # 100 EA
            1.0e25,         # Astronomical scale
        ]

    I0 = float(baseline_current_a)
    h0 = float(baseline_h00)
    E0 = float(baseline_energy_j)

    # Threshold condition: h ~ 1.0  =>  (I_ctc / I0)^2 * h0 = 1.0
    # I_ctc = I0 / sqrt(h0)
    i_ctc_threshold = I0 / math.sqrt(h0)
    e_ctc_threshold = E0 * (i_ctc_threshold / I0)**2

    log10_current_gap = math.log10(i_ctc_threshold / I0)
    log10_energy_gap = math.log10(e_ctc_threshold / E0)

    scan_records = []
    for current_a in current_scan_levels:
        ratio = current_a / I0
        energy_j = E0 * (ratio**2)
        h_scale = h0 * (ratio**2)

        # In weak-field linearized GR, g_00 = -1 + h_00.
        # Spacetime inversion / CTC formation requires h_scale >= 1.0
        has_ctc = bool(h_scale >= 1.0)

        # Light-cone tilt estimate: theta_cone ~ arctan(sqrt(h))
        cone_tilt_deg = math.degrees(math.atan(min(math.sqrt(h_scale), 1e6)))

        scan_records.append({
            "current_a": float(current_a),
            "stored_energy_j": float(energy_j),
            "h_scale": float(h_scale),
            "light_cone_tilt_deg": float(cone_tilt_deg),
            "ctc_formed": has_ctc,
            "verdict": "CTC FORMED" if has_ctc else "NO CTC (weak field)",
        })

    return {
        "baseline": {
            "current_a": I0,
            "stored_energy_j": E0,
            "h00": h0,
        },
        "threshold": {
            "I_CTC_amperes": float(i_ctc_threshold),
            "E_CTC_joules": float(e_ctc_threshold),
            "orders_of_magnitude_current_gap": float(log10_current_gap),
            "orders_of_magnitude_energy_gap": float(log10_energy_gap),
        },
        "astrophysical_comparison": {
            "solar_luminosity_per_second_j": 3.828e26,
            "supernova_total_energy_j": 1.0e44,
            "interpretation": (
                f"I_CTC = {i_ctc_threshold:.2e} A requires E_CTC = {e_ctc_threshold:.2e} J, "
                f"which is comparable to the gravitational binding energy of an entire supernova (~1e44 J). "
                f"Laboratory currents (100 A) are {log10_energy_gap:.1f} orders of magnitude below spacetime inversion."
            ),
        },
        "scan_records": scan_records,
    }
