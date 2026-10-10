# Physics gap: current device vs Terminator-5 lab conditions

This document consolidates the quantitative gap between the current
driven-EM device and the "laboratory conditions" implied by the
Terminator-5 TDE (spacetime displacement equipment).

## Consolidated numbers

| Quantity | Current | Target | Gap (orders) |
|---|---|---|---|
| Metric perturbation h | 1e-46 | 1e-3 | 43 |
| Energy density (1 m^3) | ~1e-30 J | 1.21e41 J | 71 |
| Nonlinearity eps | 1e-46 | ~1 | 46 |
| Negative energy density | Casimir 1e-3 J/m^3 | ~1e20 J/m^3 | 23 |

## Physics required (beyond current laws)

1. Macroscopic negative-energy-density control (not just Casimir-scale).
2. Non-linear GR regime with h >~ 1e-3 (back-reaction dominates).
3. Stable CTC topology (NEC violation at macroscopic scale).
4. Energy budget at stellar-mass scale in a lab-sized volume.

## What this project delivers

- A **numerical bound**, not a device claim.
- An honest FAIL register for every convergence check.
- A modular pipeline: EM source -> T_munu -> linearized GR -> BSSN.
- A public CI suite (341 fast + 15 slow, all passing).

## What it does NOT deliver

- No time machine.
- No closed timelike curve.
- No claim that Terminator-5 conditions are achievable with current physics.

## Reproduce

    python -m zt005.energy_gap_report
    pytest -m "not slow"
    pytest
