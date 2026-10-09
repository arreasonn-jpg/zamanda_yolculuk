# ZT-010 — Exotic-matter integration

**Status:** open, exploratory.

## Deliverables
- src/zt005/exotic_matter.py
- src/zt005/energy_gap_report.py
- tests/test_exotic_matter.py (7 checks)
- docs/engineering/energy_scaling.md
- results/validated/energy_gap_report.json

## Findings
| Model | log10(rho_ex/rho_EM) |
|---|---|
| Casimir (d=1um) | +26.6 |
| Morris-Thorne (b0=1m) | +38.8 |
| Tipler (w=1e6 rad/s) | +50.8 |

h-gap: 43 orders of magnitude.

Not a device. Numerical bound only.

## Nonlinear GR skeleton added

- `src/zt005/nonlinear_gr.py` — scalar source, second-order back-reaction,
  linear→nonlinear crossover.
- `tests/test_nonlinear_gr.py` — 5 checks.
- Crossover: h = 1 (order-unity), current h = 1e-46 → **46 orders**.
- Linear vs nonlinear ratio at h=1e-3: ~1e-6 (safe linear regime).
- Status: skeleton only; full BSSN solver not implemented.

## BSSN skeleton added

- `src/zt005/bssn.py` — ADM/BSSN containers, constraint residuals,
  Euler step, linear-regime check.
- `tests/test_bssn.py` — 6 checks.
- Status: skeleton only. No full BSSN evolution, no constraint damping.
- Purpose: lay out the 3+1 variable structure so future work can plug
  real RHS terms in `euler_step`.
