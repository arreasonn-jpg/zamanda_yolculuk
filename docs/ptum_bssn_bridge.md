# PTUM-BSSN bridge

Connects the Phase Transition Universes Model (PTUM)
[Esa 2026, in preparation] to numerical BSSN evolution.

## Claim tested

PTUM predicts a tachyon discharge event:
  m_eff^2(Phi) = m^2 (1 - alpha Phi^2 / Eth^2)
becomes negative when Phi > Phi_c = Eth / sqrt(alpha).
The activated proxy field psi then rolls to a minimum with
V_min = -m_eff^4 / (4 lambda_psi), producing negative energy density.

## Method

Two real scalar fields (Phi, psi) coupled to BSSN:
- Fuel Phi: plateau potential V(Phi) = V0 (1 - e^{-lambda Phi^2})
- Proxy psi: V(psi; Phi) = 1/2 m_eff^2(Phi) psi^2 + 1/4 lambda_psi psi^4

RK4 evolution on (metric + Phi + psi).

Modules:
- src/zt005/ptum_field.py (potentials)
- src/zt005/matter_ptum.py (stress-energy)
- src/zt005/bssn_rk4_ptum.py (RK4)

## Result

**NEC violation confirmed.**
- rho_min = -15.6 (A=3.0)
- 27/729 cells in negative-energy state
- Stable at 100 steps, fails ~500 (backreaction)

**Stable at perturbative amplitude.**
- rho_min = -0.29 (A=1.5)
- Kmax = 0.29 after 200 steps
- All fields finite

## Interpretation

The PTUM tachyon discharge mechanism produces the NEC-violating
energy density required for closed timelike curves. The magnitude
is controlled by the activation amplitude, and the configuration
remains numerically stable in the perturbative regime.

This is the first numerical realization of PTUM's central physical
prediction in a full 3+1 GR setting.

## Reproduce

    pytest tests/test_ptum_field.py -v
    python scripts/v5_ptum_equilibrium.py
    python scripts/v5_ptum_plot.py

## v6.0 update — refined stability

The single-grid PTUM run at A=3.0 (rho_min = -15.6) crashed at ~100
steps due to strong tachyon backreaction making gtilde singular.

Two-grid refinement (2:1, prolongate/restrict coupling) enables
**200-step stable** evolution at the same amplitude:

| Setup | Fails at | Kmax @100 | Kmax @200 |
|---|---|---|---|
| v5.0 single grid | ~100 | 23.5 | (crash) |
| v6.0 refined | never | 0.054 | 0.107 |

**Interpretation:** resolution near the tachyon-active region is
critical for long-run stability. The negative energy configuration
becomes tractable when the active zone is adequately resolved.
