# Cross-scale phase transition bridge (speculative)

This note documents a *speculative* structural comparison between
the Schwarzschild Condensation (SC) and Phase Transition Universes
Model (PTUM) phase transitions.

**This is not a unification claim.** The two models operate at a
78-order-of-magnitude scale gap and no derivation connects them.

## Common mathematical structure

Both models use a scalar order parameter `psi` with an effective
potential:

    V(psi; C) = 1/2 m_eff^2(C) psi^2 + 1/4 lambda psi^4

    m_eff^2(C) = m^2 (1 - C / C_c)

where:
- SC:   C = QCD energy density / pairing strength
        C_c = deconfinement threshold
- PTUM: C = fuel field Phi^2
        C_c = tachyon activation threshold

## Common dynamical pattern

| Stage | SC | PTUM |
|---|---|---|
| Control parameter C < C_c | hadronic phase | dormant phase |
| C crosses C_c | deconfinement | tachyon activation |
| psi forms | CFL condensate | proxy condensate |
| Metric backreaction | regular core | NEC violation |

## Numerical test

`phase_transition_common.py` implements the shared formula:

    V(psi; C) = 1/2 m^2 (1 - C/C_c) psi^2 + 1/4 lambda psi^4
    psi_0(C) = sqrt(-m_eff^2 / lambda)   for C > C_c
    growth rate = sqrt(|m_eff^2|)        for C > C_c
    timescale   = 1 / growth rate

These are identical for SC and PTUM because the *formula* is shared.
The physical difference is:

- SC: V_min is shifted to be >= 0 (regular core)
- PTUM: V_min is negative (NEC violation)

## What this shows

The two models share:
1. Mathematical form of the phase transition potential
2. Threshold-crossing mechanism
3. Condensate formation dynamics
4. BSSN infrastructure

## What this does NOT show

1. That SC and PTUM are the same physics
2. That any mechanism bridges the 78-order scale gap
3. That PTUM's negative energy is compatible with SC's regular core
4. That either model is correct

## Status

- **Mathematical parallel:** verified (this repo)
- **Physical bridge:** NOT established
- **Unification:** NOT claimed

## Reproduce

    pytest tests/test_phase_transition_common.py -v

## Reference

- Esa (2026a) SC, in preparation
- Esa (2026b) PTUM, in preparation
