# Session notes — pickup point

Last stable release: **v7.0.0** (2026-10-12)
Tests: 579 (fast) + 21 (slow), all passing.

## What works (matter modules)

- Scalar, complex scalar (U(1))
- PTUM (2-field tachyon discharge)
- SC condensate (QCD-inspired)
- Full stack: Kerr + SC + condensate + PTUM
- Refined (2-level) versions of all above

## Latest key results

| Configuration | Steps | rho_min | Notes |
|---|---|---|---|
| Refined complex+PTUM A=3.0 | 100 | -15.67 | NEC persistent |
| Refined SC+PTUM A=3.0 | 100 | -0.099 | NEC stable |
| Full stack single grid | 200 | -0.284 | stable |
| Kerr+SC 500-step | 500 | — | H growth 0.9823 |

## Critical fixes in v7.0

1. `psi_equilibrium_smooth` — smoothstep regularization at Phi_c
2. `restrict_matter=False` — one-way coupling (avoid 35% peak loss)

## Next options

**A) Paper update** — add v7.0 refined complex+PTUM to bilingual paper
**B) arXiv upload** — package ready
**C) Two-way matter coupling without peak loss**
**D) Kerr + complex + PTUM**

## Resume

    cd ~/zamanda_yolculuk
    source .venv/bin/activate
    pytest -q -m "not slow"
