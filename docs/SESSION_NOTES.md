# Session notes — pickup point

Last stable release: **v6.6.0** (2026-10-12)
Tests: 555 fast + 21 slow = 576 total, all passing.
Tags: v1.0.0, v2.0.0, v3.0.0, v4.0.0, v5.0.0, v6.0.0-6.6.0 (+ alpha tags).

## What works

**Core BSSN solver:**
- Schwarzschild (isotropic + puncture, moving puncture)
- Kerr (Kerr-Schild, a <= M, slowly-rotating and near-extremal)
- 2-level AMR-lite refinement (prolongate/restrict)

**Matter modules:**
- Scalar field, complex scalar (U(1)), PTUM (2-field)
- SC condensate (QCD-inspired)
- Full stack: Kerr + SC metric + SC condensate + PTUM

**Diagnostics:**
- Hamiltonian / momentum constraint
- Ford-Roman quantum inequality
- Kmax, alpha_min monitoring
- Refined full stack stability

## Latest key results (v6.6)

| Configuration | Steps | Result |
|---|---|---|
| Kerr + SC + condensate | 500 | H growth 0.9823 |
| Kerr + PTUM a=0.9 | 1000 | stable |
| Full stack single grid a=0.9 | 200 | stable |
| Refined SC+PTUM A=3.0 | 100 | NEC violation, rho_min = -0.099 |
| Complex + PTUM A=3.0 (eq seed) | 100 | NEC violation, rho_min = -15.05 |

## Known limitations (documented)

1. **Refined complex + PTUM loses equilibrium psi** through
   prolongate/restrict cycle → rho_min → 0. Deferred to v7.0.
2. Long-time tachyon growth from small seed not simulated
   (would require t > 2.4 for exponential activation).
3. No full 3D AMR; refinement is 2-level static.
4. No MPI / radiation extraction.

## Physics papers — status

**SC** — numerically supported (v4.0-v6.5): regular core stable,
rotating + condensate 500-step.

**PTUM** — numerically supported (v5.0-v6.6): NEC violation,
transient discharge, QFT-compatible (massless).

**SC ↔ PTUM bridge** — mathematical parallel verified; physical
bridge NOT established (78-order scale gap).

## Bilingual paper

- `paper/main_en.tex` — 10 pages, English
- `paper/main_tr.tex` — 10 pages, Turkish
- `paper/bssn_section_en.tex` / `_tr.tex` — v6.5 results
- `paper/figures/fig7_fullstack.png` — 3-panel diagnostic
- `paper/arxiv_submission_v3.zip` (~1 MB) — ready for upload

## Cihaz vs CTC

43-order gap to CTC-relevant scale: **UNCHANGED** (physics, not bug).

## What's next

**A) Paper complex scalar section** — add to bssn_section EN + TR
**B) arXiv upload** — package ready
**C) v7.0 — consistent psi init across refinement levels**
**D) v7.0 — full 3D + complex scalar + refinement + PTUM**

## How to resume

    cd ~/zamanda_yolculuk
    source .venv/bin/activate
    pytest -q -m "not slow"        # ~2 min
    git log --oneline -20
