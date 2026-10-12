# Session notes — pickup point

Last stable release: **v6.5.0** (2026-10-12)
Tests: 549 fast + 21 slow = 570 total, all passing.
Tags: v1.0.0, v2.0.0, v3.0.0, v4.0.0, v5.0.0, v6.0.0-6.5.0 + 30+ alpha/beta.

## What works

- BSSN solver (Schwarzschild, Kerr a<=0.5, moving puncture)
- SC regular core (metric only, and metric+matter)
- PTUM tachyon discharge (single + refined + Kerr)
- Complex scalar matter
- Scalar field matter coupling
- 2-level AMR-lite refinement (prolongate/restrict)
- Ford-Roman quantum inequality diagnostic
- Full stack: Kerr + SC metric + SC condensate + PTUM (200-step stable)
- Refined full stack A=3.0 (100-step, NEC violation)

## Latest key results (v6.5)

| Configuration | Steps | Result |
|---|---|---|
| Kerr + SC + condensate | 500 | H growth 0.9823 |
| Kerr + PTUM a=0.9 | 1000 | stable |
| Full stack single grid a=0.9 | 200 | stable, Kmax bounded |
| Refined full stack A=3.0 | 100 | rho_min = -0.099 (NEC violation) |

## Physics papers — status

- **SC** — numerically supported (v4.0-v6.5): regular core stable,
  rotating + condensate 500-step
- **PTUM** — numerically supported (v5.0-v6.5): NEC violation,
  transient discharge, QFT-compatible (massless), full stack refined
- **SC ↔ PTUM bridge** — mathematical parallel verified; physical
  bridge NOT established (78-order scale gap)

## Bilingual paper — status

- `paper/main_en.tex` — 10 pages, English, ready
- `paper/main_tr.tex` — 10 pages, Turkish, ready
- `paper/bssn_section_en.tex` — BSSN section (v6.5 results)
- `paper/bssn_section_tr.tex` — Turkish translation
- Both compile with 0 overfull hbox
- `paper/arxiv_submission_v3.zip` (982 KB) — ready for upload

## What's next (when returning)

Options on the table:

**A) Paper figure** — add full-stack diagnostics figure to paper
   (currently only text)

**B) Complex scalar + full stack** — extend matter to U(1)
   complex scalar combined with SC and PTUM

**C) arXiv upload** — package ready, user decision pending

**D) v7.0 — full 3D + Kerr + complex scalar + refinement + PTUM**
   (weeks of work)

## How to resume

    cd ~/zamanda_yolculuk
    source .venv/bin/activate
    pytest -q -m "not slow"        # ~2 min
    git log --oneline -20
    git tag --sort=-creatordate | head -20

## Current status

- 6 major releases (v1.0 → v6.5)
- 570 tests passing
- 43-order gap to CTC regime: UNCHANGED (physics, not a bug)
- Numerical laboratory: full stack stable + NEC violation proven
- arXiv package ready for upload
