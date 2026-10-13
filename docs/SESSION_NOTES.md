# Session notes — pickup point

Last stable release: **v7.3.0** (2026-10-13)
Tests: 595 (fast + slow), all passing.

## Milestones reached (v1.0 → v7.3)

| Version | Highlight |
|---|---|
| v1.0 | Linearized GR + exotic matter + 43-order gap |
| v2.0 | Full BSSN (Schwarzschild + Kerr) |
| v3.0 | Moving puncture + SC regular core |
| v4.0 | Matter-coupled BSSN |
| v5.0 | PTUM tachyon + NEC violation |
| v6.0 | Complex scalar + AMR-lite |
| v6.1 | SC + PTUM framework |
| v6.2 | Kerr + PTUM refined |
| v6.3 | Kerr + SC (500-step) |
| v6.4 | Full stack SC + PTUM |
| v6.5 | Refined full stack (NEC) |
| v6.6 | Complex + PTUM |
| v7.0 | Refined complex+PTUM NEC (root cause fixed) |
| v7.1 | Kerr + complex + PTUM |
| v7.2 | Transient discharge (negative result) |
| v7.3 | QFT threshold law (linear) |

## Latest key results

| Config | rho_min | Ford-Roman |
|---|---|---|
| Weak PTUM (A=1.5) | -0.28 | Compatible |
| Strong PTUM (A=3.0) | -15.66 | Violated (tau>=1) |
| Kerr+complex+PTUM A=3.0 | -15.66 | Violated |
| QFT threshold | \|rho\| = 1.258 | ratio=1 |

## QFT threshold law

**ratio(tau=1) = 0.795 * |rho_min|**

QFT-compatible: |rho_min| < 1.258 (A_fuel < 1.87)
QFT-incompatible: |rho_min| > 1.258

## Bilingual paper

- paper/main_en.tex — 11 pages, English
- paper/main_tr.tex — 11 pages, Turkish
- 8 figures (fig1-fig8)
- paper/arxiv_submission_v3.zip — READY (upload = user decision)

## Cihaz vs CTC

43-order gap: UNCHANGED.

## What's next

**A) QFT-compatible sub-regime** — A < 1.87 detailed study
**B) arXiv upload** — ready
**C) Full 3D + all mechanisms**
**D) Stop**

## Resume

    cd ~/zamanda_yolculuk
    source .venv/bin/activate
    pytest -q -m "not slow"
