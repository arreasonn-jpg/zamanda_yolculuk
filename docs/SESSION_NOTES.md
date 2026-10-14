# Session notes — pickup point

Last stable release: **v10.2.0** (2026-10-14)
Tests: 697 (fast + slow), all passing.

## Milestones (v1.0 → v8.3)

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
| v7.0 | Refined complex+PTUM NEC (root cause fix) |
| v7.1 | Kerr + complex + PTUM |
| v7.2 | Transient discharge (negative result) |
| v7.3 | QFT threshold law (linear) |
| v7.4 | 4-regime structure (3 thresholds) |
| v7.5 | Discharge invariance + Kerr-independence |
| v8.0 | Brans-Dicke scalar-tensor extension |
| v8.1 | Loop Quantum Cosmology bounce |
| v8.2 | Starobinsky f(R) + LQC invariance |
| v8.3 | Hu-Sawicki f(R) + universal invariance |
| v9.0 | Newman-Penrose GW extraction formulas |
| v9.1 | Horndeski skeleton + universal invariance extension |
| v9.2 | BSSN -> NP pipeline (null result for symmetric source) |
| v9.3 | Quadrupole GW (scalar — noise-level) |
| v10.0 | Binary orbit scalar GW (null result — physics) |
| v10.1 | GW pipeline validated (TT injection) |
| v10.2 | Proca vector GW (h+ = 4.0e-6, first real signal) |

## Key scientific results

**QFT threshold law (v7.3):**
- ratio(tau=1) = 0.795 * |rho_min|
- Threshold: |rho_min| = 1.258

**4-regime structure (v7.4):**
| Regime | A range | NEC | QFT |
|---|---|---|---|
| Deactivated | < 1.2 | No | OK |
| Activated, no NEC | 1.2-1.35 | No | OK |
| NEC + QFT-OK | 1.35-1.87 | Yes | OK |
| NEC + QFT-violating | > 1.87 | Yes | No |

**Universal invariance (v8.3):**
PTUM QFT-OK window (A in [1.35, 1.87]) invariant under:
- Kerr spin (v7.5)
- Brans-Dicke phi (v8.0)
- LQC Friedmann (v8.2)
- Any scalar-tensor F(phi) (v8.3)

**Alternative gravity (v8.x):**
- Brans-Dicke: Mercury 43 arcsec, ω > 40000
- LQC: Big Bang bounce at ρ = ρ_c
- Starobinsky: n_s = 0.965, r = 0.004
- Hu-Sawicki: chameleon screening

## Bilingual paper

- `paper/main_en.tex` — 13 pages, English
- `paper/main_tr.tex` — 13 pages, Turkish
- 15 figures (fig1-fig15)
- `paper/arxiv_submission_v3.zip` (1.2 MB) — READY

## GW extraction (v9.x)

- Newman-Penrose module (psi_4, TT, flux)
- BSSN -> TT pipeline
- Null result: spherical source -> no GW (correct)
- Nonzero: quadrupole source -> h+ = 4.46e-12 (correct)

## Cihaz vs CTC

**43-order gap: UNCHANGED** (physics, not bug).

## Physics papers supported

- **SC** — numerically supported (v3-v6)
- **PTUM** — numerically supported (v5-v8)
  - QFT-OK window: A in [1.35, 1.87]
  - Background-independent (universal)

## What's next

**A) arXiv upload** — package ready, user decision
**B) New physics** — beyond scalar-tensor (Horndeski, DHOST, etc.)
**C) 3D + GW extraction** — Newman-Penrose scalars
**D) Stop**

## Resume

    cd ~/zamanda_yolculuk
    source .venv/bin/activate
    pytest -q -m "not slow"
