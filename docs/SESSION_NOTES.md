# Session notes — pickup point

Last stable release: **v13.0.0** (2026-10-20)
Tests: fast suite all passing (see CHANGELOG per version).

## Milestones (v1.0 -> v13.0)

| Version | Highlight |
|---|---|
| v1.0 | Linearized GR + 43-order gap |
| v3.0 | Moving puncture + SC regular core |
| v5.0 | PTUM tachyon + NEC violation |
| v7.0 | Refined complex+PTUM NEC |
| v7.3 | QFT threshold law (linear) |
| v7.4 | 4-regime structure |
| v8.0 | Brans-Dicke |
| v8.1 | LQC bounce |
| v8.3 | Universal invariance |
| v9.3 | NP + GW extraction |
| v10.2 | Proca vector GW |
| v11.0 | 3D binary BH |
| v12.0 | Bowen-York orbit |
| v12.2 | Vectorized BSSN (25x) |
| v12.3 | SC+PTUM+BH triple toy |
| v12.5 | Full stress-tensor + retarded solver chain |
| v12.9 | N=17 diagnosis + 2x AMR proof of concept |
| v12.10-* v12.15 | AMR limit chain: cubic/finer-fine/Gamma rederiv; injection physics |
| v12.14 | ZT-006.2.24/.25 source chain closed |
| v13.0 | Moving patch (COM) no-op; restructure amr/ + matter/ |

## Key results

**QFT threshold (v7.3):**
- ratio(tau=1) = 0.795 * |rho_min|
- Threshold: |rho_min| = 1.258

**Universal invariance (v8.3, ext. v12.3):**
PTUM QFT window A in [1.35, 1.87] invariant under:
- Kerr spin, Brans-Dicke phi, LQC, F(phi), Horndeski, SC+BH

**GW extraction (v9.0-v10.2):**
- Scalar: h+ = 0 (no TT content)
- TT injection: validated
- Proca: h+ = 4.0e-6 (first real signal)

**3D binary BH (v11.0-v12.2):**
- Moving puncture 3D: growth 1.0001
- Bowen-York: h+ = 1.0e-4, hx = 4.4e-3 (N=13)
- Vectorized RHS: 25x speedup

## 43-order gap: UNCHANGED

## Paper
- main_en.tex / main_tr.tex (arxiv_v3 active): validation + bssn
  sections updated through v13.2 (moving patch no-op, 4x negative,
  constraint projection positive); arxiv_submission_v4.zip
  regenerated (`scripts/dev/make_zip_v4.py`)
- 15 figures
- paper/archive/ holds historical submission zips

## Current state (v13.2.0, 2026-10-23)
- v13.0.0: moving patch = no-op (COM drift < 1 cell) + subpackage
  restructure (zt005.amr/, zt005.matter/); v13.0.1: docs sync.
- v13.1.0: 4x + subcycle=4 closed NEGATIVE (H@500 1.1418 vs 1.0166).
- v13.2.0: band-local conformal projection FIRST POSITIVE lever
  (H_f/H0 @ 1000 = 0.9319 vs baseline 1.0756, -13.4%); solver core
  `src/zt005/amr/amr_projection.py` + 5 unit tests.
- v13.3.0: exact gamma-tilde covariant operator NEUTRAL (0.9382 vs
  0.9319) — gain is operator-robust; iso-reference metric delivered
  (~10% of the sub-1.0 effect is definitional). 10 projection tests.
- AMR lever chain fully mapped; the remaining defect is genuinely
  the injected band's second-derivative mismatch.

## Next
- A) v13.4: momentum-constraint projection analog (Vector-Poisson
   on Atilde divergence, plan doc increment #2)
- B) Full inspiral (weeks)
- C) arXiv upload v4 (user decision)
- D) New physics (Horndeski/DHOST)

## Resume
    cd ~/zamanda_yolculuk
    source .venv/bin/activate
    pytest -q -m "not slow"
