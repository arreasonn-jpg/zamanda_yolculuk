# Session notes — pickup point

Last stable release: **v12.3.0** (2026-10-15)
Tests: 698 fast + 21 slow, all passing.

## Milestones (v1.0 -> v12.3)

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
- main_en.tex / main_tr.tex: 15 pages each
- 15 figures
- arxiv_submission_v3.zip (~1.5 MB) READY

## Next
- A) Full inspiral (weeks)
- B) Cython for N>=15
- C) arXiv upload (user decision)
- D) New physics (Horndeski/DHOST)
- E) Stop

## Resume
    cd ~/zamanda_yolculuk
    source .venv/bin/activate
    pytest -q -m "not slow"
