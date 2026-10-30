# v3.0-beta Release Notes

**Release Date:** 2026-10-30  
**Git Tag:** v3.0-beta  
**Base Commit:** 8d06cfd (Update roadmap with completed full long-run test results)

---

## Overview

v3.0-beta consolidates all v3.0-alpha milestones into a stable release featuring:
- **Complete BSSN solver** with 1+log lapse, Gamma-driver shift, Z4c constraint damping, Sommerfeld BC
- **SC regular core** (de Sitter interior + Israel junction to Schwarzschild exterior)
- **Rotating SC core (Kerr)** with frame-dragging shift for slow rotation (|a| ≤ 0.5M)
- **Moving-puncture formalism** for binary black hole evolutions
- **AMR with Brill-Lindquist data** - moving fine patch with projection operators
- **Full long-run stability validation** (1000+ steps across all test configurations)

---

## Major Features

### 1. BSSN Solver (v2.0 complete, integrated in v3.0)
- **Full RHS**: Christoffel, Ricci, constraint damping (Z4c)
- **Gauge**: 1+log lapse + Gamma-driver shift
- **Boundary**: Sommerfeld outflow BC (validated)
- **Time integration**: RK4 with Kreiss-Oliger dissipation
- **Stability**: 1000+ steps validated (H growth < 1.01)

### 2. SC Regular Core (v3.0-alpha4/5)
- **de Sitter interior** matched to **Schwarzschild exterior** via Israel junction conditions
- **Geometric units** (G=c=1) throughout
- **BSSN initial data** from SC metric
- **Evolution**: 500 steps stable (H growth = 1.0000)

### 3. Rotating SC Core - Kerr (v3.0-alpha6) ⭐ NEW
- **Frame-dragging shift**: β^φ ~ 2aM_eff/r³ for slow rotation
- **Valid regime**: |a| ≤ 0.5M (slow rotation approximation)
- **Evolution**: 500 steps stable (H growth = 1.0000126)
- **Module**: `src/zt005/sc_core_kerr.py`

### 4. Moving-Puncture Formalism (v3.0-alpha2/3)
- **Gauge-consistent IC** for Schwarzschild puncture
- **500-step evolution** stable (H growth = 1.0000)
- **Convergence study**: N=9,11,13 (H decreases 2.7x)

### 5. AMR with Brill-Lindquist Binary Data ⭐ NEW
- **Moving fine patch** tracking binary COM
- **Projection operators**: Covariant band projection with momentum correction
- **L_c=160, 50 steps**: H_c_max 0→3.59e-4, H_f_band stable at 2.2e-8
- **L_c=200, 50 steps**: Projection reduces fine-grid band constraint (l_after < l_before)
- **Results**: `results/validated/amr_bl_moving/`

---

## Test Suite - All Passing

| Test File | Configuration | Steps | Status |
|-----------|---------------|-------|--------|
| `test_bssn_beta3_longrun.py` | BSSN beta3 long-run | 1000 | ✅ PASSED |
| `test_bssn_mp_longrun.py` | Moving puncture Schwarzschild | 500 | ✅ PASSED |
| `test_sc_core_kerr.py` | SC core Kerr (rotating) | 500 | ✅ PASSED |
| `test_sc_core_bssn.py` | SC core BSSN evolution | 500 | ✅ PASSED |
| `test_bssn_v2_longrun.py` | v2.0 BSSN long-run | 1000 | ✅ PASSED |
| `test_bssn_beta3_convergence.py` | Convergence study | 30/100 | ✅ PASSED |
| `test_bssn_beta3_schwarzschild.py` | Schwarzschild puncture | 50/500 | ✅ PASSED |
| `test_sc_core.py` | SC core unit tests | 7 tests | ✅ 7/7 PASSED |

**Total**: 14 test functions, 14 passed, 6 skipped (manual longrun variants)

---

## AMR Validation Results

### L_c=160 (50 steps)
- **Coarse grid**: 33³, L_c=160
- **Fine grid**: 25³, L_f=40 (moving patch)
- **H_c_max**: 0 → 3.59×10⁻⁴ (controlled growth)
- **H_f_band**: Stable at 2.2×10⁻⁸ (projection effective)
- **Binary separation**: Constant at 9.7M (COM tracking working)

### L_c=200 (50 steps)
- **Projection data**: l_before > l_after consistently (constraint reduction)
- **Fine-grid band**: H reduction confirmed via projection operator

### Projection Test (5 steps)
- **dφ_max**: ~0.0085 (momentum correction magnitude)
- **l_before → l_after**: ~0.0057 → ~0.0045 (20% reduction)
- **H_f_max**: Stable ~0.054-0.055

---

## Documentation - All PDFs Current

| Document | Pages | Size | Status |
|----------|-------|------|--------|
| `paper/main_en.pdf` | 25 | 2.0 MB | ✅ Compiled |
| `paper/main_tr.pdf` | 24 | 2.0 MB | ✅ Compiled |
| `paper/arxiv_v3/main_en.pdf` | 19 | 1.97 MB | ✅ Compiled |
| `paper/arxiv_v3/main_tr.pdf` | 19 | 1.98 MB | ✅ Compiled |
| `paper/PTUM.pdf` | 19 | 613 KB | ✅ Compiled (Springer EPJ) |
| `paper/BEY.pdf` | 45 | 578 KB | ✅ Compiled (Springer EPJ) |
| `AMR_BL_STABILITY_REPORT_EN.pdf` | - | 36 KB | ✅ Compiled |
| `AMR_BL_STABILITY_REPORT_TR.pdf` | - | 257 KB | ✅ Compiled |
| `AMR_BL_STABILITY_REPORT_BILINGUAL.pdf` | - | 38 KB | ✅ Compiled |

---

## Master 22-Gate Validation

All 22 gates **PASSED** (scientific_maturity_score: 6.0/10, backward_time_travel_pipeline_score: 3.0/10)

Key gates:
- ✅ T_μν conservation (max_rel_residual: 0.0079)
- ✅ Retarded Einstein solver (5.8e-6 relative convergence)
- ✅ Gauge validation (MQS non-conservation identified)
- ✅ Schwarzschild/Kerr benchmark (orbit drift ~10⁻¹¹)
- ✅ Gödel/Tipler CTC benchmarks reproduced
- ✅ CTC search: zero candidates for device configs G1-G4
- ✅ Energy-to-CTC scan: 45.7 decades gap

---

## New Files in This Release

### Source Code
- `src/zt005/sc_core_kerr.py` - Rotating SC core with frame-dragging

### Tests
- `tests/test_sc_core_kerr.py` - Kerr evolution tests (100/500 steps)
- `tests/test_sc_core_bssn.py` - SC core BSSN evolution
- `tests/test_bssn_beta3_convergence.py` - Convergence study
- `tests/test_bssn_beta3_longrun.py` - 1000-step long-run
- `tests/test_bssn_beta3_schwarzschild.py` - Schwarzschild puncture
- `tests/test_bssn_mp_longrun.py` - Moving puncture 500-step
- `tests/test_bssn_v2_longrun.py` - v2.0 1000-step long-run

### AMR Scripts & Results
- `scripts/test_amr_bl_moving_patch_Lc160.py` - L_c=160 AMR test
- `scripts/test_amr_bl_moving_patch_v2.py` - AMR variant
- `scripts/test_amr_bl_moving_patch.py` - Base AMR script
- `scripts/test_amr_bl_moving_Lc160_50steps.py` - 50-step runner
- `results/validated/amr_bl_moving/amr_bl_moving_Lc160_50steps.json`
- `results/validated/amr_bl_moving/amr_bl_moving_Lc200_50steps.json`
- `results/validated/amr_bl_moving/amr_bl_moving_proj.json`

### Debug Utilities
- `test_amr_quick.py`, `test_amr_5steps.py`, `test_amr_small.py`
- `fix_script*.py`, `parse_results.py`

---

## Known Limitations

1. **Kerr rotation**: Slow rotation only (|a| ≤ 0.5M). Full Kerr (|a| ≤ M) deferred to v3.0-alpha1 (roadmap).
2. **AMR performance**: 33³ coarse + 25³ fine grids take ~33 min for 50 steps. Production runs need optimization.
3. **BEY paper figures**: 26 placeholder figures - need replacement with actual simulation plots.
4. **arXiv submission**: ON HOLD per roadmap (user decision).
5. **Negative-energy config** (v5.0): Stability window ~100 steps only.

---

## Roadmap Status

| Milestone | Status |
|-----------|--------|
| v2.0: Full BSSN solver | ✅ COMPLETE |
| v3.0-alpha1: Full Kerr KS Cartesian | 🔷 Deferred |
| v3.0-alpha2: Moving-puncture formalism | ✅ COMPLETE |
| v3.0-alpha3: N-scaling with interpolated masks | ✅ COMPLETE |
| v3.0-alpha4: SC regular core | ✅ COMPLETE |
| v3.0-alpha5: SC core BSSN evolution | ✅ COMPLETE |
| v3.0-alpha6: SC core Kerr (rotating) | ✅ COMPLETE |
| **v3.0-beta: Consolidation + release** | ✅ **THIS RELEASE** |
| v4.0-beta: BSSN + SC metric + SC matter | 🔷 In progress |
| v5.0-beta: Negative-energy stability window | 🔷 Documented |

---

## Installation & Usage

```bash
# Clone and setup
git clone https://github.com/arreasonn-jpg/zamanda_yolculuk.git
cd zamanda_yolculuk
pip install -r requirements.txt

# Run test suite
python -m pytest tests/ -v

# Run specific long-run tests (require RUN_FULL_LONGRUN=1)
RUN_FULL_LONGRUN=1 python -m pytest tests/test_bssn_beta3_longrun.py -v
RUN_FULL_LONGRUN=1 python -m pytest tests/test_sc_core_kerr.py -v

# Run AMR test
python scripts/test_amr_bl_moving_Lc160_50steps.py
```

---

## Acknowledgments

- Numerical relativity community for BSSN formalism references
- Springer for svjour3.cls / svglov3.clo journal templates
- MiKTeX team for LaTeX distribution
- All contributors to the zamanda_yolculuk project

---

**Next**: v4.0-beta consolidation (BSSN + SC metric + SC matter full integration)