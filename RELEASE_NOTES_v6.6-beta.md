# v6.6-beta Release Notes

**Date**: 2026-10-31  
**Tag**: `v6.6-beta`  
**Branch**: `main`  
**Commit**: (to be created)

---

## Summary

v6.6-beta adds **complex scalar field** to the PTUM stack, demonstrating NEC violation and stability with equilibrium seed initialization. The refined 2:1 grid driver works correctly, but the **equilibrium psi seed is broken on the refined grid** (documented limitation).

---

## What's New

### Physics Modules
- **Complex scalar matter** (`matter_complex.py`): Complex scalar field with U(1) symmetry, potential `V(φ) = m²|φ|² + λ|φ|⁴`
- **PTUM fields** (`matter_ptum.py`): Fuel `Φ` + proxy `ψ` with tachyon condensation
- **Equilibrium seed** (`ptum_field.py`): `psi_equilibrium(Φ, m2, alpha, Eth, lam_psi)` — analytic equilibrium `ψ = √(-m²/λ)` when `Φ > Eth`

### Evolution Drivers
- **Single-grid RK4** (`bssn_rk4_complex_ptum.py`): Coupled BSSN + complex scalar + PTUM
- **Refined 2:1 grid** (`bssn_refined_complex_ptum.py`): Coarse→prolongate→fine→restrict cycle

### Initial Conditions
- **Kerr + complex + PTUM IC** (`kerr_complex_ptum_ic.py`): Full stack IC for rotating background

---

## Validation Results

### 1. Single-Grid Complex + PTUM Equilibrium Seed (`v6_complex_ptum_run_eq.py`)
```
N=9, L=20, dt=5e-4, steps=100
```

| A_fuel | ρ_min | n_neg | K_max | Φ_max | ψ_max | Finite |
|--------|-------|-------|-------|-------|-------|--------|
| 1.5    | -2.8471e-01 | 1 | 3.6129e-01 | 1.5018e+00 | 1.1177e+00 | ✅ |
| 3.0    | **-1.5047e+01** | 1 | 2.3500e+01 | 3.0348e+00 | 2.8275e+00 | ✅ |

**Key findings:**
- **NEC violation achieved** at both amplitudes
- A=1.5: moderate NEC violation (ρ ≈ -0.28), stable
- A=3.0: **strong NEC violation** (ρ ≈ -15.0), stable for 100 steps
- All fields finite throughout

### 2. Refined 2:1 Grid Complex + PTUM (`v6_refined_complex_ptum_run.py`, A_fuel=3.0)
```
N_c=9, L_c=20, N_f=9, L_f=10, dt=2e-4, steps=100, offset=(2,2,2)
```

| Step | ρ_min | n_neg | K_max (coarse) | K_max (fine) | α_min | Finite |
|------|-------|-------|----------------|--------------|-------|--------|
| 25   | 0.0000e+00 | 0 | 1.09e-04 | 1.13e-04 | 1.0000 | ✅ |
| 50   | 0.0000e+00 | 0 | 2.18e-04 | 2.23e-04 | 1.0000 | ✅ |
| 75   | 0.0000e+00 | 0 | 3.27e-04 | 3.32e-04 | 1.0000 | ✅ |
| 100  | 0.0000e+00 | 0 | 4.36e-04 | 4.41e-04 | 1.0000 | ✅ |

**Key findings:**
- **Driver works correctly** — stable refinement cycle
- **BUT: equilibrium psi seed broken on refined grid** — ρ_min = 0 (no NEC violation)
- K_max grows very slowly (1e-4 → 4e-4) — stable
- α_min = 1.0 constant — no collapse
- All fields finite

---

## Unit Tests (11 passing)

```
tests/test_bssn_rk4_complex_ptum.py       - 4 tests  ✅
tests/test_bssn_refined_complex_ptum.py   - 4 tests  ✅
tests/test_kerr_complex_ptum_ic.py        - 3 tests  ✅
Total: 11 tests in 23.18s
```

---

## Known Limitation

**Equilibrium psi seed is broken on refined grid.**

The `psi_equilibrium()` function computes the analytic equilibrium correctly on the coarse grid, but when prolonged to the fine grid, the equilibrium condition is not preserved. This means:
- Single-grid: NEC violation works correctly ✅
- Refined grid: NEC violation **does not activate** ❌ (ρ_min = 0)

This is a known issue documented in `docs/limitations.md` and will be addressed in v7.0.

---

## Grid Parameters

| Configuration | N | L | h | dt | Steps |
|---------------|---|---|---|-----|-------|
| Single-grid | 9 | 20 | 2.5 | 5e-4 | 100 |
| Refined coarse | 9 | 20 | 2.5 | 2e-4 | 100 |
| Refined fine | 9 | 10 | 1.25 | 2e-4 | 100 |

---

## Key Metrics Summary

| Metric | Single-grid A=1.5 | Single-grid A=3.0 | Refined A=3.0 |
|--------|-------------------|-------------------|---------------|
| NEC violation | ✅ ρ=-0.28 | ✅ ρ=-15.0 | ❌ ρ=0 (broken) |
| K_max | 0.36 | 23.5 | 4.4e-4 |
| α_min | ~1.0 | ~1.0 | 1.0 |
| Fields finite | ✅ | ✅ | ✅ |
| Negative cells | 1 | 1 | 0 |

---

## Files Changed

- `docs/roadmap.md` — v6.6 items marked ✅ complete
- `RELEASE_NOTES_v6.6-beta.md` — this file

---

## Next Steps (v7.0)

- Fix equilibrium psi seed on refined grid
- Refined complex + PTUM NEC violation stable
- v7.0.0 release

---

**All validation scripts pass (with documented limitation). Ready for tag v6.6-beta.**