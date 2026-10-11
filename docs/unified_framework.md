# Unified numerical framework: SC + PTUM

This document describes the shared numerical infrastructure behind two
companion physics papers:

- **SC** (Esa 2026a): Schwarzschild Condensation — QCD-inspired
  regular-core model for black-hole interiors.
- **PTUM** (Esa 2026b): Phase Transition Universes Model —
  discharge-based cosmological creation scenario.

Both papers use the same BSSN solver developed in this repository.
This document explains what is shared, what is model-specific, and
what is explicitly *not* claimed.

## What is shared (code)

| Module | Purpose | Used by |
|---|---|---|
| `bssn_state.py` | BSSN variable container | Both |
| `bssn_ricci_full.py` | Physical 3-Ricci scalar | Both |
| `bssn_full.py` | Covariant Hessian + D^2 alpha | Both |
| `bssn_rk4_beta1.py` | RK4 integrator | Both |
| `matter_scalar.py` | Real scalar stress-energy | Both |
| `matter_complex.py` | Complex scalar (U(1)) | Extensions |
| `matter_ptum.py` | Two-field PTUM stress-energy | PTUM |
| `bssn_rk4_ptum.py` | RK4 for (metric + 2 PTUM fields) | PTUM |
| `sc_like_ic.py` | SC smoothed core IC | SC |
| `sc_full_ic.py` | SC metric + SC condensate | SC |
| `refinement.py` | Two-level AMR-lite operators | Both |
| `bssn_refined_ptum.py` | Refined PTUM driver | PTUM |
| `quantum_inequality.py` | Ford-Roman diagnostic | PTUM |

## What is model-specific (physics)

### SC

- Interior: de Sitter-like with QCD Cooper pair condensate
- Order parameter: composite boson (diquark) from CFL pairing
- Potential: `V(ψ) = ½m_pair²|ψ|² + ¼g_int|ψ|⁴`
- Scales: `m_pair ~ 10 MeV`, `ρ_core ~ 10¹⁹ kg/m³`
- Physical claim: regular non-singular interior for compact objects

### PTUM

- Pre-discharge: metastable fuel field `Φ` accumulating
- Discharge: tachyonic activation of proxy field `ψ`
- Potential: `V = V_fuel(Φ) + V_proxy(ψ; Φ)`
- Scales: `V₀^{1/4} ~ 10¹⁵-10¹⁶ GeV`
- Physical claim: hot Big Bang as a finite-energy discharge event

## What is NOT claimed

**The two models are NOT unified into a single theory.**

The relationship between SC and PTUM is:

1. **Shared methodology** — both use effective field theory with
   scalar order parameters, evolved on the same BSSN infrastructure.
2. **Shared mathematical structure** — both involve tachyonic
   instabilities, condensate formation, and BSSN-evolved dynamics.
3. **Shared diagnostic toolkit** — Ford-Roman quantum inequality,
   refinement operators, PTUM/SC IC constructors.

The two models operate at different physical scales:

| | SC | PTUM |
|---|---|---|
| Length scale | ~km (stellar BH) | ~Gpc (cosmological) |
| Density scale | 10¹⁹ kg/m³ | 10⁹⁷ kg/m³ |
| Energy scale | ~10 MeV | ~10¹⁵ GeV |
| Temperature | ~10¹¹ K | ~10²⁷ K |
| Order parameter | CFL diquark | Tachyon condensate |

A 78-order-of-magnitude scale gap separates them. No derivation
currently connects the two.

## Speculative bridge (future work)

The two models share a **tachyon-like activation mechanism**:

- SC: QCD deconfinement → Cooper pairing
- PTUM: fuel field crosses Φ_c → tachyon activates

Whether these are *the same physics* at different scales is an
**open question**, not a claim of this repository. Any unification
attempt must:

1. Provide a bridge mechanism between the scales
2. Survive the 78-order scale gap
3. Be compatible with both numerical results

None of these is currently satisfied.

## Numerical results summary

### SC (v4.0)

| Setup | Result |
|---|---|
| SC regular core (metric only) | 50-step stable, H growth 1.0000 |
| SC metric + SC condensate | 50-step stable, alpha_min 0.8612 |

### PTUM (v5.0-v6.0)

| Setup | Result |
|---|---|
| PTUM A=3.0 single grid | crashed ~100 steps |
| PTUM A=3.0 refined (2:1) | 1000-step stable, rho_min = -0.097 |
| Ford-Roman massless | passes (ratio max 0.30) |
| Ford-Roman massive mu=1 | borderline (ratio 1.21) |

## How to use this framework

For a new model with scalar order parameter `Ψ`:

1. Define `T_μν(Ψ)` — copy `matter_scalar.py` or `matter_ptum.py`
2. Build initial data — copy `sc_like_ic.py` or `sc_condensate.py`
3. Choose driver:
   - `bssn_rk4_matter.py` for one field
   - `bssn_rk4_ptum.py` for two coupled fields
   - `bssn_refined_ptum.py` for two-level refinement
4. Run and check constraint preservation
5. If negative energy appears, run `quantum_inequality.py`

## References

- Esa, E. (2026a). Schwarzschild Condensation. In preparation.
- Esa, E. (2026b). Phase Transition Universes Model. In preparation.
- Ford, L. H. & Roman, T. A. (1997). Phys. Rev. D 55, 2082.
- Bruegmann, B. (2005). Phys. Rev. Lett. 96, 111101.
- Baumgarte, T. W. & Shapiro, S. L. (2010). Numerical Relativity.

## Reproduce

    pytest tests/ -m "not slow"    # ~2 min
    python scripts/v4_sc_full_run.py
    python scripts/v6_refined_long.py
