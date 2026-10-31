# Independent validation record

Three independently prompted subagents reviewed the work. They had no shared context beyond the repository files. Their full working files are in the session scratchpad; the findings and resolutions are summarised below.

## 1. Blind first-order re-derivation

**Method.** The agent derived the first-order radiation-era solution and Θ⁽¹⁾ *before* reading `approach/` or `verification/`. It used Newtonian gauge, with its own SymPy Einstein-tensor check (all 16 components of G^μ_ν − 8πG T^μ_ν, at zeroth and first order).

Its independent results:

- Φ = Ψ = −2ζ (sin x − x cos x)/x³, with x = kη/√3.
- δ = ζ[−4 cos x + 8 sin x/x + 8 cos x/x² − 8 sin x/x³].
- Velocity amplitude U = √3 ζ (sin x + 2 cos x/x − 2 sin x/x²).
- Transfer function: T_l = A(x\*) j_l(z\*) − U(x\*) j_l′(z\*) + ∫ 2 Φ_x j_l + δ_{l1} U(x_o)/3, where A = δ/4 + Φ and z\* = √3 (x_o − x\*).
- Spectrum: C_l⁽¹¹⁾ = 4π ∫ (dk/k) Δ² T_l², which gives 4π A_ζ T_l(k₀)² for the monochromatic spectrum.

**Agreement proof (exact).** Three pointwise identities, with s = 2 sin x/x − cos x:

- U = √3 ζ s′
- δ/4 + Φ = −ζ s″
- s‴ + s′ = −2 ∂ₓ(Φ/ζ)

Using these, ∂_{x\*}(T_l^{blind} − T_l^{ours}) = 0 for symbolic l(l+1), by the spherical Bessel ODE alone. The coincidence limit then gives zero for every l ≥ 1. The l = 1 observer-Doppler term cancels exactly against the emitter term.

**Verdict.** Exact agreement for every l ≥ 1, including the sign and the C_l normalisation. The l = 0 term differs only by monopole terms, which the normalisation removes. As a supplementary check, mpmath at 30 digits gave differences of ≈10⁻³¹ for l = 1, 2, 3, 7, 12. The closed-form T₂ spot values matched to all 18 printed digits.

**Minor remarks, resolved.** The agent noted an unused variable and an awkward β = 0 branch in an earlier version of `transfer_closed_form.py`. Both were removed when the partial-fraction routine was rewritten.

## 2. Physics referee

**Verdict.** Accept as an approach, with major revisions. The referee independently confirmed the following:

- the Kelvin/irrotationality lemma (by hand);
- the zero-shift metric form with N = 𝒯/T;
- the exact redshift formula, and the cancellation of T_fluid(o) in the normalised map;
- the Fermi–Walker ⇔ parallel-transport lemma. Its own SymPy check used a non-diagonal γ, unconstrained det γ and arbitrary N.
- the nonlinear O(t²) initial data, and the identification of ζ with δN;
- the normalisation expansion Θ⁽¹⁻³⁾, via an independent 4-point discrete-sky model;
- the Wick factors 2 and 3;
- the Legendre image (1/3 − μ²) → w_l, for l = 0–5;
- the C⁽¹¹⁾ normalisation, the Green's functions, and the tensor resonance at μ = −1/3;
- the completeness of the I₂/I₃ term classes.

The referee also extended the field-equation check to det γ = e^{6Z(x,z)} (diagonal); all 10 equations pass.

### Major points and resolution

| # | Point | Resolution in `approach.tex` |
|---|---|---|
| M1 | The resonance and regularity analysis is incomplete. It misses scalar collinear resonances (μ = ±1), possible secular (∝ t) growth from non-decaying temperature-frame sources, removable ln\|α\| terms in Si/Ci closed forms, and the K̂K̂ discontinuity at K → 0. "Continuous, bounded" was asserted, not proved. | The "Exchange of limits" paragraph is rewritten as items (a)–(e): what is claimed, and what must be proved. Open items now include the resonances and secular terms. New check V20 covers resonance recombination. |
| M2 | Definition 1 assumed uniqueness. The homogeneous-solution list omitted the E → E + E₀(X) relabelling. The residual conformal group should be fixed, not merely used as a check. | The definition now defers to Prop. 1, which is still Stage 2. Relabellings are added to Prop. 1. Dilations are fixed by the zero-mean Gaussian ζ and SCTs by homogeneity; their responses remain as checks V9 and V10. |
| M3 | Late-time observer-local acoustic terms cancel at first order only through ∫(1/√3 − μ)P_l dμ = 0 (l ≥ 2). The second-order analogue is not guaranteed. | Added as an open risk, with new check V19. |

### Minor points and resolution

| # | Point | Resolution |
|---|---|---|
| M4 | The "t_i → 0 commutes" claim was vacuous, because t_i drops out identically. The referee also pointed out a missed strength: ζ = (1/6) ln det γ is exactly conserved. | t_i logic rewritten. Exact conservation of ζ added as its own paragraph. |
| M5 | Only the combined C_l is orientation-invariant, not C⁽²²⁾ and C⁽¹³⁾ separately. | Caveat added after the orientation proposition; V12 now tests the combination. |
| M6 | Wording issues: "aberration by the shear"; the isotropic-emission frame was not declared; the β list was wrong. | Fixed. Emission isotropy in the fluid rest frame is declared. The "boost to dipole frame" reading of dipole subtraction is explicitly excluded. |

### Overclaims and resolution

| Overclaim | Resolution |
|---|---|
| The field-equation identities were checked only on restricted families. | The status now states the restriction and cites the covariant derivation; V1 is scoped accordingly. |
| Eq. (ray) is printed, not machine-checked. | V6 now says it was hand-derived (independently by the referee). |
| V9 was listed as "r = 1 done", but no script implements it. | Changed to "planned"; V14 is the r = 1 analogue. |
| "Constraints hold at every order" was stated. | Scoped to r = 1 and nonlinear O(t²). |
| The closed-form T_l exists only per fixed l. | The scope is stated: l = 2, 3 done; symbolic l remains the integral form. |
| §8 cited VALIDATION.md before it existed. | This file now exists. |

Gaps noted by the referee — the main deliverable not yet executed, no cheap second-order anchor, and no cost estimate — are now stated in §8 and §9. A cost estimate and an induced-GW / Poisson-gauge anchor have been added.


## 3. Adversarial code audit (mutation testing)

*(pending — filled in below when received)*
