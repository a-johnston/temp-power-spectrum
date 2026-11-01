# Unified plan: analytical temperature power spectrum through fourth order

This plan merges `orange/analysis/temperature-spectrum-plan` (`APPROACH.md`, `derivations/`,
`validation/`) and `purple/approach-proposal` (`approach/approach.tex`, `verification/`). Both
branches are kept intact in this tree's history. Where they disagree, this document is
authoritative. The execution of the plan is in
[`derivations/UNIFIED_DERIVATION.md`](derivations/UNIFIED_DERIVATION.md), with exact checks in
`solution/`.

## 0. Why the merge is natural

The two branches use **the same gauge**. Orange's comoving-orthogonal, zero-shift chart with
`N ρ^{1/4}` spatially constant is purple's temperature frame `g = N²(−dt² + γ_ij dX^i dX^j)` with
`N = 𝒯/T_phys`. The dictionary is `χ = ν`, `q_ij = e^{2ν} γ_ij`. Their linear solutions and their
nonlinear `O(t²)` data agree term by term. The merge therefore keeps one coordinate system and
takes, element by element, whichever treatment is stronger or simpler.

| Element | Taken from | Reason |
|---|---|---|
| Temperature frame; exact redshift `T_obs = T_fluid(o) exp[−∫σ_ij e^i e^j dt]`; observer temperature cancels in Θ | purple | Replaces the SW/ISW/Doppler/endpoint bookkeeping with one integral. Exact (lemma L2–L3). |
| `ζ = ⅙ ln det γ` exactly conserved in fluid labels; continuity solved by `γ = e^{2ζ} exp H`, `tr H = 0` | purple | The Gaussian variable is readable at every time, not only as `t → 0`. |
| Marked conformally Cartesian material chart; conformal flatness as a *physical* no-tensor-seed condition | orange | Fixes the residual relabelling freedom physically rather than by a gauge argument. |
| Zero-total-momentum sector solved **without projectors** | orange | Purple's SVT split needs `∇⁻²`. That fails on the whole `C^(13)` domain, because `K^(3)(k,q,−q)` contains `K^(2)(q,−q)` at exactly `K = 0`. |
| Constraint propagation (Bianchi + exact conservation) to justify using a subset of equations | orange | Turns "solve some equations" into a proof that all ten hold. |
| Implicit-function coefficients for the emission root; fixed proper time from a `T_phys = T_I` event | orange | Exact through third order; finite `T_I` is kept (the `T_I → ∞` limit is optional). |
| Fermi–Walker ⇔ parallel transport in the temperature frame (L5); realisation-wise rotation invariance of `Σ_m |a_lm|²` | purple | Reduces the tetrad to a rotation question. |
| Wick factors 2 and 3; the internal `(q,−q)` pairing of the 13 term is **kept** (no normal ordering) | both (identical) | |
| Monochromatic spectrum `P(k) = 2π²σ²/k0² · δ(k − k0)`, holding `σ² = ⟨g²⟩` fixed | orange | Purple's `approach.tex` had `k0³`, which is dimensionally wrong. Corrected there too. |
| Integral-legality and shell-restriction conditions; separation of evaluated and unevaluated integrals | orange | |
| Resonance and secular-term inventory (scalar `μ = ±1`, tensor/vector `μ = −1/3`); dilation Ward identity | purple | |

**New in the unified plan.**

1. **Mixed trace-free evolution form.** Raising an index with `γ` makes the `σσ` terms of the
   spatial Einstein equation cancel exactly:
   `TF[S′ + 2n′S + R^k_j − 2D^kD_j n + 2D^k n D_j n] = 0`, with `S = γ⁻¹σ`. Proved on generic jets.
2. **Triangular solution scheme** at every order `r`, with the same linear operators at each order:
   - `ν^(r)` from the `R_00` wave equation `3ν″ + K²ν − 6ν/t² = −c_N`. Its free regular
     coefficient of `t²` is fixed *exactly* by the Hamiltonian constraint, `ν = −t² R[e^{2ζ}δ]/36 + O(t⁴)`.
   - The scalar and vector parts of `H^(r)` follow from the momentum constraint by **single time
     integrals**.
   - The tensor part follows from a TT wave equation.
   - At `K = 0` there is no projector: `½(A″ + 2A′/t) = −c_T`.
   - All elementary Green functions: the lapse kernel uses `x j₁(x)`, `x y₁(x)`, and the tensor
     kernel is `(t′/t) sin K(t−t′)/K`.
3. **Frame lemma, strengthened.** The Fermi–Walker frame with initial symmetric-root orientation
   equals the symmetric square root `γ^{−1/2}` at the observation event **through first order**.
   So `Θ^(1)` and `Θ^(2)` are unchanged, and `Θ^(3)` changes only by an infinitesimal rotation of
   `Θ^(1)`. That leaves `C^(22)` and `Re C^(13)` *individually* unchanged, not only their sum. All
   kernels can therefore be built with `γ_o^{−1/2}`, and no transport integral is needed.
4. **Observation time by composition.** Build the coordinate-event map `Θ_c(n̂; t_o)` and substitute
   the random `t_o = t̄_o + δt_o` afterwards: `Θ^(2) = Θ_c^(2) + 2D₁∂Θ_c^(1)`, and so on.
5. **Sky-space spectrum formulas** for the monochromatic case. These avoid derivatives of
   spherical Bessel functions entirely, e.g.
   `C_ℓ^(13) = (3σ⁴/2) i^ℓ T_ℓ(k0) ∫dμ [∫dΩ P_ℓ(ẑ·n̂) 𝒦₃(k0ẑ, k0q̂, −k0q̂; n̂)]*`.

## 1. The experiment (frozen conventions)

- Signature `(−,+,+,+)`, `c = 1`. Background `a = a₁t`, `T̄ = 𝒯/(a₁t)`, `8πGα𝒯⁴ = 3a₁²`. The
  coordinate `t` is the fluid velocity potential (`T u_μ = −𝒯 ∂_μ t`, `t → 0` as `T → ∞`).
- **Gaussian variable:** `ζ(X) = ⅙ ln det γ_ij(t, X)` in the marked conformally Cartesian fluid
  labels. Exactly `t`-independent. It equals the uniform-density (δN) curvature perturbation at
  the primordial boundary. `ζ = εg` *exactly*, with `g` a homogeneous, isotropic, zero-mean
  Gaussian field.
- **No independent higher-order modes** means `ζ^(r≥2) ≡ 0`, together with:
  - **(T)** the limiting 3-metric is conformally flat, so there is no tensor seed;
  - **(V)** `d(Tu) = 0`, so there is no vorticity seed (the constant-vorticity mode is regular and
    must be excluded explicitly);
  - **(R)** `ν`, `H` are even analytic in `t` at fixed momenta, with `H(0) = 0` and `ν(0) = 0`.
    This excludes every decaying mode.
- The forced higher-order metric and fluid data are **not** zero. Examples: `γ^(2)(0) = 4g²δ`, and
  `ν^(2) = (t²/9)[2(k₁²+k₂²) − k₁·k₂] + O(t⁴)`.
- **Emission:** isotropic Planck emission in the fluid rest frame on `T_phys = T⋆`. The intersection
  branch is the one continuously connected to the background.
- **Observer:** fluid line `X = 0`. Initial event `T_phys = T_I > T⋆` on that line. Observation
  after fixed proper time `τ_o`. The FW tetrad is initially `γ(t_I, 0)^{−1/2}` (symmetric root,
  aligned with the marked labels).
- **Fourier and factorials:**
  - `g(X) = ∫_k g_k e^{ik·X}`, with `∫_k = ∫ d³k/(2π)³` and
    `⟨g_k g_k′⟩ = (2π)³δ³(k+k′)P(k)`.
  - Fields are `Y = Σ ε^r Y^(r)/r!`, with `Y^(r) = ∫_{k₁..k_r} K_r(k₁..k_r) g..g e^{iK·X}`. The
    symmetric kernels carry **no** extra `1/r!`. Equivalently, `K_r` is the multilinear coefficient of
    `z₁…z_r` when `g = Σ z_a e^{ik_a·X}`.
  - Complex `Y_ℓm` (Condon–Shortley).

## 2. Stages and gates

| Stage | Content | Gate (exact, symbolic) | Status |
|---|---|---|---|
| S0 | Exact framework: lemmas L1–L5, field equations, mixed-TF form | generic 2-jets | **done** |
| S1 | First order; `T_ℓ`; `C^(11)` | all 10 equations; Newtonian-gauge equivalence for all ℓ; Si/Ci closed form (ℓ = 2, 3) | **done** |
| S2 | Second-order fields: sources, triangular solution, `K = 0` sector, finite-`t_i` data | operators confirmed; all 10 equations as exact Taylor series (general `K` and `K = 0`); `c_M ≡ 0` at `K = 0`; dilation Ward identity in closed form | **done** (integral representation; Si/Ci reduction mechanical, not generated) |
| S3 | Third-order fields at `(k,q,−q)` | sources extracted exactly; operators identical to second order; all-order consistency proof (linearised Bianchi + regularity, derivation §4.4) | **done**; an explicit third-order Taylor residual test is **open** |
| S4 | Observable: ray, emission root, clock, frame, endpoint terms through third order; normalisation; composition with `δt_o` | every formula against direct expansion for arbitrary fields; normalisation on exact discrete skies | **done** (operator level) |
| S5 | Harmonic projection, Wick contraction, monochromatic reduction | Wick pairings and Fourier/shell measures (orange `check_algebra.py`) | **done** (formulas) |
| S6 | Evaluation of the resulting nested integrals; asymptotics; convergence proofs | dominated convergence, resonance recombination | **open** |

## 3. Validation matrix (merged)

Every row is an exact symbolic identity. Numerical spot checks may appear only as labelled
supplementary output.

| Gate | Where |
|---|---|
| Lemmas L1–L5, generic field equations | `verification/check_exact_lemmas.py`, `check_fw_conformal.py`, `check_field_equations_generic.py` |
| Mixed trace-free evolution form | `solution/check_evolution_form.py` |
| First order, two gauges, transfer-function equivalence | `verification/check_first_order.py`; engine re-check `solution/first_order_engine_check.py` |
| Quadratic/cubic data in orange's variables (cross-check of the same gauge) | `validation/check_initial_data.py`, `check_cubic.py`, `check_offshell.py` |
| Second-order operators and sources | `solution/second_order.py` |
| All ten equations at second order, general `K` and `K = 0` | `solution/second_order_consistency.py` |
| Third-order sources and operators | `solution/third_order.py` |
| Dilation Ward identity (relabelling = residual gauge) | `solution/check_ward_dilation.py` (with negative controls) |
| Ray, emission, clock, endpoint, composition through third order | `solution/check_observable_expansion.py` |
| Normalisation `Θ^(1..3)` | `verification/normalization_expansion.py`, `validation/check_algebra.py` |
| Wick pairings, angular and shell measures | `validation/check_algebra.py` |

## 4. Open items (carried forward explicitly)

> **Update.** See `report/temperature_spectrum_report.tex` for how items 1–6 were resolved, and
> for what is still open:
> - evaluation of `C^(22)` and `C^(13)` at general `k0`;
> - the late-time limit at fourth order;
> - a rigorous ultraviolet bound.
>
> The original list follows.

1. **Third-order implementation test.** Consistency is proved (derivation §4.4). An explicit
   Taylor-series residual test at third order, like the second-order one, is still to be run;
   it is expensive.
2. **Si/Ci closed forms.** Reduce the second-order kernels to closed form, as was done for
   `T_ℓ`. Third-order kernels remain one-dimensional integrals over the second-order closed forms.
3. **Evaluate `C^(22)` and `C^(13)`.** Evaluate the monochromatic nested integrals analytically
   where possible. Resonances: removable `ln|α|` at scalar `μ = ±1` and tensor/vector `μ = −1/3`,
   which must be recombined before the angular integral. Possible secular growth must be shown to
   cancel or be kept.
4. **Asymptotics.** Late time `x_o → ∞` is non-uniform at the resonances. Only the
   angle-first order of limits is justified. `x⋆ → 0` and `x⋆ → ∞` are not done.
5. **Non-monochromatic spectra.** Infrared behaviour of `C^(22)` and `C^(13)` separately; the
   dilation identity predicts the cancellation pattern. Shell-restriction existence for smooth
   spectra.
6. **Perturbative domain.** Uniqueness of the emission intersection and absence of caustics hold
   for the asymptotic series, not as a convergence statement in `ε`.
