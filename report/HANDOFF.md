# Handoff: state of the cubic/fourth-order work and the open follow-ups

Branch: `claude/magical-fermat-ufand5`. Commits from this session, after `ae83458`:

| Commit | Content |
|---|---|
| `d87b201` | Review of 5cdf048..ae83458; new report section §"Stationary reduction of the complete cubic response" (`report/stationary_cubic_reduction.tex`); constraint rows remove hyperbolic Green kernels (`solution/check_constraint_rows.py`) |
| `93a7565` | Third-order constraint compatibility; Prop. "Combined cubic response" |
| `3281705` | Finite-wavelength forward/assembled validation code (`solution/numerics/`) |
| `33b64d5` | Report text for the validation, observer-centred split, shell C13 recomputation (appendix), κ⁴ numerical check |

The PDF in the repo was **not** regenerated. The TeX compiles with `pdflatex` (installed ad hoc in the
session; three passes), with no errors or undefined references and one negligible 0.3 pt overfull line.

## 1. What a final reviewer of the whole report should know

### 1.1 Review of the earlier agent's work (5cdf048..ae83458)

- **Checks:** all 14 new exact checks pass; they are listed in `README.md` and `solution/run_all.sh`.
- **Identities re-derived by hand:**
  - the polarised stationary identity (the 2n+1 rule);
  - the normalisation adjoint;
  - the leading cubic coefficient −12τ²(g²A + gB) and the shell value 4(o²−e²)/9;
  - the spatial Green operator polynomial (all three eigenvalues);
  - the clock primitive, `W_c = (o⁴+3ι⁴)/4`, and the observer trace;
  - the first-order cancellation of hyperbolic terms;
  - the factors in the dilation channel;
  - the initial kernels `c₂`, `c₃`.
- **Outcome:** no errors found. Their markdown notes (`report/*.md`) are condensed into the new report
  section.

### 1.2 Results added in this session (all in `report/stationary_cubic_reduction.tex`)

1. **Source compatibility.** This is Prop. `prop:compat`, consisting of eqs. `compat-1`, `compat-2`. It
   follows from the contracted Bianchi identity:
   - `c_Ham' − 4c_N/t = 2iK·c_M`;
   - `c_M' + 2c_M/t = iK·c_T − (iK/6)(c_Ham+4c_N)`.

   It is verified in two ways:
   - exactly, at second order, for unequal |k₁|≠|k₂| at arbitrary angle (the identity fails off shell);
   - by exact series through t¹⁰, at third order with physical second-order fields, for noncoplanar,
     zero-pair and near-zero-pair triples (with a negative control).
2. **Constraint-row response** (Prop. `prop:constraint-response`). For K≠0,
   `K·A = V = ∫₀ᵗ[4K(ν'−ν/s) + 2i c_M]`, with non-TT(A) = 𝒫_K[V]. Only acoustic lapse kernels, the TT
   wave kernel and one cumulative integral appear; there are **no hyperbolic (−K²/3) kernels**. The
   initial Hamiltonian row is essential. The trade-off: a removable 1/K⁴ projector near K=0, where the
   Cartesian form is used instead.
3. **Constraint-row transpose** (eq. `constraint-transpose`), with explicit observer and emission point
   terms. It applies to the single-pair stationary terms at the probe momentum.
4. **Acoustic shell scalar amplitude** (Prop. `prop:shell-acoustic`):
   - `S₂ = 18c_r[T₂−2T₃/t] − ∫𝒢 C_N + (1−3μ)[3𝖻/2+3𝖺²−6ℐ] + (3μ²−1)𝖻²/8 − 3𝖺²`;
   - `ℐ = ½Cin(2x) + sin2x/(2x) − sin²x/(2x²) − ½`.

   This supersedes the hyperbolic `G_S` and the nested 𝒱_m in eqs. `moment-S`/`moment-V`, which remain
   correct. It is checked against the hyperbolic solution through t¹⁶.
5. **Combined cubic response** (Prop. `prop:combined`):
   - **Formula:** `J_abc = J₁₁₁ − Σ_{a|bc} S_meas,a[y_bc] + J_lin[y_abc]`.
   - **What combines:** the four background blocks (the direct third residual and three field-variation
     parts) collapse, by compatibility, into one pairing with the full cubic source.
   - **J_lin in constraint rows:** acoustic Wronskian + local momentum term + TT wave term.
   - **What it removes:** first field adjoints, the second ray, the third-order field and hyperbolic
     kernels. Closed-form second-order fields enter instead.
   - **Time arity:** at most 2 except J₁₁₁, which has 3.
6. **Finite-wavelength validation** (§`sec:stationary-validation`; numerical, labelled not evidence).
   - **Forward code:** `solution/numerics/forward_cubic.py` and `cubic_fields.py`. Fields come from (N)
     and the Cartesian (T) rows only; it includes the Fermi–Walker triad, geodesic, roots, clock and
     normalisation.
   - **Assembled side:** `assembled_cubic.py` evaluates eq. `combined` and uses the repository's own
     `direct_cubic_measurement.first_only_exponent_triple` for J₁₁₁.
   - **Agreement:** ≤3e-9 for three triples (noncoplanar unequal, zero pair, |k_b+k_c|=0.05) at |k|t̄_o≈3.
   - **Pieces checked separately:** J₁₁₁ (3e-8); constraint vs Cartesian J_lin (4e-14); the dilation
     channel (4e-8); the Hessian alone (1e-9).
   - **Controls:** the FW triad changes the value by 9e-5 and the clock by 0.24.
7. **κ⁴ shell coefficient** (eq. `stationary-four-gradient`). It is confirmed numerically by
   `check_kappa4_shell.py`: κ² to 1.3e-9, κ⁴ to 1.6e-7 at ι=0.5 and 1. The observer-local −567ι⁴/3402
   term is confirmed to 1e-6.
8. **Observer-centred split** (exact identity of §`sec:obs-dilation`):
   `K₃ = K^c₃ − Σ_pairs 𝔇K₂ − Σ_a 𝔇²K₁`.
   - Generic noncoplanar triple: K^c₃ is 0.9% of K₃.
   - (k,q,−q): 33%.
   - Shell average ⟨K^c₃⟩/K₁: −0.195 at x_o=0.3, against a total of 2.8955; −2.96 at x_o=1, against
     −6.156.
9. **Appendix.** An independent recomputation of (1/3)C₂⁽¹³⁾/σ²C₂⁽¹¹⁾ gives 2.8955, −6.156 and −67.29
   at x_o = 0.3, 1 and 3. The old exploratory table gives 2.895, −6.13 and −67.5, so they agree to the
   table's precision. C₂⁽¹¹⁾ = (4π/25)K₁² reproduces 0.0258 and 0.116.

### 1.3 Things to scrutinise in a final check

- **Prop. `prop:compat`, proof.** "The background Bianchi operator annihilates 𝓛_K[y] for every y" is
  argued, not separately checked off shell. The identities themselves are checked at orders 2 and 3.
- **Prop. `prop:combined`, step "the combined source is the order-three source with its initial
  row".** This relies on the stationary-identity duality. It is numerically confirmed (item 6) but not
  re-derived symbolically end to end.
- **Sign conventions.** S_meas is evaluated with averaged sky weights in `assembled_cubic.py`; the
  report uses full-dΩ weights (W = w/4π). Both conventions agree; see the comment in the code.
- **Convergence of the numerics.** Spectral time grids (44–64 nodes), a sky of 18×36 or 26×52 Gauss
  nodes, and five-point stencils for the clock composition (h=0.02). x_o=3 used 10 μ nodes, not
  convergence-tested separately; the κ⁴ fit uses κ = 0.04…0.1.
- **README and table.** `README.md` (section "Further cancellation research") and the report's status
  table were updated; check that they are consistent with the abstract.

### 1.4 External questions raised during the session

- **"A 2.2% disagreement in the k₀⁶ coefficient of the 13 term, observer-local" (another agents'
  report).** On our side the coefficient is checked by three exact symbolic routes and now numerically
  to 1.6e-7, including the ι⁴ term. The ι⁴ term is about 0.03% of the κ⁴ coefficient at typical times,
  so it cannot produce 2.2%. The likely causes are a convention difference: the physical initial event
  versus coordinate time, the clock composition, monopole normalisation, or the observer dilation. Our
  reference expression, from `check_stationary_quadrupole_four_gradients.py`, is
  `(4π/25)[J₁^[2]J₃^[4] + J₁^[4]J₃^[2]] = −2π(o²−e²)(2835ι⁴ + 9973o⁴ − 25440o²e² + 33920oe³ − 21288e⁴)/1913625`.
  Compare the time dependence term by term.
- **"Can a subhorizon process generate superhorizon modes falling off slower than k³?"** At finite
  times, no: kernels are entire in the Cartesian K (finite propagation speed), so P(K→0) is a constant
  and Δ² ∝ K³ for inputs with integrable correlations. There are two apparent exceptions, neither a
  counterexample:
  - an exactly monochromatic input, whose ξ ∝ sin(k₀r)/(k₀r), gives quadratic P ∝ 1/K (shell overlap),
    so Δ² ∝ K²; in this setup the homogeneous lapse kernel is nonzero;
  - the non-commuting late-time limit (k³ ln²k induced-tensor tails).

  This is not yet written into the report; it could become a remark in §"Limits" or §"General spectra".

## 2. Open follow-ups, in recommended order

### Done after the first version of this handoff

- **Exact first nonzero observer-centred shell coefficient** (eq. `centred-kappa4`, checked by
  `solution/check_observer_centred_kappa4.py`):
  `⟨K^c₃(k,q,−q)⟩ = k⁴[4(e−o)(12e³−8e²o+7eo²+7o³) + 5(o⁴−ι⁴)]/270 + O(k⁶)`.
  The κ² term cancels identically. The result agrees with the numerical value at x_o=0.3. It splits
  into a propagation part, which vanishes as e→o, and a pure observer-clock part, (o⁴−ι⁴)/54.
- **Remark `rem:k3`** (causal generation of superhorizon power; the answer to the author's question)
  is now in §"General spectra".

### 2b. Simplify the observer-centred cubic kernel (analytic)

**Status:** the κ⁴ step below is done (see above). The finite-wavelength step remains open.


The target is ⟨K^c₃(k,q,−q)⟩_q̂, the only genuinely cubic part of C⁽¹³⁾ after
Prop. `prop:obs-dilation`. The dilation part is already `½σ²D²C⁽¹¹⁾` plus a total derivative.

- **Leading gradient.** K^c₃ = 0 (checked by hand in this session from K^c₂ = −τ²TF(kq+qk)_nn: the
  terms are (8/3 − 16/3 + 8/3)τ² = 0).
- **Next step (cheap, exact).** Compute ⟨K^c₃⟩ at κ⁴ analytically:
  `J̄₃^[4] + 4Σ_±⟨K₂^[4](k,±q)⟩ + 16K₁^[4](k)`. Use eq. `stationary-four-gradient` for J̄₃^[4]; take J14
  for K₁^[4] from `check_stationary_quadrupole_four_gradients.py`; take the four-gradient quadratic map
  K₂^[4] from `curvature_kernels.py` / `check_cubic_curvature_hand.py`. Check it against
  `observer_centred.py` at small κ. If it vanishes or is a pure boundary term, that is a real
  cancellation result.
- **Finite wavelength.** Insert eq. `combined` into the K^c₃ identity. The dilation terms
  𝔇K₂ and 𝔇²K₁ are resolvent derivatives (∂_aG = −G ∂_aL G) of the linear and quadratic maps, so they
  may cancel analytically against parts of J_lin and S_meas. The natural decomposition is by soft legs:
  write y_bc = y_bc^(0) + (y_bc − y_bc^(0)), where y_bc^(0) is the response with one leg at zero momentum.
  Use the dilation source identities `C_{p0} = −D̂_pF₁` (direct cubic phase note) to identify which
  blocks are exactly dilations.

### 4. Evaluate the remaining shell spectrum pairings (analytic)

- **Already reduced.**
  - The top-order shell part is eq. `adjoint-nine-moment` with acoustic S moments
    (Prop. `prop:shell-acoustic`): two-time integrals of entire functions.
  - The S_meas[y_bc] shell part is the propagation adjoint of §`sec:prop-adjoint`, with its emission
    part evaluated in eq. `prop-emission-moments`.
  - J₁₁₁ on the shell has no temporal commutators and no frame term.
- **Suggested concrete targets.**
  - **(a) The shell-averaged J₁₁₁.** The q̂ average of e^{iq·n(o−v)}e^{−iq·n(o−w)} is
    sin(k(w−v))/(k(w−v)), which is elementary, so ⟨J₁₁₁⟩ should reduce to explicit low-dimensional
    integrals of elementary functions, possibly closed form.
  - **(b) The next superhorizon order.** C₂ at O(k₀⁸) needs the six-gradient quadrupole map and the
    square of the four-gradient map. The stationary route plus `quadrupole_adjoint_expansion.py`,
    extended one order, is the most direct path. Validate against `shell_c13.py` at small x_o.
  - **(c) Analytic values to check against.** Any analytic shell result can be tested against the
    numerical table in the appendix: 2.8955, −6.156 and −67.29, and the observer-centred parts.

### Also open

- Write the answer to the k³ question into the report, if desired.
- Regenerate the PDF.
- Run the full `solution/run_all.sh` (long; never run end to end in this session). The new
  `check_constraint_rows.py` (about 2.5 min) is included in it.

## 3. How to reproduce

```bash
uv sync
.venv/bin/python -B solution/check_constraint_rows.py                   # ~2.5 min, exact
.venv/bin/python -B solution/numerics/compare_cubic_quadrupole.py        # ~10 min, numerical validation
.venv/bin/python -B solution/numerics/observer_centred.py                # ~20 min
.venv/bin/python -B solution/numerics/shell_c13.py 1.0 12                # ~25 min per x_o (with split)
.venv/bin/python -B solution/numerics/check_kappa4_shell.py 0.5          # ~15 min
```

The 14 earlier checks are listed in `README.md`; all passed in this session, the slowest being
`check_quadrupole_green_source.py` at about 3 min.
