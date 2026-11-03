# The observed CMB temperature power spectrum through fourth order (radiation-only universe)

An analytical integral representation and evaluated long-wavelength results for the problem in
[`problem_statement.tex`](problem_statement.tex). The
problem is a flat, radiation-only universe in general relativity, with Planck test photons
emitted on a physical temperature surface and an observer on a fluid worldline. It asks for the
observed temperature map through third order, and its angular power spectrum through fourth
order, in the primordial amplitude.

The report is [`report/temperature_spectrum_report.pdf`](report/temperature_spectrum_report.pdf)
(source `.tex`). The general finite-wavelength fourth-order spectrum remains an **unevaluated
exact integral**, which the problem permits; it is not a closed-form evaluated spectrum.
Analytical proofs and exact SymPy checks in `solution/` have different scopes. In particular,
finite Taylor tests are not arbitrary-time substitutions of the complete cubic observable.
Supporting derivations are in `report/` (see the layout table below).
The earlier approaches, intermediate notes and reviews are kept in the git history and are not
part of the current tree.

## Main results

- **Exact redshift.** Use the fluid's velocity potential and comoving labels as coordinates
  (the *temperature frame*). Then the observed temperature is an exact line integral of the
  shear of a conformally rescaled metric, `T_obs = T_fluid(o) e^{-I}`, and the observer's own
  temperature cancels exactly from the normalised map.
- **Field solution.** A triangular solution with elementary Green functions satisfies all ten
  field equations at every order. The second-order kernels are in closed form in entire
  functions. The third-order kernels are explicit nested integrals.
- **Map and spectrum.** The complete map Θ⁽¹⁻³⁾ includes the ray, emission surface, observer
  clock, tetrad and normalisation. In the spectrum, all sky integrals are done in closed form.
- **Linear transfer function.** `T_ℓ` is given in closed form for every ℓ, in the same entire
  functions as the second-order kernels. Its superhorizon series alone determines it (rebuilt
  exactly for ℓ = 2). Inside the horizon it reduces to the Sachs–Wolfe + Doppler form
  `T_ℓ ≈ −[cos x⋆ j_ℓ(y) + √3 sin x⋆ j_ℓ′(y)]`, with y = k₀(t_o − t⋆).
- **Fourth order inside the horizon.** Not evaluated; the report sets out the proposed
  longitudinal-gauge route and its checks.
- **Structure of the fourth-order spectrum.**
  - an exact log-map decomposition: the log-temperature spectrum, plus a lognormal term, plus a
    collapsed bispectrum;
  - an exact observer-dilation decomposition;
  - at zeros of `T_ℓ`, a reduction of the fourth-order spectrum to its quadratic term alone (still an unevaluated integral);
  - on the monochromatic shell, two exact symmetry reductions: the second-order metric is two
    amplitudes with no vector part, and the angle-averaged cubic response is a single scalar;
  - the tensor amplitude equals the known radiation-era induced-GW kernel (Kohri–Terada) plus an
    elementary gauge correction, verified independently.
- **Superhorizon limit.** The quadrupole and the octupole are local curvature multipoles, exact
  in ζ. Their fourth-order factors are `1 + 29ε²σ²/3` and `1 + 46ε²σ²`. For any nonnegative
  spectrum the coefficients lie in `[8, 29/3]·m₀` and `[18, 46]·m₀`, and the monochromatic shell
  attains both upper bounds. The complete fourth-order quadrupole is also given at the next
  gradient order, O(k₀⁶), with symbolic emission, observation and initial times.
- **Fast curvature route.** `solution/curvature_kernels.py` regenerates the **full map** through
  four gradients at general `(k,q,-q)` in about 3 s, versus the legacy 13 h cubic run.
  `--compare --spectrum` checks every component against the raw saved generator output and
  derives the general-moment quadrupole and leading octupole spectra in about 12 s.
  After subtracting observer dilation, the degree-six quadrupole correction is proportional
  to just `m₂m₄(122β_H − 13β_a)`. These are gradient coefficients, not finite-wavelength results.
- **Cubic products.** An explicit orthonormal cubic integrand and combined emission/clock
  identity retain all lower-order products, including the mixed clock–emission term.
  This removes redundant local products, not the remaining propagation integrals.
- **Contracted cubic source reduction (§5.3).** The shell-averaged cubic field sources need
  only nine second-order angular moments and the forced homogeneous lapse. The direct `111`
  longitudinal momentum source vanishes pointwise; the homogeneous tensor cancels after
  averaging. A lapse normal form removes second-order response derivatives. All required
  relative-angle moments are evaluated in the existing entire functions, with at most `J₃`;
  their time integrals remain. The top-order adjoint includes explicit boundary terms.
- **Propagation adjoint (§6.1).** All cubic products linear in second-order fields can be
  written without the ray response driven by `H₂`. Their bulk contains no time derivatives
  of `H₂`; emission, physical clock and random-monopole normalization remain explicit. The
  shell emission contraction uses the same nine moments. Mixed-time phases prevent assuming
  that finite moment closure for the remaining propagation terms. No full cubic cancellation
  is claimed.
- **Open.** What remains unevaluated is listed in the report's section "What remains
  unevaluated". An appendix records exploratory finite-k₀ numerics, which are not evidence.

## Walkthrough: report sections and their checks

| Report section | Checks (in `solution/`) |
|---|---|
| Exact framework (temperature frame, exact redshift, field equations) | `foundations/check_exact_lemmas.py`, `foundations/check_fw_conformal.py`, `foundations/check_exact_field_equations.py`, `foundations/check_field_equations_generic.py`, `check_evolution_form.py` |
| Primordial state and initial data | `foundations/check_initial_data_gradient.py` |
| Field solution at every order | `check_shell_cubic_sources.py --compare-engine`, `check_tensor_newtonian.py`, `check_shell_reduction.py`, `check_green.py`, `second_order.py`, `second_order_consistency.py`, `closed_forms.py`, `closed_second_order.py`, `third_order.py`, `fast_series.py`, `third_order_consistency.py`, `check_report_identities.py` |
| The observable through third order | `check_constraint_rows.py`, `check_propagation_adjoint.py`, `check_observable_expansion.py`, `observable_generator.py`, `observable_fast.py`, `foundations/normalization_expansion.py`, `foundations/check_algebra.py` |
| Gauge invariance and Ward identities | `check_ward_dilation.py`, `check_ward_dilation3.py`, `observable_generator.py W2\|W3`, `foundations/check_first_order.py`, `first_order_engine_check.py`, `check_gradient_modes.py` |
| The spectrum (transfer function, primary form `C_ℓ = 4πvT_ℓ² + v²[Q_ℓ + 4πT_ℓ Re R_ℓ]`, shell two-amplitude reduction, log-map and observer-dilation decompositions) | `transfer_entire.py`, `check_shell_reduction.py`, `check_logmap_decomposition.py`, `check_dilation_factorisation.py`, `check_shell_locality.py` |
| Limits (superhorizon quadrupole and octupole, next gradient order, infrared) | `superhorizon.py`, `observable_generator.py SH`, `check_dilation_factorisation.py`, `check_gradient_octupole.py`, `superhorizon_next.py`, `check_quadrupole_next_hand.py`, `check_cubic_curvature_hand.py` |
| Appendix: exploratory numerics | `numerics/explore_fourth_order.py`, `numerics/explore_third_order.py`, `numerics/explore_adjoint.py` |

The new derivations are included directly in the main report from
`report/cubic_source_reduction.tex`, `report/cubic_adjoint_reduction.tex`, and
`report/propagation_adjoint.tex`. The initial review and fuller source-reduction discussion
are in `report/critical_review_and_cubic_reduction.md`.

## Further cancellation research

The four proposed directions are preserved in
[`report/future_cancellation_directions.md`](report/future_cancellation_directions.md).
Work on the first, a unified constrained adjoint for the whole experiment, is in
[`report/unified_adjoint_working_note.md`](report/unified_adjoint_working_note.md).
The polarized stationary cubic identity replaces explicit second-/third-order state responses
with first state and first adjoint responses. The note also gives the background field Green
transpose and its regular-origin covector, the exact normalization adjoint, and the endpoint,
clock and ray blocks. `solution/check_unified_adjoint.py` checks the identities, the physical
leading-gradient map and a separate algebraic boundary/measurement fixture.

The complete first field-adjoint source is now in
[`report/first_field_adjoint_source.md`](report/first_field_adjoint_source.md), implemented by
`solution/unified_adjoint_source.py`. It combines the off-shell field-operator variation with
the full measurement Hessian, eliminates the independent probe ray, and retains emission,
observer, clock and normalization terms. Its regular Green reconstruction has a generally
nonzero observer trace and an emission jump. The frame's first adjoint and cubic rotation
pairing are kept separately. `solution/check_first_field_adjoint_source.py` checks the field
coefficients against the independent Christoffel engine and the complete source/Green algebra.
The complete finite-wavelength reconstruction and cubic spectrum are **not** yet evaluated
by this route.

The continuation [`report/stationary_cubic_assembly.md`](report/stationary_cubic_assembly.md)
assembles the arbitrary-triple stationary functional, including its direct first-field-only
measurement term, with no second physical fields or second ray responses. It retains the
frame boundary, physical roots, second clock, mixed clock–emission products and exact
normalization. The direct cubic temporal field commutators vanish pointwise on the
monochromatic shell. `solution/stationary_cubic.py` implements the assembly;
`solution/check_stationary_cubic.py` checks arbitrary local jets and independently compares
the full cubic measurement with conserved-photon-momentum endpoint energies in an off-shell
anisotropic metric.

[`report/quadrupole_adjoint_gradient_expansion.md`](report/quadrupole_adjoint_gradient_expansion.md)
now evaluates the leading quadrupole field-adjoint response, including its physical boundary
traces and singular-origin covector. Its integrated first-adjoint contribution at fourth
gradient order is a polynomial: the intermediate logarithms cancel, while the homogeneous
lapse remains nonzero. `solution/quadrupole_adjoint_expansion.py` implements these formulas.
The complete cubic quadrupole through **four spatial gradients** agrees symbolically with
an independent Einstein–radiation curvature calculation at arbitrary input triples, before
Wick contraction; the shell check also reproduces the known weighted C13 coefficient
through sixth gradient order. This is **not** an unrestricted finite-wavelength or
all-multipole result.

[`report/exact_quadrupole_angular_source.md`](report/exact_quadrupole_angular_source.md)
retains the exact finite-wavelength phases in the quadrupole first source. Its measurement
part now has no sky integrals or independently solved first rays: finite combinations of
eight entire radial function families give the mixed-time kernels, including zero phases,
both random-mean normalization products, and all clock/emission/observer coefficients.
The first frame covector reduces exactly to a commutator with the known linear quadrupole.
Rotational covariance also removes the forced homogeneous tensor from the **entire**
isotropically contracted cubic response, while its homogeneous lapse/initial row remains.
`solution/entire_sphere.py` and `solution/quadrupole_angular_source.py` implement the closure.
The remaining generic source/time pairings are **not** evaluated; the follow-up below
now evaluates the internal lapse/spatial and clock Green convolutions.

[`report/direct_cubic_phase_reduction.md`](report/direct_cubic_phase_reduction.md)
assembles the complete direct cubic measurement as exact phase sectors, with physical
endpoints, clock, frame and normalization. Its sole four-time ray term reduces exactly
to six ordered three-time integrals. For one finite-momentum input and two exactly constant
primordial inputs, spatial dilation covariance fixes the **complete physical** cubic response
as `(p.partial_p)^2` of the linear response; the reduced first field adjoint is a derivative
of the background adjoint as well. This is a special channel, not the nonzero-radius shell
contraction. Dense generic standalone radial expansion is intentionally deferred so these
sectors can be combined with adjoint terms first.

[`report/regular_green_time_reduction.md`](report/regular_green_time_reduction.md)
evaluates the coupled lapse-to-spatial Green kernel and the full proper-clock Green
pairing in entire, projector-free form. The exact angular source is inserted into this
transpose, retaining the emission jump, observer trace and origin covector. Background
adjoints now need one time integral, first adjoints at most two, and the complete stationary
cubic quadrupole at most **three independent time integrations per distributed term**.
The bound is not an evaluation of the surviving kernels or spectrum. Implementation is
`solution/regular_field_green.py`; construction is guarded against accidental deep symbolic
integration before physical terms are combined.

**Report integration and constraint rows.** The results of these notes have been reviewed
and are now in the main report as §"Stationary reduction of the complete cubic response"
(`report/stationary_cubic_reduction.tex`). Its final subsection is new: the contracted
Bianchi identities give exact source compatibilities for every physical pair. With them, the
scalar/vector spatial response follows from the momentum constraint,
`K.A = ∫[4K(ν'-ν/s) + 2i c_M]`, using only acoustic and wave kernels and one cumulative integral.
No hyperbolic `-K²/3` Green function is needed. The same replacement applies to
first-adjoint pairings with pair sources, with explicit observer and emission terms. On the
shell the scalar amplitude `S₂` becomes a single acoustic Green integral plus elementary terms
and one `Cin`. This bypasses the nested hyperbolic moment kernel. The check is
`solution/check_constraint_rows.py`. The background adjoint inside the first source is
paired with an incompatible mixed residual. Its hyperbolic part cancels only when all cubic
blocks are combined. Prop. "Combined cubic response" does this: at third order the
same compatibility holds (checked through `t¹⁰` for noncoplanar, zero- and near-zero-pair
triples). The four background blocks therefore collapse to one acoustic/local/wave pairing
with the full cubic source. The complete cubic response becomes
`J₁₁₁ − Σ S_meas[y_bc] + J_lin[y₃]`, with closed-form second-order fields, no first field
adjoint and no hyperbolic kernel. Only `J₁₁₁` still has three-time terms.

**Finite-wavelength validation (numerical, not evidence).**
`solution/numerics/compare_cubic_quadrupole.py` compares this assembled cubic quadrupole with an
independent forward calculation. The forward side solves the fields from the evolution equations
alone and transports the photon and Fermi–Walker triad, with physical roots, clock and exact
normalisation. Agreement is `≤3e-9` for noncoplanar, zero-pair and near-zero-pair triples at
`|k|t̄_o ≈ 3` (about 10 minutes). `observer_centred.py` shows that at finite wavelength the dilation
of the linear and quadratic maps carries most of the cubic response (99% for a generic triple, 67%
for `(k,q,−q)`). `shell_c13.py` evaluates the shell ratio `(1/3)C₂⁽¹³⁾/σ²C₂⁽¹¹⁾` (App. A). The PDF
in the repository has not been regenerated from the updated TeX.

## Setup and reproduction

```
uv sync                     # .venv with sympy, mpmath, numpy, and the dev tools (pyproject.toml, uv.lock)
.venv/bin/python -B solution/check_shell_cubic_sources.py --compare-engine
.venv/bin/python -B solution/check_propagation_adjoint.py
.venv/bin/python -B solution/check_unified_adjoint.py
.venv/bin/python -B solution/check_first_field_adjoint_source.py
.venv/bin/python -B solution/check_stationary_cubic.py
.venv/bin/python -B solution/check_quadrupole_adjoint_expansion.py
.venv/bin/python -B solution/check_stationary_quadrupole_four_gradients.py --symbolic
.venv/bin/python -B solution/check_entire_sphere.py
.venv/bin/python -B solution/check_quadrupole_angular_source.py
.venv/bin/python -B solution/check_phase_integral_contraction.py
.venv/bin/python -B solution/check_direct_cubic_measurement.py
.venv/bin/python -B solution/check_dilation_ward.py
.venv/bin/python -B solution/check_regular_field_green.py
.venv/bin/python -B solution/check_quadrupole_green_source.py
.venv/bin/python -B solution/check_constraint_rows.py
solution/run_all.sh         # the exact checks; uses .venv by default (override with PYTHON=...)
.venv/bin/python -B solution/curvature_kernels.py --spectrum  # cache-free map + spectra, seconds
.venv/bin/python -B solution/curvature_kernels.py --compare   # optional: needs saved legacy kernels
SLOW=1 solution/run_all.sh  # also runs the independent legacy cubic generator (~13 h)
cd report && pdflatex temperature_spectrum_report.tex   # three times, or latexmk -pdf
```

Development: `uv sync` also installs the dev tools (prek, ruff, ty). `prek install` sets up the
pre-commit hook, which runs `ruff check --fix`, `ruff format` and `ty check` from the project venv
(`prek.toml`); run them on demand with `prek run --all-files`. The ty configuration in
`pyproject.toml` relaxes three rules for the symbolic modules only, where SymPy's partly typed
dynamic API causes false positives; the numerics keep the full rule set.

The slow steps of `run_all.sh` are listed in its header. Individual checks run as
`.venv/bin/python solution/<script>.py`. `.devbox/` defines an equivalent container with LaTeX
and an independent CAS (Maxima).

## Layout

| Path | Contents |
|---|---|
| `problem_statement.tex` | the task |
| `report/` | the report (`.tex` and built `.pdf`) and four supporting derivations: `tensor_newtonian_crosscheck.md`, `quadrupole_next_hand_check.md`, `longitudinal_redshift.md`, `curvature_and_adjoint_derivations.md` |
| `solution/` | exact checks for the field solution, observable, spectrum and limits |
| `solution/foundations/` | exact checks of the framework: lemmas, generic-jet field equations, first order in two gauges, normalisation and Gaussian algebra (`gr.py` holds the shared tensor utilities) |
| `solution/results/` | exact expanded sources and kernels referenced by the report |
| `solution/numerics/` | exploratory finite-k₀ numerics (not evidence) |
| `pyproject.toml`, `uv.lock`, `prek.toml` | environment, lint and type-check configuration |

## Exact checks in `solution/`

| `solution/` script | What it establishes |
|---|---|
| `check_evolution_form.py` | mixed trace-free evolution form (the σσ terms cancel), on generic 2-jets |
| `first_order_engine_check.py` | the multilinear engine reproduces the first-order solution for a generic direction |
| `second_order.py` | exact second-order sources and linear operators (written to `results/`) |
| `second_order_consistency.py N` | all ten equations at second order, as exact Taylor series through `O(t^{N−2})`, for general `K` and `K = 0` |
| `check_green.py` | lapse and tensor Green functions, `K → 0` limits, the regular homogeneous mode |
| `check_ward_dilation.py` | dilation Ward identity for second-order fields, in closed form, with negative controls |
| `third_order.py` | exact third-order sources at `(k, q, −q)` (written to `results/`); operators identical to second order |
| `check_observable_expansion.py` | ray, emission root, moving endpoint, clock and composition formulas through third order, for arbitrary fields |
| `check_report_identities.py` | exact residual Bianchi identities and homogeneous constraint counterexample; factored conformal driver and cubic matrix shear; FW rotation; centred/polarised map; two-norm spectrum; general-moment and shell leading-gradient resummation. Fast, read-only, no numerical checks. |
| `check_gradient_modes.py` | exact special-conformal null seed and quadratic completion; comparison of the primordial profile with arXiv:1409.4777v2; nonzero uncompensated-gradient quadrupole; soft generator and leading-gradient Ward check, with a negative control. The arbitrary-time Ward identity is proved by covariance in the report, not by full-kernel substitution. |
| `closed_forms.py`, `closed_second_order.py` | second-order kernels in closed form in entire functions (Si, Cin, 𝔢_m), regular on every resonance; every antiderivative verified by exact differentiation; `--check N` substitutes into the lapse equation (passes); `--check all` (hours) also substitutes into (M) and (T) |
| `fast_series.py N` | third-order implementation test, fast and exact: the full nonlinear equations at `(k, q, −q)`, including the exact `K = 0` pair, vanish through `O(t^{N−2})`; `N = 8` takes about 70 s; negative controls `mutate-c3`, `mutate-K2` |
| `third_order_consistency.py N` | slower SymPy reference implementation of the same test |
| `check_ward_dilation3.py` | third-order field dilation identity `K₃(k,0,0) = (k·∂_k)²K₁` |
| `observable_generator.py W2\|W3\|SH\|SH3` | independent generator of Θ⁽¹⁻³⁾ from the exact definitions; observable Ward identities through `O(κ⁴)`; superhorizon closed form at orders 1–2 (`SH`). `SH3` (third order) did not complete within 50 min: the SymPy assembly of the observable is the bottleneck |
| `observable_fast.py SH\|SH3\|W2\|W3 [KMAX] [mutation]` | fast exact port of the independent generator: same definitions, sparse-polynomial arithmetic over ℚ(i); SH in 11 s, W2/W3 through O(κ⁴), and SH3 (third-order superhorizon map at (k, q, −q), general momenta) in about 12 min; mutations show that the dilation Ward tests are blind to the ray, endpoint and normalisation terms but fail without e^{−ζ} in the initial direction |
| `check_tensor_newtonian.py` | independent cross-check of the shell tensor amplitude against the finite-time Kohri–Terada induced-GW kernel through the nonlinear Newtonian → temperature-frame map (derivation: `report/tensor_newtonian_crosscheck.md`); `check_shell_reduction.py` ties it to this repo's `c_T` |
| `superhorizon_next.py kernels2\|kernels3\|spectrum` | next superhorizon order: map kernels through O(κ⁴) from the fast generator; C₂⁽¹¹⁾, ¼C₂⁽²²⁾ and ⅓Re C₂⁽¹³⁾ at O(k₀⁶), and the complete ε⁴k₀⁶ coefficient. `kernels3` takes about 13 h on a laptop; `spectrum` takes about 25 min the first time and 4 s once cached |
| `check_quadrupole_next_hand.py [--compare]` | independent hand-map calculation of C₂⁽¹¹⁾ and ¼C₂⁽²²⁾ at O(k₀⁶) from its own O(t⁴) field data and the explicit ray, endpoint and clock formulas; endpoint and clock pieces separately nonzero; `--compare` checks every projected-kernel component against the generator output (derivation: `report/quadrupole_next_hand_check.md`) |
| `check_bessel_translation.py` | separated (translation) form of the mixed-time sky integrals: single-time Bessel functions, 3j weights and P_L(μ), checked exactly against the Funk–Hecke form through degree 10 |
| `check_cartesian_green.py` | projector-free Cartesian form of the spatial field equation: eigenvalues of 𝔅_K and exact equality of its polynomial-in-K series with the projector-based kernels through O(t⁷) (pairs, the K = 0 pair, third order); hence the recombined kernels are analytic in K, continuous at K = 0 |
| `check_orthonormal_driver.py` | in the orthonormal photon direction v = e^ζ e^{H/2} e the redshift driver is ½H′ + [H,[H,H′]]/48 + O(ε⁴) for non-commuting H, so the conformal and (H₁²)′ terms drop out; on the shell it is ½H′ through third order |
| `check_orthonormal_map.py` | the map in orthonormal directions: the first-order equation for v is transverse, and the second-order integrand becomes ½H₂′ₙₙ + X₁·∂H₁′ₙₙ − 2H₁′(n, v₁), because the local shear products cancel identically against e − v |
| `check_adjoint_exact.py` | exact ingredients of the adjoint form of the top-order cubic response: lapse Green function and its jump, the factor 14, the source compatibility c_Ham′ − 4c_N/t = 2iK·c_M that removes the contact term, and the absence of d₃, D₃ |
| `check_cubic_curvature_hand.py` | independent covariant-curvature calculation of the complete quadrupole at O(k₀⁶) (C⁽¹¹⁾, ¼C⁽²²⁾ and ⅓Re C⁽¹³⁾); no generator, field engine or cached kernel; verifies the three-tensor form and has negative controls; ~6 s |
| `reconstruct_transfer.py [l N]` | rebuilds T_ℓ from its superhorizon series alone: an ansatz read off the integrand (252 unknowns for ℓ = 2), matched exactly through O(k₀⁶⁴); unique, predicts the next orders exactly, and agrees with the direct integral deep inside the horizon (~7 min, cached) |
| `transfer_subhorizon.py [l order]` | subhorizon limit of T_ℓ (k₀ → ∞ at fixed times), derived exactly from Prop. Tl: T_ℓ = −[cos x⋆ j_ℓ(y) + √3 sin x⋆ j_ℓ′(y)] + O(x_o⁻²), the Sachs–Wolfe + Doppler form; all growing and logarithmic terms cancel; the observer first enters at O(x_o⁻²) |
| `curvature_kernels.py --spectrum [--compare]` | cache-free complete map kernels through gradient degree four for general `(k,q,-q)`; all three general-moment quadrupole coefficients at degree six, separate leading octupole contributions, and the pure `m₂m₄` remainder after observer dilation. Optional exact comparison clears denominators in the **raw** saved generator output, not its simplified cache. Does not compute the six-/five-gradient map needed for the next quadrupole/octupole coefficient |
| `check_cubic_product_identity.py` | generic non-commuting cubic orthonormal driver, ray pullback, and combined emission/clock identity, with omission controls; no numerical checks |
| `check_longitudinal_redshift.py` | exact longitudinal-gauge form of the redshift (shift, tensor shear, Doppler factors, displaced endpoints) and its equality with the temperature-frame exponent (derivation: `report/longitudinal_redshift.md`) |
| `superhorizon.py` | exact superhorizon `C₂⁽¹¹⁾`, `C₂⁽²²⁾`, `C₂⁽¹³⁾` |
| `transfer_entire.py L...` | closed form of the linear transfer function `T_ℓ` for every ℓ in the entire kernels `J_m` (no logs or Ci); coefficients are Laurent coefficients of two rational functions; proved exactly for each given ℓ (x*-derivative identity plus the value at x* = x_o); a printed high-precision comparison is supplementary |
| `check_dilation_factorisation.py` | the leading-gradient map is exactly covariant under ζ → ζ + c (= momentum dilation); the mean of ζ factors out as `⟨e^{−4ζ_o}⟩`; the invariant ratio `1 + (5/3)ε² m₂²/m₄`; Gaussian average of the dilation operator |
| `check_logmap_decomposition.py` | exact check (~2 min) of the log-map decomposition `¼C⁽²²⁾ + ⅓Re C⁽¹³⁾ = C^{L,4} + N_ℓ[C] + M_ℓ` and of the explicit 3j lognormal term, in a band-limited Gaussian sky with non-Gaussian `L₂`, `L₃` |
| `check_shell_reduction.py` | on the monochromatic shell the second-order vector source vanishes identically in t and μ, and reflection leaves one TT polarisation, so A₂ is two amplitudes (scalar and tensor); off-shell negative control |
| `check_shell_locality.py` | on the monochromatic shell the first-order fields and θ are local in g; σ⁽²⁾ = ½A₂′ + ½(H₁²)′ + 2gH₁′, so A₂ is the only non-local second-order input to the map |
| `check_gradient_octupole.py` | leading octupole: conformal identity ∇⟨aR_bc⟩ = (f/2)[∂³f²]_STF, leading T_ℓ coefficient (ℓ = 2..7), β₃, exact Wick moments of the resummed octupole spectrum and the shell ratios 148σ², 27σ²; ⟨θθL₂⟩ = 0 at leading gradient order (so M_ℓ = O(κ⁸)) |

`foundations/` additionally contains the framework checks listed in the walkthrough table.
`foundations/check_algebra.py` is a `unittest` suite: run it with
`python -m unittest discover -s solution/foundations -p check_algebra.py`.

## Exploratory numerics in `solution/numerics/` (not evidence)

The problem forbids numerics in place of derivation. These scripts only look for structure at
finite k₀, to be proved analytically; the report records them in an appendix.

| Script | What it does |
|---|---|
| `explore_fourth_order.py ward\|harden\|scan` | second-order part of the fourth-order spectrum on the monochromatic shell from the exact sources and the report's Green and map formulas; validated by the closed-form T₂, a finite-k dilation Ward identity with negative controls, an independent emission-endpoint test, analytic A₂′, cutoff and grid convergence; exploratory fourth-order C_ℓ at nodes of T_ℓ; amplitude-level recombination diagnostics (`recombination`) and the observer-dilation cross term (`cross_term`) |
| `explore_third_order.py check2\|ward3\|c13` | numerical multilinear engine that solves the map from its definitions through third order at (k, q, −q); reproduces the second-order map of the other script to 1e-6; third-order Ward identity with O(q²) residual and a negative control; superhorizon limit 4; node check; gives ⅓Re C[θ, L₃] and hence the complete fourth-order C₂ away from nodes |
| `explore_adjoint.py kernel\|map` | the adjoint form of the top-order cubic response (Prop. adjoint): a fixed kernel built from the growing- and decaying-mode partial transfer functions, integrated against the shell-averaged sources, with no third-order solve; checked against the direct Green solution and the full numerical map |

Run from `solution/numerics/`, e.g. `../../.venv/bin/python explore_fourth_order.py scan`. The
numeric source caches `solution/results/*_numeric.pkl` are rebuilt automatically from
`second_order.py` and `third_order.py` on first use and are not tracked. Large-x_o runs are
memory-hungry (several GB at x_o = 16).
