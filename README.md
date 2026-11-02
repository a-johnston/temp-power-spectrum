# The observed CMB temperature power spectrum through fourth order (radiation-only universe)

Analytical solution of the problem in `problem_statement.tex`: the observed temperature anisotropy
through third order and its angular power spectrum through fourth order in the primordial
amplitude, for a flat radiation-only universe in general relativity, with Planck test photons
emitted on a physical temperature surface and an observer on a fluid worldline.

## Start here

- **`report/temperature_spectrum_report.tex`** (and `.pdf`): the final report. It contains every
  result, proof and evaluated limit, the verification ledger, what remains unevaluated, and an
  appendix of exploratory (not-evidence) numerics. Build with `pdflatex` run three times in
  `report/` (or `latexmk -pdf`).
- `problem_statement.tex`: the task.

Main results, in report order: the exact temperature-frame redshift formula; the triangular
field solution with all ten equations proved; second-order kernels in closed form and
third-order kernels as nested integrals; the complete map Θ⁽¹⁻³⁾; the entire-function closed form
of the linear transfer function T_ℓ for every ℓ; the log-map and observer-dilation
decompositions of the fourth-order spectrum; complete fourth-order spectra at zeros of T_ℓ; and
the superhorizon quadrupole and octupole exact in ζ (fourth-order factors 1 + 29ε²σ²/3 and
1 + 46ε²σ²).

## Setup

```
uv sync            # creates .venv with sympy, mpmath, numpy (pyproject.toml, uv.lock)
```

## Reproduce

- Exact checks (SymPy, no numerical evidence): `solution/run_all.sh` uses `.venv` by default
  (override with `PYTHON=...`). The slow steps are listed in its header. Individual checks run
  as `.venv/bin/python solution/<script>.py`; the table below says what each establishes.
- Older layers, kept for provenance: `verification/run_all.sh` and
  `python -B -m unittest discover -s validation -p 'check_*.py'`.
- Exploratory numerics (report appendix; **not evidence**): see below.

## Layout

| Path | Contents |
|---|---|
| `report/` | the final report |
| `solution/` | the exact SymPy checks behind the report (table below); `solution/results/` holds the exact expanded sources and kernels referenced by the report |
| `solution/numerics/` | exploratory finite-k₀ numerics (below) |
| `verification/` | exact lemmas, generic-jet field equations, first order in two gauges, the earlier Ci/Si transfer closed form |
| `reviews/` | independent reviews and review notes (the report incorporates their accepted points) |
| `arXiv-1409.4777v2/` | Mirbabayi & Zaldarriaga, *CMB Anisotropies from a Gradient Mode*, used for comparison |
| `PLAN.md`, `derivations/UNIFIED_DERIVATION.md` | the unified plan and its execution notes; the report supersedes them where they differ |
| `APPROACH.md`, `approach/`, `derivations/{QUADRATIC,CUBIC}_INITIAL_DATA.md`, `derivations/OFFSHELL_VALIDATION.md`, `validation/` | earlier (orange and purple) branches, kept for provenance |

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
| `superhorizon.py` | exact superhorizon `C₂⁽¹¹⁾`, `C₂⁽²²⁾`, `C₂⁽¹³⁾` |
| `transfer_entire.py L...` | closed form of the linear transfer function `T_ℓ` for every ℓ in the entire kernels `J_m` (no logs or Ci); coefficients are Laurent coefficients of two rational functions; proved exactly for each given ℓ (x*-derivative identity plus the value at x* = x_o); a printed high-precision comparison is supplementary |
| `check_dilation_factorisation.py` | the leading-gradient map is exactly covariant under ζ → ζ + c (= momentum dilation); the mean of ζ factors out as `⟨e^{−4ζ_o}⟩`; the invariant ratio `1 + (5/3)ε² m₂²/m₄`; Gaussian average of the dilation operator |
| `check_logmap_decomposition.py` | exact check (~2 min) of the log-map decomposition `¼C⁽²²⁾ + ⅓Re C⁽¹³⁾ = C^{L,4} + N_ℓ[C] + M_ℓ` and of the explicit 3j lognormal term, in a band-limited Gaussian sky with non-Gaussian `L₂`, `L₃` |
| `check_shell_locality.py` | on the monochromatic shell the first-order fields and θ are local in g; σ⁽²⁾ = ½A₂′ + ½(H₁²)′ + 2gH₁′, so A₂ is the only non-local second-order input to the map |
| `check_gradient_octupole.py` | leading octupole: conformal identity ∇⟨aR_bc⟩ = (f/2)[∂³f²]_STF, leading T_ℓ coefficient (ℓ = 2..7), β₃, exact Wick moments of the resummed octupole spectrum and the shell ratios 148σ², 27σ²; ⟨θθL₂⟩ = 0 at leading gradient order (so M_ℓ = O(κ⁸)) |

## Exploratory numerics in `solution/numerics/` (not evidence)

The problem forbids numerics in place of derivation. These scripts only look for structure at
finite k₀, to be proved analytically; the report records them in an appendix.

| Script | What it does |
|---|---|
| `explore_fourth_order.py ward\|harden\|scan` | second-order part of the fourth-order spectrum on the monochromatic shell from the exact sources and the report's Green and map formulas; validated by the closed-form T₂, a finite-k dilation Ward identity with negative controls, an independent emission-endpoint test, analytic A₂′, cutoff and grid convergence; complete fourth-order C_ℓ at nodes of T_ℓ; amplitude-level recombination diagnostics (`recombination`) and the observer-dilation cross term (`cross_term`) |
| `explore_third_order.py check2\|ward3\|c13` | numerical multilinear engine that solves the map from its definitions through third order at (k, q, −q); reproduces the second-order map of the other script to 1e-6; third-order Ward identity with O(q²) residual and a negative control; superhorizon limit 4; node check; gives ⅓Re C[θ, L₃] and hence the complete fourth-order C₂ away from nodes |

Run from `solution/numerics/`, e.g. `../../.venv/bin/python explore_fourth_order.py scan`. The
numeric source caches `solution/results/*_numeric.pkl` are rebuilt automatically from
`second_order.py` and `third_order.py` on first use and are not tracked. Large-x_o runs are
memory-hungry (several GB at x_o = 16).
