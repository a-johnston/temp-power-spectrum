# Analytical CMB temperature power spectrum through fourth order (radiation-only universe)

- `problem_statement.tex`: the task.
- **`PLAN.md`**: the unified plan, merging the orange and purple branches. It is authoritative
  where they differ.
- **`derivations/UNIFIED_DERIVATION.md`**: the execution of the plan, through third order in the
  map and fourth order in `C_ℓ`, with a status ledger. Each result is marked as evaluated, an
  exact integral representation, or open.
- **`report/temperature_spectrum_report.tex`**: the final report for the problem. It covers the
  results, the proofs, the evaluated limits, the verification ledger, and what remains
  unevaluated. The strengthened report includes a constraint-proof repair, an explicit
  Fermi–Walker map correction, compact finite-time source/map/spectrum representations, and
  an all-amplitude Gaussian average of the leading-gradient quadrupole (not a finite-wavelength
  resummation). A comparison with Mirbabayi & Zaldarriaga, *CMB Anisotropies from a Gradient
  Mode* (arXiv:1409.4777v2), adds an exact adiabatic-gradient null test and a finite-time
  soft-gradient identity. It distinguishes a nonlinear coordinate mode from the physical
  gradient-squared signal for the stated Gaussian seed. Matter-era transfer coefficients
  are not used. The linear transfer function is now given for every ℓ in the same entire
  kernels as the second-order fields, and the superhorizon spectrum is factorised into an
  observer dilation (convention dependent) and an invariant curvature norm. The fourth-order
  spectrum is split exactly via the centred log-temperature map into its own spectrum, a
  closed-form lognormal 3j term and a collapsed bispectrum; on the monochromatic shell the
  non-local (resonant) content is confined to A₂, A₃ and ν₂. The large expanded
  kernels remain in the referenced ancillary files.
  The report's repaired constraint proof and explicit frame completion supersede the
  corresponding older derivation notes.
  Build with `latexmk -pdf` in the devbox image, or run `pdflatex` twice in `report/`.
- `solution/`: exact SymPy checks for the unified derivation. Run them with `solution/run_all.sh`.
  The exact sources are in `solution/results/`.

Earlier material, kept for provenance (the unified documents supersede it):

- orange: `APPROACH.md`, `derivations/{QUADRATIC,CUBIC}_INITIAL_DATA.md`, `OFFSHELL_VALIDATION.md`,
  `reviews/`, and `validation/` (run with `python -B -m unittest discover -s validation -p 'check_*.py'`);
- purple: `approach/approach.tex`, `approach/VALIDATION.md`, and `verification/` (run with
  `verification/run_all.sh`).

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
| `check_logmap_decomposition.py` | exact check (~6 min) of the log-map decomposition `¼C⁽²²⁾ + ⅓Re C⁽¹³⁾ = C^{L,4} + N_ℓ[C] + M_ℓ` and of the closed-form 3j lognormal term, in a band-limited Gaussian sky with non-Gaussian `L₂`, `L₃` |
| `check_shell_locality.py` | on the monochromatic shell the first-order fields and θ are local in g; σ⁽²⁾ = ½A₂′ + ½(H₁²)′ + 2gH₁′, so A₂ is the only non-local second-order input to the map |
