# Analytical CMB temperature power spectrum through fourth order (radiation-only universe)

- `problem_statement.tex`: the task.
- **`PLAN.md`**: the unified plan, merging the orange and purple branches. It is authoritative
  where they differ.
- **`derivations/UNIFIED_DERIVATION.md`**: the execution of the plan, through third order in the
  map and fourth order in `C_ℓ`, with a status ledger. Each result is marked as evaluated, an
  exact integral representation, or open.
- **`report/temperature_spectrum_report.tex`**: the final report for the problem. It covers the
  results, the proofs, the evaluated limits, the verification ledger, and what remains
  unevaluated. Build it with `latexmk -pdf` in the devbox image.
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
| `closed_forms.py`, `closed_second_order.py` | second-order kernels in closed form in entire functions (Si, Cin, 𝔢_m), regular on every resonance; every antiderivative verified by exact differentiation; `--check N` substitutes into the lapse equation (passes); `--check all` (hours) also substitutes into (M) and (T) |
| `fast_series.py N` | third-order implementation test, fast and exact: the full nonlinear equations at `(k, q, −q)`, including the exact `K = 0` pair, vanish through `O(t^{N−2})`; `N = 8` takes about 70 s; negative controls `mutate-c3`, `mutate-K2` |
| `third_order_consistency.py N` | slower SymPy reference implementation of the same test |
| `check_ward_dilation3.py` | third-order field dilation identity `K₃(k,0,0) = (k·∂_k)²K₁` |
| `observable_generator.py W2\|W3\|SH\|SH3` | independent generator of Θ⁽¹⁻³⁾ from the exact definitions; observable Ward identities through `O(κ⁴)`; superhorizon closed form at orders 1–2 (`SH`). `SH3` (third order) did not complete within 50 min: the SymPy assembly of the observable is the bottleneck |
| `superhorizon.py` | exact superhorizon `C₂⁽¹¹⁾`, `C₂⁽²²⁾`, `C₂⁽¹³⁾` |
