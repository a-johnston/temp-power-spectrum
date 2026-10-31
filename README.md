# Analytical CMB temperature power spectrum through fourth order (radiation-only universe)

- `problem_statement.tex`: the task.
- **`PLAN.md`**: the unified plan, merging the orange and purple branches. It is authoritative
  where they differ.
- **`derivations/UNIFIED_DERIVATION.md`**: the execution of the plan, through third order in the
  map and fourth order in `C_ℓ`, with a status ledger. Each result is marked as evaluated, an
  exact integral representation, or open.
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
