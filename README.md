# Analytical CMB temperature power spectrum through fourth order (radiation-only universe)

- `problem_statement.tex`: the task.
- `approach/approach.tex`: the proposed approach. It covers the exact "temperature-frame" reformulation, the primordial state and initial data, the construction of the observable, the second- and third-order plan, the validation matrix, open items and risks. Build it with `latexmk -pdf` or `tectonic`.
- `approach/VALIDATION.md`: findings from three independent subagents (a blind re-derivation, a physics referee and a mutation-testing code audit), and how each was resolved.
- `verification/`: exact SymPy checks. Run them all with `verification/run_all.sh` (about 10 minutes; the closed-form step is the slow one).

| Script | What it establishes |
|---|---|
| `check_exact_lemmas.py` | L1 (continuity ⇔ expansion-free; Euler is automatic), L2 (exact ray and redshift equations), L3 (conformal map, T_obs = T_fluid(o) Ẽ_o/Ẽ_e), L4 (Hamiltonian constraint) |
| `check_fw_conformal.py` | L5: Fermi–Walker transport in g ⇔ parallel transport in the temperature frame |
| `check_field_equations_generic.py` | exact field equations for generic γ with t-independent det (proof via symbolic 2-jets) |
| `check_exact_field_equations.py` | the same identities, fully symbolic, on a restricted det-1 family |
| `check_first_order.py` | first-order solution in two gauges (all 10 Einstein equations); exact equivalence of the transfer functions for all l |
| `check_initial_data_gradient.py` | nonlinear initial data at O(t²), exact in ζ, with the momentum constraint satisfied |
| `normalization_expansion.py` | Θ⁽¹⁾, Θ⁽²⁾, Θ⁽³⁾ in terms of the line-of-sight exponent, checked on exact discrete skies with negative controls |
| `transfer_closed_form.py L...` | closed-form T_l in Si/Ci, verified by exact differentiation; output in `results/` |

**Status:** Stages 0–1 are done and verified: the exact framework, the full first order and C_l⁽¹¹⁾. Stages 2–5 are specified but not yet executed: the second- and third-order kernels, Θ⁽²⁾ and Θ⁽³⁾, C_l⁽²²⁾ and C_l⁽¹³⁾, and the monochromatic specialisation. See §8–9 of `approach.tex`.
