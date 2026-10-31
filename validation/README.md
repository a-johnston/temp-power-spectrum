# Exact validation and evidence ledger

## Scope

The repository now contains three analytical layers:

1. `../APPROACH.md`: the full-spectrum plan and verified composition/statistical algebra.
2. `../derivations/QUADRATIC_INITIAL_DATA.md`: independently accepted quadratic
   data within the stated local, fixed-momentum regular perturbative class.
3. `../derivations/CUBIC_INITIAL_DATA.md`: a new explicit cubic source and nested
   retarded construction, **not yet independently reviewed**.

`../derivations/OFFSHELL_VALIDATION.md` preserves the upstream ADM/Bianchi
argument recommended by the quadratic reviewer. None of these layers supplies
the complete observed temperature map or requested spectrum, or proves
unrestricted spectral convergence. No full-spectrum validation is claimed.

## Reproduce

With the proposed Devbox image, from the repository root:

```sh
python -B -m unittest discover -s validation -p 'check_*.py' -v
python -m pip check
```

Without rebuilding the image:

```sh
python3 -m venv /tmp/david-analysis-venv
/tmp/david-analysis-venv/bin/pip install -r .devbox/requirements.txt
/tmp/david-analysis-venv/bin/python -B -m unittest discover -s validation -p 'check_*.py' -v
/tmp/david-analysis-venv/bin/python -m pip check
```

Dependencies are SymPy 1.14.0 and its pinned dependency mpmath 1.3.0.
The tests use exact symbolic expressions only; they do not invoke mpmath,
`evalf`, numerical integration/evolution, floating-point comparison, or sampling.
Symbolic definite integration of polynomial Green-function terms is analytical,
not numerical quadrature.

## Current local results

Python 3.12.3, SymPy 1.14.0: **40 tests passed**, exit code 0.
The suite contains 16 composition/statistics tests, 16 quadratic initial-data
tests, 3 offshell checks, and 5 new cubic checks. The full suite currently takes
about 90–100 seconds in this environment; the cubic expressions dominate runtime.

### Composition, statistics and elementary functions: 16 tests

| Test | What it establishes |
|---|---|
| `test_radiation_background` | Friedmann, continuity, temperature scaling, and background proper-time differential |
| `test_sky_division_and_zero_mean` | Three sky-normalization coefficients and their zero sky means |
| `test_log_to_temperature` | Cubic log-redshift-to-temperature product coefficients |
| `test_common_temperature_factor_cancels` | A genuinely common sky-independent factor cancels through cubic |
| `test_cubic_emission_root` | Implicit root coefficients annihilate the composed surface equation |
| `test_taylor_power_factors` | Taylor factors in `11`, `22`, `13`, and odd-Gaussian-degree bookkeeping |
| `test_wick_pairings` | One linear contraction, two connected quadratic contractions, three linear–cubic contractions |
| `test_scalar_fundamental_solutions` | Radiation scalar homogeneous equation, regular limit and Wronskian |
| `test_scalar_retarded_kernel` | Scalar Green equation, coincidence value, unit derivative jump and zero-momentum limit |
| `test_tensor_fundamental_solutions` | Tensor homogeneous solutions at nonzero and zero momentum |
| `test_shell_measure_and_linear_angular_factor` | Delta-shell Jacobian, fixed-variance measure and linear angular prefactor |
| `test_endpoint_energy_inverse_and_log` | Endpoint energy to log-temperature and temperature, with the correct inverse-redshift signs |
| `test_actual_angular_normalization_and_dipole_projection` | Exact axisymmetric angular averaging, surviving nonlinear multipoles, and final dipole projection |
| `test_tensor_retarded_green` | Explicit tensor Green equation and both retarded boundary conditions |
| `test_shifted_scalar_endpoint_cubic` | Cubic evaluation of an arbitrary scalar on a shifted event |
| `test_gaussian_moments_without_pair_enumeration` | Independent Gaussian polynomial-moment checks of `22` and the essential cubic traces in `13` |

The angular-factor test uses the harmonic addition theorem rather than proving
it. The Green-function tests do not establish convergence against arbitrary
nonlinear sources. The surface-root test assumes a transverse intersection.

### Quadratic initial data: 16 tests

`initial_data.py` implements the component source and elementary matrix
operators. `check_initial_data.py` supplies these checks:

| Test | What it establishes |
|---|---|
| `test_exact_fluid_reduction` | Integrated continuity and Euler identities in the chosen gauge |
| `test_seed_ricci_from_connection` | General conformal seed curvature derived from its Christoffel symbols |
| `test_leading_constraints_and_evolution` | Leading nonlinear time coefficient solves evolution, Hamiltonian and momentum conditions, using contracted Bianchi |
| `test_linear_full_equations` | All linear metric components, lapse definition, Hamiltonian and momentum constraints for arbitrary magnitude; rotation covariance supplies general direction |
| `test_linear_zero_mode_and_boundary` | Removable zero-input-momentum limit and leading time series |
| `test_quadratic_source_against_direct_metric_expansion` | A separately coded nilpotent-amplitude expansion of the exact metric RHS agrees with the bilinear source on noncommuting time-dependent symmetric probes with general symbolic momenta |
| `test_regular_source_boundary_general_momenta` | Cancellation of the apparent time singularity and the correct leading coefficient for arbitrary vector inputs |
| `test_source_symmetry` | Symmetry of spatial indices and input momenta on the boundary jets |
| `test_zero_total_momentum_forced_tracefree_response` | Nonzero forced homogeneous trace/trace-free response, and vanishing spatial operator at exactly K=0 |
| `test_lapse_density_boundary_kernel` | Quadratic lapse/density coefficients, including the antipodal input pair |
| `test_matrix_green_boundary_and_equation` | Both trace-sector Green equations and retarded boundary conditions for an arbitrary symmetric matrix |
| `test_indicial_inverse_and_exact_monomial_integral` | Indicial inverse and exact time integration of constant, quadratic and quartic sources |
| `test_quadratic_constraints_through_next_time_coefficient` | Independently assembled Hamiltonian residual through time power η² and all momentum residuals through η³, using the metric through η⁴ and general planar momentum pair |
| `test_constraint_propagation_energy_identity` | Hermitian constraint-energy identity used in the analytic all-time uniqueness argument |
| `test_quadratic_source_fully_symbolic_matrices` | Source comparison with independent symbolic symmetric matrix functions and first derivatives, with all six momentum components free |
| `test_constant_dilation_identity` | Exact finite-time mixed kernel with one zero-input-momentum dilation |

The original direct source test uses explicit noncommuting probes. The added
fully symbolic test removes that probe restriction, but both start from the
exact metric equation. Its upstream ADM reduction is checked separately below. The quadratic constraint test
uses `(a,0,0)` and `(b,c,0)` with symbolic independent `a,b,c`, which covers general
pairs modulo rotations; no collinear or equal-magnitude assumption is made.
It verifies a finite number of time coefficients, **not** all-time constraints
by brute-force substitution. The all-time constraint argument is the separate
ADM/Bianchi proof in the derivation. The second user-supplied review accepts
that proof and the exact reduction; reproducible upstream checks are now retained.

The retarded series bound proves finite-time convergence at fixed input
momenta. No numerical evaluation of that series is used. No test here certifies
unbounded spectral integration, the eventual temperature-kernel contraction
restriction, or gauge invariance of a temperature kernel not yet constructed.

### Offshell reduction: 3 tests

`check_offshell.py` checks the mixed ADM-to-metric matrix identity, the offshell
trace-reversed spatial relation `E_ij=γ_ij H_c/2`, and all four coordinate Bianchi
divergences at a spatial normal point with arbitrary lapse derivatives,
extrinsic curvature and constraint residuals. It does not assume the constraints
vanish. These checks are upstream of the constraint-energy identity.

### Cubic dynamics: 5 tests

| Test | What it establishes |
|---|---|
| `test_cubic_determinant_and_fluid_composition` | Generic symbolic matrix log-determinant coefficient and density/lapse composition |
| `test_cubic_source_against_exact_metric_rhs` | The explicit source matches a separately organized full inverse-metric/connection/Ricci expansion on noncommuting time profiles with all nine momentum components symbolic |
| `test_cubic_regular_boundary_general_momenta` | Singular and constant terms cancel after seed subtraction for arbitrary input vectors |
| `test_cubic_lapse_boundary_and_internal_zero_pair` | Leading cubic lapse coefficient, including the exact `(k,q,−q)` restriction |
| `test_two_constant_dilations_exact_finite_time` | Full finite-time cubic source and lapse identity with two zero-input-momentum dilation legs |

The cubic comparison uses the same finite amplitude ring and elementary
connection/Ricci contraction helpers in its two constructions, but not the same
cubic source formula. It is a noncommuting-probe regression, not a fully symbolic
arbitrary-matrix cubic proof or an independent tensor implementation. The new
cubic derivation's weighted time-nesting bounds and constraint induction require
independent review. Boundary jets are used only in explicitly labelled tests,
not as finite-time field solutions. No cubic all-time temperature observable is
implemented or validated.

## Independent review history

### Native execution failures — historical record

Two read-only lanes were requested for the original proposal:

- Workflow `c52527ab-ba2f-4d03-8d66-584d0ba4c22c`: both failed before review with
  `No API key found for openai.`
- One retry selecting registered model `ninfer/qwen3.8-27b`, workflow
  `f1ba8cc6-15fb-4abd-a9b0-fae3743cb745`: both failed with `Request timed out.`

No review conclusions are attributed to these runs. No external CLI fallback
was launched.

### User-supplied independent review — received

The user subsequently supplied a complete independent review of the **original
plan and 11-test suite**. Its verdict was **“Sound approach with stated gaps.”**
It reported no principal algebraic discrepancy and identified initial-data
construction, explicit kernels, exact restrictions, and gauge verification as
blocking obligations. Attribution and disposition are recorded in
`../reviews/INDEPENDENT_REVIEW_01.md`.

The first review recommended the general-momentum quadratic boundary problem.
The user then supplied a second independent review of `ea6b053`, with verdict
**“Quadratic milestone established—within the explicitly declared,
fixed-momentum perturbative class.”** See `../reviews/INDEPENDENT_REVIEW_02.md`.
It found no blocker or major defect within that scope and recommended retaining
the stronger offshell and fully symbolic regressions now added here.

The new cubic derivation is not covered by either favorable review; its
independent acceptance is pending.

## Container validation

`.devbox/Containerfile` installs Python/venv and the pinned requirements inside
`/opt/analysis-venv`, placed on `PATH`. The build context is `.devbox`, so
`COPY requirements.txt` refers to `.devbox/requirements.txt`.
The dependency install and checks were exercised in the local venv. No further
dependencies were needed for this milestone.
**The container image itself has not been built in this environment.**

## Remaining acceptance gates

- Independent review of the new cubic source, polarization, nested retarded
  representation, regular endpoint, and constraint induction.
- Observer/clock/tetrad/ray/endpoint construction and residuals through cubic.
- Component-level observable gauge verification at fixed physical labels.
- Final harmonic kernels and analytical spectral convergence/restriction proofs.

Passing the current tests does not waive any of these gates.
