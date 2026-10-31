# Exact validation and evidence ledger

## Scope

These tests validate selected algebraic identities used in `../APPROACH.md`.
They are not a third-order Einstein–fluid solver, a complete temperature map,
or a computed fourth-order angular spectrum. No full-physics validation or
independent subagent approval is claimed.

## Reproduce

With the proposed Devbox image, from the repository root:

```sh
python -m unittest discover -s validation -p 'check_*.py' -v
python -m pip check
```

Without rebuilding the image:

```sh
python3 -m venv /tmp/david-analysis-venv
/tmp/david-analysis-venv/bin/pip install -r .devbox/requirements.txt
/tmp/david-analysis-venv/bin/python -m unittest discover -s validation -p 'check_*.py' -v
/tmp/david-analysis-venv/bin/python -m pip check
```

Dependencies are SymPy 1.14.0 and its pinned dependency mpmath 1.3.0.
The tests use exact symbolic expressions only; they do not invoke mpmath,
`evalf`, numerical integration/evolution, floating-point comparison, or sampling.

## Results in this session

Python 3.12.3, SymPy 1.14.0: **11 tests passed**, exit code 0.
`pip check`: **No broken requirements found**, exit code 0.

| Test | What it establishes |
|---|---|
| `test_radiation_background` | Friedmann, continuity, temperature scaling, and background proper-time differential |
| `test_sky_division_and_zero_mean` | All three sky-normalization coefficients and their zero sky means |
| `test_log_to_temperature` | All cubic log-redshift-to-temperature product coefficients |
| `test_common_temperature_factor_cancels` | A common sky-independent factor cancels from normalized temperature through cubic |
| `test_cubic_emission_root` | The three implicit-function root coefficients annihilate the composed surface equation |
| `test_taylor_power_factors` | Taylor factors in `11`, `22`, and `13`, and the odd-Gaussian-degree bookkeeping |
| `test_wick_pairings` | One linear contraction, two connected quadratic contractions, and all three linear–cubic contractions |
| `test_scalar_fundamental_solutions` | Radiation scalar homogeneous equation, regular limit, and Wronskian |
| `test_scalar_retarded_kernel` | Scalar Green function homogeneous equation, coincidence value, unit derivative jump, and zero-momentum limit |
| `test_tensor_fundamental_solutions` | Tensor homogeneous solutions for nonzero and zero momentum |
| `test_shell_measure_and_linear_angular_factor` | Delta-shell Jacobian, fixed-variance measure, and linear angular prefactor |

The harmonic addition theorem is used as a standard mathematical identity in
the angular-factor test; the script checks its consequences for the chosen
Fourier conventions rather than proving the theorem. The Green function tests
check the elementary retarded kernel, not convergence of integrals against
unknown nonlinear sources. The root test assumes a nonzero background surface
derivative. These scope restrictions are intentional.

## Independent review: BLOCKED, not passed

Two read-only review lanes were requested: (1) nonlinear initial data,
constraints, physical observer and gauge definition; (2) normalization,
Gaussian contractions, integral legality and validation criteria.

- Workflow `c52527ab-ba2f-4d03-8d66-584d0ba4c22c`: both reviewers failed
  before producing a review with `No API key found for openai.`
- One retry through the same subagent interface, explicitly selecting the
  registered `ninfer/qwen3.8-27b` model:
  workflow `f1ba8cc6-15fb-4abd-a9b0-fae3743cb745`.
  Both reviewers failed with `Request timed out.` No review was produced.

No conclusions are attributed to those reviewers. No external CLI fallback
was launched. Resume independent review after repairing native provider access,
or after the user approves a different reviewer execution mechanism. A future
review should read both the original statement and the proposed approach, and
independently reproduce the algebra rather than treating passing tests as proof
of the unimplemented physics.

## Container validation

`.devbox/Containerfile` installs Python/venv and the same pinned requirements
inside `/opt/analysis-venv`, placed on `PATH`. The manifest's build context is
`.devbox`, so `COPY requirements.txt` refers to `.devbox/requirements.txt`.
The package install and checks were successfully exercised in the local venv.
**The container image itself has not been built:** neither Docker nor Podman
is available in this container. Host-side image build/restart is a remaining
infrastructure check, not evidence for the analytical calculation.

## Completion criteria

The full analytical validation matrix and unresolved derivations are in
`../APPROACH.md`, sections 8–10. In particular the nonlinear constraints,
third-order field/observer/ray residuals, gauge invariance of explicit kernels,
and convergence proofs remain unimplemented and unvalidated.
