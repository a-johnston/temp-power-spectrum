# User-supplied independent review 02: quadratic milestone accepted

## Provenance

The user supplied a full independent review targeting commit `ea6b053` and the
seven requested source/derivation/evidence documents. This file is a **summary
and disposition**, not a verbatim transcript or a native subagent artifact.

The reviewer reported running all 30 tests successfully, performing additional
exact symbolic calculations, and leaving the repository unchanged. Those
execution claims are attributed to the reviewer. No numerical evolution,
quadrature, fitting or momentum sampling was reported.

## Verdict and boundaries

> **Quadratic milestone established—within the explicitly declared,
> fixed-momentum perturbative class.**

No blockers or major defects were found for that milestone. The acceptance
covers the fixed material chart, background normalization, conformally flat
primordial seed and regular analytic class. It is local in the smooth coordinate
domain and modewise at finite input momenta.

It does not establish global nonlinear existence for arbitrary Gaussian
realizations, unrestricted momentum-integral convergence, cubic dynamics,
observable kernels, or the temperature spectrum.

## Independently checked claims reported by the reviewer

- Exact fluid integration and the treatment of its spatial and temporal
  integration functions; agreement with the primordial uniform-density limit.
- The exact mixed ADM reduction and constraints, including lapse-gradient and
  matrix-ordering signs, without using gravitational constraints in the
  evolution reduction.
- The nonlinear leading boundary coefficient and its covariant constraints.
- Full finite-time linear fields, including the anisotropic metric equation.
- Quadratic determinant/connection/Ricci/source factorials and Fourier
  polarization, using fully symbolic symmetric matrices and independent first
  time derivatives rather than fixed matrix probes.
- General-momentum cancellation of the primordial singularity, seed subtraction,
  and the mixed constant-dilation solution `J(η;0,p)=(2−p∂p)H(η;p)`.
- Retarded-series order, operator norms, finite-time convergence, differentiation,
  regular endpoint and uniqueness in the specified class.
- Actual values at zero total momentum, retaining forced homogeneous trace and
  trace-free sectors.
- Offshell ADM/Bianchi constraint propagation and its energy argument; finite-time
  Cauchy data and matching-time independence.

## Minor recommendation and response

The reviewer noted that the direct-source regression starts from the proposed
exact metric equation and cannot validate that equation against ADM by itself.
Likewise, the energy regression assumes the constraint-propagation equations.
The recommendation was to preserve the upstream offshell evidence reproducibly.

Response:

1. Added `derivations/OFFSHELL_VALIDATION.md`, including the important offshell
   relation `E_ij=γ_ij H_c/2`, rather than incorrectly imposing `E_ij=0`.
2. Added `validation/check_offshell.py`, checking the mixed ADM-to-metric algebra,
   trace-reversed residual relation, and explicit coordinate Bianchi calculation
   with arbitrary lapse derivatives, extrinsic curvature and residuals.
3. Added a fully symbolic quadratic matrix/source comparison and a finite-time
   constant-dilation regression to `check_initial_data.py`.

These retained regressions run locally. They do not waive the review's scope
qualifications.

## Next gate

The reviewer explicitly found no material correction necessary before cubic
extension. The new `derivations/CUBIC_INITIAL_DATA.md` proceeds on that basis.
**This review does not approve the new cubic derivation.** Its source,
polarization, lower-order integral nesting, boundary cancellation, and constraint
induction require a separate independent review before that milestone is accepted.
