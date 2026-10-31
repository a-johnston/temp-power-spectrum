# User-supplied independent review: disposition

## Provenance and scope

The user supplied a full independent review in the conversation after the
initial proposal (`9b61870`). This file is a **summary and response ledger**, not
a verbatim transcript and not output from our failed native subagent runs.
The reviewer reported reading the four original project documents, passing all
11 original symbolic tests, and independently checking Gaussian moments, an
actual angular average, and a tensor Green function. Those execution claims
are attributed to the reviewer. The additional identities have now also been
implemented and run locally as exact symbolic regressions.

**Reviewer verdict: “Sound approach with stated gaps.”** No discrepancy was
found in the displayed normalization, logarithmic-temperature, root-shift,
Wick, Taylor, Fourier/angular, shell-measure, background or elementary
homogeneous-solution formulas. This is not certification of the missing
nonlinear kernels or the requested spectrum.

## Findings and response

| Finding | Accepted response | Current status |
|---|---|---|
| Blocker 1: primordial asymptotic class, labels, constraints and uniqueness not constructed | Fix a comoving-orthogonal gauge; distinguish the conformally flat physical seed restriction from gauge fixing; keep forced homogeneous trace/trace-free response | New quadratic construction in `derivations/QUADRATIC_INITIAL_DATA.md`; awaiting independent review; cubic still absent |
| Blocker 2: abstract operator/source recursions are not explicit analytical integrands | Give exact component equations, bilinear source, regular seed subtraction and completely specified retarded integrands | Supplied for the quadratic metric; not for cubic dynamics or the observable |
| Blocker 3: smooth-spectrum integrability does not prove contraction or shell restriction | Require values on exact zero-partial-sum surfaces as well as integrability near them; prove shell restriction separately | Explicit quadratic K=0 sector now available; full cubic/observable restriction and spectral proofs unresolved |
| Major 1: exact covariance is not verification of cubic kernels | Retain fixed physical/primordial labels in a component-level transformation or second-gauge derivation | Unresolved for the observable |
| Major 2: observer/ray prescriptions are not their cubic realization | Derive every clock, trajectory, tetrad, ray and endpoint coefficient before normalization | Unresolved; endpoint composition regressions added |
| Major 3: algebra tests are not physical closure | Add tests of independently expanded equations and retain exact scope statements | Expanded tests and proof ledger; no full-spectrum validation claim |
| Minor: freeze the experiment and formal remainder meaning | Require observation after emission after clock start; keep the Gaussian spectrum remainder explicitly formal | Updated `APPROACH.md` |

## Additional identities adopted

For endpoint-evaluated emission energy
`E_e/E_e0=1+εe1+ε²e2/2+ε³e3/6`, with `E_o=1`:

- `q1=−e1`, `q2=−e2+e1²`, `q3=−e3+3e1e2−2e1³`;
- `t2=2e1²−e2`, `t3=−e3+6e1e2−6e1³`.

For an axisymmetric log-temperature `q=εAμ`, exact angular normalization gives
`Θ1=Aμ`, `Θ2=A²(μ²−1/3)`, `Θ3=A³(μ³−μ)` before final dipole removal. The quadratic
quadrupole remains; deleting the logarithmic dipole first changes the map.

The tensor retarded Green function is
`G_T(η,s;K)=1_(η≥s)(s/η) sin[K(η−s)]/K`, with zero-momentum limit
`s(1−s/η)`. Tests also check cubic composition of an arbitrary endpoint scalar
and Gaussian polynomial moments independently of pairing enumeration.

## Integral-legality distinction retained

An almost-everywhere-defined cubic kernel is insufficient for
`K^(3)(k,q,−q)`. The restriction itself must exist. As the reviewer illustrates,
a local `|k+q|^(−p)` singularity can have the generic three-dimensional
integrability threshold `p<3` but the equal-shell antipodal threshold `p<2`.
Neither threshold is asserted to describe the actual final observable.

## Next review

The reviewer recommended the general-momentum quadratic boundary problem as
the next milestone. Review the new derivation independently, especially:

1. The exact fluid integration and ADM reduction in the fixed gauge.
2. The quadratic source coefficients and Fourier symmetrization.
3. The regular-boundary subtraction and all retarded series bounds.
4. Constraint propagation from the chosen regular class.
5. The exact K=0 sector and retention of forced trace-free response.

The favorable review of the **earlier plan** does not automatically approve
these newly derived results.
