# Four directions toward a minimal finite-wavelength nonlinear response

Recorded following the review and the cubic-source/propagation-adjoint reductions.

## Objective and limits

The target is an elegant reduction of the complete third-order temperature observable and its fourth-order power spectrum, not a literal vanishing of all nonlinear corrections. The established leading-gradient shell results already give

\[
 \frac{\tfrac13 C_2^{(13)}}{C_2^{(11)}}=4\sigma^2,
 \qquad C_2=\epsilon^2 C_2^{(11)}
 \left(1+\frac{29}{3}\epsilon^2\sigma^2\right)+\cdots.
\]

A promising stronger organization would express the full response through lower-order responses, explicit physical boundaries and dilation operators, plus a small number of irreducible tidal integrals. No such complete reduction has yet been established.

The four proposals below are directions to investigate, not established results or claims that nobody else has explored them. The existing source and propagation reductions are in the main report, §5.3.1 and §6.1.

## 1. One constrained adjoint for the whole experiment — first priority

Collect the metric, fluid, photon trajectory, emission intersection, observer clock, tetrad and sky normalization into a single constrained system. For a real angular test functional of the observable, use

\[
 \mathscr L[y,\lambda;\epsilon]
 =J[y;\epsilon]+\lambda\cdot E[y;\epsilon],
 \qquad E=0,\quad \mathscr L_y=0.
\]

Include physical boundary and initial-data constraints, not just the bulk equations. Gauge-fix the state description and remove redundant Einstein equations before interpreting the linearized inverse.

**Hypothesis to test.** Stationarity of the augmented functional allows its cubic response to be expressed using only first variations of the physical state and adjoint. This could eliminate explicit second- and third-order state responses while combining field, ray and measurement cancellations before expanding their separate kernels.

**Important limitation.** First-order adjoint variations still require Green operators and may carry the same comparable-momentum dependence as the second-order physical fields. Eliminating higher-order state variables is not automatically an evaluation of their integrals.

**First milestones.** Prove and symbolically check the polarized stationary-response identity; construct a boundary-safe residual system for this experiment; identify the background adjoint and its first variation; test a complete constrained measurement model before attempting a full angular spectrum.

## 2. Total-derivative search modulo the field and transport equations

Build a controlled basis of bilinear first-/second-order expressions and search for

\[
 \mathcal I_{\rm combined}
 =\frac{d\mathcal B}{dt}+\sum_A F_A E_A
 +\mathcal I_{\rm irreducible}.
\]

Here the equations include field evolution and constraints, null transport, emission and clock conditions. Apply the search to the combined observable, not an isolated final kernel, and retain the boundary functional explicitly.

**Goal.** Find a minimal remainder, ideally a few invariant structures, rather than a large collection of component expressions.

**Limitation.** Failure of a finite ansatz is not an impossibility proof. Record the derivative order, basis, boundary assumptions and equivalences used in every search.

## 3. Canonical photon transport, Lie series and commutators

Use the conformal photon Hamiltonian

\[
 \mathcal H(t,X,p)=\sqrt{\gamma^{ij}p_i p_j}
\]

to reorganize ordered transport into products of lower-order evolution, Poisson-bracket commutators and nested commutators. Include the emission, observer-frame and clock transformations in the composition.

**Hypothesis to test.** Symmetric product pieces cancel against exponentiation, normalization and endpoint composition, while commutators isolate genuinely ordered propagation effects. This extends the proposed scalar shuffle identities to the full transport operator.

**Limitations.** Equal acoustic frequencies do not make spatially displaced transport operators commute. Also, the naive on-shell massless canonical action is degenerate; a generating-function/eikonal construction must fix the appropriate boundary data and energy normalization rather than treating a vanishing null action as the observable.

## 4. A finite-wavelength tidal normal form

Generalize the successful local-curvature gradient construction to curvature and optical transport along the physical light ray. Separate

\[
 \text{dilation/frame effects}
 +\text{physical endpoint effects}
 +\text{integrated tidal response}.
\]

**Goal.** Replace component metric/ray redundancies with a small invariant basis. At finite wavelength, bilocal or ordered curvature integrals are more plausible than a purely local observer-curvature formula.

**Limitations.** Coordinate changes do not remove physical lensing or redshift evolution. Any observational/light-cone coordinate system must respect the accelerated radiation-fluid observer and the specified physical labels; familiar formulas with a geodesic observer cannot simply be imported.

## How to assess a possible obstruction

A surviving mixed-time phase is not a no-go theorem: it invalidates an automatic finite polynomial-moment closure, but the complete sum could still simplify. A more meaningful obstruction would be an analytically isolated, gauge-invariant response with a characteristic Green-function dependence that survives in the complete observable and cannot be reproduced by boundary or dilation terms. Resonance arguments must preserve the finite-time prescriptions; a late-time limit cannot be exchanged with angular integration without proof.

## Working order

Start with direction 1, then use direction 2 to simplify its remaining functional. Directions 3 and 4 are alternative organizations and possible ways to identify the physical remainder. Record both successful identities and failed restricted ansätze, without turning the latter into unsupported impossibility claims.

Initial work on direction 1 is recorded in `report/unified_adjoint_working_note.md` and checked by `solution/check_unified_adjoint.py`.
