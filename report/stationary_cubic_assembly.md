# Stationary cubic assembly without second physical fields or second rays

## Checkpoint and scope

The preceding source-reduction/unified-adjoint work was committed as `76389f6`.
This note records the continuation after that checkpoint.

The general stationary cubic functional is now assembled with an explicit **direct first-field-only measurement term**. Its second ray response is eliminated, including the physical frame boundary. The field correction and nonlinear initial-data terms are then added with the correct arbitrary-triple partitions.

Implementation: `solution/stationary_cubic.py`.
Checks: `solution/check_stationary_cubic.py`.

This is an organization of the exact finite-order response, not a claim that its remaining Green/sky integrals cancel or are evaluated. The measurement sector is independently cross-checked against conserved canonical photon momentum in a homogeneous anisotropic off-shell metric. The unrestricted finite-wavelength Einstein–radiation arbitrary-triple observable has **not** yet been evaluated by the new assembly. A subsequent physical comparison of its complete quadrupole through **four spatial gradients**, at arbitrary symbolic input triples, is now in `quadrupole_adjoint_gradient_expansion.md`; this is stronger than the off-shell fixture but remains a gradient expansion. The main report PDF is unchanged by these continuations.

## 1. The assembled functional and its polarization

Use a repeated first direction initially. Let

\[
 \zeta=\epsilon g,\quad \nu=\epsilon a,\quad H=\epsilon h
\]

be the **first-field-only family**, and let `J_111` be the third derivative of the complete physical normalized measurement on this family. It is deliberately off shell beyond first order. Nonlinear Einstein evolution is restored by the stationary correction, not by pretending that higher physical fields vanish.

This is the **field-reduced** stationary identity. The geometric/measurement adjoint terms have been absorbed into `J_111`; adding those full augmented-system terms separately would double-count them.

Let `C_N,r,C_H,r` be the second/third residual derivatives of the independent field equations along that family. The spatial residual is `2 exp(H/2) E_T exp(-H/2)`, as in the preceding notes. Let `c_r` be the corresponding nonlinear initial Hamiltonian coefficients. Then

\[
\boxed{
 J^{(3)}=J_{111}
 +\int\!dt\,[\Lambda_{N,0}C_{N,3}+\Lambda_{H,0}:C_{H,3}
 +3\Lambda_{N,1}C_{N,2}+3\Lambda_{H,1}:C_{H,2}]
 -\Lambda_{c,0}c_3-3\Lambda_{c,1}c_2.}
\tag{1}
\]

The spatial/Fourier dual pairings in this expression are understood; the initial terms also have that pairing. The adjoints are precisely those defined in the earlier notes, with first-field source `S_a` from `first_field_adjoint_source.md`.

For three distinct primordial directions the actual kernel is

\[
\boxed{
 J_{abc}=J_{111,abc}
 +\int_0^o\!dt\,[\Lambda_{N,0}C_{N,abc}+\Lambda_{H,0}:C_{H,abc}
 +\sum_{a|bc}(\Lambda_{N,a}C_{N,bc}+\Lambda_{H,a}:C_{H,bc})]
 -\Lambda_{c,0}c_{abc}-\sum_{a|bc}\Lambda_{c,a}c_{bc}.}
\tag{2}
\]

There are three partitions, each with coefficient **one**, not three. For each partition the background adjoint uses total momentum `k_a+k_b+k_c`, while the first adjoint is reconstructed at probe momentum `k_b+k_c`. The latter can be zero. Neither its lapse sector nor its initial covector is deleted.

`assemble_polarized_cubic` implements (2). It accepts the already polarized source coefficients and exactly three single–pair partitions, and returns the specified definite time integral plus measurement and initial terms.

### Two different coefficient conventions

To obtain the arbitrary-triple kernel from a **diagonal third derivative** `J^(3)[g,g,g]`, substitute `g=alpha g_a+beta g_b+gamma g_c` and take

\[
 J_{abc}=\tfrac16[\alpha\beta\gamma]\,J^{(3)}[g,g,g].
\tag{3}
\]

By contrast, taking three distinct labels directly from the **full Taylor observable/residual** already cancels the `1/3!`; no division remains. The helpers `polarize_diagonal_cubic` and `direct_field_residuals` deliberately use these different conventions and the check script tests their agreement.

## 2. First ray and local notation

All quantities below are evaluated on the background ray `X_0=n(o-t)` after taking their field derivatives. The first ray responses `X,v` are at fixed coordinate observation time, with `X(o)=v(o)=0`, before physical-clock composition. The vector `v` is transverse to `n`.

Set

\[
 G=2g\,1+h,\quad P=1-nn^T,\quad
 r=\tfrac12\nabla h'_{nn},\quad s=-P h'n,
\]
\[
 R(t)=\int_e^t r,\qquad K(t)=\int_e^t(s-R).
\]

Here `e=bar t_star`, `o=bar t_o`, and the physical initial time is denoted by `iota=bar t_I`; Fourier imaginary units are unrelated to `iota`. Spatial derivatives act only on fields, holding `n,X,v` fixed. Primed fields are partial time derivatives, while ray/adjoint primes are total pulled-back derivatives.

The first direction equation is

\[
 v'=\mathcal G_1=\tfrac12P h'n
 +\tfrac12[\nabla G_{nn}-(n\cdot\nabla G)n].
\tag{4}
\]

This includes the primordial spatial-gradient contribution; the first lapse does not enter conformal null propagation directly.

## 3. Explicit second-ray forcing, with no second ray left to solve

Write the constrained second direction as

\[
 v_2=u_2+n|v|^2,\qquad n\cdot u_2=0.
\]

The radial part is fixed by unit normalization and cannot be omitted. Define

\[
 Q_X=n|v|^2-Gv-\tfrac14G^2n+(X\cdot\nabla G)n.
\tag{5}
\]

Then `X_2'=u_2+Q_X`. The transverse second equation is `u_2'=Q_v`, with

\[
 Q_v=P\Big\{
 2(X\cdot\nabla)\mathcal G_1-h'v+h'_{nn}v-\tfrac14[G,h']n
\]
\[
 \quad+(v\cdot\nabla G)n+(n\cdot\nabla G)v
 -2\nabla(v^TGn)
 +\tfrac12((Gn)\cdot\nabla G)n
\]
\[
 \quad+\tfrac14[G,n\cdot\nabla G]n
 -\tfrac12G\nabla G_{nn}\Big\}.
\tag{6}
\]

All terms in (5)–(6) use only first fields and first ray responses. The commutator in (6) with a spatial derivative is distinct from the temporal frame commutator; the former need not vanish on the equal-frequency shell.

The code checks (5)–(6) directly against the exact coordinate geodesic transformed by `v=exp(G/2)e`, using independent symmetric metric and spatial derivative jets at a point. This is an arbitrary-jet identity, not a short-time test.

## 4. Direct cubic redshift and the frame boundary

Define the known local remainder

\[
 L_{111}=\tfrac18[h,[h,h']]_{nn}
 +\tfrac32X^iX^j\partial_i\partial_jh'_{nn}
 -6v^T(X\cdot\nabla h')n
 +3v^Th'v-3|v|^2h'_{nn}.
\tag{7}
\]

The exact first-field-only third driver is

\[
 F_3=3(r\cdot X_2+s\cdot u_2)+L_{111}.
\]

The same adjoints `R,K` give

\[
 \frac d{dt}(R\cdot X_2+K\cdot u_2)
 =r\cdot X_2+s\cdot u_2+R\cdot Q_X+K\cdot Q_v.
\]

At emission `R=K=0`; at observation `X_2(o)=0`, but the frame gives

\[
 u_2(o)=-R_{{\rm FW},2}n,\qquad
 R_{{\rm FW},2}=\tfrac14\int_\iota^o[h(t,0),h'(t,0)]dt.
\tag{8}
\]

The integral in (8) is on the observer worldline, not on the null ray. Its background endpoints suffice: their first shifts change the third frame response, which acts on the vanishing background sky and does not enter the cubic observable. The second frame response itself must be kept. Hence the direct cubic line-of-sight term is

\[
\boxed{
 U_3=\int_e^o[L_{111}-3R\cdot Q_X-3K\cdot Q_v]dt
 -3K(o)\cdot R_{{\rm FW},2}n.}
\tag{9}
\]

No second ray survives. The last term is retained, not set to zero using the vanishing of the coordinate ray at the observer. It is the same frame sector as (18) of `first_field_adjoint_source.md`, written as a ray-adjoint boundary, so the two representations must **not** be added a second time.

The lower-order fixed-coordinate quantities needed next are

\[
 F_1=\tfrac12h'_{nn},\quad U_1=\int_e^oF_1dt,
 \qquad F_2=(X\cdot\nabla)h'_{nn}-2v^Th'n,\quad U_2=\int_e^oF_2dt.
\]

## 5. Physical roots, second clock and all mixed products

For the first-field-only family the physical emission root gives

\[
 d_1=-e a_e,\qquad
 d_2=e[a_e^2-2X_e\cdot\nabla a_e-2d_1(D_ta)_e],
 \quad D_t=\partial_t-n\cdot\nabla.
\tag{10}
\]

At the observer the complete first and second clock shifts simplify to

\[
 D_1=-\frac{\iota^2a(\iota,0)+\int_\iota^ot a(t,0)dt}{o},
\]
\[
\boxed{
 D_2=\frac{2\iota^3a(\iota,0)a'(\iota,0)
 -\int_\iota^ot a(t,0)^2dt-2oa(o,0)D_1-D_1^2}{o}.}
\tag{11}
\]

This includes the perturbation of the physical initial-temperature event; the initial quadratic terms have been simplified using its root equation, not deleted.

The physical exponent derivatives are

\[
 I_1=U_1,\qquad I_2=U_2-2F_{1,e}d_1+2D_1\partial_oU_1,
\]
\[
 I_3=U_3-3F_{1,e}d_2-3F_{2,e}d_1-3(D_tF_1)_e d_1^2
\]
\[
 \quad+3D_1\partial_o(U_2-2F_{1,e}d_1)
 +3D_2\partial_oU_1+3D_1^2\partial_o^2U_1.
\tag{12}
\]

All observer derivatives in (12) are taken at fixed coordinate observer before composing with the clock. They differentiate the known first ray integrals and their boundary data. In particular,

\[
 -6D_1[(\partial_oF_{1,e})d_1+F_{1,e}\partial_od_1]
\]

is present. It is not inferred from the homogeneous verification example, where it happens to vanish; it is checked separately by an arbitrary endpoint/observer-jet identity. Third endpoint and clock shifts are unnecessary because the background redshift driver vanishes.

For `u=I_1-<I_1>` and `b=I_2-<I_2>`, exact sky normalization gives

\[
 \Theta_{111}^{(3)}=P_0[-I_3+3ub-u^3+3\langle u^2\rangle u].
\tag{13}
\]

Thus `J_111=int w Theta_111^(3) dOmega`, with only the prescribed final harmonic/dipole operation. Equations (4)–(13) specify the direct measurement term in (1), without second physical fields or second rays.

## 6. Direct field residuals and an additional shell simplification

The direct field coefficients in (1)–(2) are finite local operations, not new differential equations to solve. The helper `direct_field_residuals` passes first-label-only fields through the exact Christoffel engine, converts the spatial residual to the orthonormal frame, and returns the residual and nonlinear initial-data fields. Querying two or three labels gives the physical-input coefficients directly.

Every primordial input label is included here. These are **not** the state-probe derivatives with `delta zeta=0` used to construct the first adjoint source.

There is a useful derivative reduction. Along the first-field-only family the direct cubic lapse kinetic invariant vanishes:

\[
 D^3\operatorname{tr}(S^2)=0.
\]

The direct temporal part of the cubic spatial residual is exactly

\[
\boxed{
 C_{H,3}^{\rm temporal}
 =\tfrac14[h,[h,h'']]+[h',[h,h']]
 +\frac1{2t}[h,[h,h']].}
\tag{14}
\]

The remaining direct lapse contribution consists of spatial geometry and `-24a^3/t^2`; the remaining spatial contribution is the local curvature/lapse-gradient expression. Equation (14) is checked in the free noncommuting algebra, not by assuming the linear field equation.

On the monochromatic shell `h(t,X)=f(t)B(X)`, **every temporal commutator in (14) vanishes pointwise**, before angular averaging. The temporal commutators in (7), (8) and (6) vanish as well. Thus the direct three-first-field Einstein source has no time derivatives of first fields on this shell. Spatial nonlinearities and transport terms remain; this is not a full cancellation.

The nonlinear initial coefficients used in (1) are

\[
 c_2[g,g]=\tfrac19[|\nabla g|^2-4g\nabla^2g],\qquad
 c_3[g,g,g]=\tfrac13[4g^2\nabla^2g-2g|\nabla g|^2].
\tag{15}
\]

Their arbitrary-label kernels agree with the report after the appropriate diagonal polarization factors. Both negative signs in (1)–(2) come from the residual `lim(nu/t^2)-c[g]`.

## 7. Independent cubic observable comparison

A useful off-shell verification metric is spatially homogeneous but anisotropic:

\[
 \gamma=\exp(H(t,\epsilon)),\qquad N=a_1t\exp(\nu(t,\epsilon)),
\]

with noncommuting symmetric `H_1(t)`, prescribed nonzero `H_2,H_3,nu_2`, and exact polynomial time profiles. The polynomials are the specified test functions at arbitrary finite time, **not** a time-series approximation to an Einstein solution.

Spatial homogeneity conserves canonical photon momentum. With the physical observer normalization the exact endpoint expression is

\[
 e^{-I(n)}=
 \left[n^TR^Te^{H_o/2}e^{-H_e}e^{H_o/2}Rn\right]^{-1/2},
\tag{16}
\]

where the endpoint matrices are evaluated at the physical roots/clock times, and `R` is the Fermi–Walker rotation. This provides an independent measurement calculation without solving the ray equations or invoking the adjoint formula.

The test script expands (16) in amplitude through cubic order and evaluates the sky normalization/harmonic weight with exact sphere monomial moments. It verifies:

1. With higher fields zero, (9)–(13) reproduce the independent cubic measurement.
2. With nonzero `H_2,H_3,nu_2`, the independent full measurement satisfies

\[
 J^{(3)}_{\rm meas}=J_{111}
 -3S_{{\rm meas},1}[\nu_2,H_2]+J_{{\rm linear}}[H_3],
\tag{17}
\]

using the **actual source implementation** from `unified_adjoint_source.py` and the exact random-mean normalization. This tests compatibility of the new cubic measurement term with the previously assembled first source.
3. The frame contribution is nonzero for the selected test weight, so deleting it fails; the clock and both physical endpoint shifts are included.

This metric is not an Einstein–radiation solution. The test verifies the complete measurement algebra and its assembly, not a physical spectrum or a finite-wavelength Einstein cancellation. General spatial dependence in (5)–(9) is checked separately using arbitrary metric/spatial jets; the homogeneous example is not used to infer those identities.

## 8. Verification and remaining work

```bash
.venv/bin/python -B solution/check_stationary_cubic.py
```

The script checks the arbitrary local ray identities, radial normalization, free-algebra temporal field identity, exact endpoint/clock composition, physical-input spatial field sources in the conformal sector at three arbitrary momenta, nonlinear initial data, both polarization conventions, the executable arbitrary-triple assembly, and the independent full measurement comparison above.

The previous generic stationary-response proof establishes the field-response elimination in (1)–(2). The subsequent `quadrupole_adjoint_gradient_expansion.md` evaluates the leading first adjoint and establishes the complete physical cubic quadrupole comparison through four gradients for arbitrary triples. A full finite-wavelength/all-multipole evaluation is still required before making a new physical completion claim.

The new script and nine existing targeted checks passed: first field-adjoint source, unified adjoint, cubic product identities, orthonormal driver, propagation adjoint, Cartesian Green, observable expansion, report identities and shell cubic sources with `--compare-engine`. Repository-wide Ruff lint/format and type checking also passed, as did shell syntax and whitespace checks. The full long-running test suite was not rerun.

Next steps:

- insert the now-closed exact quadrupole angular source from `exact_quadrupole_angular_source.md` into the first field-adjoint Green pairing for the shell partitions;
- extend the physical comparison beyond the quadrupole/four-gradient sector, preserving the homogeneous lapse, initial covectors and physical endpoint traces;
- assemble the direct cubic multi-phase measurement contractions and seek cancellation between direct residuals, adjoint products and measurement terms before evaluating the remaining time integrals.

All results here are finite-time response identities. No late-time exchange or arbitrary-spectrum UV/IR convergence is claimed. The full finite-wavelength spectrum remains unevaluated by this route.
