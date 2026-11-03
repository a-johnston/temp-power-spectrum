# Direct cubic phase reduction and the exact dilation Ward channel

## Status

The previous exact angular-source checkpoint is committed as **`d5a18c3`**. This continuation supplies the direct first-field-only measurement block needed beside that source in the stationary functional.

There are two analytical reductions beyond simply rearranging angular kernels:

1. The sole four-time ray contribution has one time integral evaluated exactly. The complete direct measurement can then be organized in products involving at most three independent times.
2. With one finite-momentum input and two **exactly constant** primordial inputs, the complete physical cubic response is a radial Euler derivative squared of the known linear response. In the field-reduced formulation, its first field adjoint is likewise a derivative of the background adjoint, not another unknown response.

A phase-sector implementation assembles the complete direct measurement, including physical roots, clock, frame and normalization, without expanding dense components prematurely. Exact finite-phase contraction is implemented. Dense generic noncoplanar standalone `F_r` materialization remains expensive and was not completed in the exploratory timed runs; the unexpanded phase plan is intended for combining with the other stationary blocks first.

This is **not** an evaluated generic Einstein cubic map, shell spectrum, or convergence result. The main report PDF is unchanged.

## 1. Contract spatial jets before assembling the cubic response

For three first Fourier inputs write

\[
 G_a(t,X)=[2z_a I+h_a(t)]e^{ip_a\cdot X},\qquad
 \nu_a(t,X)=a_a(t)e^{ip_a\cdot X}.
\]

The first profiles are arbitrary smooth input amplitudes for the measurement identity. To obtain the physical radiation response they must be the prescribed regular linear radiation solutions. Polynomial off-shell fixtures below are exact functions, not truncated radiation transfer functions.

All spatial jets of a single input reduce immediately to momentum contractions:

\[
 (X_b\cdot\nabla)G_a=i(p_a\cdot X_b)G_a,\quad
 (X_b\cdot\nabla)h'_a=i(p_a\cdot X_b)h'_a,
\]
\[
 X_b^iX_c^j\partial_i\partial_j h'_a
 =-(p_a\cdot X_b)(p_a\cdot X_c)h'_a.
\]

The same replacement applies to the linear direction driver with `n` held fixed. Thus the complete cubic driver can be polarized directly as six trilinear terms and three single–pair terms, rather than obtained by expanding a large diagonal expression and differentiating its components.

The polarized second coordinate forcing, excluding its radial contribution, is

\[
 Q^{\mathrm{rest}}_{X,bc}=-\tfrac12(G_bv_c+G_cv_b)
 -\tfrac18(G_bG_c+G_cG_b)n
 +\tfrac i2[(p_c\cdot X_b)G_c+(p_b\cdot X_c)G_b]n.
\]

The excluded term is `n(v_b.v_c)` and is restored by the evaluated kernel in the next section. The complete transverse `Q_V,bc` is the Fourier specialization of the preceding generic formula, with no transverse or connection term discarded.

Implementation is in `bilinear_ray_forcing` and `polarized_bulk` in `solution/direct_cubic_measurement.py`. Their contracted identity is checked against the previous generic spatial-jet formula at **nine independent symbolic momentum components** and arbitrary first-ray probe vectors. Every single–pair term has coefficient one; there is no extra factor three per partition.

## 2. Evaluate the common fourth ray time exactly

Let the phase-stripped first direction and gradient amplitudes be `V_a(u)` and `r_a(u)` from the preceding angular-source note:

\[
 v_a(t)=-\int_t^o V_a(v)e^{ip_a\cdot n(o-v)}dv,\qquad
 R_a(t)=\int_e^t r_a(u)e^{ip_a\cdot n(o-u)}du.
\]

The diagonal redshift third derivative contains

\[
 U_{3,R}=-3\int_e^o(n\cdot R(t))|v(t)|^2dt.
\]

Its polarized kernel is `-sum_a integral (n.R_a)(v_b.v_c) dt`. For each triple of Volterra times, the common original time lies between the `R` time and the earlier of the two `v` times. It appears nowhere else in this integrand. Its integral is therefore the **length of that interval**, not an unevaluated response:

\[
\boxed{\begin{aligned}
 U_{3,R;abc}={}&-\sum_{(a,b,c)\in S_3}
 \int_e^o dw\int_e^w dv\int_e^v du\,(v-u)\\
 &\times(n\cdot r_a(u))[V_b(v)\cdot V_c(w)]
 e^{iK_{abc}(u,v,w)\cdot n},\\
 K_{abc}={}&p_a(o-u)+p_b(o-v)+p_c(o-w).
\end{aligned}}
\]

The six orderings include the three choices of `R` and both orderings of the two `v` times. No further factorial remains. The two minus signs in the first `v` responses multiply to plus, leaving the displayed overall minus sign.

This is an exact finite-time Fubini identity for the specified smooth first profiles. It does not assume equal frequencies, small wavelengths or a late-time limit. The code verifies its sign and factor against an independently integrated nonzero exponential fixture.

**This is not removal of the unit-direction constraint.** The separate local contribution `-3 |v|^2 h'_nn` in the direct cubic driver stays present. Only the `n|v|^2` part of the second coordinate forcing has had its common integration time evaluated.

All other bulk terms have at most two first Volterra responses plus their common local time. Physical root, clock, frame and normalization products have at most three independent times as well. The three-time bound here concerns **the direct measurement block**, not the subsequent field Green pairings of the complete stationary functional.

## 3. Assemble all physical measurement operations before contraction

`first_only_exponent_triple` returns the first exponents `I_a`, polarized second exponents `I_bc` and the complete first-only third exponent `I_abc`.

It includes:

- both physical emission-root orders and spatial pullback of the emission lapse;
- total background-ray derivatives rather than partial field-time derivatives at emission;
- both proper-clock orders and the physical initial-temperature event;
- the full noncommuting observer-frame boundary, on the observer worldline;
- both mixed observer–emission derivatives.

`compose_mixed_exponents` is checked against the earlier diagonal composition formula for arbitrary functions of observer time, not only homogeneous metrics where some mixed derivatives vanish.

For `u_a=I_a-average(I_a)` and `b_bc=I_bc-average(I_bc)`, the third normalized map before its **final** monopole centering is

\[
 -I_{abc}+\sum_{a|bc}u_ab_{bc}-u_au_bu_c
 +\sum_{a|bc}\langle u_bu_c\rangle u_a.
\]

There is no intermediate dipole or higher-multipole projection. A quadrupole weight annihilates the final monopole term, so omitting only that final centering from the weighted objective is exact. The frame is already included; adding the separate frame-adjoint representation would double count it.

## 4. Phase plans, bound variables and entire angular moments

`solution/phase_integral_contraction.py` expands finite products into sectors

\[
 \int\cdots\int P(n;\mathbf t)e^{iK(\mathbf t)\cdot n}\,d\mathbf t.
\]

Each sector contains an explicit sky polynomial, phase vector and ordered definite limits. Repeated copies of an integral acquire **independent bound times**, including a squared integral written with one source variable. Nested causal bounds retain their dependencies. Already contracted means remain scalar factors: their products do not become a single combined-phase average.

The only additional time simplifications are ordinary polynomial integrations in variables absent from the phase, and polynomial triangular convolutions such as

\[
 \int_e^o dt\int_e^t du\,P(t,u)e^{iK(u)\cdot n}
 =\int_e^o du\left[\int_u^o P(t,u)dt\right]e^{iK(u)\cdot n}.
\]

No oscillatory, acoustic or field Green integral is silently solved by this operation.

Angular contraction uses the preceding exact Cartesian moment formula and entire `F_r(K.K)`. A conservative polynomial degree count gives `F_0,...,F_13` for the direct cubic quadrupole before optional coefficient cancellations. This is a **safe upper bound**, not a claim of a minimal basis or thirteen fixed moments. Time-independent sky coefficients are protected from deep expansion when their internal functions or contracted integrals are irrelevant to collecting powers of `n`.

Use:

- `measurement_phase_plan(data)` to retain the unexpanded polynomial/phase sectors for combination with field-adjoint blocks;
- `contract_direct_measurement(data,C)` to materialize the quadrupole objective in entire radial functions.

The generic dense expansion should not be done simply to produce a longer final kernel. Exploratory full noncoplanar materializations timed out; they are not counted as completed tests. An exact noncoplanar fixture assembles **41 sectors**, with no remaining fourth ray time, but leaves the dense final radial expansion intentionally deferred. Forty-one is a fixture count, not a universal closure count.

## 5. A complete finite-phase measurement calibration: two constant inputs

Take the two other primordial directions to be constant, with `z=1`, `a=h=0`. The finite-momentum direction has arbitrary time profiles, and its field profile is held fixed in this subsection.

A constant spatial conformal factor `s` rescales the background ray to `X=e^{-s}(o-t)n`. Therefore the direct measurement with one finite-momentum input and two constant inputs is the **phase-only** Euler derivative squared of its linear measurement.

Writing `x=(p.n)(o-u)`, its one-time sky kernel is

\[
\boxed{
 \Theta^{(3)}_{111}(p,0,0;n)
 =P_0\,\frac12\int_e^o h'_{p,nn}(u)
 [x^2-ix]e^{ix}du,}
\]

where the profile `h_p` is held fixed. The two-time `R` contribution collapses by the triangular identity and combines with the direct Hessian term to give this expression.

The constant–constant **field** residual is zero, but the second coordinate forcing is **not**: `Q_X,00=-n`. Its pairing with `R_p` supplies the `-ix` term above after the common time is integrated. Dropping the entire constant pair would fail this identity.

The executable check establishes this identity at arbitrary symbolic momentum, arbitrary endpoint times and an arbitrary first time profile. Its entire angular contraction agrees with `(alpha partial_alpha)^2` of the linear phase kernel at `alpha=1`. This is not a gradient or finite-time-series check.

## 6. The complete physical dilation Ward channel

For the **physical** response, field profiles must also transform. Set

\[
 g(X)=s+\epsilon e^{ip\cdot X},\qquad Y=e^sX,\qquad q=e^{-s}p.
\]

The coordinate change removes the constant part of the spatial conformal factor and turns the remaining primordial wave into `epsilon exp(iq.Y)`. The regular radiation solution has the same form with its linear profiles evaluated at `q`.

This transformation preserves the full experiment:

- the observer stays at the coordinate origin;
- emission and initial-clock surfaces are fixed physical temperatures;
- proper time is a scalar and the observer-frame orientation is unchanged by a uniform spatial scale;
- the nonlinear regular matching transforms covariantly, including its `exp(-2g)` factor;
- sky normalization and the prescribed final multipole projection commute with the dilation.

Uniqueness of the specified regular branch then gives, to first order in the nonconstant input but to all orders in the constant `s`, the linear sky evaluated at `e^{-s}p`. With `D_p=p.partial_p`,

\[
\boxed{
 \Theta^{(3)}(p,0,0;n)=D_p^2\Theta^{(1)}(p;n).}
\]

For each retained multipole, this reduces the cubic transfer to

\[
 \boxed{J^{(3)}_{\ell m}(p,0,0)=D_p^2J^{(1)}_{\ell m}(p),\qquad
 T^{(3)}_\ell(k,0,0)=(k\partial_k)^2T_\ell(k).}
\]

The derivative acts on **both** the field transfer profiles and the ray phase. It is not merely the phase-only identity in the previous section. The existing finite-time entire representation of the linear transfer therefore determines this special cubic channel without a new nonlinear time integration.

These are two exactly zero input momenta, not a proof about exchanging soft limits with observation time or spectral integrals. In particular, the channel is not the nonzero-radius monochromatic-shell Wick contraction and cannot replace an arbitrary-triple map.

### Source and initial-data identities

Let `D_hat` act on explicit spatial momenta while holding the first field profile fixed. Every spatial term in the frozen linear field residual is homogeneous of degree two. Consequently

\[
 C_{p0}=-\widehat D_pF_1,\qquad
 C_{p00}=\widehat D_p^2F_1=-2C_{p0},\qquad C_{00}=0.
\]

The initial rows are explicitly

\[
 c_{p0}=2p^2/9,\qquad c_{p00}=-4p^2/9=-2c_{p0}.
\]

These are checked with the exact Einstein residual engine at symbolic `p`, arbitrary lapse profile and five independent STF field profiles. No time or gradient expansion is used, and the initial rows are not discarded.

### No new first field adjoint for the constant direction

In the **field-reduced** system, the constant direction has zero first lapse/spatial field response; its nonzero coordinate-ray displacement has already been included in the reduced measurement. Its field operator and linear measurement kernel are the background ones evaluated at `q=e^{-s}p`. Differentiating the background adjoint equation, including the regular initial row, gives

\[
 \boxed{
 \Lambda_{N,\mathrm{const}}=-D_p\Lambda_{N,0},\qquad
 \Lambda_{H,\mathrm{const}}=-D_p\Lambda_{H,0},\qquad
 \Lambda_{c,\mathrm{const}}=-D_p\Lambda_{c,0}.}
\]

This assertion is about the reduced field/initial covectors. It is **not** an assertion that the coordinate ray response vanishes or that every unreduced geometric adjoint transforms without scaling weights.

Thus only two nonzero **field/initial single–pair** terms remain for `(p,0,0)`; both reuse derivatives of the known background adjoint. The nonzero constant–constant coordinate forcing is already included in the direct measurement and is not deleted by `C_00=0`. The direct source is `-2 C_p0`, and the initial source has the same relation. Their complete stationary sum, with the direct phase contribution, reassembles into `D_p^2` of the linear observable. It is generally nonzero.

## 7. Verification and remaining work

```bash
.venv/bin/python -B solution/check_phase_integral_contraction.py
.venv/bin/python -B solution/check_direct_cubic_measurement.py
.venv/bin/python -B solution/check_dilation_ward.py
```

All three new checks and fourteen existing targeted regressions passed, including the symbolic arbitrary-triple physical four-gradient comparison, exact angular source, stationary assembly, first/unified/propagation adjoints, Cartesian Green, observable expansion, report identities, shell sources against the independent engine, Fermi–Walker conformal transport and Bessel translation. Repository-wide Ruff lint/format and type checks, runner shell syntax and whitespace checks passed. The full long-running suite was not rerun.

The checks distinguish:

- exact phase moments, bound-variable independence, causal domains and products of means;
- the exact evaluated common-time kernel and arbitrary-function endpoint/clock polarization;
- exact Fourier-specialized local cubic algebra at nine symbolic momentum components;
- the exact finite-wave measurement dilation identity at arbitrary times/profile;
- a complete homogeneous anisotropic **off-shell** measurement comparison using conserved photon momentum and independent scalar physical roots, with noncommuting Fermi–Walker rotation;
- the exact Einstein-source dilation identities;
- a noncommuting four-state **algebraic fixture**, not a radiation solution, which independently solves the first adjoint source equation and reproduces the complete stationary Ward identity with a nonzero initial-data row.

The physical Ward result is proved by coordinate covariance of the regular physical experiment, not inferred from the algebraic fixture or the noncoplanar phase-plan example.

For generic triples, the next step remains to merge the direct phase sectors with the first-source Green pairings and direct field/initial terms, then seek cancellations before materializing or evaluating the remaining kernels. General mixed-time/Green integration, the unrestricted physical comparison, all-multipole cubic evaluation and the fourth-order spectrum remain open. No arbitrary-spectrum UV/IR convergence or late-time integral exchange is asserted.
