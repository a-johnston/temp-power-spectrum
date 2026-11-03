# Evaluated regular Green blocks and a three-time cubic quadrupole bound

## Status

The direct-phase/dilation checkpoint is committed as **`ba28b26`**. This continuation evaluates previously nested Green operations and inserts the exact angular source into the resulting transpose.

**New analytical results:**

- An explicit, entire, projector-free Cartesian spatial Green operator.
- An evaluated lapse-to-spatial Green convolution, including its regular-origin column.
- An evaluated complete proper-clock Green pairing, including the physical initial-event term; also an evaluated first clock shift for the regular physical scalar mode.
- The background lapse, spatial and origin adjoints now each require only **one** source-time integration. The complete first field adjoints require at most **two**.
- With the preceding direct-measurement reduction, the complete stationary cubic **quadrupole** can be written with at most **three independent time integrations per distributed term**. This is a bound on total time arity, not merely nesting depth.

The remaining two-/three-time source pairings have **not** been evaluated, nor has a new unrestricted finite-wavelength spectrum. No old nine-moment closure or complete nonlinear cancellation follows from the bound. The main report TeX/PDF is unchanged.

Implementation: `solution/regular_field_green.py`, checked by `check_regular_field_green.py` and `check_quadrupole_green_source.py`.

## 1. One entire factorial-tail family

Use

\[
 T_j(a,d)=\sum_{m\ge0}\frac{(-a)^m d^{2m+j}}{(2m+j)!}
 =\frac{d^j}{j!}{}_1F_2\left(1;\frac{j+1}{2},\frac{j+2}{2};-\frac{ad^2}{4}\right).
\]

It is entire in both arguments. In particular,

\[
 T_j(0,d)=d^j/j!,\quad
 \partial_dT_j=T_{j-1}\ (j\ge1),\quad
 \partial_dT_0=-aT_1,
\]
\[
 aT_{j+2}=d^j/j!-T_j,\qquad
 \partial_aT_j=-\tfrac12(dT_{j+1}-jT_{j+2}).
\]

`T_0` and `T_1` are the entire cosine and sine-over-frequency functions. Their series/hypergeometric definitions, **not formulas dividing by a square root**, define the zero and negative arguments as well. Indices through five suffice for the Green kernels and required first time derivatives. Optional momentum differentiation remains entire but may introduce higher indices.

These are known universal functions, not new state responses or unknown transfers. Away from `a=0` they are simply elementary trigonometric tails, for example `T_2=(1-cos(sqrt(a)d))/a` and `T_3=(d-sin(sqrt(a)d)/sqrt(a))/a`. The entire definitions supply their removable limits; the implementation does not evaluate them by these singular-looking quotients.

## 2. The complete spatial inverse is a quadratic operator polynomial

For momentum `q`, set `z=q.q`, `Q=TF(qq^T)` and

\[
 \mathfrak B_q[B]=\operatorname{TF}
 [zB-q(Bq)^T-(Bq)q^T].
\]

On STF matrices it obeys the exact identity

\[
 \mathfrak B_q^3=\frac{2z}{3}\mathfrak B_q^2+
 \frac{z^2}{3}\mathfrak B_q.
\]

Its tensor, vector and scalar eigenvalues are `z,0,-z/3`. Put `d=t-s`. The retarded spatial inverse of `partial_t^2+2 partial_t/t+B_q` is

\[
\boxed{\begin{aligned}
 \mathfrak G_q(t,s)={s\over t}\bigg[&d\,\mathbf1
 -{T_3(z,d)+3T_3(-z/3,d)\over4}\mathfrak B_q\\
 &+{3T_5(z,d)+T_5(-z/3,d)\over4}\mathfrak B_q^2\bigg].
\end{aligned}}
\]

No inverse power of `q.q`, momentum-direction projector or special zero-momentum prescription occurs. The three eigenvalue evaluations reproduce `s T_1(lambda,d)/t` exactly, and the kernel has zero diagonal value and unit diagonal first derivative. Self-adjointness follows from that of `B_q` under the Cartesian trace pairing.

At `q=0`, this is simply `s(t-s) B/t` on every STF component. It does **not** delete a homogeneous tensor response from the uncontracted map.

## 3. Evaluate the triangular lapse-to-spatial convolution

Write `a=z/3`. The inverse of `3 partial_t^2+z-6/t^2`, with positive forcing convention, becomes

\[
\boxed{
 G_N(t,s)=\frac13\left[T_1(a,d)+
 \frac{dT_2(a,d)-T_3(a,d)}{ts}\right].}
\]

It has diagonal derivative `1/3`, as required by the leading coefficient three. At zero momentum it is `(t^3-s^3)/(9ts)`, retaining the singular source-origin behavior.

Since `Q` is a scalar eigenmatrix of the spatial operator, the previously nested convolution needs no general matrix integration. Define

\[
 \mathcal A(a,d)=\tfrac12[T_3(a,d)+T_3(-a,d)],
\]
\[
 \mathcal B(a,d)=\tfrac12[dT_3(a,d)+T_4(-a,d)-T_4(a,d)],
\]
\[
 \mathcal D(a,d)=\tfrac12[dT_4(a,d)-2T_5(a,d)+T_5(-a,d)].
\]

Then its **evaluated** value is

\[
\boxed{
 J_q(t,s):=\int_s^t\frac utT_1(-a,t-u)G_N(u,s)du
 ={s\mathcal A(a,d)+\mathcal B(a,d)+\mathcal D(a,d)/s\over3t}.}
\]

One derivation is the convolution table for `S_o=T_1(a,.)`, `S_h=T_1(-a,.)`, `C_o=T_0(a,.)`:

\[
 S_h*S_o=\mathcal A,\qquad
 S_h*(uS_o)=\mathcal B,\qquad
 S_h*(S_o-uC_o)=a\mathcal D.
\]

The last equality is entire at `a=0`; the implementation never divides by `a`. The same result is certified by the exact forced equation

\[
 (\partial_t^2+2\partial_t/t-a)J_q=G_N,
\]

and its zero causal initial data. Thus this is an evaluated definite Green integral, **not** a renamed unevaluated response.

The complete regular probes are now

\[
 \Phi_H=(0,\mathfrak G_q[B]),\qquad
 \Phi_N=(G_N,-4QJ_q),
\]
\[
\boxed{
 \Phi_c=(v_q,H_c),\quad
 v_q(t)=3[T_2(a,t)-T_3(a,t)/t],\quad
 H_c(t)=-12Q\mathcal D(a,t)/t.}
\]

They obey `9s Phi_N(t;s) -> Phi_c(t)`, including first time derivatives. The origin column has not disappeared. At zero momentum,

\[
 J_0(t,s)=\frac{s(t-s)^3}{18t}+
 \frac{(t-s)^4}{36t}+\frac{(t-s)^5}{180st}.
\]

The coupled spatial response itself has the additional factor `Q`, so it vanishes at zero momentum while the homogeneous lapse/origin response does not.

### A genuine physical sum cancellation, not deletion of a Green family

For the physical first scalar input, `c_1=-z/9=-a/3`. Its direct primordial spatial forcing contributes

\[
 H_{\zeta}=-2Q T_3(-a,t)/t.
\]

The regular initial column contributes `c_1 H_c`. Their sum is

\[
 \boxed{H_1=-2Q[T_2(a,t)-2T_3(a,t)/t],\qquad
 \nu_1=-zv_q/9.}
\]

All hyperbolic/negative-argument terms cancel **between the physical source and initial matching**. These are the known regular oscillatory radiation profiles, without inverse momentum projectors. Dropping the initial row would miss this cancellation.

This first-order example does not justify deleting hyperbolic functions from the off-shell inverse or adjoints. Analogous cancellations in the complete generic cubic sum remain to be established.

## 4. Evaluate the entire proper-clock Green pairing

For the already sky-contracted coefficient `A_a` of the preceding source, the clock functional is

\[
 -\frac{A_a}{o}\left[\int_i^o t\,m(t)dt+i^2m(i)\right].
\]

A primitive of `t G_N(t,s)` is

\[
\boxed{\begin{aligned}
 P_N(T;s)=\frac13\bigg[&T T_2(a,T-s)-T_3(a,T-s)\\
 &+\frac{(T-s)T_3(a,T-s)-2T_4(a,T-s)}s\bigg].
\end{aligned}}
\]

It vanishes at `T=s`. Thus the complete clock contribution to the lapse adjoint is

\[
 \boxed{\Lambda_{N,a}^{\rm clock}(s)=-\frac{A_a}{o}
 [P_N(o;s)-P_N(\max(i,s);s)+i^2\mathbf1_{s<i}G_N(i,s)].}
\]

For the origin column set

\[
 P_c(T)=3[T T_3(a,T)-2T_4(a,T)],\quad
 W_c=P_c(o)-P_c(i)+i^2v_q(i).
\]

Then `Lambda_c,a^clock=-A_a W_c/o`. At zero momentum,

\[
 W_c=(o^4+3i^4)/4,
\]

which is generally nonzero. For a physical first input `p`, its first proper-clock shift is itself evaluated:

\[
 \boxed{D_p=\frac{p^2}{9o}W_c(p^2/3;i,o).}
\]

`QuadrupoleAngularSource.assemble` accepts this as `clock_shift=`; the original general off-shell clock functional remains the default. The coefficient `A_a` is the known linear-quadrupole observer derivative from the preceding note. The factory still displays its linear integral rather than claiming to reevaluate the linear transfer here. No new Green-time integral remains in its clock pairing.

## 5. Insert the complete angular source without time integration by parts

`WeakGreenAdjoint` accepts `ContractedMeasurement` and the already Fourier-dual `FieldVariationCoefficients`. There is no sky callback: field terms and common clock coefficients cannot accidentally receive a second sky average.

For the spatial covector the construction is explicitly

\[
\begin{aligned}
 \Lambda_H(s)={}&-\int_s^o[\mathfrak G_q V_H+
 \partial_t\mathfrak G_q D_H]dt\\
 &+\int_{\max(e,s)}^o[\mathfrak G_q M_H+
 \partial_t\mathfrak G_q M_{H'}]dt\\
 &+\mathbf1_{s<e}\partial_e\mathfrak G_q(e,s)[M_{H',e}]
 +\partial_o\mathfrak G_q(o,s)[M_{H',o}].
\end{aligned}
\]

Here the derivatives act **only on the Green kernel**, never on the supplied weak-source coefficients. `RegularFieldGreen.spatial(B,derivative=1)` enforces this distinction even when `B` depends on time.

The lapse and origin covectors use the same complete weak action on the closed `Phi_N` and `Phi_c`, plus the evaluated clock expressions. The emission lapse evaluation, both emission/observer derivative evaluations and the initial-temperature point event stay explicit. Consequently

\[
 \Lambda_H(o^-)=M_{H',o},\qquad
 \Lambda_H(e^-)-\Lambda_H(e^+)=M_{H',e},\qquad
 \Lambda_c=9\lim_{s\to0}s\Lambda_N(s).
\]

One-sided branches are explicit. In particular, the value exactly at emission denotes the right trace, not a spurious half-weight.

For the background source `M_H'=C_P(t)/2` alone, the lapse covector reduces to

\[
 \boxed{\Lambda_{N,0}(s)=-2\int_{\max(e,s)}^o
 (Q_P:C_P(t))\,\partial_tJ_P(t,s)dt,}
\]

with the origin expression obtained from `H_c` in the same way. Both now have one integral rather than a nested lapse/spatial Green pair.

`physical_quadrupole_adjoints(p,q,C,t,u,e,o,i)` supplies:

- the background adjoint at **total** momentum `P=p+q`;
- the complete physical first source from the entire radiation profiles and angular coefficients;
- its reconstruction at **probe/pair** momentum `q`, including `q=0`;
- the stationary pairing `integral Lambda_N C_N + Lambda_H:C_H - Lambda_c c`.

Every reconstruction uses a fresh time dummy. In particular, querying the background covector at the very symbol `t` used by the outer source does not set both Green times equal.

## 6. Why the complete cubic quadrupole has time arity at most three

Count **independent integration variables after distributing finite sums/products**. A product of three separate one-time means counts as three, not one. Exact angular moment contraction adds no time variable.

| Block | Maximum time arity per term |
|---|---:|
| Background `Lambda_N,0`, `Lambda_H,0`, `Lambda_c,0` | 1 |
| Field-variation weak coefficients, linear in the background covector | 1 |
| Measurement weak coefficients, after first-ray Volterra substitution | 1 |
| Complete first field/origin adjoints with the closed probes | 2 |
| Background adjoint paired with the direct local cubic field source | 2 |
| First adjoint paired with a local second field source | 3 |
| Nonlinear initial rows | at most 2 |
| Complete direct cubic measurement, by the preceding common-time reduction | 3 |

A separate frame-adjoint representation has at most two times; the direct measurement already includes its equivalent boundary, so it is not added again.

The direct second/third Einstein residuals and initial kernels are local expressions in known first profiles, not hidden time responses. These observations establish the bound for the complete stationary quadrupole functional. It does **not** depend on equal frequencies, a shell contraction, a gradient expansion or the special dilation channel.

The code checks a real Einstein single–pair source at a noncoplanar unequal-frequency triple, not merely an abstract forcing symbol. It independently checks the complete direct normalized measurement at noncoplanar unequal-frequency physical first inputs, retaining clock and frame terms. These fixtures corroborate the structural bound; they are not an independently evaluated generic physical cubic sky.

All reorganizations here are finite-time, fixed-momentum identities on the prescribed regular branch. The regular-origin source powers make the required pairings integrable; the singular covector is retained. None of these bounds is a uniform estimate in loop momentum, a late-time exchange or an arbitrary-spectrum UV/IR convergence proof.

## 7. Symbolic evaluation must not happen accidentally

Insertion initially stalled because `sympy.trace(Matrix(...))` invokes deep `Matrix.doit()`, which attempts every definite integral stored inside the matrix coefficients. A nominal tensor contraction was therefore triggering unspecified integrations before physical terms had been combined.

`inner`, `stf` and the linear Ricci trace in `unified_adjoint_source.py` now use explicit finite entry sums. Their algebra is unchanged. A regression guard raises if `Integral.doit` is invoked while constructing the generic source, Green covectors or stationary pairing. Intentional elementary/gradient fixture integrations are performed outside that guard.

This is important for maintaining the intended analytic organization, not just a timing optimization.

## 8. Checks and remaining calculation

```bash
.venv/bin/python -B solution/check_regular_field_green.py
.venv/bin/python -B solution/check_quadrupole_green_source.py
```

The exact certificates cover the Cartesian operator polynomial, spectral evaluations, differential equations, causal jumps, origin limits, physical scalar source/initial cancellation, complete weak duality, observer/emission traces and evaluated clock primitives. Independent finite-frequency quadrature validates the coupled/clock integrals and parameter-derivative resolvent integrals, including exact and near-zero momenta; it is labeled numerical validation, not the proof of the identities.

The source-insertion check verifies all three background and first covectors, their time arities, the actual Einstein single–pair functional and its nonzero initial row. At zero probe momentum its nonzero origin covector reproduces the independent leading physical gradient coefficient. This last check is a gradient coefficient, not a new unrestricted physical comparison.

Both new checks and seventeen prior targeted regressions passed, including the symbolic arbitrary-triple physical four-gradient comparison, exact angular source, dilation identities, stationary assembly, first/unified/propagation adjoints, Cartesian Green, observable and report identities, shell-source engine comparison, conformal frame transport and Bessel translation. Repository-wide Ruff lint/format and type checks, runner shell syntax and whitespace checks passed. The full long-running suite was not rerun.

Next combine the remaining three-time phase/source terms before radial expansion, looking especially for constraint-correlated cancellations like the scalar source/initial example above. The exact pure-parameter resolvent identity `partial_a G=-G*(partial_a L)*G` explains additional collapses for dilation/operator-variation weights; it does not apply to arbitrary time-dependent source weights. Evaluation of the general surviving two-/three-time kernels, an independent unrestricted physical comparison, and the complete fourth-order spectrum remain open.
