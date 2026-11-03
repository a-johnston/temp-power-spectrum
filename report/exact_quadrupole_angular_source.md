# Exact mixed-time angular closure of the quadrupole first source

## Status and scope

The evaluated gradient-adjoint work was committed as **`611f92a`**. This continuation retains the exact finite-wavelength phases rather than expanding them.

The complete **measurement part of the first quadrupole field-adjoint source** now has no remaining sky integrals or independently solved first-ray responses. Its mixed-time kernels are finite combinations of eight universal entire radial functions, with the physical emission, observer, proper-clock and random-mean normalization terms retained. The remaining field-operator variation is already sky contracted; its background angular tensor is supplied by the same closure.

A separate exact rotation Ward identity reduces the first frame covector to a commutator with the already-known linear quadrupole. The frame contribution is retained once, not deleted or double counted.

This is an exact source reduction at finite wavelength, **not** an evaluation of the remaining time/Green integrals or the cubic spectrum. It does not establish full nonlinear cancellation, arbitrary-spectrum convergence, or a complete all-multipole cubic map. The main report PDF is unchanged.

Implementation:

- `solution/entire_sphere.py`
- `solution/quadrupole_angular_source.py`

Checks:

- `solution/check_entire_sphere.py`
- `solution/check_quadrupole_angular_source.py`

## 1. An entire Cartesian phase-moment basis

For a nonnegative integer `r`, define

\[
 F_r(z)=\frac{{}_0F_1(;r+3/2;-z/4)}{(2r+1)!!}
 =\sum_{m\ge0}\frac{(-1)^mz^m}{2^m m!(2m+2r+1)!!}.
\tag{1}
\]

For real nonzero `x`, this equals `j_r(x)/x^r` at `z=x^2`, but (1), not a quotient, defines the function. In particular

\[
 F_r(0)=1/(2r+1)!!,\qquad F_r'(z)=-F_{r+1}(z)/2,
\]
\[
 F_r(z)=(2r+3)F_{r+1}(z)-zF_{r+2}(z).
\tag{2}
\]

There is no division by the phase length, no square-root branch and no separate zero-phase numerical prescription.

For a Cartesian monomial with powers `m=(m_1,m_2,m_3)` and degree `d=sum m_i`,

\[
\boxed{
 \langle n_1^{m_1}n_2^{m_2}n_3^{m_3}e^{iK\cdot n}\rangle
 =i^d\!\sum_{0\le r_i\le m_i/2}
 (-1)^rF_{d-r}(K\cdot K)
 \prod_{i=1}^3
 \frac{m_i!\,K_i^{m_i-2r_i}}{2^{r_i}r_i!(m_i-2r_i)!},
 \quad r=\sum_i r_i.}
\tag{3}
\]

This is the explicit evaluation of `i^(-d) partial_K^m j_0(|K|)`. It is valid for all phase vectors by entire continuation, not only small phases. At `K=0` it reduces to the standard exact sphere moments.

The implementation first reduces polynomials modulo `n.n=1`. Different polynomial representatives can give different-looking radial expressions; (2) reconciles them. Such representational differences are not failures of the angular identity.

For a trace-free symmetric quadrupole covector `C`, write

\[
 w=\tfrac{15}{2}n^TCn,\quad N_n=nn^T-I/3,\quad
 \mathcal U[P;K]=\langle Pe^{iK\cdot n}\rangle,
\]
\[
 \mathcal A_C[P;K]=\langle wPe^{iK\cdot n}\rangle.
\]

Thus `average(w)=0` and `A_C[N;0]=C`. Formula (3) evaluates both operators componentwise. The first-source polynomials below have weighted degree at most seven after unit-sphere reduction, so **only the families `F_0,...,F_7` are needed**. This is not eight fixed angular moments: their arguments vary with the two times and the momenta. The bound refers to the weak source with the probe derivative `j'` retained; optional time differentiation can introduce higher radial indices, which remain entire.

## 2. Eliminate the first-ray unknowns before contracting the sky

One first physical direction is

\[
 (\zeta,\nu,H)=(z,a(u),h(u))e^{ip\cdot X},\qquad G(u)=2zI+h(u),
\]

with `z` time independent. The independent state probe has momentum `q` and `delta zeta=0`. All field amplitudes here are **before** ray pullback, and `h'` is a partial field-time derivative.

Let `P_n=I-nn^T`, `phi_p(u)=exp(ip.n(o-u))`, and define phase-stripped vectors

\[
 v_*(u)=\tfrac12P_nh'(u)n
 +\tfrac i2[pG_{nn}(u)-(p\cdot n)G(u)n],
\]
\[
 x_*(t,u)=(u-t)v_*(u)-G(u)n/2,
\]
\[
 r_*(u)=ip\,h'_{nn}(u)/2,\qquad
 k_*(t,u)=-P_nh'(u)n-(t-u)r_*(u).
\tag{4}
\]

The complete first-ray and propagation-adjoint responses are then

\[
 X(t)=\int_t^o\phi_p(u)x_*(t,u)du,
 \qquad v(t)=-\int_t^o\phi_p(u)v_*(u)du,
\]
\[
 R(t)=\int_e^t\phi_p(u)r_*(u)du,
 \qquad K(t)=\int_e^t\phi_p(u)k_*(t,u)du.
\tag{5}
\]

These satisfy the original first equations and `X(o)=v(o)=0`, `R(e)=K(e)=0` by the fundamental theorem of calculus. The upper-interval minus sign in `v` is essential. These are fixed-coordinate-observer responses: **do not add the clock displacement `D_a n` to `X`**.

The one mixed phase entering the nonseparable ray terms is

\[
\boxed{
 K_{tu}=q(o-t)+p(o-u)=(p+q)(o-t)+p(t-u).}
\tag{6}
\]

There is no equal-time or monochromatic restriction in (4)–(6).

## 3. Three finite mixed-time ray polynomials

With `sym` followed by STF projection as in the preceding source note, define

\[
 B_1^>(t,u)=\tfrac i2(q\cdot x_*)N_n
 +\operatorname{TF}\operatorname{sym}(v_*n^T),
\]
\[
 B_1^<(t,u)=-\tfrac12\operatorname{TF}\operatorname{sym}((P_nk_*)n^T),
\]
\[
 B_0^<(t,u)=-\tfrac12\operatorname{TF}\operatorname{sym}(r_*n^T)
 -\tfrac i2(q\cdot k_*)N_n
 +\tfrac i2(q\cdot n)\operatorname{TF}\operatorname{sym}(k_*n^T).
\tag{7}
\]

Their weighted unit-sphere-reduced ranks are respectively `7,7,6`. They reproduce exactly the old `E_1,E_0` coefficients after (5) is inserted. No partial integration in the independent state probe is performed, so its endpoint and origin prescriptions are unchanged.

## 4. Exact normalization: retain both separated phase products

Set `K_p=p(o-u)`, `K_q=q(o-t)`. The normalization contribution to the coefficient of `j'(t)` has the time density

\[
\boxed{\begin{aligned}
 \mathcal N(t,u)={}&-\tfrac14\mathcal A_C[h'_{nn}(u)N_n;K_p+K_q]\\
 &+\tfrac14\mathcal U[h'_{nn}(u);K_p]\,
                 \mathcal A_C[N_n;K_q]\\
 &+\tfrac14\mathcal A_C[h'_{nn}(u);K_p]\,
                 \mathcal U[N_n;K_q].
\end{aligned}}
\tag{8}
\]

Both product terms are required. To see their signs directly, use

\[
 U_a(n)=\tfrac12\int_e^o\phi_p(u)h'_{nn}(u)du,\quad
 \Theta_a=-U_a+\langle U_a\rangle,\quad J_a=-\langle wU_a\rangle,
\]

and contract `(w Theta_a-J_a) N_n phi_q(t)/2`. This gives (8) by linearity. Thus not every term can be represented by a single combined phase vector.

The product terms are generically nonzero at finite wavelength. The leading-gradient comparison alone cannot test their retention: they start above its leading first-adjoint coefficient. The new tests include an exact finite-phase identity at six symbolic momentum components and numerical negative controls for each product.

No dipole or higher-multipole projection of `Theta_a` has been made inside this normalization term.

## 5. Complete angularly contracted measurement source

The first observer-clock shift is still

\[
 D_a=-[\iota^2a(\iota)+\int_\iota^o u a(u)du]/o.
\]

The already sky-contracted common-clock coefficient is

\[
 A_a=\partial_o\langle wU_a\rangle
 =\tfrac12 C:h'(o)
 +\tfrac12\int_e^o\mathcal A_C[i(p\cdot n)h'_{nn}(u);p(o-u)]du.
\tag{9}
\]

The remaining bulk matrices are

\[
 \overline M_H(t)=\int_e^t\mathcal A_C[B_0^<(t,u);K_{tu}]du,
\]
\[
\boxed{\begin{aligned}
 \overline M_{H'}(t)={}&\int_t^o\mathcal A_C[B_1^>(t,u);K_{tu}]du
 +\int_e^t\mathcal A_C[B_1^<(t,u);K_{tu}]du\\
 &+\int_e^o\mathcal N(t,u)du
 +\tfrac{iD_a}{2}\mathcal A_C[(q\cdot n)N_n;q(o-t)].
\end{aligned}}
\tag{10}
\]

For `K_e=(p+q)(o-e)`, the endpoint and clock coefficients are

\[
 \overline M_{\nu,e}=\tfrac e2\mathcal A_C[h'_{nn}(e);K_e],\quad
 \overline M_{H',e}=\tfrac{ea(e)}2\mathcal A_C[N_n;K_e],\quad
 \overline M_{H',o}=D_aC/2,
\]
\[
 M_{\nu,\mathrm{clock}}(t)=-A_at/o,\qquad
 M_{\nu,\mathrm{initial}}=-A_a\iota^2/o.
\tag{11}
\]

Equations (3)–(11) specify the **entire measurement contribution** to the first field-adjoint source without angular integrals, independent first rays, a gradient expansion, or momentum divisions. The definite time integrals have explicit limits and known first-field integrands; they are not claimed to be evaluated.

The background angular tensor is likewise just

\[
 C_P(t)=\mathcal A_C[N_n;P(o-t)].
\]

This is the angular tensor **entering** the background Green transpose, not the background field adjoint itself. The Green reconstruction and regular-origin covector from the earlier note must still be applied.

### API and normalization contract

`QuadrupoleAngularSource.assemble(z,h(u),a(u),t,u)` returns `ContractedMeasurement`. Its coefficients are already normalized-sky contracted. The method `.with_field(field_coefficients(...))` combines them with the background-dual field variation using an **identity** sky operation.

Do not pass these coefficients through a full `dOmega` integral again. In particular, the common clock and initial-clock terms are not second sky integrals. Homogeneous probes `q=0`, opposite momenta and zero/near-zero `K_tu` require no special momentum projector.

## 6. The frame adjoint is a commutator with the known linear quadrupole

For any first sky map, define its quadrupole tensor

\[
 S_a=\tfrac{15}{2}\langle\Theta_aN_n\rangle
 =-\tfrac{15}{4}\int_e^o\mathcal U[h'_{nn}(u)N_n;p(o-u)]du.
\]

The exact rotation Ward identity is

\[
\boxed{
 \mathcal K_a=\langle\operatorname{skew}(w\nabla_{S^2}\Theta_an^T)\rangle
 =[S_a,C],\qquad \Lambda_{R,a}=-\mathcal K_a.}
\tag{12}
\]

Proof: infinitesimal rotations preserve sphere area, so the rotation generator can be moved from `Theta_a` to the quadrupole weight by integration by parts. The resulting trace pairing is the commutator in (12). The code also checks this at generic momentum and generic STF matrices using the exact Cartesian phase moments.

This reduces a sky-gradient frame functional to the **already-known linear quadrupole**. It uses only unweighted moments through rank four. For a physical scalar first mode with the report's transfer convention,

\[
 S_a=-\tfrac{15}{2}T_2(|p|)Q_{\hat p}\quad(p\ne0).
\]

The same linear tensor also gives `J_a=C:S_a` and the common-clock coefficient `A_a=-C:partial_o S_a`. Thus the frame covector and this clock scalar can both reuse the known linear quadrupole, rather than require separate angular calculations.

The entire matrix representation, not division by `p^2`, should be used at zero momentum. Rotation preserves multipole degree, so only the quadrupole is needed **in this frame sector**. This does not justify projecting the first map inside the normalization or propagation terms.

The ray-adjoint observer vector `K(o)` from (5) differs from `grad_S Theta_a` only by a radial vector. Their antisymmetric pairings are identical, as checked pointwise before angular integration. Thus (12) is also an exact reduction of that observer boundary.

The cubic frame pairing remains `sum_cyclic integral_iota^o K_a:Omega_bc dt`. It is not zero in general. If the direct measurement is evaluated with `stationary_cubic.py`, its frame-exponent boundary already represents this contribution: **do not add (12) again**. The first field source still has no extra frame row, since the rotation driver has no linear variation.

## 7. A homogeneous-tensor cancellation for the entire contracted response

At fixed loop radius `r`, rotational covariance of the exact second Einstein residual gives

\[
 C_{H,bc}(t;\alpha,-\alpha)=f_H(t,r)Q_\alpha,\qquad
 C_{N,bc}(t;\alpha,-\alpha)=f_N(t,r),\qquad
 c_{bc}=5r^2/9.
\]

Here `Q_alpha=TF(alpha alpha^T)`, with `alpha.alpha=r^2`. The check verifies the matrix and scalar source identities at symbolic loop radius and symbolic unit-sphere orientation, with arbitrary first time profiles and no time or gradient expansion. On the prescribed regular branch, the homogeneous spatial Green operator is rotationally invariant, so the forced homogeneous STF response also has the form `H_bc(t)=g_H(t,r) Q_alpha`.

At cubic order every occurrence of this second response is **linear** in `H_bc` (or its time derivatives or traces at physical endpoints). For the partition with first direction `p` and probe momentum zero, the complete first-source/Green functional is independent of the internal orientation `alpha`. Consequently

\[
 \left\langle\mathscr L_{p,0}[H_{bc}]\right\rangle_{\hat\alpha}
 =\mathscr L_{p,0}[\langle H_{bc}\rangle_{\hat\alpha}]=0,
 \qquad \langle Q_\alpha\rangle_{\hat\alpha}=0.
\]

Thus the forced homogeneous tensor can be removed from the **whole isotropically contracted cubic response**, not only its bulk field equations. This includes the measurement and physical boundary contributions. The argument is linearity and rotational covariance, not a quadrupole truncation; any fixed sky weight has this cancellation for this one forced block. It assumes the prescribed regular branch, with no independent primordial homogeneous tensor data.

This operation is allowed **only after the isotropic loop contraction**: the uncontracted arbitrary-triple map still contains the tensor. The homogeneous lapse and its nonzero scalar initial row are not removed. Neither are the direct `111` terms or the nonhomogeneous pair partitions. This is a redundant-tensor reduction, not cancellation of the nonlinear spectrum. In particular, it is not a claim that a squared second-order response in `C^(22)` vanishes.

## 8. Verification and limits

```bash
.venv/bin/python -B solution/check_entire_sphere.py
.venv/bin/python -B solution/check_quadrupole_angular_source.py
```

Exact checks cover:

- all 120 Cartesian moments through rank seven, against derivatives of the entire characteristic function;
- exact zero-phase values and raw unit-sphere trace identities;
- generic-momentum first-ray/Volterra signs, observer/emission data and source polynomial identities;
- the full finite-phase random-mean normalization identity and frame Ward identity;
- every previously evaluated leading-gradient measurement coefficient at symbolic endpoint times;
- the exact constant-seed measurement dilation identity `M_H'=-(q.partial_q) C_q(t)/2`;
- exact rotational covariance of the homogeneous Einstein pair residual, with a nonzero lapse/initial source retained;
- a complete homogeneous off-shell clock/emission/observer action and its no-double-sky adapter.

Independent finite-phase sphere quadrature additionally validates the closed kernels, including exact/near-zero combined phases and homogeneous probes. These numerical checks are labeled separately: the analytical source identities do not follow merely from quadrature agreement. The homogeneous fixture is not an Einstein–radiation solution.

Both new scripts and twelve existing targeted regressions passed: the quadrupole adjoint, symbolic physical four-gradient comparison, stationary cubic, first source, unified adjoint, propagation adjoint, Cartesian Green, observable expansion, report identities, shell sources with the independent engine, Fermi–Walker conformal transport and Bessel translation checks. Repository-wide Ruff lint/format and type checking, shell syntax and whitespace checks passed. The full long-running suite was not rerun.

The exact sky closure is now available for insertion into the regular field-adjoint Green reconstruction. The next open calculation is the remaining mixed-time/Green pairing in the stationary cubic functional, looking for physical sum/total-derivative cancellations before evaluating it. The direct cubic measurement's multi-phase sky contractions also remain to be assembled with these blocks; the universal moment operator applies to them, but their rank is not restricted to seven. Neither the eight function families nor (12) imply automatic closure on the earlier nine equal-time moments or vanishing of the nonlinear response.

No late-time integral exchange or arbitrary-spectrum UV/IR convergence is asserted.
