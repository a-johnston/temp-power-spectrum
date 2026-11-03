# Derivations: adjoint cubic response, orthonormal map, and the curvature form of the next gradient order

Supporting derivations for Prop. adjoint, the orthonormal-direction map (Prop. orthonormal) and its
cubic product identity (eqs. cubic-orthonormal, cubic-composed), the curvature form of the map through
four gradients (eq. curvature-full-map), and the next-gradient-order spectra (eqs. three-tensor,
general-next). The checks are

```sh
.venv/bin/python -B solution/check_adjoint_exact.py
.venv/bin/python -B solution/check_cubic_curvature_hand.py   # ~6 s; imports no generator or cached kernel
```

These derivations come from two external correctness reviews; the reviews themselves are in the git history.

## 1. Adjoint form of the top-order cubic response

Write `k=k_0`, `c=k/sqrt(3)`, `e=bar t_star`, `o=bar t_o`, and

\[
w(t)=w_\ell(k(o-t)).
\]

Let `u_1,u_2` be the report's growing and decaying lapse modes and define

\[
G_N(t,s)=-\frac{u_2(ct)u_1(cs)-u_1(ct)u_2(cs)}{3c},\qquad t\ge s.
\]

It obeys

\[
L_NG_N=0\quad(t>s),\qquad G_N(s,s)=0,\qquad
\partial_tG_N(s,s)=-\frac13.
\]

The minus sign is important: the particular solution is for `L_N nu = -c_N`.

### 1.1 Why the third-order shifts disappear

The background conformal shear is identically zero. Consequently the redshift integrand has no zeroth-order term:

\[
F=\epsilon F_1+\frac{\epsilon^2}{2}F_2+\frac{\epsilon^3}{6}F_3+\cdots.
\]

Expansion of the moving lower endpoint gives

\[
I_3=\int_e^{t_o}F_3\,dt
-3F_1d_2-3F_2d_1-3F_1'd_1^2.
\]

There is no `d_3` term: it would multiply `F_0`. Equivalently, moving a boundary by an order-three displacement changes an order-one integrand first at order four. Likewise the background normalized map is zero and independent of observation time, so

\[
\Theta_3=\Theta_{c,3}+3D_2\partial_o\Theta_{c,1}
+3D_1\partial_o\Theta_{c,2}+3D_1^2\partial_o^2\Theta_{c,1}
\]

contains no `D_3`. The same reasoning excludes `X_3,e_3`: they multiply zero background shear or its derivatives.

This proves the intended claim, but not the stronger and false statement that endpoint and clock corrections disappear from the cubic observable. The terms displayed above remain, and Section 3 gives their nonzero contributions to the requested spectrum coefficient.

The genuinely new cubic field enters only as

\[
-\frac12\int_e^o A'_{3,nn}(t)e^{ik\cdot n(o-t)}dt.
\]

Centering removes only its sky monopole. Lower-order normalization products cannot involve a cubic field until a higher amplitude order. After the internal-direction average only its scalar metric profile survives.

### 1.2 Homogeneous term and contact term

With

\[
\mathcal H[\nu]=6\nu+\frac{18}{k^2}\left(\frac{\nu'}t+\frac\nu{t^2}\right),
\]

one has exactly

\[
\mathcal H[u_1]=6-h_{S1}.
\]

The constant six disappears on differentiation in the transfer integral. On the shell the regular-lapse coefficients are

\[
c_1=-\frac{k^2}{9},\qquad c_3=-\frac{14k^2}{9}.
\]

Hence the homogeneous cubic transfer is `14 T_ell`, as stated.

When differentiating `H[nu_part]`, the `nu_part''` term sees the diagonal jump of the Green function. Its contribution to the half-derivative transfer is

\[
\frac12\frac{18}{k^2t}\left(-\frac13\right)c_N(t)
=-\frac{3}{k^2t}c_N(t).
\]

This is the local term of the contact form of the adjoint. Omitting it from the Hamiltonian-source version would be an error.

For `s<e`, the partial transfer functions have fixed lower limit `e>0`. They are finite constants as functions of `s`; `u_1(cs)=O(s^2)` and `u_2(cs)=O(s^{-1})`. The cubic lapse source is `c_N(s)=O(s^2)`, because `nu,H=O(s^2)` and all nonlinear terms in the lapse equation start at that order. The Fubini integrand near `s=0` is therefore `O(s)`, hence absolutely integrable. The remaining domain is compact and regular. This supplies the missing endpoint estimate.

### 1.3 Source compatibility removes the explicit contact term

For solved lower-order fields, the scalar Hamiltonian and momentum constraints imply a particularly useful source identity:

\[
\boxed{c_{\rm Ham}'-\frac4t c_N=2iK\cdot c_M.}
\]

To derive it, differentiate `h_S=H[nu]+3 c_Ham/(2K^2)` and use the lapse equation. This gives

\[
h_S'=6\left(\nu'-\frac\nu t\right)
+\frac{3}{2K^2}\left(c_{\rm Ham}'-\frac4t c_N\right).
\]

The longitudinal momentum constraint gives independently

\[
h_S'=6\left(\nu'-\frac\nu t\right)+\frac{3i}{K^2}K\cdot c_M.
\]

Comparing proves the identity. It is a compatibility relation for physical lower-order solutions, not an identity for arbitrary off-shell placeholders in an expanded source file.

Define the contact-free kernel

\[
\mathcal W^0_\ell(s)=-\frac1{3c}\left[
 u_1(cs)\mathcal T_\ell[\mathcal H u_2](s)
-u_2(cs)\mathcal T_\ell[\mathcal H u_1](s)\right].
\]

The contact-form kernel is `W = W^0 - 3 w(s)/(k^2 s)` on `s>e`. Its contact term and Hamiltonian term combine **pointwise**:

\[
-\frac{3}{k^2t}w\,c_N+\frac{3}{4k^2}w\,c_{\rm Ham}'
=\frac{3i}{2k^2}w\,K\cdot c_M.
\]

Thus a cleaner adjoint transfer is

\[
\boxed{
\mathcal A_\ell^{\rm top}
=14T_\ell+\int_0^o\mathcal W^0_\ell(s)\bar c_N(s)ds
+\frac{3}{2k^2}\int_e^o w(t)\,\overline{iK\cdot c_M}(t)dt.
}
\]

The bar denotes the same real shell average as in the report. In particular, take the real part of `iK dot c_M` as a whole; `c_M` itself is imaginary in the usual Fourier representation.

The unnormalized spectrum identity is

\[
\boxed{
\left.\frac13\operatorname{Re}C_\ell^{(13)}\right|_{\rm top}
=4\pi\sigma^4 T_\ell\,\mathcal A_\ell^{\rm top}.
}
\]

It remains meaningful at `T_ell=0`, whereas the displayed ratio in a ratio form does not. The new form uses the same partial transfer functions but no differentiated Hamiltonian source and no separate contact term. It does not assert that all remaining finite-wavelength cancellation has disappeared.

`check_adjoint_exact.py` verifies the Green identities, homogeneous factor, contact coefficient, arbitrary-function operator identity, source combination, and moving-endpoint/clock combinatorics exactly. The Fubini argument is the analytical estimate above, not a numerical test.

## 2. Orthonormal-direction F_2 identity

Set `M=exp(H)` and `v=exp(zeta) exp(H/2)e`. The exact driver is

\[
e^T\sigma e=v^T\mathcal Bv,\qquad
\mathcal B=\frac12e^{-H/2}(e^H)'e^{-H/2}.
\]

The differential of the exponential gives

\[
\mathcal B=\frac12\int_{-1/2}^{1/2}e^{s\,\mathrm{ad}_H}H'\,ds
=\frac12H'+\frac1{48}[H,[H,H']]+\cdots.
\]

Odd commutator powers integrate to zero. Thus, in the factorial amplitude convention, `B_1=H_1'/2`, `B_2=H_2'/2`; no quadratic matrix product is present. At first order along the background ray,

\[
v_1=e_1-\left(\zeta_1\mathbf1+\frac12H_1\right)n.
\]

Expanding the exact driver at the displaced point gives

\[
F_2=\mathcal B_{2,nn}+2X_1\cdot\partial\mathcal B_{1,nn}
-4\mathcal B_1(n,v_1),
\]

hence

\[
\boxed{F_2=\frac12H_{2,nn}'+X_1\cdot\partial H_{1,nn}'-2H_1'(n,v_1).}
\]

For a direct check in the original variables,

\[
\sigma_2=\frac12H_2'+\frac12(H_1^2)'+2\zeta_1H_1'.
\]

Substituting `e_1=v_1+(zeta_1+H_1/2)n` into `-4 sigma_1(n,e_1)` cancels `2 zeta_1 H_1'` exactly. The remaining matrix product cancels because

\[
n^T(H_1H_1'+H_1'H_1)n=2n^TH_1'H_1n
\]

for symmetric `H_1,H_1'`, even when the matrices do not commute. No shell assumption is needed for `F_2`.

The first-order transverse direction equation in the report is also correct. Transversality follows from `P_perp` in its first term and

\[
n\cdot\left[\nabla H_{1,nn}-(n\cdot\partial)(H_1n)\right]=0.
\]

On the shell, `H_1(t,X)=f(t)B(X)` additionally kills the cubic double commutator. Off shell that commutator must be kept. The existing two orthonormal checks were rerun and passed.

## 3. Independent curvature route to the k_0^6 quadrupole

This route is independent of both `observable_fast.py` and the earlier quadratic hand-map script. It starts with the exact equations and expands in spatial gradients, keeping the amplitude dependence until the final kernel extraction.

### 3.1 Local field data, exact in the primordial curvature

Let

\[
\gamma_0=e^{2\zeta}\delta,
\qquad
\mathfrak a=-\frac{R[\gamma_0]}{36},
\qquad
\mathsf A_{ab}=-\frac13R_{\langle ab\rangle}[\gamma_0].
\]

All components and covariant derivatives below are in the seed metric, with orthonormal components at the observer. These tensors are not the Ricci curvature of the physical spacetime at observation. Write

\[
\nu=t^2\mathfrak a+t^4\mathfrak b+O(\partial^6),\qquad
H=t^2\mathsf A+t^4\mathsf B+O(\partial^6).
\]

Here `a,A` have two spatial gradients and `b,B` have four. The lapse equation gives

\[
30\mathfrak b=\nabla^2\mathfrak a-6\mathfrak a^2-\operatorname{tr}\mathsf A^2.
\]

The mixed spatial equation gives

\[
10\mathsf B=-\operatorname{TF}\left[
\delta R^a{}_b[\delta\gamma=\gamma_0\mathsf A]
-2\nabla^a\nabla_b\mathfrak a+4\mathfrak a\mathsf A^a{}_b\right].
\]

The variation of mixed Ricci, combined with the three-dimensional curvature identity, reduces this to

\[
\boxed{\mathsf B=\frac1{20}\nabla^2\mathsf A
+\frac15\mathfrak a\mathsf A+\frac35(\mathsf A^2)_{\rm TF}.}
\]

There is a further simplification because the seed is conformally flat. The vanishing Cotton tensor and contracted Bianchi identity give

\[
\nabla^b\mathsf A_{ab}=2\nabla_a\mathfrak a,
\qquad
\boxed{\nabla^2\mathsf A
=3(\nabla\nabla\mathfrak a)_{\rm TF}
-9(\mathsf A^2)_{\rm TF}-18\mathfrak a\mathsf A.}
\]

Consequently

\[
\boxed{\mathsf B=\frac3{20}(\nabla\nabla\mathfrak a)_{\rm TF}
+\frac3{20}(\mathsf A^2)_{\rm TF}
-\frac7{10}\mathfrak a\mathsf A.}
\]

One derivation of the first formula uses the Ricci variation for `delta gamma = -Ric/3 + R gamma/9`; the Hessian-of-`R` terms cancel against the lapse Hessian in the spatial equation. The script does not assume that simplification as its only field check: it independently constructs `delta Gamma` and `delta Ric` using covariant derivatives and verifies the unsimplified spatial equation through cubic amplitude. It also checks the Cotton identity directly. No inverse momentum projectors are needed.

### 3.2 Hand ray and endpoint expansion

Put

\[
\Delta=o^2-e^2,\qquad E=o^4-e^4,
\qquad L=\int_e^o t(o-t)^2dt
=\frac{o^4}{12}-\frac{o^2e^2}{2}+\frac{2oe^3}{3}-\frac{e^4}{4}.
\]

Use the radial geodesic of `gamma_0` and parallel-transport the seed frame along it. Covariant Taylor expansion supplies the spatial dependence of `A` along that geodesic. The time-dependent metric changes its position only at gradient order two; multiplying that change by `nabla A=O(partial^3)` first affects the driver at order five. The change in direction must nevertheless be kept at order four:

\[
\delta v(t)=\frac{t^2-o^2}{2}P_\perp\mathsf A n+O(\partial^3).
\]

Since the local driver is `t A + 2t^3 B` through this order, its even-sky part is

\[
F_{\rm even}=t\mathsf A_{nn}
+\frac{t(o-t)^2}{2}(\nabla_n\nabla_n\mathsf A)_{nn}
+2t^3\mathsf B_{nn}
+t(o^2-t^2)\left[(\mathsf A^2)_{nn}-\mathsf A_{nn}^2\right]
+O(\partial^6).
\]

The odd-gradient terms are odd sky multipoles and do not contribute to the quadrupole. The physical emission and clock shifts at the required gradient order are

\[
\delta t_e=-e^3\mathfrak a+O(\partial^3),\qquad
\delta t_o=-\frac{o^4+3i^4}{4o}\mathfrak a+O(\partial^4).
\]

The exponent therefore receives `+e^4 a A_nn` from the emission boundary and `-(o^4+3i^4) a A_nn/4` from the observer clock. These signs are reversed in the linear-in-exponent part of the map.

After normalization, with `P_0` denoting sky centering, the result is

\[
\boxed{\begin{aligned}
\Theta_{\rm even}=\mathsf P_0\bigg[&
-\frac\Delta2\mathsf A_{nn}
-\frac L2(\nabla_n\nabla_n\mathsf A)_{nn}
-\frac E2\mathsf B_{nn}
-\frac{\Delta^2}{4}(\mathsf A^2)_{nn}
+\frac{3\Delta^2}{8}\mathsf A_{nn}^2\\
&+\left(\frac{o^4+3i^4}{4}-e^4\right)\mathfrak a\mathsf A_{nn}
\bigg]+O(\partial^6).
\end{aligned}}
\]

In particular, the coefficient `3 Delta^2/8` includes both the direction response and nonlinear map normalization. The latter has not been discarded. The Fermi–Walker commutator starts beyond the gradient order relevant here.

### 3.3 Three-tensor quadrupole

The angular projection identities are

\[
[\mathsf A_{nn}^2]_{\ell=2}=\frac47(\mathsf A^2)_{\langle ab\rangle}n^an^b,
\]

\[
[(\nabla_n\nabla_n\mathsf A)_{nn}]_{\ell=2}
=\frac17\left[11\nabla_{\langle a}\nabla_{b\rangle}\mathfrak a
-27(\mathsf A^2)_{\langle ab\rangle}-54\mathfrak a\mathsf A_{ab}\right]n^an^b.
\]

The second follows by tracing the symmetrized rank-four derivative, commuting derivatives, and using the two seed-curvature identities above. Inserting these formulas and the simplified `B` gives

\[
\boxed{\Theta_{\ell=2}=Q_{ab}n^an^b+O(\partial^6),\qquad
Q=-\frac\Delta2\mathsf A
+\beta_H(\nabla\nabla\mathfrak a)_{\rm TF}
+\beta_S(\mathsf A^2)_{\rm TF}
+\beta_a\mathfrak a\mathsf A,}
\]

where

\[
\begin{aligned}
\beta_H&=-\frac{11L}{14}-\frac{3E}{40}
=-\frac{o-e}{420}(59o^3+59o^2e-106oe^2+114e^3),\\
\beta_S&=\frac{27L}{14}-\frac{3E}{40}-\frac{\Delta^2}{28}
=\frac{o-e}{140}(7o^3+7o^2e-118oe^2+62e^3),\\
\beta_a&=\frac{27L}{7}+\frac{7E}{20}+\frac{o^4+3i^4}{4}-e^4
=\frac3{140}(43o^4-90o^2e^2+120oe^3-108e^4+35i^4).
\end{aligned}
\]

This is an all-amplitude local gradient formula, not a finite-wavelength closed solution. It is substantially simpler than component ray kernels at this gradient order. The new script checks every quadrupole-tensor component of this reduced expression against the preceding unreduced even-sky expression through cubic amplitude.

### 3.4 Exact kernel extraction and angular averages

The script sets `k_0=1` during the homogeneous gradient calculation and uses three distinct nilpotent labels for

\[
k=(1,0,0),\quad q=(\mu,\sqrt{1-\mu^2},0),\quad -q.
\]

The coefficient of the three labels is the report's symmetric cubic kernel, with no extra factorial. Homogeneity restores `k_0^6` at the end. All angular integrations use exact polynomial sphere moments and the exact `mu` integral on `[-1,1]`. The square-root variable is treated algebraically with its quadratic relation, not sampled numerically.

For the cubic coefficient both products are retained:

\[
K_1^{[2]}K_3^{[4]}+K_1^{[4]}K_3^{[2]}.
\]

Omitting the second is another way to get a wrong `k_0^6` answer even with a correct cubic fourth-gradient map.

Let `p_{r,d}=<P_2(khat dot n) K_r^[d]>_n`, with the internal-direction average also taken for the cubic kernel. The independent calculation gives

\[
p_{1,2}=\frac\Delta{45},\quad p_{1,4}=\frac{2\beta_H}{135},
\qquad p_{3,2}=\frac{4\Delta}{45},\quad
p_{3,4}=\frac{196\beta_H}{405}-\frac{2\beta_a}{45}.
\]

Thus the cubic contribution follows from the exact Wick/angular normalization

\[
\left[\frac13\operatorname{Re}C_2^{(13)}\right]_{k_0^6}
=4\pi\sigma^4k_0^6(p_{1,2}p_{3,4}+p_{1,4}p_{3,2}).
\]

### 3.5 Final coefficient and a cleaner spectral presentation

The resulting cubic coefficient is

\[
\boxed{\begin{aligned}
\frac13\operatorname{Re}C_2^{(13)}={}&\frac{16\pi}{2025}\sigma^4k_0^4\Delta^2\\
&-\frac{2\pi}{1913625}\sigma^4k_0^6\Delta
(2835i^4+9973o^4-25440o^2e^2+33920oe^3-21288e^4)
+O(k_0^8).
\end{aligned}}
\]

This agrees exactly with the report.

The same independent map reproduces the reported `C_2^(11)` and `C_2^(22)/4`. A more compact presentation of all their next-gradient coefficients is

\[
\boxed{\begin{aligned}
[C_2^{(11)}]_{k_0^6}
&=\frac{16\pi}{6075}\sigma^2k_0^6\Delta\,\beta_H,\\
[\tfrac14C_2^{(22)}]_{k_0^6}
&=\frac{32\pi}{18225}\sigma^4k_0^6\Delta\,(30\beta_H-\beta_a),\\
[\tfrac13\operatorname{Re}C_2^{(13)}]_{k_0^6}
&=\frac{8\pi}{18225}\sigma^4k_0^6\Delta\,(110\beta_H-9\beta_a),\\
[C_2^{(4)}]_{k_0^6}
&=\frac{8\pi}{18225}\sigma^4k_0^6\Delta\,(230\beta_H-13\beta_a).
\end{aligned}}
\]

The `beta_S` term drops out of this fourth-order spectrum only **after** the Gaussian average. At this amplitude order its interference is proportional to the order-four part of `<tr A^3>`. Write the leading conformal tensor in terms of the Gaussian trace-free Hessian and trace-free gradient square. The scalar field value is independent of both; the Hessian is independent of the gradient; the mean trace-free gradient square is zero. The order-four contractions therefore vanish. This is a genuine structural cancellation, not permission to delete that tensor from the map.

Expanding the last line recovers exactly

\[
-\frac{2\pi}{1913625}\sigma^4k_0^6\Delta
(4095i^4+18601o^4-48480o^2e^2+64640oe^3-38856e^4).
\]

The initial-time dependence and the boundary terms are real. In particular, the independent cubic calculation isolates

\[
\boxed{\delta_{\rm emission}[\tfrac13\operatorname{Re}C_2^{(13)}]_{k_0^6}
=\frac{8\pi}{2025}\sigma^4k_0^6\Delta e^4,}
\]

\[
\boxed{\delta_{\rm clock}[\tfrac13\operatorname{Re}C_2^{(13)}]_{k_0^6}
=-\frac{2\pi}{2025}\sigma^4k_0^6\Delta(o^4+3i^4).}
\]

These are nonzero despite the absence of `d_3,D_3`. The script has negative controls for omitting the spatial Hessian, fourth-gradient field, emission term, clock term, or `K_1^[4] K_3^[2]` ordering. It separately verifies the contractions that vanish.

## 4. Explicit cubic product identity

The exact orthonormal driver is

\[
\mathcal B=\frac12e^{-H/2}(e^H)'e^{-H/2}
=\frac12H'+\frac1{48}[H,[H,H']]+O(\epsilon^4).
\]

With factorial amplitude coefficients and the pulled-back direction `v = -n + epsilon v1 + epsilon² v2/2`, expansion gives

\[
\begin{aligned}
F_3={}&\tfrac12H'_{3,nn}+\tfrac18[H_1,[H_1,H_1']]_{nn}
+\tfrac32X_1\cdot\partial H'_{2,nn}-3H_2'(n,v_1)\\
&+\tfrac32X_2\cdot\partial H'_{1,nn}
+\tfrac32X_1^iX_1^j\partial_i\partial_jH'_{1,nn}
-6(X_1\cdot\partial)H_1'(n,v_1)
-3H_1'(n,v_2)+3H_1'(v_1,v_1).
\end{aligned}
\]

The time derivative of `H` is taken before pullback; spatial derivatives act only on fields. The local conformal and quadratic matrix products have disappeared. The double commutator is zero on the shell, but generally nonzero off it. Unit norm also fixes `n·v2 = |v1|²`; its radial part is not optional.

Let `Ur(o) = integral_e^o Fr(t;o) dt`. Clock composition of the moving-emission integral gives

\[
\begin{aligned}
\widehat I_3={}&U_3-3[F_1d_2+F_2d_1+F_1'd_1^2]_e\\
&+3D_1\partial_o\{U_2-2[F_1d_1]_e\}
+3D_2\partial_o U_1+3D_1^2\partial_o^2U_1.
\end{aligned}
\]

The endpoint prime is the derivative along the background ray. In particular, the mixed clock–emission product `-6 D1 partial_o(F1 d1)` survives. Center this exponent and normalize once; do not apply a second clock correction afterwards.

`solution/check_cubic_product_identity.py` proves the matrix identity for generic symmetric noncommuting matrices, checks the pulled-back coefficient rule, and checks the endpoint/clock composition on independent jets. Omission controls detect the commutator and six boundary/clock pieces. This closes the **algebraic** identity. It does not turn all lower-order propagation products into the top-order adjoint source, or evaluate their integrals.

## 5. Full map through four gradients from local curvature

Use the seed metric `gamma0 = exp(2 zeta) delta` and its orthonormal curvature tensors

\[
\mathfrak a=-R/36,\qquad \mathsf A=-R_{\rm TF}/3,\qquad
\mathsf B=\frac3{20}(\nabla\nabla\mathfrak a)_{\rm TF}
+\frac3{20}(\mathsf A^2)_{\rm TF}-\frac7{10}\mathfrak a\mathsf A.
\]

Set `Delta=o²-e²`, `E=o⁴-e⁴`, `M=integral_e^o t(o-t)dt`, and `L=integral_e^o t(o-t)²dt`. The unprojected map through four gradients is

\[
\begin{aligned}
\Theta=\mathsf P_0\bigg[&-\frac\Delta2\mathsf A_{nn}
-M(\nabla_n\mathsf A)_{nn}-\frac L2(\nabla_n\nabla_n\mathsf A)_{nn}
-\frac E2\mathsf B_{nn}\\
&-\frac{\Delta^2}4(\mathsf A^2)_{nn}
+\frac{3\Delta^2}8\mathsf A_{nn}^2
+\left(\frac{o^4+3i^4}4-e^4\right)\mathfrak a\mathsf A_{nn}\bigg]
+O(\partial^5).
\end{aligned}
\]

Section 3 already established the even part. The odd part is its first covariant Taylor term along the seed radial geodesic; it contains both the dipole and octupole. Terms from the time-dependent displacement first affect the driver at fifth gradient order. The direction response already contributes at fourth order and is retained above. Fermi–Walker rotation is beyond this gradient order.

`solution/curvature_kernels.py` extracts nilpotent-label coefficients for general

\[
k=(p,0,0),\qquad q=(b_1,b_2,0),\qquad -q.
\]

It imports the independent local-curvature algebra, **not** the field engine or ray generator. The leading, third-gradient and fourth-gradient terms are centered and assembled into the full map. The resulting polynomial term counts are 37 for `K1`, 478 for `K2(k,q)` and 318 for `K3(k,q,-q)`.

The optional comparison uses the **raw** saved legacy generator outputs, not their simplified SymPy cache. It groups terms by formal inverse-momentum powers and clears a common polynomial denominator over exact Gaussian rationals. Equality holds component by component at general momenta and symbolic emission, observation and initial times. Polynomial equality establishes the continuation across the removable denominator loci as well.

Construction takes about 3 s; construction, raw comparison and all spectral contractions take roughly 12–15 s on a laptop. This replaces the approximately 13-hour run **for the full map through gradient degree four in this configuration**, not for arbitrary wavelengths or higher gradient orders. The legacy generator run remains an independent derivation. Its saved pickle files are ignored local artifacts, so `--compare` requires them; `--spectrum` and kernel construction are cache-free. `--output NEW_PATH` exports exact polynomial kernels without overwriting an existing file.

## 6. General spectra at the next gradient order

Define `mj = integral k^j Pcal(k) d ln k`, and let `[ ]_[6]` denote total gradient degree six at fixed times. The three-tensor quadrupole of Section 3.3 has time coefficients

\[
\beta_H=-\frac{11L}{14}-\frac{3E}{40},\qquad
\beta_a=\frac{27L}{7}+\frac{7E}{20}+\frac{o^4+3i^4}4-e^4.
\]

The general-momentum contractions give

\[
\begin{aligned}
[C_2^{(11)}]_{[6]}&=\frac{16\pi}{6075}\Delta\beta_H m_6,\\
[\tfrac14C_2^{(22)}]_{[6]}&=\frac{32\pi}{18225}\Delta
\{12\beta_Hm_0m_6+(18\beta_H-\beta_a)m_2m_4\},\\
[\tfrac13\operatorname{Re}C_2^{(13)}]_{[6]}&=\frac{8\pi}{18225}\Delta
\{60\beta_Hm_0m_6+(50\beta_H-9\beta_a)m_2m_4\}.
\end{aligned}
\]

Consequently,

\[
\boxed{
[C_2^{(4)}]_{[6]}=
18m_0[C_2^{(11)}]_{[6]}
+\frac{8\pi}{18225}\Delta(122\beta_H-13\beta_a)m_2m_4.
}
\]

**After removing observer dilation, only one spectral moment product remains.** The entire `m0 m6` term is precisely the degree-six dilation operator `D²/2 = 18`. The trace-free `A²` tensor drops out after the Gaussian contractions, not from the map. Setting `mj = sigma² k0^j` recovers the report's shell coefficients and the combination `230 betaH - 13 betaa`.

For reproducibility, all sphere averages use

\[
\langle n_x^{2a}n_y^{2b}n_z^{2c}\rangle
=\frac{(2a-1)!!(2b-1)!!(2c-1)!!}{(2a+2b+2c+1)!!},
\]

with odd moments zero, followed by `1/2 integral_-1^1 dmu` and exact radial moment substitution. Every degree-two map is pure quadrupole, so in the degree-six cross products its partner needs no explicit harmonic projection. Both `K1[2] K3[4]` and `K1[4] K3[2]` are included.

The same construction independently reproduces the leading octupole. With `B3=16 pi beta3²/1225`,

\[
C_3^{(11)}=B_3m_6,\qquad
\tfrac14C_3^{(22)}=B_3(9m_0m_6+28m_2m_4),\qquad
\tfrac13\operatorname{Re}C_3^{(13)}=9B_3m_0m_6.
\]

The dipole is subtracted from the odd map before these contractions, and complex conjugation is retained. On the shell this yields `37 + 9 = 46`, equivalently the reported separate ratios `148` and `27` before the `1/4` and `1/3` factors.

Compact spectral support justifies the gradient expansion under the integrals. Finiteness of the moments defines these coefficients more generally but does not by itself bound the remainder for arbitrary noncompact spectra.

## 7. Reproduction

```sh
.venv/bin/python -B solution/check_cubic_product_identity.py
.venv/bin/python -B solution/check_cubic_curvature_hand.py
.venv/bin/python -B solution/curvature_kernels.py --spectrum
# Optional independent comparison, when the legacy raw kernels are present:
.venv/bin/python -B solution/curvature_kernels.py --compare --spectrum
```

These checks use exact arithmetic, not fitting, numerical evolution or quadrature. Runtimes above are performance observations, not evidence for a mathematical identity. The raw comparison is stronger than agreement of final spectrum coefficients, but is still restricted to the stated momentum configurations and gradient/amplitude orders.
