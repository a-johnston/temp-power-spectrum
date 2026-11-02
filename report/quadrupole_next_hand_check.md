# Independent hand-map check of the next quadrupole coefficient

## Scope and result

The new saved calculation available during this check produces `C_2^(11)` and `C_2^(22)/4` through `O(k_0^6)`. **Both are independently reproduced here.** More strongly, every component of the `O(k_0^4)` quadratic **quadrupole map kernel** agrees with the saved output, before the relative-momentum integral.

This does **not** compute the `O(k_0^6)` part of `Re C_2^(13)/3`. Thus it must not be described as a completed check of the full fourth-order-in-amplitude correction. The corresponding cubic saved output was not present when this check was constructed.

The derivation uses the report's hand-expanded map formulas, a direct local field expansion through `t^4`, and exact sphere moments. It imports **no generator, shared field engine, or existing source cache**. The optional comparison reads only the already saved generator output, after the independent calculation has finished.

```sh
.venv/bin/python -B solution/check_quadrupole_next_hand.py
.venv/bin/python -B solution/check_quadrupole_next_hand.py --compare
```

The second command compares against the trusted local file `solution/results/superhorizon_kernels_2.pkl`; it does not run or import `observable_fast.py` or `superhorizon_next.py`.

## 1. Final coefficients

Use background times

\[
o=\bar t_o,\qquad e=\bar t_\star,\qquad i=\bar t_I,\qquad
\Delta=o^2-e^2,
\]

and hold the primordial variance `sigma^2` and physical measurement labels fixed as `k_0` varies. Then

\[
\boxed{\begin{aligned}
\frac14C_2^{(22)}={}&\frac{68\pi}{6075}\sigma^4 k_0^4\Delta^2\\
&-\frac{8\pi}{637875}\sigma^4 k_0^6\Delta
\left(105i^4+719o^4-1920o^2e^2+2560oe^3-1464e^4\right)
+O(k_0^8).
\end{aligned}}
\]

The first term equals `(17/3) sigma^2 C_2,lead^(11)`, as required by the previously checked leading-gradient result. The new coefficient includes the physical observer clock; its initial-event dependence is not an arbitrary coordinate-time convention.

The independently obtained linear result is

\[
\boxed{\begin{aligned}
C_2^{(11)}={}&\frac{4\pi}{2025}\sigma^2k_0^4\Delta^2\\
&-\frac{4\pi}{637875}\sigma^2 k_0^6(o-e)^2(o+e)
\left(59o^3+59o^2e-106oe^2+114e^3\right)+O(k_0^8).
\end{aligned}}
\]

It agrees with a separate expansion of the report's linear transfer integral. As usual, these coefficients enter the physical spectrum with prefactors `epsilon^2` and `epsilon^4`, respectively.

## 2. Independent local field data through fourth gradient order

For a shell pair choose

\[
p=(b,0,c),\qquad q=(-b,0,c),\qquad
k_0^2=b^2+c^2,\quad p\cdot q=c^2-b^2,
\]

\[
b^2=\tfrac12 k_0^2(1-\mu),\qquad c^2=\tfrac12 k_0^2(1+\mu).
\]

This treats every non-collinear shell pair by rotation, and the polynomial data extend to the degenerate pairs. Define

\[
P=pp^T-\tfrac13 k_0^2\mathbf1,\qquad
Q=qq^T-\tfrac13 k_0^2\mathbf1.
\]

At first order,

\[
H_{1,p}=\left(-\frac{t^2}{3}+\frac{k_0^2t^4}{60}\right)P+O(k_0^6),
\qquad
\nu_{1,p}=-\frac{k_0^2t^2}{9}+\frac{k_0^4t^4}{270}+O(k_0^6),
\]

and similarly for `q`. For the report's symmetric second-order kernel write

\[
A_2(p,q;t)=t^2U+t^4V+O(k_0^6).
\]

Direct expansion of the spatial equation, with no scalar/vector/tensor inverse projectors, gives

\[
U=\frac29\operatorname{diag}
(2b^2-3c^2,\,-b^2-3c^2,\,-b^2+6c^2),
\]

\[
V=\frac1{810}\operatorname{diag}\left(
-22b^4+97b^2c^2+101c^4,\quad
11b^4+184b^2c^2+101c^4,\quad
11b^4-281b^2c^2-202c^4\right).
\]

The lapse data are

\[
\nu_2=\frac{3b^2+5c^2}{9}t^2
-\frac{(b^2+c^2)(13b^2+25c^2)}{405}t^4+O(k_0^6).
\]

The latter fourth-order coefficient is not needed for this quadratic observable, but is independently checked. The script verifies the Hamiltonian constraint through `t^2` and the momentum constraints through `t^3` using these data.

### How the field coefficients are obtained

Let `h=-t^2/3+...` be the first-order Hessian coefficient. With the unknown `A_2` temporarily zero, the coefficients in `gamma=exp(2g) exp(H)` are

\[
\gamma_p=2\mathbf1+hP,\qquad\gamma_q=2\mathbf1+hQ,
\]

\[
\gamma_{pq}=4\mathbf1+2h(P+Q)+\tfrac12h^2(PQ+QP).
\]

Compute Christoffels, spatial Ricci curvature, and the raised lapse Hessian directly from these matrices and Fourier derivatives. If the resulting mixed trace-free source is `c_T=c_0+t^2 c_2+...`, then

\[
U=-c_0/3,\qquad
V=-\frac1{10}\left(c_2+\operatorname{TF}R_1[U]
+2\operatorname{TF}(KK^T)a_2\right),
\quad a_2=(4k_0^2-p\cdot q)/9.
\]

The factors 3 and 10 come from the time operator `(A''+2A'/t)/2`. This local construction explains why no `1/|p+q|^2` terms are needed in the gradient-expanded map or field data here.

## 3. Hand ray expansion

Let `lambda=o-t`, `p_n=p.n`, and `V_p=2p_n n-p`. A superscript in square brackets in this section denotes **momentum degree**, not primordial-amplitude order. The first-order ray displacement and direction through the necessary degrees are

\[
X_p^{[0]}=-\lambda n,\qquad e_p^{[0]}=n,
\]

\[
X_p^{[1]}=-\frac{i\lambda^2}{2}V_p,\qquad
e_p^{[1]}=i\lambda V_p,
\]

\[
\begin{aligned}
X_p^{[2]}={}&\frac{\lambda o^2}{6}Pn
-\frac{\lambda^2(t+2o)}9Pn
+\frac{\lambda^2(t+2o)}{18}P_{nn}n
+\frac{\lambda^3}{6}p_nV_p,\\
e_p^{[2]}={}&\left(\frac{o^2}{6}-\frac{t^2}{3}\right)Pn
-\frac{o^2-t^2}{6}P_{nn}n
-\frac{\lambda^2}{2}p_nV_p.
\end{aligned}
\]

These follow by integrating the report's first-order ray source and the observer direction `e_p(o)=(1+H_{1,p}(o)/2)n`. The identity `dX_p^[2]/dt=e_p^[2]` is checked symbolically. There is no fixed-point ray generator.

The first-order shear coefficients are

\[
\sigma_p^{[2]}=-\frac t3P,\qquad
\sigma_p^{[4]}=\frac{k_0^2t^3}{30}P,
\]

and the quadratic shear coefficients are

\[
\sigma_{pq}^{[2]}=tU-\frac{2t}{3}(P+Q),
\]

\[
\sigma_{pq}^{[4]}=2t^3V+\frac{t^3}{9}(PQ+QP)
+\frac{k_0^2t^3}{15}(P+Q).
\]

They are inserted directly into the report's polarized expression

\[
I_{pq}=\int_e^o dt\left[
\sigma_{pq,nn}e^{iK_n\lambda}
+i(q\cdot X_p)\sigma_{q,nn}e^{iq_n\lambda}
+i(p\cdot X_q)\sigma_{p,nn}e^{ip_n\lambda}
-2\sigma_q(n,e_p)e^{iq_n\lambda}
-2\sigma_p(n,e_q)e^{ip_n\lambda}\right]
-[F_p d_q+F_q d_p]_e.
\]

Only the even-in-`n` parts of degrees two and four can contribute to the quadrupole. Odd-degree pieces are odd sky polynomials and have exactly zero `ell=2` projection; they are not needed for the `k_0^6` quadrupole power.

For auditability, the script separates the degree-four integrand into:

- field: `sigma_pq^[4] - (K_n lambda)^2 sigma_pq^[2]/2`;
- position: the two explicit `i k.X` terms with `X^[0,1]`;
- direction: the explicit `sigma(n,e)` terms with `e^[0,1,2]`;
- the moving emission endpoint;
- clock composition;
- the normalized-map product.

All time integrals are integrals of finite polynomials over `[e,o]`.

## 4. Endpoint, clock, and normalization

At the degree required here,

\[
d_p=\frac{k_0^2e^3}{9}+O(k_0^3),\qquad
D_p=\frac{k_0^2(o^4+3i^4)}{36o}+O(k_0^4).
\]

`d_p` is the first-order displacement of the physical emission surface; `D_p` follows from the fixed elapsed-proper-time condition starting at `T=T_I`, not from a fixed coordinate observation time.

Their contributions to the degree-four quadratic map kernel are

\[
\mathcal K_{2,\mathrm{end}}^{[4]}
=-\frac{k_0^2e^4}{27}(P+Q)_{nn},
\qquad
\mathcal K_{2,\mathrm{clock}}^{[4]}
=\frac{k_0^2(o^4+3i^4)}{108}(P+Q)_{nn}.
\]

The normalization term is

\[
\mathsf P_0\left[\tau^4P_{nn}Q_{nn}\right],\qquad
\tau^2=(o^2-e^2)/6.
\]

It is retained in the kernel. Its contribution to the `k_0^6` *spectrum* vanishes only after the exact sky and relative-momentum averages. This is consistent with the report's general cancellation of the mixed log-map term at this gradient order, and is not permission to omit normalization from the map.

The leading quadratic map is the pure quadrupole

\[
\mathcal K_2^{[2]}=-\frac{o^2-e^2}{2}U_{nn}.
\]

## 5. Exact angular averages and independent comparison

For a normalized sky mean,

\[
\langle n_x^a n_y^b n_z^c\rangle_\Omega=
\begin{cases}
\dfrac{(a-1)!!(b-1)!!(c-1)!!}{(a+b+c+1)!!},&a,b,c\text{ even},\\
0,&\text{otherwise}.
\end{cases}
\]

Since the leading kernel is a pure quadrupole, its inner product with any higher-degree sky polynomial automatically selects that polynomial's quadrupole. The shell measures and the two connected Wick pairings give

\[
\left[\frac14C_2^{(22)}\right]_{k^4}
=\frac{\pi\sigma^4}{5}\int_{-1}^1d\mu\,
\left\langle(\mathcal K_2^{[2]})^2\right\rangle_\Omega,
\]

\[
\left[\frac14C_2^{(22)}\right]_{k^6}
=\frac{2\pi\sigma^4}{5}\int_{-1}^1d\mu\,
\left\langle\mathcal K_2^{[2]}\mathcal K_2^{[4]}\right\rangle_\Omega.
\]

Every integral is evaluated exactly. For the optional output comparison, the complete quadrupole tensor is extracted as

\[
Q_{ij}[f]=\frac{15}{2}\left\langle f(n)
\left(n_i n_j-\frac13\delta_{ij}\right)\right\rangle_\Omega.
\]

The saved output uses `p=(k_0,0,0)`, `q=k_0(mu,sqrt(1-mu^2),0)`. After an explicit rotation into the symmetric-pair basis, **all components** of its degree-four quadrupole tensor equal the hand result, for symbolic `mu,e,o,i`. This is stronger than merely matching one integrated spectrum coefficient.

### Nonzero endpoint and clock tests

Separately, the endpoint and clock contributions to the new spectrum coefficient are

\[
\boxed{
\delta_{\rm end}\!\left(\tfrac14C_2^{(22)}\right)
=\frac{32\pi}{18225}\sigma^4 k_0^6\Delta e^4,
\qquad
\delta_{\rm clock}\!\left(\tfrac14C_2^{(22)}\right)
=-\frac{8\pi}{18225}\sigma^4 k_0^6\Delta(o^4+3i^4).
}
\]

For example, at fixed `o,e`, the coefficient of `i^4` is

\[
-\frac{8\pi}{6075}\sigma^4 k_0^6\Delta.
\]

The position and direction contributions are also individually nonzero. Their exact expressions are printed by the script. Dropping any of these four pieces fails the spectrum comparison. This is the desired check of terms not isolated by the soft Ward tests. The separated pieces are bookkeeping contributions, not individually gauge-invariant observables.

## 6. Remaining work

- The next cubic contribution `Re C_2^(13)/3` is not included in this check. It is needed before stating the full `epsilon^4 k_0^6` temperature spectrum.
- No late-time limit is taken: this is a simultaneous small-momentum expansion at fixed finite physical labels.
- No main-report or shared-generator files are changed by this check. The author can incorporate the coefficient and ledger entry independently of the ongoing cubic calculation.
