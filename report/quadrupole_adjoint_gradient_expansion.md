# Evaluated quadrupole adjoints and a physical cubic comparison through four gradients

## Status

The preceding cubic assembly was committed as **`1b80932`**. This note records the continuation.

Three additional results are established:

1. The first nonzero field-adjoint response for an arbitrary quadrupole objective is evaluated explicitly, including the spatial observer trace, emission jump, singular lapse origin and initial-clock covector.
2. Its contribution to the fourth-gradient cubic response has no remaining time or angular integrals. The homogeneous-lapse partition is nonzero.
3. The **complete cubic quadrupole through four spatial gradients** agrees with an independent Einstein–radiation curvature calculation for arbitrary input triples, before Wick contraction. The comparison is a symbolic polynomial identity in six independent momentum components and all five quadrupole components, with arbitrary endpoint times. The shell contraction also recovers the known weighted `C_2^(13)/3` coefficient through sixth gradient order, including both linear–cubic orderings.

These are exact **spatial-gradient coefficients**, not a full finite-wavelength identity. They are not merely a short-time or numerical test: emission, observation and initial-clock times remain symbolic. The general cubic sky map, higher multipoles and remaining finite-wavelength integrals have not been evaluated by this continuation. The main report PDF is unchanged.

Implementation: `solution/quadrupole_adjoint_expansion.py`.
Checks: `solution/check_quadrupole_adjoint_expansion.py` and `solution/check_stationary_quadrupole_four_gradients.py`.

## 1. Objective, scaling and power counting

Scale all momenta by a common parameter `kappa`. A unit primordial input at momentum `kappa p` has the regular first fields

\[
 a_p=-\kappa^2p^2t^2/9+O(\kappa^4),\qquad
 h_p=-\kappa^2t^2Q_p/3+O(\kappa^4),\quad
 Q_p=pp^T-p^2 1/3.
\]

Use the objective

\[
 J=\langle w\Theta\rangle_n,\qquad
 w=\tfrac{15}{2}n^T\mathsf Cn,\qquad
 \mathsf C=\mathsf C^T,\quad\operatorname{tr}\mathsf C=0.
\]

Thus `average(w)=0` and `average(w STF(nn))=C`. In the full-`dOmega` source API, the weight is `w/(4 pi)`; the common clock is still contracted only once. The matrix `C` here is a quadrupole covector, not an angular power spectrum.

For physical first fields, the relevant orders are

| Quantity | First nonzero gradient order |
|---|---:|
| direct physical-input `C_H,2`, `C_H,3` | 2 |
| direct physical-input `C_N,2`, `C_N,3` | 4 |
| background spatial adjoint | 0 |
| first spatial/lapse adjoint and background lapse adjoint | 2 |
| background/first regular-origin covector | 2 |
| nonlinear initial coefficients `c_2,c_3` | 2 |

Consequently the bulk lapse pairings start at order six in the cubic quadrupole, but **both initial-data pairings contribute at order four**. Their omission is not justified by the higher order of the bulk lapse term.

The direct first-field-only cubic measurement has no emission/clock/normalization correction at order four, by power counting. This does **not** remove these effects from the full cubic measurement: they survive in the first-adjoint source paired with the physical quadratic field/initial sources. The full finite-wavelength code continues to retain the frame sector; its physical contribution starts above the order tested here.

Throughout, `0<iota<e<o`, with `iota=bar t_I`. Set

\[
 \tau^2=(o^2-e^2)/6,\qquad D=o^4-e^4,\qquad
 L=\int_e^o t(o-t)^2dt
 =o^4/12-o^2e^2/2+2oe^3/3-e^4/4.
\]

## 2. Exact angular reduction and scalar-seed curvature dual

Define

\[
 \mathsf A(u,v)=\frac17\operatorname{TF}\big[
 (u\cdot v)\mathsf C+2\operatorname{sym}(u(\mathsf Cv)^T)
 +2\operatorname{sym}(v(\mathsf Cu)^T)\big],
\]
\[
 \mathsf B(u,v)=\operatorname{TF}\operatorname{sym}(u(\mathsf Cv)^T).
\]

These are exactly

\[
 \mathsf A(u,v)=\langle w(u\cdot n)(v\cdot n)N_n\rangle,
 \qquad \mathsf B(u,v)=\langle w(v\cdot n)\operatorname{TF}\operatorname{sym}(un^T)\rangle.
\]

The finite curvature construction in the first-source note reduces, for the leading scalar primordial metric `G=2I`, to

\[
 \rho_{p,q}=\operatorname{TF}\big[
 (-q^2+\tfrac12p\cdot q)\mathsf C
 +2\operatorname{sym}(q(\mathsf Cq)^T)
 -\operatorname{sym}(p(\mathsf Cp)^T)
 -\operatorname{sym}(p(\mathsf Cq)^T)\big].
\tag{1}
\]

It satisfies `rho:j=C:ricci_mixed(p,q,2I,j)` for every STF probe `j`. The check proves this identity at arbitrary symbolic vectors and covector, independently of the closed formula. In particular `rho_(p,0)=-B(p,p)`; there is no zero-momentum singularity.

For one physical input `p` and independent **state** probe momentum `q`, put

\[
 U=-\tfrac14(p\cdot q)\mathsf C+\tfrac32\mathsf A(p,q)
 +\tfrac12\mathsf A(q,q)-\mathsf B(p,q)
 +\tfrac12[\mathsf A(p,p)-\mathsf B(p,p)],
\]
\[
 V=\tfrac32\mathsf A(p,p)-\mathsf B(p,p)-p^2\mathsf C/6.
\tag{2}
\]

The sky-integrated order-two measurement source multiplying `j'(t)` is simply

\[
 M(t)=(o-t)^2U+\tau^2V.
\tag{3}
\]

Its part `tau^2 [A(p,p)-p^2 C/3]/2` comes from exact sky normalization and is generally nonzero. The bulk coefficient of `j(t)` starts above this gradient order.

## 3. Fully evaluated first field-adjoint response

The leading background spatial adjoint is `W_0(s) C`, where

\[
 W_0(s)=\begin{cases}
 s^2(1/e-1/o)/2,&s<e,\\
 s(o-s)/(2o),&e<s<o.
 \end{cases}
\]

For `r=-1,0,2`, define the already evaluated moments

\[
 J_r(s)=\frac{o^{r+2}-b^{r+2}}{2(r+2)}
 -\frac{o^{r+3}-b^{r+3}}{2o(r+3)}
 +\mathbf1_{s<e}\frac{(1/e-1/o)(e^{r+3}-s^{r+3})}{2(r+3)},
\quad b=\max(s,e).
\tag{4}
\]

Equivalently `J_r=integral_s^o t^r W_0(t)dt`; (4) contains no unevaluated integral. In particular `J_2(0)=D/40`.

Write `Lambda_H,a=kappa^2 Lambda_H,a^[2]+O(kappa^4)`. For `b=max(s,e)`, set

\[
 F(b)=o^2/b-b-2o\log(o/b),\qquad G(b)=1/b-1/o.
\]

The first spatial adjoint is

\[
\boxed{\begin{aligned}
 \Lambda_{H,a}^{[2]}(s)={}&s^2[F(b)U+\tau^2G(b)V]
 -2s[J_0(s)-sJ_{-1}(s)]\rho_{p,q}\\
 &+\frac{4p^2s^2}{9}J_{-1}(s)\mathsf C
 -\mathbf1_{s<e}\frac{p^2es^2}{18}\mathsf C
 +\frac{p^2s^2(o^4+3\iota^4)}{72o^3}\mathsf C.
\end{aligned}}
\tag{5}
\]

The last two terms are respectively the physical emission and observer boundary contributions. They give

\[
 \Lambda_{H,a}^{[2]}(o^-)=\frac{p^2(o^4+3\iota^4)}{72o}\mathsf C,
\qquad
 \Lambda_{H,a}^{[2]}(e^-)-\Lambda_{H,a}^{[2]}(e^+)
 =-p^2e^3\mathsf C/18.
\]

Neither can generically be set to zero.

For the lapse, let `C_p=C:Q_p`, `B_pq=(p+q)^T C q`, and `b_I=max(s,iota)`. Its order-two coefficient is

\[
\boxed{\begin{aligned}
 \Lambda_{N,a}^{[2]}(s)={}&\frac{C_p}{27}\left[
 \frac{o^4-b_I^4}{4s}-s^2(o-b_I)
 +\mathbf1_{s<\iota}(\iota^4/s-s^2\iota)
 -\mathbf1_{s<e}(e^4/s-s^2e)\right]\\
 &+\frac{8B_{pq}}9[J_2(s)/s-s^2J_{-1}(s)]
 +\frac{4C_p}{27}[2J_2(s)/s+s^2J_{-1}(s)].
\end{aligned}}
\tag{6}
\]

The singular origin is physical adjoint data, not a discarded irregular mode. It gives

\[
\boxed{
 \Lambda_{c,a}^{[2]}=9\lim_{s\to0}s\Lambda_{N,a}^{[2]}(s)
 =\frac{3o^4+5\iota^4-8e^4}{20}\,\mathsf C:Q_p
 +\frac D5(p+q)^T\mathsf Cq.}
\tag{7}
\]

Equations (5)–(7) were reconstructed independently on the causal zero-gradient Green probes `G_H=s(t-s)/t` and `G_N=(t^3-s^3)/(9ts)`, with the one-sided endpoint prescription. They also agree with the coefficient extracted from the complete first-source implementation.

## 4. Logarithms cancel from the integrated first-adjoint block

The leading physical quadratic spatial source is constant in time. Its pairing requires only the integral of (5), which reduces to

\[
\boxed{
 \overline\Lambda_{H,a}^{[2]}
 =\frac L3U+(\tau^2)^2V-\frac D{120}\rho_{p,q}
 +\frac{p^2(3o^4+5\iota^4-8e^4)}{360}\mathsf C.}
\tag{8}
\]

This has **no logarithms and no integrals**. The logarithms in (5) are generally nonzero; their cancellation uses both sides of the emission surface. The check verifies the scalar integration identities for arbitrary symbolic `e,o`, not only numerical endpoints.

For physical primordial inputs `u,v`, rather than state probes,

\[
 C_{H,uv}^{[2]}=-4(Q_u+Q_v)-4\operatorname{TF}\operatorname{sym}(uv^T),
 \qquad c_{uv}^{[2]}=[2(u^2+v^2)-u\cdot v]/9.
\]

One single–pair contribution to the cubic quadrupole at order four is therefore the fully evaluated expression

\[
 \overline\Lambda_{H,a}^{[2]}:C_{H,bc}^{[2]}
 -\Lambda_{c,a}^{[2]}c_{bc}^{[2]},\qquad p=k_a,\quad q=k_b+k_c.
\tag{9}
\]

Sum the three partitions **once each**. No further factor three is attached to an individual partition. The sources in (9) include both physical primordial inputs; they are not the `delta zeta=0` state-probe derivative used in constructing the adjoint.

## 5. Shell contraction and the nonzero homogeneous lapse

For `|q|=|k|`, normalized isotropic averaging gives

\[
\boxed{
 \left\langle\sum_{a|bc}(9)\right\rangle_q
 =\frac{k^2\,\mathsf C:Q_k}{11340}
 [13908e^4-20960e^3o+15720e^2o^2-5833o^4-2835\iota^4].}
\tag{10}
\]

This is **only the first-adjoint sector**, not the full physical cubic response or spectrum.

For `(k|q,-q)`, the homogeneous tensor source is `-4Q_q`; its isotropic average is zero. The homogeneous lapse has `c_(q,-q)=5q^2/9` and contributes

\[
\boxed{
 -\frac{k^2\,\mathsf C:Q_k}{36}(3o^4+5\iota^4-8e^4).}
\tag{11}
\]

It is generically nonzero. The bulk lapse term starts too high in gradients to cancel it at this order. The physical initial-clock term remains visible.

## 6. Complete physical cubic comparison, not just a shell block

The full order-four comparison adds to (9):

- the direct first-field-only measurement, using the second-ray-free forcing from `stationary_cubic.py`;
- the direct third Einstein residual with first fields only;
- the order-two background spatial Green/phase corrections paired with the leading third residual;
- the background regular-origin covector paired with `c_abc`.

No second or third physical field response is solved in this calculation. The old fast observable module supplies **only its rational polynomial arithmetic**: its geometry, ray, root, clock and observable generators are not called.

For the independently computed reference, `check_cubic_curvature_hand.py` constructs the covariant seed-curvature tensors directly, without the field engine, ray generator or adjoint implementation. Its compact quadrupole expression supplies the complete Einstein–radiation gradient coefficient, including physical emission, clock and normalization.

The symbolic comparison uses

\[
 k_a=(p,0,0),\quad k_b=(b_1,b_2,0),\quad k_c=(d_1,d_2,d_3),
\]

with all six components independent, and tests a basis of all five STF quadrupole covectors. The actual polynomial-ring names for `d_1,d_2,d_3` are `ip,u,v`; here they are **not inverse momenta**. No inverse momentum occurs. Rotation covers generic triples, and polynomial continuity extends the identity to degenerate/zero-momentum configurations. This is an arbitrary-triple quadrupole comparison **through four gradients**, not an inference from its shell average.

For the shell, choose unit `k=(1,0,0)` and `C=Q_k`, so `w=5 P_2(n_x)`. The complete averaged cubic kernel is

\[
 \overline J_{abc}
 =\kappa^2\frac{4(o^2-e^2)}9
 +\frac{\kappa^4}{3402}
 [3984e^4-6256e^3o+4692e^2o^2-1853o^4-567\iota^4]
 +O(\kappa^6).
\tag{12}
\]

For a shell radius `k_0`, the two displayed coefficients acquire `k_0^2,k_0^4` respectively. Gaussian amplitudes have not been inserted into (12).

Combining (12) with the linear coefficient, including its order-four part, recovers the existing weighted `C_2^(13)/3` coefficient through order six. Both `J_1^[2] J_3^[4]` and `J_1^[4] J_3^[2]` are kept. This is a **reproduction** of the known low-gradient spectrum, not a newly evaluated finite-wavelength spectrum.

## 7. Reproduction and next step

```bash
.venv/bin/python -B solution/check_quadrupole_adjoint_expansion.py
.venv/bin/python -B solution/check_stationary_quadrupole_four_gradients.py
.venv/bin/python -B solution/check_stationary_quadrupole_four_gradients.py --symbolic
```

The symbolic comparison takes about 30 seconds locally. The default second check also tests two noncoplanar numerical momentum triples, but with exact arithmetic and arbitrary symbolic times. The `--symbolic` run is the general momentum-polynomial proof, not a numerical fit.

Both new check scripts, including default and symbolic modes, and eleven existing targeted regressions passed: stationary cubic assembly, first field source, unified adjoint, independent cubic curvature, propagation adjoint, Cartesian Green, observable expansion, report identities, shell cubic sources with `--compare-engine`, cubic product identities and the orthonormal driver. Repository-wide Ruff lint/format and type checking passed, along with shell syntax and whitespace checks. The independent curvature reference changed only by a type annotation permitting its already-configurable momentum fixture to contain arbitrary triples. The full long-running suite was not rerun.

The subsequent `exact_quadrupole_angular_source.md` now retains exact finite-wavelength phases in the quadrupole first source and removes its sky integrals. The mixed-time angular factors have phase vector

\[
 K=(p+q)(o-t)+p(t-u),
\]

not just one equal-time relative-angle polynomial. The exact closure uses Cartesian moments of the entire spherical phase average `j_0(|K|)`, with regular zero-phase limits and both separated normalization products retained. The frame covector also reduces to a commutator with the linear quadrupole. This is not automatic nine-moment closure; the remaining time/Green pairings, direct cubic multi-phase contractions and their physical cancellations still need to be analyzed.

No arbitrary-spectrum UV/IR convergence, late-time integral exchange, complete finite-wavelength cancellation or all-multipole cubic evaluation follows from the results above.
