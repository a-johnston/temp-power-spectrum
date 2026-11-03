# Critical review and an exact reduction of the cubic field sources

## Executive assessment

The report's principal weakness is **not an identified major correctness error**. It is the distance between a constructive answer and a substantially evaluated finite-wavelength answer. The problem explicitly permits completely specified definite integrals, so an unevaluated integral is not, by itself, a failure. Nevertheless, the finite-wavelength fourth-order response still has much more structure than its presentation in §5.3 reveals.

A useful new reduction is available **before constructing the cubic kernels**:

1. On the monochromatic shell, the longitudinal cubic momentum source contains **no term made of three first-order fields alone**, pointwise in relative angle and at arbitrary finite time.
2. The homogeneous second-order tensor cancels out of both angle-averaged cubic sources. The homogeneous lapse does not cancel.
3. The two averaged sources require only the degree-zero, degree-one and degree-two moments of the second-order lapse and the two second-order metric amplitudes: **nine angular moments**, plus the separately retained homogeneous lapse.
4. A normal-form change of the cubic lapse removes every derivative of a second-order metric response from its lapse source. This change has a general-momentum version, not just a shell version.
5. The nine moments have explicit time-integral representations with **all primordial relative-angle integrals evaluated**. They can be expressed in entire hypergeometric kernels, or, more economically, in finite sums of the report's existing `Ephi` and `J` functions, with at most `J_3`.
6. Combining this normal form with the momentum adjoint gives a top-order cubic transfer formula involving no third-order kernel and no derivatives of second-order response functions. The integration-by-parts boundary terms are explicit.

These are genuine finite-time reductions of the intermediate field problem. They are **not** a cancellation of the full cubic observable, and do not evaluate all of `Q_l` or `R_l`. In particular, the lower-order propagation products in the cubic map remain to be combined with the reduced field contribution.

These findings have now been integrated into the main report through `cubic_source_reduction.tex` and `cubic_adjoint_reduction.tex`, supported by `solution/check_shell_cubic_sources.py`. The problem statement remains unchanged. Follow-up work on the propagation target is included in `propagation_adjoint.tex` and checked by `solution/check_propagation_adjoint.py`; the main report records the current combined status.

## 1. Critical review of the completion claim

### Reading the problem statement correctly

The statement deliberately lets the solver choose and declare the nonlinear Gaussian variable and the physical initial-event label. Those are required parts of an answer, not missing inputs that invalidate the problem. It also explicitly allows exact definite integrals and explicitly defined special functions: an elementary closed form is a desirable improvement, not a literal requirement for every remaining integral. Conversely, a mere instruction to solve a differential equation would not meet its analytical requirement; complete source prescriptions, Green kernels, limits and boundary data matter.

No primordial spectral shape is fixed before the final shell specialization. Therefore claims about arbitrary spectra must carry explicit existence conditions or acknowledge divergent cases; finite-time compact-support results do not settle arbitrary UV/IR tails. Finally, the full map through third order and the much smaller Gaussian contraction needed for its fourth-order spectrum are different deliverables. A contracted simplification cannot substitute for the former.

The shell specialization in this note is applied to the report's already constructed general observable and field equations. It is not an assumption used to omit terms while constructing that observable.

### What is already strong

- The physical measurement is specified before perturbation theory. The cancellation of the observer temperature is an exact cancellation in the normalized map, not an illicit removal of observer terms.
- The Gaussian primordial variable, regularity prescription and nonzero constrained higher-order initial data are explicit.
- The triangular field scheme has a constraint-propagation proof. The report correctly distinguishes that proof from finite-depth implementation tests.
- The log-map and observer-dilation decompositions retain interference and random-monopole terms. Neither is incorrectly advertised as evaluating the finite-wavelength answer.
- The finite-time entire-function formulation is preferable to taking late-time resonant limits under an angular integral.
- The curvature-based gradient results are substantive evaluated results, not merely a restatement of perturbation theory.

### What still prevents an elegant, fully evaluated answer

**First priority: intermediate cubic sources.** In §5.3 the expanded component sources are an audit trail, not an economical mathematical object. The source reductions below are a better next step than trying to simplify an already assembled final kernel.

**Second priority: combine field and propagation terms.** Reducing only the top-order field cannot establish the cancellation of the full `13` contribution. The report's orthonormal direction, centered exponent and physical clock should remain the working variables. Any proposed cancellation must include the mixed clock–emission term and all boundary terms.

**Third priority: distinguish the complete map from the spectrum's needs.** The arbitrary-triple map is constructively specified, whereas the implemented cubic configuration is `(k,q,-q)`. The new reduction further exploits the angular contraction needed by `13`; it does not replace the arbitrary-triple map demanded by the problem.

**Fourth priority: do not overextend the spectrum class or limits.** Compact spectral support and finite physical times are the established setting. Noncompact-spectrum UV control, uniform IR estimates over unbounded external momentum, and the resonance-sensitive late-time fourth-order limit remain separate tasks.

There is no reason to delete the report's cautionary language in §11. Conversely, it would be too strong to declare the whole report incorrect merely because its permitted exact integrals are unevaluated. The appropriate verdict is: **a carefully qualified analytical construction, with strong evaluated limits, but not yet the elegant finite-wavelength reduction one would hope for**.

## 2. Notation and perturbative factors

Use exactly the report's factorial and multilinear conventions. In this note `nu_3` is a cubic kernel, not its coefficient divided by `3!`.

Set

\[
 k= k_0>0,\qquad c=k/\sqrt3,\qquad x=ct,
\]
\[
 a(t)=\cos x-\frac{\sin x}{x},\qquad
 b(t)=-6\left(\frac{2\sin x}{x}-\cos x-1\right).
\]

For any shell mode `p`,

\[
 \nu_p=a,\qquad H_p=bQ_p,\qquad
 Q_p=\hat p\hat p^T-\frac13\mathbf1.
\]

Two useful exact identities are

\[
 L_k a=0,\qquad L_k=3\partial_t^2+k^2-\frac6{t^2},\qquad
 b'=6(a'-a/t).
\tag{1}
\]

Let `mu = khat.qhat`, and write the pair response at `(k,q)` as

\[
 n(t,\mu)=\nu_2(k,q),\qquad
 H_2=S(t,\mu)(uu^T-\mathbf1/3)+T(t,\mu)(vv^T-ww^T),
\]

with the same `u,v,w` as the report. Here `S=F_S`, `T=F_T`; `T` in this note is not a temperature. Define

\[
 \langle f\rangle_\mu=\frac12\int_{-1}^{1} f(\mu)\,d\mu,
 \qquad q_S=\frac{1+3\mu}{6},\qquad q_T=\frac{1-\mu}{2}.
\]

The pair `(k,-q)` contributes an identical integral after `mu -> -mu`. This accounts for the factors **2** below; no additional Wick or perturbative factorial is included in the sources.

Write `n_h(t)=nu_2(q,-q)`. The exact homogeneous tensor is `Z_h(t) Q_q`. It is not set to zero.

## 3. Exact averaged cubic sources

Let

\[
 \bar c_N=\langle c_N^{(3)}(k,q,-q)\rangle_{\hat q},\qquad
 \bar m=\langle i k\cdot c_M^{(3)}(k,q,-q)\rangle_{\hat q}.
\]

The full angular average includes the azimuth. For scalars it reduces to the `mu` mean above. Define

\[
 D_0(t,\mu)=\frac{12a}{t^2}
 -k^2(1+\mu)\left[5+4a+b(\mu-1/3)\right],
\]
\[
 B_0(t)=ak^2\left(16+8a+\frac29b^2\right)-\frac{24a^3}{t^2}.
\]

Then

\[
\boxed{
\begin{aligned}
 \bar c_N={}&B_0+\frac{12a}{t^2}n_h\\
 &+2\left\langle
 D_0 n+ak^2\left(\frac{3+\mu}{6}S-q_T T\right)
 +\frac{b'}2\left(q_S S'+q_T T'\right)
 \right\rangle_\mu .
\end{aligned}}
\tag{2}
\]

For momentum, introduce elementary coefficients

\[
\begin{aligned}
 A_S&=(3/2+a)\frac{3+\mu}{6}-\frac{\mu b}{4}q_S,&
 A_T&=-(3/2+a)q_T-\frac{\mu b}{4}q_T,\\
 B_S&=\frac{(1+\mu)b'}4q_S,&
 B_T&=\frac{(1+\mu)b'}4q_T,\\
 B_n&=-(1+\mu)\left[b'(\mu-1/3)+2a'\right].
\end{aligned}
\]

The result is

\[
\boxed{
 \bar m=-2ak^2n_h'
 +2k^2\left\langle
 A_S S'+A_T T'+2\mu a n'+B_S S+B_T T+B_n n
 \right\rangle_\mu .}
\tag{3}
\]

Every weight multiplying a response or its derivative has degree at most two in `mu`. Thus only

\[
 N_r=\langle\mu^r n\rangle_\mu,\quad
 S_r=\langle\mu^r S\rangle_\mu,\quad
 T_r=\langle\mu^r T\rangle_\mu,\qquad r=0,1,2,
\tag{4}
\]

and their time derivatives are needed.

### Why the momentum cancellation occurs

The exact mixed-matrix identities are

\[
 S_{\rm mix}=\frac12M^{-1}M',\qquad C_i=M^{-1}\partial_iM,\qquad M=e^H,
\]
\[
 D_jS_{\rm mix}{}^j{}_i
 =\partial_jS_{\rm mix}{}^j{}_i
 +3\zeta_{,j}S_{\rm mix}{}^j{}_i
 -\frac12\operatorname{tr}(C_iS_{\rm mix}).
\tag{5}
\]

On the shell, `H_1(t,X)=b(t) B(X)`, so `[H_1,H_1']=0` pointwise. Consequently:

- the quadratic commutator in the mixed shear vanishes;
- its purely first-order cubic double commutator vanishes;
- the remaining first–second-order commutator in its cubic part is antisymmetric, so its contraction with `k_i k_j` is zero;
- the quadratic commutator in `C_i` has zero trace against the symmetric first-order shear.

All direct `111` terms in the longitudinal momentum source therefore vanish **before angular averaging**. This does not require `H_2` to commute with `H_1`.

For the homogeneous pair, the tensor-dependent terms in both scalar sources are proportional to

\[
 \hat k^TQ_q\hat k=\mu^2-\frac13,
\]

whose mean is zero. The lapse contributions in (2) and (3) survive. Discarding the entire forced homogeneous pair would be wrong.

For comparison, the unaveraged pure-first-order lapse source is

\[
 c_{N,111}=ak^2\left[
 16+8a+4ab(\mu^2-1/3)
 +b^2(\mu^2/3+1/9)+10b(\mu^2-1/3)
 \right]-\frac{24a^3}{t^2}.
\]

Its mean is precisely `B_0`.

## 4. A derivative-free cubic lapse normal form

Define

\[
 \Phi(t)=a\left\langle q_S S+q_T T\right\rangle_\mu,
 \qquad \eta=\bar\nu_3+\Phi.
\tag{6}
\]

Since `a=O(t^2)` and `S,T=O(t^2)`, `Phi=O(t^4)`. Thus this change does **not** alter the constrained regular `t^2` initial coefficient.

The second-order equations are

\[
 L_K n=-C_N,\qquad K^2=2k^2(1+\mu),
\]
\[
 S''+\frac2tS'-\frac{K^2}{3}S=-4K^2n+J_S,
 \qquad T''+\frac2tT'+K^2T=J_T.
\tag{7}
\]

Using (1) and (7), the equation for `eta` becomes

\[
\boxed{L_k\eta=-\bar r_N,}
\tag{8}
\]
\[
\boxed{
\begin{aligned}
 \bar r_N={}&B_R+\frac{12a}{t^2}n_h\\
 &+2\left\langle
 D_R n+\frac{ak^2}{6}(2-3\mu-3\mu^2)S
 +\frac{ak^2}{2}(2+\mu-3\mu^2)T
 \right\rangle_\mu,
\end{aligned}}
\tag{9}
\]

where

\[
 D_R=\frac{12a}{t^2}
 -k^2(1+\mu)\left[5+2a-6a\mu+b(\mu-1/3)\right],
\]
\[
 B_R=k^2\left(8a-8a^2-\frac83a^2b+\frac49ab^2\right)
 +4aa'b'-\frac{24a^3}{t^2}.
\tag{10}
\]

There are **no derivatives of any second-order response** in (9). In obtaining `B_R`, the elementary forcing moment is evaluated exactly:

\[
 \left\langle q_SJ_S+q_TJ_T\right\rangle_\mu
 =\frac2{27}\left[k^2(12ab+72a-b^2+36)-18a'b'\right].
\tag{11}
\]

### General-momentum version

The cancellation of the metric-response derivatives is not restricted to the shell. For a first-order scalar mode of momentum `p`, let

\[
 Q_p^{\rm dim}=p p^T-\mathbf1 p^2/3,
 \qquad H_p'=\frac6{p^2}(\nu_p'-\nu_p/t)Q_p^{\rm dim}.
\]

For each partition into a single mode `p` and a pair response `J=H_{qr}` of momentum `P=q+r`, with total `K=p+P`, introduce

\[
 \phi_{p|qr}=\frac{\nu_p}{2p^2}Q_p^{\rm dim}:J.
\]

The metric-dependent part of the cubic lapse source is

\[
 -\nu_p K_iJ_{ij}p_j+\frac12H_p':J'.
\]

It equals `L_K phi_{p|qr}` plus

\[
 -\nu_p K_iJ_{ij}p_j
 -\frac{\nu_p}{2p^2}\left[
 (K^2-p^2)Q_p^{\rm dim}:J
 +3Q_p^{\rm dim}:(J''+2J'/t)
 \right].
\tag{12}
\]

Replace `J''+2J'/t` by its second-order spatial equation. The result has no derivatives of the second-order metric. Summing the three partitions gives the general cubic change of variable. At `p=0`, use the regular analytic continuation of `nu_p/p^2`; the physical growing scalar has `nu_p=O(p^2)` at finite time. Independent vector/tensor primordial seeds are not assumed.

This normal form is not permission to omit `Phi` when reconstructing `nu_3` or the observable.

## 5. Elementary second-order sources

For reproducibility, the complete source polynomials in (7) are

\[
\begin{aligned}
 C_N={}&k^2\left[
 \frac23ab-4a+\mu(4a^2-\tfrac43ab+2a)-2\mu^2ab
 \right]\\
 &+\frac{b'^2}{2}(\mu^2-1/3)+\frac{12a^2}{t^2},
\end{aligned}
\tag{13}
\]
\[
\begin{aligned}
 J_S={}&k^2\left[
 2(3+\mu)a^2+\frac{15\mu^2+13\mu-8}{3}ab
 +16(1+\mu)a\right.\\
 &\left.\qquad+\frac{5-5\mu-15\mu^2-9\mu^3}{36}b^2
 +\frac{-10+7\mu+15\mu^2}{6}b+5+7\mu
 \right]\\
 &-(1+3\mu)a'b',
\end{aligned}
\tag{14}
\]
\[
 J_T=\frac{1-\mu}{36}\left\{
 k^2\left[b^2(9\mu^2+12\mu+1)
 -6b(6\mu a+3\mu+2)-72a^2+36\right]-36a'b'
 \right\}.
\tag{15}
\]

The initial lapse coefficient is `c_2(mu)=k^2(4-mu)/9`. The scalar and tensor metric amplitudes have zero constant initial data and are the regular solutions of (7). No scalar/vector/tensor sector has been removed from the general construction; vanishing shell vectors follow from the report's exchange symmetry.

## 6. Evaluating all nine angular moments

This is an evaluation of the relative-angle integral, not a claim to evaluate the remaining time integrals.

Put

\[
 y=(1+\mu)/2,\qquad 0\le y\le1,\qquad F_\beta(z)={}_0F_1(;\beta;z).
\]

Let `d=t-s`. The exact regular Green kernels are

\[
 G_N(t,s;y)=-\frac13\left[
 d F_{3/2}(-c^2d^2y)+\frac{d^3}{3ts}F_{5/2}(-c^2d^2y)
 \right],
\tag{16}
\]
\[
 G_S(t,s;y)=\frac{s}{t}dF_{3/2}(+c^2d^2y),\qquad
 G_T(t,s;y)=\frac{s}{t}dF_{3/2}(-k^2d^2y).
\tag{17}
\]

In particular `G_N(t,s;0)=-(t^3-s^3)/(9ts)`, the exact homogeneous-momentum lapse kernel. The scalar `G_S` has hyperbolic rather than acoustic behavior; this is the scalar eigenvalue `-K^2/3` of the Cartesian spatial operator. It is not a physical instability claim, and the lapse-forced terms must be retained for its cancellations.

Expand the explicit polynomials (13)–(15), after `mu=2y-1`, as

\[
 C_N(t,y)=\sum_{j=0}^{2}c_j(t)y^j,\quad
 J_S(t,y)=\sum_{j=0}^{3}s_j(t)y^j,\quad
 J_T(t,y)=\sum_{j=0}^{3}t_j(t)y^j,
\]
\[
 c_2^{\rm init}(y)=d_0+d_1y,\qquad d_0=5k^2/9,\quad d_1=-2k^2/9.
\]

Here `c_j` are source-polynomial coefficients, not the initial coefficient `c_2` used elsewhere in the report.

Define entire universal angular kernels

\[
 \mathsf A_m^\beta(z)=\frac1{m+1}{}_1F_2(m+1;\beta,m+2;z),
\tag{18}
\]
\[
 \mathsf B_m^{\beta,\delta}(z,w)
 =\sum_{r,h\ge0}
 \frac{z^r w^h}{(\beta)_r(\delta)_h r!h!(m+r+h+1)}.
\tag{19}
\]

They are the evaluated moments of one or two `F` factors. The series in (19) is absolutely and locally uniformly convergent in both complex arguments. Both functions are universal, not definitions of the original cosmological spectrum.

For clarity, the complete angle-free moment representations follow. Write

\[
 \mathcal N_m(t)=\int_0^1 y^m n(t,y)\,dy,
\]

and analogously `mathcal S_m, mathcal T_m`. For `m=0,1,2`,

\[
\begin{aligned}
 \mathcal N_m(t)={}&t^2\sum_{j=0}^{1}d_j\mathsf A_{m+j}^{5/2}(-c^2t^2)\\
 &+\int_0^t ds\sum_{j=0}^{2}c_j(s)
 \left[-\frac{t-s}{3}\mathsf A_{m+j}^{3/2}(-c^2(t-s)^2)
 -\frac{(t-s)^3}{9ts}\mathsf A_{m+j}^{5/2}(-c^2(t-s)^2)\right],
\end{aligned}
\tag{20}
\]
\[
 \mathcal T_m(t)=\int_0^t ds\,\frac{s(t-s)}{t}
 \sum_{j=0}^{3}t_j(s)\mathsf A_{m+j}^{3/2}(-k^2(t-s)^2),
\tag{21}
\]
\[
\begin{aligned}
 \mathcal S_m(t)={}&\int_0^t ds\,\frac{s(t-s)}{t}
 \sum_{j=0}^{3}s_j(s)\mathsf A_{m+j}^{3/2}(c^2(t-s)^2)\\
 &-16k^2\int_0^t ds\,\frac{s(t-s)}{t}\,\mathcal V_m(t,s),
\end{aligned}
\tag{22}
\]

where `z=c^2(t-s)^2` and

\[
\begin{aligned}
 \mathcal V_m(t,s)={}&s^2\sum_{j=0}^{1}d_j
 \mathsf B_{m+1+j}^{3/2,5/2}(z,-c^2s^2)\\
 &+\int_0^s dr\sum_{j=0}^{2}c_j(r)
 \left[-\frac{s-r}{3}\mathsf B_{m+1+j}^{3/2,3/2}(z,-c^2(s-r)^2)\right.\\
 &\left.\hspace{42mm}-\frac{(s-r)^3}{9sr}
 \mathsf B_{m+1+j}^{3/2,5/2}(z,-c^2(s-r)^2)\right].
\end{aligned}
\tag{23}
\]

The original `mu` moments (4) are finite combinations:

\[
 N_r=\sum_{m=0}^{r}\binom{r}{m}2^m(-1)^{r-m}\mathcal N_m,
\tag{24}
\]

and the same relation holds for `S_r,T_r`. All universal kernel indices needed in (20)–(23) are at most five.

The homogeneous lapse retained in the sources is specified without a projector or limiting ambiguity:

\[
 n_h(t)=\frac{5k^2t^2}{9}
 -\int_0^t\frac{t^3-s^3}{9ts}\,C_N(s,\mu=-1)\,ds.
\tag{25}
\]

### No new special functions are actually necessary

The hypergeometric form is convenient for manifest regularity. All the required kernels also reduce to finite combinations of the existing entire functions in §5.3.

Set `alpha^2=-4z`, `gamma^2=-4w`. Then

\[
 F_{3/2}(-\alpha^2u^2/4)=\frac{\sin\alpha u}{\alpha u},\qquad
 F_{5/2}(-\alpha^2u^2/4)=
 \frac{3(\sin\alpha u-\alpha u\cos\alpha u)}{\alpha^3u^3}.
\tag{26}
\]

To evaluate (18) or (19):

1. Set `y=u^2`, so the measure times the monomial is `2u^(2m+1) du`.
2. Insert (26), multiply out, and express sines/cosines as exponentials. This gives a **finite** Laurent-exponential sum `sum C_(omega,p) u^p exp(i omega u)`.
3. For `p>=0`, integrate it as `C Ephi_p(i omega)`.
4. For `p=-M<0`, integrate it as `C J_M(omega,1)`, using precisely the Taylor subtraction of the report. The discarded principal parts cancel in the sum because the original integrand is regular at `u=0`.

The largest negative power is `u^-3` (for a `3/2,5/2` product at `m=0`), so **only `J_1,J_2,J_3` are required**. Frequencies are `+/- alpha`, or `+/- alpha +/- gamma`. Hyperbolic kernels correspond to imaginary frequencies; the entire definitions apply unchanged.

For example,

\[
 \mathsf A_0^{3/2}(-\alpha^2/4)=\frac{2(1-\cos\alpha)}{\alpha^2},\qquad
 \mathsf A_0^{5/2}(-\alpha^2/4)=\frac{6(\alpha-\sin\alpha)}{\alpha^3},
\]
\[
 \mathsf B_0^{3/2,3/2}(-\alpha^2/4,-\gamma^2/4)
 =\frac{\operatorname{Cin}(\alpha+\gamma)-\operatorname{Cin}(\alpha-\gamma)}{\alpha\gamma}.
\tag{27}
\]

Every apparent zero-frequency or zero-argument division is prescribed by the entire series (18)–(19). No singular resonance boundary prescription is needed at finite time.

### Legitimacy of the exchanges

At finite `t`, the arguments of the entire functions range over compact sets. Their power series converge absolutely and uniformly there. The elementary source `C_N(s,mu)=O(s^2)` uniformly in `mu`, while `J_S,J_T=O(1)`. Thus the apparent `1/s` and `1/r` factors in (20), (23), and (25) are integrable against the sources. This justifies the angular/time exchanges and termwise summations. Nothing in this argument permits a late-time limit under those integrals.

## 7. A derivative-free top-order cubic adjoint

There is a useful further combination with the existing adjoint in the report. Define `P` and `Q` (scalar auxiliary functions, not momenta or spectra) by

\[
 P=-2an_h+2\langle A_SS+A_TT+2\mu an\rangle_\mu,
\]
\[
\begin{aligned}
 Q={}&2a'n_h+2\left\langle
 \left[\frac{(1+2\mu)b'}4q_S-\frac{a'(3+\mu)}6\right]S\right.\\
 &\quad+\left[\frac{(1+2\mu)b'}4q_T+a'q_T\right]T\\
 &\quad\left.+\left[-(1+\mu)b'(\mu-1/3)-2a'(1+2\mu)\right]n
 \right\rangle_\mu .
\end{aligned}
\tag{28}
\]

Equation (3) is exactly `bar m = k^2(P'+Q)`. Again the coefficients require moments of degree at most two only.

Let `e=bar t_star`, `o=bar t_o`,

\[
 w(t)=w_\ell(k(o-t)),\qquad Z=P-2\Phi,
\]

and retain the report's adjoint kernel `mathcal W_l^0`. The top-order cubic transfer appearing in its Proposition on the adjoint can be written

\[
\boxed{
\begin{aligned}
 \mathcal A_\ell={}&14T_\ell
 +\int_0^o\mathcal W_\ell^0(s)\bar r_N(s)\,ds\\
 &+\frac32\,[wZ]_{e}^{o}
 +\frac32\int_e^o\left[wQ-w'Z+\frac{2w}{t}\Phi\right]dt.
\end{aligned}}
\tag{29}
\]

Here `w'` is its total derivative with respect to `t`. The corresponding spectrum contribution is still

\[
 \frac13\operatorname{Re}C_\ell^{(13)}\big|_{\rm top}
 =4\pi\sigma^4 T_\ell\mathcal A_\ell.
\]

To prove (29), use the scalar momentum constraint

\[
 \bar h_{S3}'=6(\bar\nu_3'-\bar\nu_3/t)+3\bar m/k^2,
\]

solve (8) with the same regular homogeneous coefficient `14 a`, recover `bar nu_3=eta-Phi`, and integrate the `Phi'` and `P'` terms by parts. The boundary term in (29) is mandatory, including the observer boundary. The factor 14 is not a new approximation or renormalization.

This formula is a reduction of the **top-order field** response. It does not include the lower-order field products in the ray, endpoints, clock, tetrad or normalization. Those must still be added exactly as in the original report.

## 8. Verification and recommended next work

The new exact check is

```
.venv/bin/python -B solution/check_shell_cubic_sources.py --compare-engine
```

It verifies:

- (2)–(3), starting with arbitrary functions `a(t), b(t), n_+/-, S_+/-, T_+/-, n_h, Z_h`;
- an independent comparison with the original Christoffel-based `engine.equations`, at arbitrary finite times, not a truncated series;
- pointwise absence of the pure-first-order longitudinal momentum source;
- the homogeneous tensor cancellation and the elementary lapse-source mean;
- the general-momentum normal-form identity (12) and the shell coefficients in (9);
- the complete quadratic source polynomials (13)–(15) against the original engine;
- the elementary forcing average (11), the momentum integration by parts, and all boundary coefficients in (29);
- the exact Green equations and diagonal jumps, including zero momentum;
- all-order coefficient identities defining the universal angular kernels;
- cancellation of all Laurent principal parts for every angular kernel needed here, and the maximum pole order three;
- omission controls for mixed tensor terms and the normal-form derivative terms.

The existing `check_adjoint_exact.py`, `check_cubic_product_identity.py`, and `check_report_identities.py` also pass. The new script passes Ruff lint/format checks and the project's type checker. The full repository suite, including the long legacy generator runs, was not rerun for this review.

The moment formulas (20)–(23) are derived by substitution into the verified Green representations and the analytically justified angular integration. They have not been turned into a new end-to-end implementation of the full spectrum. Nor does this check independently rederive the report's entire observable.

### Next milestone

Do **not** start by writing another expanded cubic kernel. Insert (29) into the report's top-order contribution and organize the remaining second-order ray/endpoint/clock products in the same nine-moment language wherever their angular weights allow it. Mixed line-of-sight phases will require additional moment kernels; the finite-polynomial reduction proved here must not be assumed for those phases. Integrate their time derivatives by parts only in combined expressions, retaining the physical endpoints.

A successful next cancellation would be an identity between those propagation products and the explicit `Phi`, `Z` and boundary pieces in (29). This note does not claim that identity has been established. It supplies a much smaller, independently checked object against which to look for it.

### Follow-up now included in the report

The new propagation adjoint eliminates the parts of `X_2,v_2` driven by `H_2`, and removes partial time derivatives of `H_2` from the bulk, at general momenta and finite times. It retains the observer and emission boundaries, and explicitly includes every term in this sector from emission, proper clock and centered normalization. The shared-time shell emission contraction reduces to `n_h+3N_2-N_0` and the derivative of a single metric moment. Mixed-time phases still survive in the remaining bulk and clock/normalization terms, so a nine-moment closure for those terms is not assumed. No cancellation of the full cubic response has been established. See §6.1 of the rebuilt main report for the complete formulas and proof.
