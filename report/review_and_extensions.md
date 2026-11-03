# Critical review and extensions

## Verdict and scope

The report is a substantial **constructive analytical integral solution**, with genuinely evaluated long-wavelength results. It is **not an evaluated general-wavelength fourth-order spectrum**. That distinction matters, but an unevaluated integral is not automatically a failure: `problem_statement.tex` expressly allows exact definite-integral representations with specified integrands, limits and boundary prescriptions.

I found no discrepancy in the quadrupole coefficients through total gradient degree six or in the leading octupole coefficients. The new independent calculation below reproduces considerably more than those numbers: the full normalized map kernels through four gradients at general `(k,q,-q)`, before angular projection or Gaussian averaging. The review does not constitute an independent arbitrary-time substitution of the complete cubic observable, a proof of ultraviolet convergence, or an evaluation of the remaining finite-wavelength integrals.

This review covers the pre-existing analytical framework and the curvature/product-identity extensions described here. The author's concurrent `reconstruct_transfer.py` work is separate and is not certified by this review.

## Is the report honest about its results?

**Mostly yes in its detailed qualifications, but its headline presentation was stronger than its actual level of evaluation.** There is no basis to infer an intent to conceal. The report repeatedly distinguishes proofs, finite exact series tests, unevaluated integrals, and exploratory numerics. It explicitly warns that dilation tests are blind to several propagation terms, that the adjoint identity covers only the top-order field contribution, and that resummation applies to the leading-gradient expression rather than the full nonlinear solution. Those are unusually useful disclosures.

Nevertheless, several claims needed correction or greater prominence:

1. **“Final report”, “we compute”, and “full solution” can suggest a closed evaluated spectrum.** The general fourth-order answer is a nested-integral representation. The title, opening abstract, README and a new completion-scope paragraph now make this visible immediately. The problem permits this form of answer; the issue is presentation, not the mere presence of integrals.
2. **“Every analytical step … backed by an exact SymPy check” overstated the verification coverage.** Constraint propagation, covariance, compact-domain convergence and exchange of integrations have analytical arguments. Computer-algebra checks cover particular identities, configurations and series depths. Neither is interchangeable with an independent finite-time verification of the assembled cubic map. The README now says so.
3. **The claim that the finite-wavelength difficulty was localized to mixed-time sky factors was too strong.** Single-time factors still leave angular moments of field amplitudes and nested time integrations. Entire dependence establishes regularity, not an elementary antiderivative. That passage is corrected.
4. **The open-items section was stale in the opposite direction:** it said the cubic degree-six coefficient came from the generator alone, although the independent curvature check already reproduced it. It also understated the octupole comparison now available. Both are updated.
5. **“UV convergence … established by power counting only” was internally inconsistent.** Power counting suggests a condition; it does not establish convergence. The wording now reflects that.
6. **A transfer node eliminates the cubic interference, not the unevaluated quadratic integral.** The node identity is exact, but the numerical node positions and coefficients in the exploratory appendix are not evaluated analytical node spectra. Readers should not interpret the phrase “complete fourth order at nodes” as providing those numbers analytically.
7. **The long-wavelength resummations must not be marketed as finite-wavelength results.** The report already makes this qualification explicitly; it should remain attached to every headline coefficient such as `29/3` or `46`.

The most accurate description is: **a complete prescription at the requested amplitude orders within the stated domain, with selected closed-form evaluations and substantial but non-exhaustive independent verification.** The stronger claim that all remaining finite-time cancellations have been found, or that every assembled kernel has been independently verified at arbitrary times, is not justified.

## Requirements audit

| Problem requirement | Assessment |
|---|---|
| Specify the Gaussian nonlinear variable and absence of independent higher modes | The marked-chart Gaussian `zeta = epsilon g`, conformally flat seed, no vorticity seed and regular branch specify an ensemble. Nonlinear metric/fluid initial data are not incorrectly set to zero. |
| Satisfy all Einstein constraints | The triangular construction and residual-Bianchi argument include the essential homogeneous Hamiltonian matching. Regularity alone would not suffice; the report correctly provides a counterexample. Exact series checks supplement, rather than replace, the proof. |
| Exact physical measurement, fixed clock, transported tetrad | The conformal redshift construction retains the emission root, clock and frame. The general-spectrum Fermi–Walker rotation must still be added to the symmetric-frame map; its real interference vanishes in the spectrum, not necessarily in the map. |
| Normalize before dipole subtraction | The centered-exponent construction preserves random-monopole and dipole products. The new curvature implementation retains the dipole until the final projection as well. |
| Complete map through cubic order, including vectors/tensors | General momenta are constructively specified by coefficient extraction and Green operators. Only `(k,q,-q)` is expanded in the cubic ancillary implementation. Shell symmetries do not authorize dropping lower-order tensor products. No arbitrary-triple finite-time substitution was added in this review. |
| Correct amplitude, Fourier and Wick normalization | The factorial convention and connected pairings give the stated `1/4 C22 + 1/3 Re C13`. The new independent raw-kernel and spectral comparisons agree with these factors, including the second `K1[4] K3[2]` ordering at degree six. |
| Analytical rather than numerical answer | The framework is analytical. The remaining finite-dimensional integrals are permitted by the problem. Exploratory numerical evidence cannot certify their analytical evaluation or asymptotics. |
| Finite times, general momenta, declared limits | The field/map construction precedes shell restriction. The general-spectrum convergence claim is restricted to finite measures with compact support away from zero. UV/IR extensions are not established generally. The new general-moment formulas are gradient coefficients, not a finite-wavelength replacement. |
| Verify gauge invariance at physical labels | Covariance and Ward arguments are meaningful. Soft/dilation tests alone are not enough to test emission, clock and propagation products. The raw full-map comparison now provides an additional independent check through four gradients. |

These conclusions do not prove that a smaller finite-time closed form is impossible. “Derive closed forms wherever possible” remains a research obligation, not a theorem discharged by failing to find another simplification.

## 1. Explicit cubic product identity

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

## 2. Replace the expensive kernel construction with local curvature

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

The earlier supporting derivation already established the even part. The odd part is its first covariant Taylor term along the seed radial geodesic; it contains both the dipole and octupole. Terms from the time-dependent displacement first affect the driver at fifth gradient order. The direction response already contributes at fourth order and is retained above. Fermi–Walker rotation is beyond this gradient order.

`solution/curvature_kernels.py` extracts nilpotent-label coefficients for general

\[
k=(p,0,0),\qquad q=(b_1,b_2,0),\qquad -q.
\]

It imports the independent local-curvature algebra, **not** the field engine or ray generator. The leading, third-gradient and fourth-gradient terms are centered and assembled into the full map. The resulting polynomial term counts are 37 for `K1`, 478 for `K2(k,q)` and 318 for `K3(k,q,-q)`.

The optional comparison uses the **raw** saved legacy generator outputs, not their simplified SymPy cache. It groups terms by formal inverse-momentum powers and clears a common polynomial denominator over exact Gaussian rationals. Equality holds component by component at general momenta and symbolic emission, observation and initial times. Polynomial equality establishes the continuation across the removable denominator loci as well.

Construction takes about 3 s; construction, raw comparison and all spectral contractions take roughly 12–15 s in the review runs. This replaces the approximately 13-hour run **for the full map through gradient degree four in this configuration**, not for arbitrary wavelengths or higher gradient orders. The old run remains a useful independent derivation. Its saved pickle files are ignored local artifacts, so `--compare` requires them; `--spectrum` and kernel construction are cache-free. `--output NEW_PATH` exports exact polynomial kernels without overwriting an existing file.

## 3. A new general-spectrum simplification

Define `mj = integral k^j Pcal(k) d ln k`, and let `[ ]_[6]` denote total gradient degree six at fixed times. The report's three-tensor quadrupole has time coefficients

\[
\beta_H=-\frac{11L}{14}-\frac{3E}{40},\qquad
\beta_a=\frac{27L}{7}+\frac{7E}{20}+\frac{o^4+3i^4}4-e^4.
\]

The new general-momentum contractions give

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

## 4. What this pass does not solve

- **Closed finite-time averaged sources:** no closed evaluation of `bar cN` or `overline(i K·cM)` is supplied. The adjoint still contains those source integrals and does not include all propagation products.
- **Quadrupole through spectrum degree eight:** one needs the map through six gradients. `|Q[4]|²` alone misses `2 Q[2]:Q[6]`, including its amplitude-order decompositions.
- **Octupole through spectrum degree eight:** one needs the five-gradient map in its interference with the three-gradient map. The leading octupole is now independently checked, not extended to that next order.
- **General-wavelength fourth order:** the nested integrals are not evaluated. Local curvature polynomials at fixed gradient order do not establish a finite-wavelength closure.
- **Late-time fourth order:** no uniform treatment of the resonant angular regions is supplied. A pointwise limit of a regular finite-time kernel is insufficient.
- **Unbounded spectra and fully nonlinear evolution:** UV/IR remainder estimates, possible caustics and convergence of the full amplitude expansion remain outside the established result.
- **Deferred expensive checks:** replacing the cubic gradient-kernel bottleneck does not complete the unrelated arbitrary-time closed-form substitutions into all field equations.

The practical next step is to extend the covariant field/ray expansion by one gradient order for the octupole and two for the quadrupole, retaining the clock and emission terms before projection. This is more promising than repeating large component-kernel assembly. At finite wavelength, the scalar adjoint and separated angular moments remain useful reductions, but they have not yet supplied the requested smaller evaluated answer.

## Reproduction

```sh
.venv/bin/python -B solution/check_cubic_product_identity.py
.venv/bin/python -B solution/check_cubic_curvature_hand.py
.venv/bin/python -B solution/curvature_kernels.py --spectrum
# Optional independent comparison, when the legacy raw kernels are present:
.venv/bin/python -B solution/curvature_kernels.py --compare --spectrum
```

These checks use exact arithmetic, not fitting, numerical evolution or quadrature. Runtimes above are performance observations, not evidence for a mathematical identity. The raw comparison is stronger than agreement of final spectrum coefficients, but is still restricted to the stated momentum configurations and gradient/amplitude orders.
