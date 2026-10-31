# Cubic constrained radiation data: explicit modewise construction

## Status and dependency

This extends the independently accepted **quadratic** construction at `ea6b053` in `QUADRATIC_INITIAL_DATA.md` (abbreviated Q below). The acceptance is recorded in `../reviews/INDEPENDENT_REVIEW_02.md` and is restricted to the declared local, fixed-momentum, regular perturbative class.

The cubic construction below is **new and not independently reviewed**. Its source is explicit; its finite-time solution is an unevaluated convergent integral representation. Parent-authored exact checks compare its source to a separate expansion of the nonlinear metric equation, verify general-momentum boundary cancellations and fluid composition, and test a finite-time constant-dilation identity. No numerical evolution or quadrature is used.

This document does not construct an observed temperature kernel. Spatial convolution convergence, stochastic expectations over unbounded spectra, the observer/ray calculation, and observable gauge verification remain separate obligations.

## 1. Same experiment, gauge and seed

Keep all definitions in Q §§1–2: marked material coordinates, irrotational comoving-orthogonal flow on its smooth coordinate domain, zero shift, `a=Aη`, Gaussian `ζ=εg`, and

\[
q=I+\epsilon h+\frac{\epsilon^2}{2}j+\frac{\epsilon^3}{6}l+O(\epsilon^4),
\quad
\chi=\frac16\ln\det q-\epsilon g,
\quad \rho=\bar\rho e^{-4\chi},\quad N=a e^\chi.
\tag{C1}
\]

The full spatial metric is still general; no scalar-only projection is introduced. The seed and regular class fix

\[
h(0)=2gI,\qquad j(0)=4g^2I,\qquad l(0)=8g^3I,
\quad q-e^{2\epsilon g}I=O(\eta^2),\quad q'=O(\eta).
\]

These are **forced metric Taylor coefficients**, not independent higher-order primordial modes. In particular the cubic seed must not be set to zero, and its homogeneous components must not be removed.

Let

\[
c=\chi^{(1)}=\operatorname{tr}h/6-g,
\quad d=\chi^{(2)}=\operatorname{tr}j/6-\operatorname{tr}h^2/6,
\]
\[
w=\tfrac12\operatorname{tr}(hj)-\tfrac13\operatorname{tr}h^3,
\qquad \chi^{(3)}=\operatorname{tr}l/6-w.
\tag{C2}
\]

The trace cyclic property, not commutation of `h` and `j`, produces the displayed `tr(hj)` coefficient. All subsequent matrix products retain their order.

## 2. Explicit connection and curvature coefficients

Use the linear connection and Ricci operators of Q (12). Define

\[
\Gamma_1=\Gamma_1[h],\qquad
\Gamma_2=\Gamma_1[j]-2h\Gamma_1[h],
\]
\[
\Gamma_{3N}=-3h\Gamma_1[j]+(6h^2-3j)\Gamma_1[h].
\tag{C3}
\]

For example `(h Γ_1)^a_ij=Σ_b h_ab Γ_1^b_ij`, with Euclidean index contractions; the subscript `N` denotes only the part not containing the new unknown `l`. The full cubic connection is `Γ^(3)=Γ_1[l]+Γ_3N`.

To fix all contractions, for any connection-shaped arrays `A,B` write

\[
(\mathscr D A)_{ij}=\partial_a A^a{}_{ij}-\partial_j A^a{}_{ia},
\quad
\mathscr B(A,B)_{ij}=A^a{}_{ab}B^b{}_{ij}-A^a{}_{jb}B^b{}_{ia}.
\]

All repeated spatial indices run from 1 to 3. Then

\[
R_1=\mathscr D\Gamma_1,\qquad
R_2=\mathscr D\Gamma_2+2\mathscr B(\Gamma_1,\Gamma_1),
\]
\[
R_{3N}=\mathscr D\Gamma_{3N}
+3\mathscr B(\Gamma_1,\Gamma_2)+3\mathscr B(\Gamma_2,\Gamma_1),
\qquad R^{(3)}=R_1[l]+R_{3N}.
\tag{C4}
\]

Here `R_2` is the full quadratic Ricci coefficient, including `R_1[j]`; omitting that term would corrupt the cubic source.

For compactness define two already known covariant two-tensors

\[
A_{1,ij}=\partial_i\partial_j c-R_{1,ij},
\]
\[
A_{2,ij}=\partial_i\partial_j d-2\Gamma_1^a{}_{ij}\partial_a c
+2\partial_i c\partial_j c-R_{2,ij}.
\tag{C5}
\]

These are finite expressions in the specified lower-order fields, not unsolved source functions.

## 3. Cubic source with all factorials retained

Put `H=1/η`. Expanding the **exact** metric equation Q (3) gives

\[
\mathcal L l=\mathcal S_3[h,j],
\tag{C6}
\]

with exactly the same linear operator Q (15). Its complete cubic source is

\[
\begin{split}
\mathcal S_{3,ij}={}&3(h'j'+j'h')_{ij}-6(h'hh')_{ij}
-6c'j'_{ij}-6d'h'_{ij}\\
&+\{4(H^2w+Hw')+24H^2cd-16H^2c^3\}\delta_{ij}\\
&+[-12H^2d+24H^2c^2-12Hd']h_{ij}
-12(H^2c+Hc')j_{ij}\\
&-2\partial_i\partial_j w
-6\Gamma_1^a{}_{ij}\partial_a d
-6\Gamma_2^a{}_{ij}\partial_a c\\
&+6(\partial_i c\partial_j d+\partial_i d\partial_j c)
-2R_{3N,ij}\\
&+12c A_{2,ij}+(12d+24c^2)A_{1,ij}.
\end{split}
\tag{C7}
\]

This contains every term of amplitude degree three, including the induced components of `j`. No constraints were used to eliminate terms from (C7). Some checks on its origin are:

- `[q' q^(-1) q']^(3)=3(h'j'+j'h')−6h'hh'`.
- `[−2χ'q']^(3)=−6(c'j'+d'h')`.
- `[e^(−2χ)]^(3)=−2χ^(3)+12cd−8c³`.
- `[2e^(2χ) A]^(3)=2A^(3)+12c A_2+(12d+24c²)A_1` because the background spatial bracket is zero.

After writing `χ^(3)=tr(l)/6−w`, the terms linear in `l` are precisely moved into `L`. These identities also explain why replacing the two matrix products `h'j'` and `j'h'` by a single unsymmetrized product is invalid.

## 4. Fully specified Fourier source and lower-order nesting

Let `H_ij(η;k)` be the exact linear kernel Q (11) and `J_ij(η;k,p)` the exact quadratic kernel Q (20). Distinguish the matrix kernel `H_ij` from the scalar `H=1/η`. Define

\[
l_{ij}(\eta,x)=\int_{k_1}\int_{k_2}\int_{k_3}
L_{ij}(\eta;k_1,k_2,k_3)\prod_{a=1}^3g(k_a)
 e^{i(k_1+k_2+k_3)\cdot x},\quad
\int_k=\int_{\mathbb R^3}\frac{d^3k}{(2\pi)^3}.
\tag{C8}
\]

There is no extra `1/3!` in this convolution. The symmetric source kernel `S_3(η;k_1,k_2,k_3)` is this explicit finite polarization of (C2)–(C7):

1. Introduce three distinct formal amplitudes `z_1,z_2,z_3`.
2. Substitute `h=Σ_a z_a H(η;k_a)` and `j=2Σ_(a<b) z_a z_b J(η;k_a,k_b)`.
3. Substitute `g=Σ_a z_a`; this also fixes `c,d,w` by (C2).
4. Apply `∂_i z_a=i k_ai z_a`; on products this uses the sum of their momenta. Time derivatives act only on `H,J` and their displayed integrals.
5. Take **one sixth** of the coefficient of `z_1 z_2 z_3`.

Equivalently work in the finite ring `z_a²=0`; that cannot change the coefficient being extracted. The labels remain distinct even if some input vectors coincide or are opposite. This prescription is a finite sum of indexed products and permutations, not a differential equation or unspecified functional inverse. It fixes every spatial derivative and convolution coefficient.

Every pair kernel occurring at source time `s` means exactly

\[
J(s;k_a,k_b)=4I+\int_0^s dt\,G_{k_a+k_b}(s,t)F_2(t;k_a,k_b),
\]
\[
\partial_s J(s;k_a,k_b)=\int_0^s dt\,\partial_sG_{k_a+k_b}(s,t)F_2(t;k_a,k_b).
\tag{C9}
\]

Here `F_2` is the explicit regular quadratic source Q (12)–(17), and `G` is the explicitly nested kernel Q (18)–(19). There is no upper-end term because `G(s,s)=0`. Equations (C2)–(C9), together with those exact Q definitions, specify the complete integrands and inner limits. Each monomial in (C7) contains at most one quadratic field or its derivative; no unlisted product of two quadratic solutions is required at cubic order.

`../validation/cubic_source.py` implements (C7) with the finite ring in `amplitude_jets.py`. The separately organized test expansion forms the full inverse metric, log determinant, connection and Ricci tensor directly before taking its triple-amplitude coefficient. It does not call (C3)–(C7), but shares the finite ring and elementary connection/Ricci contraction helpers; it is not a wholly independent tensor engine.

## 5. Seed subtraction and the regular primordial endpoint

For `K=k_1+k_2+k_3`, the seed kernel is `L_0=8I`, and

\[
\mathcal L_K(8I)=\frac{16}{\eta^2}I+8K^2I+16KK^T.
\]

Define the combined regular source

\[
F_3(\eta;k_1,k_2,k_3)=S_3(\eta;k_1,k_2,k_3)
-\frac{16}{\eta^2}I-8K^2I-16KK^T.
\tag{C10}
\]

The separated singular terms must **not** be integrated from zero individually. From `h→2gI`, `j→4g²I`, one obtains `w→4g³`, whereas `c,d=O(η²)` and `h',j'=O(η)`. Thus `4H²wI` supplies exactly the subtracted singularity, and all remaining singularities are removable.

In fact the constant term cancels too:

\[
F_3(0;k_1,k_2,k_3)=0,\qquad F_3=O(\eta^2).
\tag{C11}
\]

This follows by direct substitution of the general-momentum boundary coefficients Q (7)–(10) into (C7), and is symbolically checked with all nine input-momentum components independent. A separate geometric explanation is that the exact nonlinear `η²` metric coefficient Q (7) is only quadratic in `ζ=εg`; it has no cubic amplitude coefficient. Therefore

\[
L(\eta;k_1,k_2,k_3)=8I+O(\eta^4).
\tag{C12}
\]

This does **not** mean the cubic lapse or density begins at `η⁴`. The determinant products in (C2) give

\[
\chi^{(3)}(\eta,x)=\eta^2\left[
\frac43g^2\Delta g-\frac23g|\nabla g|^2\right]+O(\eta^4),
\]
\[
\chi^{(3)}(\eta;k_1,k_2,k_3)=\frac{\eta^2}{9}
\left[-4\sum_a k_a^2+2\sum_{a<b}k_a\cdot k_b\right]+O(\eta^4).
\tag{C13}
\]

The existence of a nonzero leading fluid coefficient despite a vanishing leading metric evolution coefficient is an important normalization/product check.

## 6. Exact finite-time representation and proof domain

With the same explicit retarded kernel `G_K` as Q (19), the cubic solution is

\[
\boxed{L(\eta;k_1,k_2,k_3)=8I+
\int_0^\eta ds\,G_K(\eta,s)F_3(s;k_1,k_2,k_3).}
\tag{C14}
\]

This is an **unevaluated analytical integral representation**. Further closed-form reductions have not been exhausted. No finite-time conclusion is obtained by truncating (C12).

For every fixed finite momentum triple and finite `T`, the accepted quadratic solution and its first derivative are bounded on `[0,T]`; their even regular series hold near zero. Spatial derivatives in the source multiply fixed momenta. Hence (C7) and the cancellation (C10)–(C11) define a bounded continuous `F_3` on `[0,T]`, analytic near zero and smooth for positive time. If `M_3=sup||F_3||` and `d=||D_K||`, the `n`th term obtained by expanding the outer Green function is bounded by

\[
M_3d^n\frac{T^{2n+2}}{(2n+2)!}.
\]

The first derivative obeys the corresponding differentiated bound. This proves convergence of the outer retarded representation and permits the same Volterra/indicial uniqueness argument as at order two.

For the inner representations in (C9), the quadratic solution and derivative series converge absolutely with bounds uniform on `[0,T]`. Substitution into (C7) involves only a finite sum of products with fixed linear kernels. Where a coefficient carries `s^(−1)` or `s^(−2)`, first group the seed-subtracted even Taylor remainder: `J−4I=O(s²)`, `J'=O(s)`, and their series terms have these weighted bounds. The quadratic `n`th solution term is bounded by `M_2 d_ab^n s^(2n+2)/(2n+2)!`, and its first derivative by `M_2 d_ab^n s^(2n+1)/(2n+1)!`; dividing by at most `s²` or `s`, respectively, preserves summability. All time nesting and series exchanges used in (C9)–(C14) are therefore legitimate **at fixed momenta**. Singular seed pieces are canceled before these bounds are applied.

If `L−8I=Σ_(n≥2)η^n Z_n` and `F_3=Σ_mη^m F_(3,m)`,

\[
Z_{m+2}=\left[\frac{P_f}{(m+2)(m+3)}+
\frac{P_t}{(m+3)(m+4)}\right]
(F_{(3,m)}-D_K Z_m),\qquad Z_0=Z_1=Z_2=0.
\]

There is no positive-order resonance. Regularity and the specified seed exclude the same free singular branches and extra constant trace-free seed as at quadratic order. These arguments do not supply bounds for unbounded primordial integrations or global nonlinear existence for arbitrary random fields.

## 7. Zero total momentum and exact zero partial sums

If `K=0`, the outer kernel is exactly the elementary trace/trace-free kernel Q (18); keep the entire source in (C14). If an internal pair has `k_a+k_b=0`, its `J` and `J'` in (C9) are evaluated using Q (22) and its derivative, including the forced homogeneous trace-free response. Neither operation requires choosing a momentum direction at zero.

Consequently (C14) gives an actual **modewise metric-kernel value** at `(k,q,−q)`. Its leading lapse restriction is

\[
\chi^{(3)}(\eta;k,q,-q)
=-\frac{\eta^2}{9}(4k^2+10q^2)+O(\eta^4).
\]

This is not yet the harmonic temperature kernel in `C_l^(13)`. Propagation, physical endpoints, observer terms, normalization and angular projection remain to be constructed. Their combined restriction and stochastic integrability must be established separately; the existence of the metric-kernel restriction alone cannot certify the spectrum.

## 8. Constraints, fluid coefficients and finite-time data

The accepted offshell propagation equations Q (23) apply because the integrated fluid relations remain exact for any admissible `q`. At cubic order all lower-order constraint residuals vanish. The order-three residual therefore obeys exactly the homogeneous background propagation system Q (24), without an unaccounted quadratic source.

The constructed solution obeys `q=b+O(η²)`, `q'=O(η)` and `χ=O(η²)` through cubic order. At fixed momenta the same boundedness estimates hold:
`a² H_c^(3)=O(1)`, `a M_i^(3)=O(η)`. The accepted energy argument then forces both residuals to vanish for every finite positive time, also at zero total momentum. This is an analytical induction on perturbative order, not a claim that a finite time-jet test establishes all-time constraints. `OFFSHELL_VALIDATION.md` and `check_offshell.py` now preserve the offshell reduction used in this argument.

Let `e=χ^(3)`. The cubic physical field coefficients include

\[
(N/a)^{(3)}=e+3cd+c^3,\quad
(\rho/\bar\rho)^{(3)}=-4e+48cd-64c^3,
\]
\[
(T/\bar T)^{(3)}=(a u^0)^{(3)}=-e+3cd-c^3,
\qquad u^i=0.
\tag{C15}
\]

At the specified finite `η_i`, evaluate the exact kernels `H,J,L` and their first derivatives. In particular

\[
L'(\eta_i)=\int_0^{\eta_i}ds\,\partial_{\eta_i}G_K(\eta_i,s)F_3(s),
\]

with no upper-end term. Insert the resulting `q_i,q'_i` into the exact extrinsic-curvature formula Q §2, retaining all cubic products, and obtain lapse/density/temperature/velocity from (C1). This supplies constrained Cauchy data through cubic order. No new arbitrary matching time enters (C14); positive-time uniqueness gives independence of a split at `η_i` when both value and derivative are retained.

## 9. Validation status and next gate

The user-supplied independent review accepted the quadratic result, not this cubic extension. Current tests check:

- the cubic determinant and fluid composition with symbolic noncommuting matrices;
- the displayed source against a separately expanded exact metric RHS on noncommuting time profiles with all nine momentum components symbolic;
- cancellation of the singular and constant source terms for general momenta;
- the cubic lapse boundary kernel and its exact `(k,q,−q)` restriction;
- a finite-time mixed kernel with two constant-dilation inputs.

The finite-time dilation check uses the independently inferred identities

\[
L(\eta;0,0,p)=(2-p\partial_p)^2H(\eta;p)
=4H-3p\partial_p H+p^2\partial_p^2H,
\quad
\chi^{(3)}(\eta;0,0,p)=(p\partial_p)^2c(\eta;p).
\]

They follow from `exp(2λ) H(η;exp(−λ)p)` for the metric and `c(η;exp(−λ)p)` for the lapse, rather than from integrating the proposed cubic source. The test checks these identities in the full sourced equation at arbitrary positive time and nonzero remaining input magnitude. It does not add a random zero-mode component to the primordial spectrum or certify observable gauge invariance.

The source comparison is not a fully symbolic arbitrary-matrix cubic identity test; its probe coverage is stated in the test. The full construction and the weighted nesting bounds require independent review. That review should focus on factorial polarization, use of the full quadratic Ricci/connection coefficients, seed subtraction, inner-integral derivative limits, and the constraint induction.

After that gate, the next physical task is the fixed-proper-time observer, transported tetrad, perturbed null ray, and temperature-surface intersection—not a premature contraction of metric kernels as if they were temperature kernels.
