# Complete first field-adjoint source

## Status and result

This continues `unified_adjoint_working_note.md`. The previously unassembled first **field**-adjoint source is now specified explicitly. It contains:

1. the off-shell variation of the independent lapse/spatial field equations;
2. the complete normalized measurement Hessian, with the independent probe's ray eliminated;
3. both emission terms, the observer endpoint, the physical clock integral and its initial-event term;
4. a precise treatment of the Fermi–Walker sector: it does not source the first field adjoint in the rotation chart, but its first adjoint and nonzero cubic pairing are retained separately.

Everything below is at arbitrary finite time and general momenta. No shell restriction or time expansion is used. All spatial coefficient operations are finite polynomials in the momenta; homogeneous probes are included without projectors. The regular Green map gives a fully specified integral reconstruction of the first field adjoint in §6.

**Continuation:** `stationary_cubic_assembly.md` now assembles the arbitrary-triple stationary functional, including the direct first-field-only measurement contribution and its explicit frame boundary.

**Limits:** the reconstruction integrals have not been evaluated into new closed angular/time kernels, and an independent full Einstein–radiation cubic comparison remains open. In particular, this does not establish a full finite-wavelength cancellation or a newly evaluated spectrum.

Implementation: `solution/unified_adjoint_source.py`.
Exact checks: `solution/check_first_field_adjoint_source.py`.

## 1. What is being differentiated

Let `y=(nu,H)` be the field state, and eliminate the null paths, physical endpoints, clock, frame and normalization only for the purpose of defining the exact reduced observable `J_red(y,g)`. This is an exact block elimination, not a fixed-endpoint or fixed-clock approximation.

Let `E(y,g)=0` be the independent field residuals (6)–(7) of the working note. The field component of the full unified adjoint satisfies

\[
 E_y^*\Lambda=-J_{{\rm red},y}.
\]

For a fixed first primordial direction `a`, the first adjoint therefore satisfies

\[
 A^*\Lambda_a=\mathcal S_a,\qquad
 \mathcal S_a[b]=-D_aJ_{{\rm red},y}[b]
 -\langle\Lambda_0,D_aE_y[b]\rangle.
\tag{1}
\]

The independent test `b` here is a **state probe**, not a second primordial input direction:

\[
 a:\quad (\zeta,\nu,H)=(z,a(t),h(t))e^{ip\cdot X},
 \qquad b:\quad(\delta\zeta,\delta\nu,\delta H)=(0,m(t),j(t))e^{iq\cdot X}.
\tag{2}
\]

Thus the primordial input is held fixed when differentiating in `b`. Using a second physical seed with `delta zeta=1` in place of this probe would compute a different derivative.

Write `P=p+q`, `d=p.q`, and `G=2z 1+h`. All tensor probes and field adjoints are symmetric and trace-free, paired by `A:B=tr(A^T B)`. Primes on field amplitudes are partial time derivatives. Fourier duals are indexed by the state coefficient they pair with, with the usual opposite-momentum convention implicit in the real-space pairing.

The background field adjoint entering the field term is `Lambda_0(t,P)`, whereas the first-adjoint reconstruction uses the **probe** operator at momentum `q`. These are different momenta. For the shell cubic contraction the probe can be homogeneous; it must not be deleted.

For a physical single plane-wave input with `p!=0`, the known first state is

\[
 a=z(\cos x-\sin x/x),\qquad
 h=-6z(2\sin x/x-\cos x-1)
 (\hat p\hat p^T-1/3),\qquad x=pt/\sqrt3.
\]

The coefficient identities below hold even off shell, and their checks use arbitrary functions instead of these special solutions.

## 2. The field-operator variation, off shell

The residual normalization is

\[
 E_N=3\nu''-\Delta_\gamma\nu-2|D\nu|^2+\operatorname{tr}S^2
 +3(e^{-2\nu}-1)/t^2,
\]
\[
 E_H=2e^{H/2}E_Te^{-H/2},
\]

where `E_T` is the report's exact mixed spatial trace-free equation. In particular, the spatial residual has been converted to five symmetric orthonormal-frame components **before** differentiation.

The mixed lapse residual is

\[
 \boxed{C_N=A_Nm-a\,p^TjP+\tfrac12h':j',}
\]
\[
 A_N=-q^ThP-z(2q^2-d)+4da+12a/t^2.
\tag{3}
\]

For the spatial residual, define

\[
 \mathcal U=\frac12\{m[P(Gq)^T+(Gq)P^T-dG]
 +a[P(jp)^T+(jp)P^T-dj]\}.
\]

Then

\[
 \boxed{C_H=2a'j'+2m'h'
 +2\operatorname{TF}[\mathcal R_{p,q}(G,j)-2\mathcal U
 -2am(pq^T+qp^T)].}
\tag{4}
\]

There are no second time derivatives of the probe. The disappearance of quadratic noncommuting kinetic products in (4) follows from the orthonormal conversion; it is not an assumption that the probe solves its linear field equation.

### Explicit finite construction of the curvature operator

The notation `R_{p,q}` is not an unevaluated field equation. It is this finite bilinear polynomial operation:

\[
 \Gamma_k(B)^l{}_{ij}=\frac i2(k_iB_{lj}+k_jB_{li}-k_lB_{ij}),
\]
\[
 R_k(B)=\tfrac12[k^2B+(\operatorname{tr}B)kk^T-k(Bk)^T-(Bk)k^T].
\]

Set `Gamma_a=Gamma_p(G)`, `Gamma_b=Gamma_q(j)` and

\[
 B_{ab}=\tfrac12(Gj+jG),\qquad
 \Gamma_{ab}^l=\Gamma_P(B_{ab})^l-G_{lr}\Gamma_b^r-j_{lr}\Gamma_a^r.
\]

Repeated spatial indices are summed. The lower-index mixed Ricci tensor is

\[
 (R_{ab})_{ij}=iP_l\Gamma_{ab}^l{}_{ij}-iP_j\Gamma_{ab}^l{}_{il}
 +\Gamma_a^l{}_{lr}\Gamma_b^r{}_{ij}
 +\Gamma_b^l{}_{lr}\Gamma_a^r{}_{ij}
 -\Gamma_a^l{}_{jr}\Gamma_b^r{}_{il}
 -\Gamma_b^l{}_{jr}\Gamma_a^r{}_{il}.
\]

The required orthonormal-frame coefficient is

\[
 \mathcal R_{p,q}(G,j)=R_{ab}
 -\tfrac12\{G,R_q(j)\}-\tfrac12\{j,R_p(G)\}.
\tag{5}
\]

This construction also makes the absence of momentum denominators transparent.

### Explicit weak transpose coefficients

Put `ell=Lambda_{N,0}(t,P)` and `L=Lambda_{H,0}(t,P)`. Define `R^*` by the finite trace pairing of (5). Then

\[
 \ell C_N+L:C_H
 =V_\nu m+D_\nu m'+V_H:j+D_H:j',
\tag{6}
\]

with

\[
 V_\nu=\ell A_N-4P^TLGq+2d(L:G)-8a p^TLq,
 \qquad D_\nu=2L:h',
\]
\[
 V_H=-\ell a\operatorname{TF}[\operatorname{sym}(pP^T)]
 +2\mathcal R^*_{p,q,G}[L]
 -2a\operatorname{TF}[p(LP)^T+(LP)p^T-dL],
\]
\[
 D_H=\tfrac12\ell h'+2a'L.
\tag{7}
\]

For completeness, `R^*` can be obtained with five evaluations and no tensor inversion ambiguity. Use basis matrices `diag(1,0,-1)`, `diag(0,1,-1)`, and the three symmetric off-diagonal unit pairs. If `c_i` are the five values of `L:R(G,E_i)`, its dual has diagonal `((2c_0-c_1)/3,(2c_1-c_0)/3,-(c_0+c_1)/3)` and off-diagonal entries `c_2/2,c_3/2,c_4/2`. This is implemented by `stf_dual`.

The contribution to `S_a` is the **negative** of (6). If integrating by parts, retain

\[
 -\int_0^o(V_\nu m+D_\nu m'+V_H:j+D_H:j')dt
 =\int_0^o[(-V_\nu+D_\nu')m+(-V_H+D_H'):j]dt
 -[D_\nu m+D_H:j]_0^o.
\tag{8}
\]

The weak form is preferable: it avoids introducing distributional endpoint shorthand. For the actual background adjoint and regular physical first direction, the boundary in (8) vanishes using their specific endpoint/asymptotic data. That fact does **not** remove the separate measurement boundary functionals below.

### Two exact structural checks

For a homogeneous scalar probe `q=0,j=0`,

\[
 C_N=12am/t^2,\qquad C_H=2m'h'.
\]

For a constant primordial shift `p=0,h=a=0`,

\[
 C_N=-2zq^2m,\qquad
 C_H=-2z\mathfrak B_qj-8zQ_qm,
 \quad Q_q=qq^T-q^2/3.
\]

The latter is precisely the finite-wavelength spatial-operator dilation `e^{-2z}`. Both are checked without a small-momentum expansion.

## 3. Measurement Hessian with no independent probe ray

Let `e=bar t_star`, `o=bar t_o`, denote the physical initial time by `iota=bar t_I`, and set `X_0(t,n)=n(o-t)`. The symbol `i` in Fourier phases is the imaginary unit, not this initial time. Write

\[
 \phi_q(t,n)=e^{iq\cdot n(o-t)},\quad N_n=nn^T-1/3,\quad
 F_a(t,n)=\tfrac12\phi_p h'_{nn},\quad U_a=\int_e^oF_a dt.
\]

The test ray is not needed in the final source. Only the known first direction's **fixed-coordinate-observer** responses enter. Their observer data are zero; do not add a clock displacement `D_a n` to `X_a`, since clock composition is included separately in (11)–(14):

\[
 X_a'=v_a+\phi_p(zn+hn/2),\qquad
 v_a'=\phi_p\{(1-nn^T)(h'n/2+izp)
 +\tfrac i2[p h_{nn}-(p\cdot n)hn]\},
\]
\[
 X_a(o)=v_a(o)=0.
\tag{9}
\]

These are total pulled-back derivatives. Set

\[
 r_a=\tfrac12\nabla(\phi_p h'_{nn}),\quad
 s_a=-(1-nn^T)\phi_p h'n,\quad
 R_a(t)=\int_e^t r_a,\quad K_a(t)=\int_e^t(s_a-R_a).
\tag{10}
\]

Here `grad` acts on the plane-wave field before evaluating on the ray. The use of the full longitudinal `R_a`, rather than its transverse projection, is essential.

The normalized first map and first exponent multiplier are

\[
 \Theta_a=-U_a+\langle U_a\rangle_\Omega,\quad
 J_a=\int w\Theta_a\,d\Omega,\quad
 \lambda_{I,a}=w\Theta_a-J_a/(4\pi).
\]

Define the known first clock shift and its coordinate-observer sensitivity

\[
 D_a=\mathcal C[a]
 =-\frac{\iota^2a(\iota,0)+\int_\iota^ot a(t,0)dt}{o},\qquad
 A_a=\int w\,\partial_oU_a\,d\Omega=-\partial_oJ_a.
\tag{11}
\]

`A_a` is already sky-integrated. It is not an additional direction-dependent integrand.

### Ray bulk coefficients

Let `P_n=1-nn^T` and `sym(uv^T)=(uv^T+vu^T)/2`. Define trace-free matrices

\[
 E_1=\tfrac i2(q\cdot X_a)N_n
 -\operatorname{TF}\operatorname{sym}(v_an^T)
 -\tfrac12\operatorname{TF}\operatorname{sym}((P_nK_a)n^T),
\]
\[
 E_0=-\tfrac12\operatorname{TF}\operatorname{sym}(R_an^T)
 -\tfrac i2(q\cdot K_a)N_n
 +\tfrac i2(q\cdot n)\operatorname{TF}\operatorname{sym}(K_an^T).
\tag{12}
\]

The complete bulk coefficients of the **negative measurement Hessian** are

\[
 M_H=\phi_q wE_0,
\]
\[
 M_{H'}=\phi_q\left[wE_1+\tfrac12\lambda_{I,a}N_n
 +\tfrac i2wD_a(q\cdot n)N_n\right].
\tag{13}
\]

The three contributions in the second line are respectively ray propagation, normalized-map response, and the bulk part of the physical observer-clock composition.

To see why no test ray remains, its original bulk terms are

\[
 \phi_q[\tfrac i2(q\cdot X_a)j'_{nn}-v_a\cdot j'n]
 +r_a\cdot X_b+s_a\cdot v_b.
\]

The exact duality identity replaces the last two terms by

\[
 -\phi_q(R_a\cdot jn/2+K_a\cdot G_j)
 +\frac d{dt}(R_a\cdot X_b+K_a\cdot v_b),
\]

where `G_j=P_n j'n/2+i(q j_nn-(q.n)jn)/2`. The total derivative vanishes only after using `R_a(e)=K_a(e)=0` and `X_b(o)=v_b(o)=0`. The remaining terms give (12), not a loss of either physical metric endpoint.

### Explicit emission, observer and clock functionals

Let `a_e=phi_p(e) a(e)` and `F_{a,e}=F_a(e,n)`. Then

\[
 S_{{\rm meas},a}[m,j]
 =\int_e^o dt\int d\Omega\,[M_H:j+M_{H'}:j']
\]
\[
 \quad+\int d\Omega\,\phi_q(e)
 \left[e wF_{a,e}m(e)+\frac e2w a_e N_n:j'(e)\right]
 +\frac{D_a}{2}\int d\Omega\,wN_n:j'(o)
\]
\[
 \quad-\frac{A_a}{o}\left[\iota^2m(\iota)+\int_\iota^ot m(t)dt\right].
\tag{14}
\]

The last line is evaluated on the observer worldline, so it has **no ray phase** `phi_q`. It is not integrated over sky again.

The observer term follows from the exact identity

\[
 \partial_oU_j=\tfrac12j'_{nn}(o)
 +\int_e^o\tfrac i2(q\cdot n)\phi_q j'_{nn}dt.
\]

Neither the observer term nor the physical initial-event term `iota^2 m(iota)` can be dropped. The two emission terms likewise have distinct origins: the lapse probe shifts the physical emission time, while the first lapse shifts the lower limit of the test metric's redshift integral.

### Independent Hessian check

At second mixed order, the complete exponent is

\[
 I_{ab}=B_{ab}+eF_{a,e}m_e+eF_{b,e}a_e
 +D_a\partial_oU_b+\mathcal C[m]\partial_oU_a.
\]

Direct expansion of the exact random-mean normalized map gives

\[
 J_{ab}=-\int w I_{ab}\,d\Omega
 +\int w(U_a-\langle U_a\rangle)(U_b-\langle U_b\rangle)\,d\Omega.
\]

The second term equals `-int lambda_{I,a} U_b dOmega`. Thus (14) is exactly `-D_a J_red,y[b]`, with the bilinear ray term represented by (12). The check script verifies the raw coordinate-driver identity and the root/clock/normalization algebra independently, not by assuming the desired Hessian formula. The finite weighted-sky normalization check is an exact algebra test, not a quadrature approximation.

## 4. The assembled first field source

Combining (6) and (14) gives the explicit result

\[
\boxed{
 \mathcal S_a[m,j]=S_{{\rm meas},a}[m,j]
 -\int_0^o[V_\nu m+D_\nu m'+V_H:j+D_H:j']dt.}
\tag{15}
\]

This is a linear functional of the independent field probe, with every coefficient determined by first physical fields, their first ray responses, the background field adjoint, and first sky averages. It contains no second physical field or independent probe-ray response.

The executable assembly is `FirstFieldAdjointSource.action(m,j)`. Construct it from `field_coefficients(...)`, `measurement_source(...)`, the three background times and a callable for the full `dOmega` integral. It returns the specified definite time/sky integrals with the sign in (15). The field coefficients and common-clock terms are already sky-contracted and are deliberately **not** passed through that sky integral again. The check script exercises the assembly API on regular probes and causal Green probes, including this distinction.

The initial residual `E_c=lim(nu/t^2)-c[g]` has a state derivative independent of `g`, so it adds no separate `D_a E_{c,y}` term to (15). This does **not** set nonlinear primordial data to zero. Their contributions to the cubic formula remain

\[
 -\langle\Lambda_{c,0},c_{abc}\rangle
 -\sum_{\rm cyc}\langle\Lambda_{c,a},c_{bc}\rangle,
\]

as in (8f) of the earlier note. The new initial covector `Lambda_{c,a}` is obtained from (15), not guessed from a homogeneous terminal condition.

## 5. Why the frame is absent here but not from the cubic map

Use the exact change of frame state

\[
 R=e^{g(0)}e^{H(t,0)/2}E_{\rm FW},\qquad R(t_I)=1.
\]

Its exact rotation equation is `R'=Omega R`. At the background the driver has no linear term, and its mixed quadratic term is

\[
 \Omega_{ab}=\tfrac18([h_a,h_b']+[h_b,h_a']).
\tag{16}
\]

Consequently `R_a=0`. A quadratic frame rotation acts on a vanishing background sky, so it does not enter the measurement Hessian or the first **field** source (15). This is an exact consequence of the state transformation and order counting, not an instruction to discard frame physics.

The first **frame** adjoint is generally nonzero. Define

\[
 \mathcal K_a=\int d\Omega\,
 \operatorname{skew}[w(\nabla_{S^2}\Theta_a)n^T].
\]

With frame residual `R'-Omega R`, its first multiplier is

\[
 \Lambda_{R,a}(t)=-\mathcal K_a,\qquad \iota\le t\le o.
\tag{17}
\]

The cubic stationary-response identity retains the frame contribution

\[
 \sum_{\rm cyc}\int_\iota^o\mathcal K_a:\Omega_{bc}\,dt.
\tag{18}
\]

For the equal-frequency shell the quadratic local rotation driver vanishes, but no such restriction was made in deriving or checking (16)–(18). The arbitrary-triple cubic map must keep this sector until the applicable restriction or final observable operation is established.

## 6. Regular Green reconstruction and nonhomogeneous adjoint boundaries

Equation (15) gives a direct reconstruction without integrating its probe derivatives by parts. At momentum `q`, let `G_N^+`, `v_q` and `mathfrak G_q` be the regular Green objects from the earlier note. Extend the following probes by zero for `t<s`:

\[
 \Phi_H(t;s,B)=\left(0,\mathfrak G_q(t,s)[B]\right),
\]
\[
 \Phi_N(t;s)=\left(G_N^+(t,s),
 -4\int_s^t\mathfrak G_q(t,u)[Q_q]G_N^+(u,s)du\right),
\]
\[
 \Phi_c(t)=\left(v_q(t),
 -4\int_0^t\mathfrak G_q(t,u)[Q_q]v_q(u)du\right).
\tag{19}
\]

Then the complete first field adjoint and its initial-data covector are specified by

\[
\boxed{
 \Lambda_{H,a}(s):B=\mathcal S_a[\Phi_H(\cdot;s,B)],\quad
 \Lambda_{N,a}(s)=\mathcal S_a[\Phi_N(\cdot;s)],\quad
 \Lambda_{c,a}=\mathcal S_a[\Phi_c].}
\tag{20}
\]

Five trace-free test matrices reconstruct `Lambda_H,a` explicitly via the same finite dual-basis operation used in §2. Derivatives in (15) act on the causal Green probes before any integral exchange; endpoint traces use the corresponding one-sided Volterra prescription. No new arbitrary homogeneous constants are allowed.

The source is no longer unspecified: substituting (7), (12)–(14) into (20) defines only integrals of known first-order objects and regular Green kernels. This does not mean those integrals have been evaluated or that their mixed-time phases have canceled.

### New boundary diagnostics

The observer point evaluation `j'(o)` in (14) forces the generally nonzero left terminal trace

\[
\boxed{\Lambda_{H,a}(o^-)=\frac{D_a}{2}\int w(n)N_n\,d\Omega.}
\tag{21}
\]

The emission derivative evaluation similarly forces the jump

\[
\boxed{
 \Lambda_{H,a}(e^-)-\Lambda_{H,a}(e^+)
 =\frac e2\int d\Omega\,w(n)a_e(n)\phi_q(e,n)N_n.}
\tag{22}
\]

These follow from the unit diagonal derivative of the spatial Green operator. Imposing zero terminal data and a continuous first spatial adjoint would lose physical clock/emission terms. For a pure harmonic weight, (21) vanishes except at `ell=2`, but no such cancellation may be assumed for (22).

The regular-origin covector remains

\[
\boxed{\Lambda_{c,a}=9\lim_{s\to0}s\Lambda_{N,a}(s).}
\tag{23}
\]

For the physical first direction, the field weak-source coefficients have `V_nu=O(1/t)`, `D_nu=O(t^2)`, `V_H=O(t)` and `D_H=O(1)` near the origin. Their action on the regular coefficient probe is integrable. The measurement endpoints lie at strictly positive times. Moreover `9s Phi_N(t;s) -> Phi_c(t)`, including first time derivatives, follows from (19) and the regular lapse Green-origin limit. The coupled spatial response therefore has the same limit too; it is not being omitted in (23). These bounds allow the weak functional to inherit the limit.

The script verifies (20)–(23), including derivative endpoint evaluations and the clock point source, in an exact homogeneous-momentum operator example. It is not a physical approximation or a new universe model.

## 7. Verification and what remains

Run

```bash
.venv/bin/python -B solution/check_first_field_adjoint_source.py
```

The script checks:

- (3)–(4) against `engine.equations` at general momenta and arbitrary time functions, including off-shell orthonormal conversion;
- all transpose coefficients (7) against contraction of the independent engine result;
- homogeneous lapse probes and the finite-wavelength constant-dilation spatial identity;
- raw coordinate-driver jets against the no-probe-ray form, including its exact dual total derivative;
- both emission terms, clock initial data, observer endpoint and normalized physical Hessian;
- the noncommuting quadratic frame generator and its first terminal covector;
- the complete source assembly API on regular and causal Green probes, without a second sky integration of the common clock;
- the Green solution, nonzero observer trace, emission jump and regular-origin covector with exact symbolic integrals.

There is no quadrature, finite time-series evidence or numerical physical model in these checks. The source construction is a formal finite-order identity with the same regular/finite-time prescriptions as the report; no arbitrary-spectrum UV/IR convergence or late-time integral exchange follows.

The new script and ten existing targeted regressions passed, including `check_unified_adjoint.py`, the propagation/Cartesian Green/observable expansion checks, `check_shell_cubic_sources.py --compare-engine`, and `foundations/check_fw_conformal.py`. Ruff lint/format, project type checking, shell syntax and whitespace checks also passed. The full long-running repository suite was not rerun. These are working research notes; the main report PDF was not changed by this step.

**Next:** specialize and evaluate the polarized stationary functional now specified in `stationary_cubic_assembly.md`, using this first field source and retaining the direct third residual and initial-data terms. Its complete measurement sector has an independent off-shell cubic check, but comparison of the full Einstein–radiation observable and any further physical cancellation or angular-kernel reduction remain open.
