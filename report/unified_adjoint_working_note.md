# Unified adjoint: first stationary-response results and a boundary-safe formulation

## Status

This is the first implementation of direction 1 in `future_cancellation_directions.md`.

Established here:

- an exact **polarized** cubic-response identity that eliminates explicit second- and third-order physical state responses in favor of first state and first adjoint responses;
- an explicit residual formulation of the complete temperature experiment on fixed parameter domains, including its nonlinear primordial data, physical endpoints, clock, Fermi–Walker frame and random sky normalization;
- the background field Green transpose, including the regular-origin covector `Lambda_c=9 lim(t Lambda_N)` and an explicit background harmonic adjoint;
- the exact normalization adjoint and its mean-free first variation;
- identification of the existing propagation adjoint as a block of the unified first adjoint response;
- a check against the established physical leading-gradient curvature map;
- an independent exact boundary/clock/redshift/normalization verification fixture.

**Follow-up:** the complete first field-adjoint source is now assembled in [`first_field_adjoint_source.md`](first_field_adjoint_source.md), with finite matrix coefficients and a specified regular Green reconstruction. Its implementation is `solution/unified_adjoint_source.py`, checked by `solution/check_first_field_adjoint_source.py`. The first spatial adjoint has an explicit observer trace and emission jump; the frame's cubic contribution is retained separately.

**Cubic continuation:** `stationary_cubic_assembly.md` now specifies the full stationary functional and its direct first-field-only measurement term, eliminating the remaining second ray and preserving the physical frame boundary.

**Evaluated gradient continuation:** `quadrupole_adjoint_gradient_expansion.md` gives the leading field-adjoint coefficient and its integrated fourth-gradient contribution in closed form. It establishes the complete cubic quadrupole comparison through four gradients at arbitrary symbolic input triples, using an independent Einstein–radiation curvature calculation.

**Exact angular continuation:** `exact_quadrupole_angular_source.md` eliminates the sky integrals and first-ray unknowns from the finite-wavelength quadrupole measurement source. Its frame covector is a commutator with the known linear quadrupole.

**Not established:** the unrestricted finite-wavelength Green/time reconstruction, all-multipole cubic comparison, or a full finite-wavelength cancellation/newly evaluated spectrum. The fixture below is not a radiation solution and is not used as physical evidence.

Checks for this initial stage: `solution/check_unified_adjoint.py`.

## 1. The stationary cubic identity

Let the input be the primordial field `g`, and write the complete independent residual system as

\[
 F(y,g)=0.
\]

Here `y` includes fields, paths, endpoints and measurement variables. Assume the linearized residual operator

\[
 A=F_y\big|_0
\]

has the inverse specified by the regular initial branch and physical boundary conditions. The inverse and its transpose must include those conditions; the inverse of a bulk differential expression alone is not sufficient.

For a real angular observable `J(y,g)`, define

\[
 L(y,\lambda,g)=J(y,g)+\langle\lambda,F(y,g)\rangle,
 \qquad A^*\lambda_0=-J_y\big|_0.
\tag{1}
\]

The pairing includes bulk and boundary residuals. For primordial directions `g_a`, let

\[
 Ay_a=-F_g[g_a],\qquad d_a=(y_a,g_a).
\]

Define `D_a` to differentiate in the **frozen** direction `d_a`, not to differentiate `y_a` again. With `L_0=J+<lambda_0,F>` and `lambda_0` held fixed, the first adjoint responses solve

\[
 A^*\lambda_a=-D_a(L_{0,y}).
\tag{2}
\]

Then the full third derivative is

\[
\boxed{
 J_{abc}
 =D_aD_bD_cL_0
 +\langle\lambda_a,D_bD_cF\rangle
 +\langle\lambda_b,D_aD_cF\rangle
 +\langle\lambda_c,D_aD_bF\rangle .}
\tag{3}
\]

There is **no `y_ab`, `y_abc` or second adjoint response** in this expression. The three single–pair terms have coefficient one. For a repeated direction, the formula is

\[
 J'''=D^3L_0+3\langle\lambda_1,D^2F\rangle.
\tag{4}
\]

If the normalized sky map is itself included as a state variable, the harmonic objective can be made linear in that state and independent of the input. Equation (3) simplifies to

\[
\boxed{
 J_{abc}=\langle\lambda_0,D_aD_bD_cF\rangle
 +\sum_{\rm cyc}\langle\lambda_a,D_bD_cF\rangle .}
\tag{5}
\]

This still includes every nonlinear measurement term: they are now derivatives of the measurement residuals, not explicit nonlinearities of `J`.

### Proof and relation to the previous construction

The second state response satisfies

\[
 Ay_{bc}=-D_bD_cF.
\]

Differentiate the objective and use (1) to eliminate `y_abc`. The remaining second-state terms are

\[
 J_{abc}=D_aD_bD_cL_0
 +\sum_{\rm cyc}\langle D_aL_{0,y},y_{bc}\rangle.
\]

Using (2) and the second-state equation gives

\[
 \langle D_aL_{0,y},y_{bc}\rangle
 =-\langle\lambda_a,Ay_{bc}\rangle
 =\langle\lambda_a,D_bD_cF\rangle,
\]

which proves (3). This is the cancellation of the explicit higher-state responses that follows from stationarity; it is **not a vanishing of `J_abc`**.

The polarized form matters for `J_3(k,q,-q)`: a repeated-direction formula alone would not fix the kernel combinatorics. The test script checks (3) for three independent input labels in a nonlinear system with a nonsymmetric Jacobian, a nonlinear objective and explicit input dependence, against recursive substitution in the actual implicit equations.

## 2. A complete residual system for this experiment

The following is a construction of the independent state system, not a claim that its adjoint has already been solved.

All physical labels `T_star,T_I,tau_o` and the initial material orientation are fixed. Write `zeta=g` as the input; directional/amplitude differentiation supplies the perturbative powers. Work on a fixed field-time interval containing the relevant perturbed events. Null paths and the observer frame are parameterized on fixed unit intervals so that varying endpoints cannot disappear from integral differentiation.

### 2.1 Independent radiation-field residuals

Eliminate the fluid variables using the exact temperature frame:

\[
 N=a_1t e^\nu,\qquad u=N^{-1}\partial_t,\qquad
 \gamma=e^{2g}e^H,\quad H=H^T,\quad\operatorname{tr}H=0,
 \quad S=\tfrac12e^{-H}(e^H)'.
\]

The equation of state, Euler equation, fluid normalization and continuity are already built into this parametrization. Induced vector/tensor metric components are retained in the five components of `H`; no extra vorticity seed is introduced.

Use the lapse equation and five independent orthonormal-frame spatial trace-free equations:

\[
 F_N=3\nu''-\Delta_\gamma\nu-2|D\nu|^2+\operatorname{tr}S^2
 +\frac3{t^2}(e^{-2\nu}-1)=0,
\]
\[
 F_H=2e^{H/2}\operatorname{TF}\left[
 S'+2(t^{-1}+\nu')S+R^i{}_j[\gamma]
 -2D^iD_j\nu+2D^i\nu D_j\nu
 \right]e^{-H/2}=0.
\tag{6}
\]

The conjugated residual is symmetric and trace-free, hence has five independent components. Do **not** add all Hamiltonian and momentum equations as additional unconstrained multiplier rows: that would introduce Bianchi redundancy into the Jacobian.

Instead select the report's regular branch, including its initial Hamiltonian matching:

\[
 H(0)=\nu(0)=0,\quad H,\nu\ \text{even and regular in }t,
\]
\[
 \lim_{t\to0}\frac{\nu}{t^2}
 =\frac1{18}e^{-2g}(2\nabla^2g+|\nabla g|^2).
\tag{7}
\]

The leading `H/t^2` coefficient is forced by (6), not imposed as an additional independent boundary row. The report's Cartesian/triangular uniqueness and constraint-propagation construction specify the constrained regular solution order by order. Nonlinear derivatives of (7) are essential: setting higher-order primordial metric/lapse data to zero would change the answer.

At the background the field block of `A` has exactly the known operators

\[
 L_K=3\partial_t^2+K^2-6/t^2,
\]
\[
 \partial_t^2H+\frac2t\partial_tH+\mathfrak B_KH
 +4\operatorname{TF}(KK)\nu,
\tag{8}
\]

with regular boundary prescriptions. Thus the transpose inverse is to be built from the already known regular Green map, not from a newly chosen homogeneous solution. Its background transpose, including the regular-singular-origin pairing, is given in §2.5 below. The complete first field-adjoint source is now assembled in `first_field_adjoint_source.md`.

### 2.2 Observer event and frame

Introduce scalar unknowns `t_I,t_o`. The two residuals are

\[
 F_{\rm init}=\ln t_I+\nu(t_I,0)-\ln\bar t_I=0,
\]
\[
 F_{\rm clock}=a_1\int_{t_I}^{t_o}t e^{\nu(t,0)}dt-\tau_o=0.
\tag{9}
\]

For the Fermi–Walker frame, set `t_F(u)=t_I+u(t_o-t_I)`, `0<=u<=1`, and impose

\[
 F_E=\frac{dE}{du}+(t_o-t_I)S(t_F(u),0)E=0,
 \qquad E(0)=\gamma(t_I,0)^{-1/2}.
\tag{10}
\]

The prescribed initial orientation fixes the rotational freedom. Orthonormality is preserved by this evolution; it is not added as redundant residual equations at every `u`.

### 2.3 Null paths on a fixed parameter domain

For each observed direction `n`, introduce `t_e(n)` and let

\[
 \Delta(n)=t_o-t_e(n),\qquad t(s,n)=t_o-s\Delta(n),\qquad 0\le s\le1.
\]

The ray starts at the observer at `s=0` and reaches emission at `s=1`. Use the coordinate direction `e^i` as a state variable and define

\[
 Q=\sigma_{ij}e^ie^j,\qquad \sigma=\tfrac12\partial_t\gamma,
\]
\[
 F_X=X_s+\Delta e=0,
\]
\[
 F_e=e_s+\Delta\left[-2Se-\Gamma(e,e)+eQ\right]=0,
\]
\[
 X(0,n)=0,\qquad e(0,n)=-E(1)n.
\tag{11}
\]

The initial null/unit-direction condition follows from the tetrad and is preserved by the ODE. It is not imposed as an extra redundant scalar equation throughout the ray.

Emission and redshift residuals are

\[
 F_{\rm em}(n)=\ln t_e+\nu(t_e,X(1,n))-\ln\bar t_\star=0,
\]
\[
 F_I(n)=I(n)-\Delta(n)\int_0^1Q(t(s,n),X(s,n),e(s,n))\,ds=0.
\tag{12}
\]

This formulation retains the moving endpoints in `Delta`, the evaluation times and the ray boundary conditions. The fixed `s` interval is not a fixed-emission approximation.

### 2.4 Exact measurement and harmonic objective

The direction-independent observer temperature has already canceled exactly. Introduce

\[
 F_Z(n)=Z(n)-e^{-I(n)}=0,\qquad
 F_M=M-\frac1{4\pi}\int Z(n)\,d\Omega=0,
\]
\[
 F_\Theta(n)=M(1+\Theta(n))-Z(n)=0.
\tag{13}
\]

For a real test weight `w(n)` with zero mean, take

\[
 J_w=\int w(n)\Theta(n)\,d\Omega.
\tag{14}
\]

Real and imaginary spherical-harmonic weights give the harmonic kernels. For `ell>=2`, the final measured-dipole projection does not affect (14), but the dipole remains inside (13) and all preceding nonlinear operations.

Equations (6)–(13), their initial/observer boundary residuals, and the regular branch define the augmented system. The generic identity (5) is applicable with the corresponding transpose inverse if all these blocks and their boundary functionals are retained. It is not permissible to substitute an unconstrained inverse of all ten Einstein equations, or to transpose a bulk operator while dropping its regular-origin/physical-endpoint terms.

### 2.5 Background field transpose and its regular-origin covector

This part of the field adjoint can be made explicit now. At a fixed Fourier momentum `K`, put `Q_K=KK-delta K^2/3` and let `B_K` be the report's Cartesian trace-free spatial operator. For residual forcing `(f_N,f_H)` and regular lapse coefficient `b_c`, the linear state equations are

\[
 L_K\nu=f_N,\qquad L_HH+4Q_K\nu=f_H,\qquad
 L_H=\partial_t^2+2\partial_t/t+B_K.
\]

Their regular inverse is

\[
 \nu(t)=b_c v_K(t)+\int_0^tG_N^+(t,s)f_N(s)ds,
 \quad v_K(t)=\frac9{K^2}u_1(Kt/\sqrt3),
\]
\[
 G_N^+(t,s)=\frac{u_2(ct)u_1(cs)-u_1(ct)u_2(cs)}{3c},\qquad c=K/\sqrt3,
\]
\[
 H(t)=\int_0^t\mathfrak G_K(t,s)[f_H(s)-4Q_K\nu(s)]ds.
\tag{8a}
\]

Here `G_N^+` is the inverse for `L_K nu=f_N`; it has the opposite sign to the report's Green coefficient multiplying a source `c_N` in `L_K nu=-c_N`. At `K=0`, `v_0=t^2` and `G_N^+=(t^3-s^3)/(9ts)`. The spatial Green operator `mathfrak G_K` is the existing entire Cartesian one, and `B_K` is self-adjoint on symmetric trace-free matrices under the trace pairing.

For ordinary bulk test sources `p_N,p_H` under the plain time/trace pairing, the transpose inverse is therefore

\[
 \Lambda_H(s)=\int_s^T\mathfrak G_K(t,s)[p_H(t)]dt,
 \qquad p_{\rm eff}(t)=p_N(t)-4Q_K:\Lambda_H(t),
\]
\[
 \Lambda_N(s)=\int_s^TG_N^+(t,s)p_{\rm eff}(t)dt,
 \qquad \Lambda_c=\int_0^Tv_K(t)p_{\rm eff}(t)dt.
\tag{8b}
\]

It satisfies the full duality identity

\[
 \int_0^T(p_N\nu+p_H:H)dt
 =\Lambda_c b_c+\int_0^T(\Lambda_N f_N+\Lambda_H:f_H)dt.
\tag{8c}
\]

In particular, the initial coefficient has **not** disappeared from the transpose. Since

\[
 \lim_{s\to0}9sG_N^+(t,s)=v_K(t),
\]

its covector is

\[
\boxed{\Lambda_c=9\lim_{s\to0}s\Lambda_N(s).}
\tag{8d}
\]

This also follows from the scalar Green boundary concomitant: `nu~b_c t^2` and `Lambda_N~L/t` contribute `9 b_c L` at the regular origin. Discarding that contribution is not justified by regularity of the physical state. The script checks (8c) with a coupled scalar/spatial field example and verifies that dropping the origin term changes the result. It also checks Cartesian self-adjointness, the Green equations and diagonal jumps. For smooth bulk sources the spatial formal transpose is `partial_t^2-2 partial_t/t+2/t^2+B_K`; using `L_H` itself as its transpose would be wrong.

For the actual harmonic objective there is an explicit background adjoint without distributional endpoint shorthand. Define

\[
 C_{K,ij}(t)=\int d\Omega\,w(n)\operatorname{TF}(n_i n_j)
 e^{iK\cdot n(o-t)}.
\]

The background observable derivative is

\[
 J_{w,y}[\delta H]=-\frac12\int_e^o C_K(t):\delta H'(t)dt.
\]

There is no direct lapse, ray or frame derivative at the background, because the background shear vanishes. The field multipliers in the sign convention (1) are

\[
\boxed{
 \Lambda_{H,0}(s)=\frac12\int_{\max(s,e)}^o
 \partial_t\mathfrak G_K(t,s)[C_K(t)]dt,}
\]
\[
 \Lambda_{N,0}(s)=-4\int_s^oG_N^+(t,s)Q_K:\Lambda_{H,0}(t)dt,
\]
\[
 \Lambda_{c,0}=-4\int_0^ov_K(t)Q_K:\Lambda_{H,0}(t)dt
 =9\lim_{s\to0}s\Lambda_{N,0}(s).
\tag{8e}
\]

Set these bulk kernels to zero above `o`. The angular tensor `C_K` is given by the already evaluated plane-wave sky integrals. The derivative on the spatial Green operator is taken **before** interchanging integrations: no observer/emission boundary is discarded. Equivalently, the source `-J_{w,y}` has bulk weight `-C_K'/2` on `[e,o]` and boundary weights `+C_K(o)/2` and `-C_K(e)/2`. This is why simply imposing homogeneous adjoint terminal data after integrating by parts would be incorrect.

The nonlinear initial matching also has explicit cubic contributions. For `F_c=lim(nu/t^2)-c[g]`, the stationary formula includes

\[
 -\langle\Lambda_{c,0},c_{abc}\rangle
 -\sum_{\rm cyc}\langle\Lambda_{c,a},c_{bc}\rangle,
\tag{8f}
\]

where the initial kernels are precisely the report's

\[
 c_{bc}=[2(k_b^2+k_c^2)-k_b\cdot k_c]/9,
\]
\[
 c_{abc}=\frac29[-2\sum_a k_a^2+\sum_{a<b}k_a\cdot k_b].
\]

These are the contributions of this particular initial residual, not all explicit primordial dependence of the bulk and measurement equations. The latter also enter (5). Equations (8a)–(8f) establish the background field transpose and its initial functional. The follow-up `first_field_adjoint_source.md` uses them to specify the full first field-adjoint weak source and its regular Green reconstruction.

## 3. A solved block: exact normalized-map adjoint

Use the residual signs in (12)–(13), pair sky residuals with `dOmega`, and suppose `int w dOmega=0`. Eliminating the algebraic measurement adjoints gives the multiplier of `F_I` exactly:

\[
\boxed{
 \lambda_I(n)=(1+\Theta(n))
 \left[w(n)-\frac{J_w}{4\pi}\right].}
\tag{15}
\]

In particular

\[
 \int\lambda_I\,d\Omega=0
\]

exactly, not only after ensemble averaging. The sign convention is such that the reduced sensitivity to an independent change of exponent is

\[
 \delta J_w=-\int\lambda_I\,\delta I\,d\Omega.
\]

At the background and first variation,

\[
 \lambda_{I,0}=w,\qquad
\boxed{\lambda_{I,a}=w\Theta_a-\frac1{4\pi}\int w\Theta_a\,d\Omega.}
\tag{16}
\]

The uncentered random monopole has canceled from the **first adjoint source itself**. This does not remove the normalization physics: it packages it into a manifestly mean-free dual weight before propagating the adjoint back through the ray and field equations.

To derive (15), stationarity in `Theta,M,Z` gives

\[
 \lambda_\Theta=-w/M,\quad \lambda_M=J_w/M,\quad
 \lambda_Z=(-w+J_w/(4\pi))/M,
\]

and stationarity in `I` then gives `lambda_I=-lambda_Z e^{-I}`. The exact finite weighted-measure analogue and its first variation are checked symbolically; this is an algebraic identity of the normalization operation, not a numerical sky quadrature.

## 4. Physical endpoint and clock sensitivities

These checks are most simply expressed after reducing the ray block at fixed coordinate observation time. They fix signs and boundary couplings that the full adjoint must reproduce.

For the logarithmic emission constraint, with `F= sigma_ij e^i e^j` and `dnu/dt` along the ray, its multiplier is

\[
 \lambda_{\rm em}
 =-\frac{\lambda_I F_e}{t_e^{-1}+(d\nu/dt)_e}.
\]

Consequently

\[
 \lambda_{{\rm em},0}=0,\qquad
 \lambda_{{\rm em},a}=-\bar t_\star\,wF_a(\bar t_\star).
\tag{17}
\]

Let `J_w^c(o)` denote the coordinate-time observable after ray and emission elimination, and before clock composition. The common proper-clock multiplier obeys

\[
 \lambda_{\rm clock}=-\frac{\partial_o J_w^c}{a_1o e^{\nu(o,0)}},\qquad
 \lambda_{{\rm clock},a}=-\frac{\partial_oJ_{w,a}^c}{a_1\bar t_o}.
\tag{18}
\]

The physical initial-event multiplier, in the logarithmic convention of (9), has first variation

\[
 \lambda_{{\rm init},a}=a_1\bar t_I^2\lambda_{{\rm clock},a}
 =-\frac{\bar t_I^2}{\bar t_o}\partial_oJ_{w,a}^c.
\tag{19}
\]

At first order the transported frame is `E_1(o)=-g(0) 1-H_1(o,0)/2`, independent of the initial frame time, and `J_{w,1}^c` is independent of `t_I`. This justifies the absence of an extra frame term in (19) at this order; higher variations of the frame are not being set to zero.

Together the first adjoints in (18)–(19) reproduce the known second-order-lapse clock functional

\[
 \mathcal C[m]=\frac{-\bar t_I^2m(\bar t_I,0)-\int_{\bar t_I}^{\bar t_o}t m(t,0)dt}{\bar t_o}.
\]

A rescaling of a residual rescales its multiplier: the test fixture uses the equivalent unlogged temperature roots, and checks the appropriately rescaled formulas rather than mixing conventions.

## 5. The existing ray adjoint is part of the unified first adjoint

At the background the ray state does not directly enter the exponent, since the background shear vanishes. Its background adjoint is zero. At first order, in the symmetric orthonormal-direction variables on a fixed `[e,o]`, use the notation of report §6.1:

\[
 r=\tfrac12\nabla h'_{nn},\quad s=-P h'n,\quad
 R(t)=\int_e^t r,\quad K(t)=\int_e^t(s-R),\quad P=1-nn^T.
\]

The first adjoint components are

\[
\boxed{\lambda_{X,a}=-wR,\qquad \lambda_{v,a}=-wPK.}
\tag{20}
\]

Indeed, the linear ray residuals are `X'-v-Hn/2-zeta n` and `v'-G`. Their adjoint equations at first order are

\[
 -wr-\lambda_{X,a}'=0,\qquad
 -ws-P\lambda_{X,a}-\lambda_{v,a}'=0,
\]

with the natural lower-end conditions. Equation (20) solves both. This is exactly the ray duality used previously to remove `X_2[H_2],v_2[H_2]`, now identified as a block of (2). The full field adjoint must additionally include the mean-free source (16), endpoint/clock multipliers, frame boundary conditions and the state-dependence of the field operator.

## 6. Physical leading-gradient check

The known leading-gradient observable is

\[
 \Theta_{\rm lead}=-\tau^2e^{-2\epsilon g}
 [\epsilon A-\epsilon^2B],
\]

where `A=(partial_i partial_j g)_TF n^i n^j` and `B=(partial_i g partial_j g)_TF n^i n^j` at the observer, and `tau^2=(bar t_o^2-bar t_star^2)/6`.

Introduce one implicit state `f` with residual `f exp(epsilon g)-1=0`, and write the objective as `-tau^2 f^2(epsilon A-epsilon^2 B)`. The stationary formula uses only `f_1=-g` and the first adjoint response, and gives

\[
 \Theta^{(3)}_{\rm lead}=-12\tau^2(g^2A+gB),
\]

in agreement with direct differentiation of the established curvature expression. Omitting the nonlinear-seed residual term loses part of the `g^2 A` contribution. This check confirms both that the method retains the physical cubic response and that nonlinear primordial boundary data cannot be treated as dispensable bookkeeping. It is a consistency check of the established limit, not a new finite-wavelength solution.

## 7. Verification fixture and scope

The check script also builds a 15-variable implicit measurement fixture containing:

- three direction-dependent temperature roots `t_e+epsilon c_i t_e^2=bar t_star`;
- a physical-form initial-temperature root and exact proper-clock integral for `N=t(1+epsilon c_o t)`;
- an integrated redshift driver with both linear and quadratic amplitude terms;
- the redshift exponential, its random mean and the normalized map.

The stationary cubic response is compared with independent exact amplitude expansions of the roots, clock equation, exponential and normalization. Clock, emission and normalization omission controls fail as required. All arithmetic is rational/symbolic; the amplitude series is the required perturbative expansion, not a time-series approximation or numerical evidence.

**This fixture is not a solution of the Einstein–radiation equations.** Its purpose is to test the full stationary-response plumbing and boundary algebra before implementing the field adjoint.

## 8. What to do next

The first response-elimination milestone is established algebraically. The follow-up `first_field_adjoint_source.md` has now assembled the field-operator variation and complete measurement Hessian, with emission/observer/clock functionals and the separate frame sector. It also specifies the regular Green reconstruction. The remaining substantive work is:

1. Insert the exact finite-wavelength quadrupole angular source from `exact_quadrupole_angular_source.md` into the first-adjoint Green pairing. Keep the direct third residual, direct cubic multi-phase measurement and nonlinear initial terms (8f).
2. Extend the independent physical comparison beyond the now-established arbitrary-triple cubic quadrupole through four spatial gradients. Distinguish exact gradient-coefficient identities from unrestricted finite-wavelength results.
3. Reduce the mixed-time angular integrals beyond the evaluated gradient sector, keeping the homogeneous lapse and the physical boundary traces from the first adjoint.
4. Determine whether the remaining Green integrals cancel or reduce. The first adjoint response can carry the same internal-momentum structure as a second-order state, so removing `y_2,y_3` does not by itself remove resonant dynamics or mixed-time phases.

The full finite-wavelength cancellation remains open. The first field source is now explicit, but its integral reconstruction and the stationary cubic functional have not yet replaced the report's evaluated/nested-integral answer.

## Reproduction and regression status

```bash
.venv/bin/python -B solution/check_unified_adjoint.py
.venv/bin/ruff check solution/check_unified_adjoint.py
.venv/bin/ruff format --check solution/check_unified_adjoint.py
.venv/bin/ty check solution/check_unified_adjoint.py
```

All pass. Targeted regressions also passed: `check_adjoint_exact.py`,
`check_propagation_adjoint.py`, `check_cartesian_green.py`,
`check_cubic_product_identity.py`, `check_observable_expansion.py`,
`check_orthonormal_map.py`, `check_report_identities.py`, and
`check_shell_cubic_sources.py --compare-engine`. The new check is included in
`solution/run_all.sh`; the full long-running suite has not been rerun for this work.
The working note is separate from the main PDF, which has not been changed by this step.
