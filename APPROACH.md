# Analytical radiation-temperature spectrum: proposed approach and validation contract

## 1. Review verdict and scope

The statement in `problem_statement.tex` calls for a nonlinear relativistic observable, not just a fluid density spectrum or a lensing correction to a linear temperature map. The appropriate route is to expand an **exact, physically labelled measurement** through cubic order, then perform Gaussian contractions. Deriving only second-order evolution cannot give the requested fourth-order power: the linear–cubic interference is equally necessary.

This document is a proposal with some exactly verified algebraic components, **not a completed derivation of the spectrum**. The subsequent [quadratic initial-data derivation](derivations/QUADRATIC_INITIAL_DATA.md) supplies explicit quadratic metric sources and a convergent finite-time retarded representation in a fixed gauge. Cubic dynamics and the resulting temperature kernels remain absent. Section 10 makes the completion gates explicit. A formal inverse differential operator, an unnamed “source”, or a numerical evolution is not an acceptable final substitute for those kernels.

The model can be specified consistently, but several choices must be made explicit before there is a unique conditional answer:

- the nonlinear Gaussian variable and its spatial coordinate/threading convention;
- the regular primordial boundary condition, including homogeneous and zero-momentum sectors;
- the physical initial event, material observer tag, initial triad, and elapsed proper time;
- the primordial covariance and its regularity, or a statement of which integrals remain formal/divergent;
- the prescription for selecting the emitting intersection if more than one exists.

These are specifications of the stochastic experiment, not licenses to discard terms. A spectrum can remain a functional of an arbitrary admissible primordial power spectrum and of the measurement labels; no unique numerical curve is determined by the statement.

## 2. Conventions and a definite experiment

Use signature `(-,+,+,+)`, `c=1`, and retain Newton's constant `G`. Let the unperturbed flat radiation solution have

\[
a(\eta)=A\eta,\quad \mathcal H=1/\eta,\quad
\bar\rho=\frac{3}{8\pi G A^2\eta^4},\quad
\bar T=(\bar\rho/\alpha)^{1/4},\qquad \eta>0.
\]

Fix the background normalization `A` and the coordinate length convention once. Distinguish:

1. the regular primordial boundary `η → 0+` used to define the growing solution;
2. a finite auxiliary Cauchy time `η_i > 0`;
3. an initial physical event on the observer's worldline;
4. the background emission time `η_e`, defined by `T̄(η_e)=T_star`;
5. the background observation time `η_o`, defined by the physical clock prescription.

Define the initial event to be the crossing of `T_phys=T_i`, with `T_i>T_star`, on a fixed material fluid tag, and observe after a specified proper time `Δτ`. Require `η_o>η_e>η(T_i)` in the background. Choose the initial triad as additional marked physical data. Transport it by the rule below. The background clock condition is

\[
\Delta\tau=\int_{\eta(T_i)}^{\eta_o}a(\eta)\,d\eta
=\frac A2\left[\eta_o^2-\eta(T_i)^2\right].
\]

The auxiliary `η_i` is not this physical clock event. The tag is held fixed across realizations; do not silently replace this ensemble by density-, volume-, or peak-weighted observers. A coordinate origin may represent the tag in a fixed gauge, but the tag must transform with the fields in a gauge comparison. The initial triad is a marked orthonormal frame, not a set of coordinate axes reconstructed in every gauge.

Choose the branch of the emission intersection continuously connected to the background one, and require a nonzero temperature derivative along the background null ray. Work perturbatively on the smooth interval before fluid shocks or breakdown of this intersection prescription. This is a domain-of-existence condition, not a short-wavelength, small-angle, or late-time approximation. Global existence for every realization of an unbounded Gaussian field is not asserted.

### Primordial random variable

A candidate precise convention is the regular primordial limit of the logarithmic local scale factor on uniform-density hypersurfaces, with comoving material spatial labels and the residual spatial coordinate freedom fixed:

\[
\zeta_{\rm p}(\mathbf x)
=\lim_{\eta\to0^+}\frac16\log\det(\gamma_{ij}/a^2),
\qquad
\gamma_{ij}/a^2=e^{2\zeta}\,[e^h]_{ij},\quad \mathrm{tr}\,h=0.
\]

The determinant is coordinate-dependent unless this spatial gauge has been fixed. The boundary construction now uses a marked Cartesian material chart and declares its limiting seed spatial metric to be `exp(2 ζ_p) δ_ij`. Conformal flatness is a physical restriction on primordial tensor data, not merely a gauge choice; see `derivations/QUADRATIC_INITIAL_DATA.md` §1. Uniform density is a slicing condition; material threading is a separate spatial condition. This is a definition at the regular primordial boundary, **not a claim that finite-time ζ is exactly conserved at finite momentum**.

Set

\[
\zeta_{\rm p}(\mathbf x;\epsilon)=\epsilon\,g(\mathbf x),
\qquad \langle g\rangle=0,
\]

with `g` Gaussian. All derivatives `∂_ε^r ζ_p|_0` for `r≥2` vanish. There are no independent primordial transverse vector, tensor, decaying, or higher-order scalar seeds. This does **not** set the nonlinear metric or fluid initial data to zero: even `exp(2εg)` has nonzero quadratic and cubic parts. Induced tensor and vector fields at finite time are mandatory.

Prove that regularity and these boundary data fix the required solution branch, rather than assuming it. Any residual homogeneous modes must be fixed by a stated background convention. Keep forced zero-total-momentum perturbations, including those sourced by `(k,-k)`; zero ensemble mean for the seed does not eliminate them.

### Fourier and angular conventions

\[
g(\mathbf x)=\int_{\mathbf k}g(\mathbf k)e^{i\mathbf k\cdot\mathbf x},
\quad \int_{\mathbf k}:=\int_{\mathbb R^3}\frac{d^3k}{(2\pi)^3},
\quad g(-\mathbf k)=g(\mathbf k)^*,
\]
\[
\langle g(\mathbf k)g(\mathbf k')\rangle
=(2\pi)^3\delta^3(\mathbf k+\mathbf k')P_g(k),
\quad \mathcal P_g(k)=\frac{k^3P_g(k)}{2\pi^2}.
\]

Use complex spherical harmonics with the Condon–Shortley phase,
`∫ dΩ Y_lm* Y_l'm' = δ_ll' δ_mm'` and `Y_00=1/sqrt(4π)`.
Every field is expanded as `X=Σ ε^r X^(r)/r!`. Symmetric kernels below represent `X^(r)` itself: there is **no additional `1/r!` in their convolution definition**. Products therefore generate binomial/multinomial coefficients.

## 3. Dynamical construction: constraints first, evolution second

Start from the full equations, not a scalar-only truncation:

\[
G_{\mu\nu}=8\pi G T_{\mu\nu},\quad
T_{\mu\nu}=\frac43\rho u_\mu u_\nu+\frac13\rho g_{\mu\nu},
\quad u_\mu u^\mu=-1,
\]
\[
u^\mu\nabla_\mu\rho+\frac43\rho\nabla_\mu u^\mu=0,
\quad
u^\mu\nabla_\mu u_\nu=-\frac14
(\delta_\nu{}^\mu+u_\nu u^\mu)\nabla_\mu\ln\rho.
\]

Equivalently the fluid acceleration is `a_ν=−(δ_ν^μ+u_ν u^μ)∇_μ ln T_phys`; radiation-fluid observers generally are **not geodesic**.

In an ADM decomposition let `n` be the future unit normal, `γ` the spatial metric, `K_ij=−(1/2) L_n γ_ij`, `E=T_μν n^μ n^ν`, and `S_i=−γ_i^μ n^ν T_μν`. The constraints are

\[
{}^{(3)}R+K^2-K_{ij}K^{ij}=16\pi G E,
\qquad D_j(K^{ij}-\gamma^{ij}K)=8\pi G S^i.
\]

**Construction task:** derive the regular scalar seed solution and its forced order-two/order-three corrections near the primordial boundary, solve both constraints coefficient by coefficient, and propagate that solution to finite `η_i`. No gradient expansion is allowed as the finite-time answer. A near-boundary expansion may select the exact solution branch only if convergence or a well-defined asymptotic boundary prescription is established. If that boundary construction fails, the primordial convention must be revised explicitly; arbitrary finite-time constraint data are not automatically the same model.

For finite-time calculation use the comoving-orthogonal, zero-shift gauge fixed in `derivations/QUADRATIC_INITIAL_DATA.md`, retaining a general six-component spatial metric. Poisson gauge, including transverse shift and transverse-traceless spatial metric, remains an option for an independent comparison, not a scalar-only replacement. Retain lapse, scalar curvature, velocity and density, transverse metric/fluid components, and tensor components. Tensor and vector constraints and propagation equations must be derived in that same convention. Barotropic flow with initially vanishing vorticity provides a useful conservation check; it does not imply vanishing vector metric perturbations or vanishing transverse parts of every coordinate velocity variable.

Let `F[X]=0` collect the exact equations and gauge conditions, expanded around the radiation background. At each order the same linearized operator acts on the new unknowns:

\[
L X^{(1)}=0,\qquad
L X^{(2)}=-F''[X^{(1)},X^{(1)}],
\]
\[
L X^{(3)}=-3F''[X^{(1)},X^{(2)}]
-F'''[X^{(1)},X^{(1)},X^{(1)}].
\]

These equations specify a generation algorithm, **not the final integral answer**. Publish every component of `L`, its constraints, its source kernels, the initial kernels at `η_i`, and the regular boundary prescription. Solve scalar/vector/tensor sectors with explicit retarded propagators. Treat total wavevector `K=0` directly instead of assigning arbitrary values to `K_i K_j/K²` projectors or inverse Laplacians.

### Exact linear building block

For the scalar metric convention
`ds²=a²[−(1+2Φ)dη²+(1−2Φ)δ_ij dx^i dx^j]` at first order,

\[
\Phi''+\frac4\eta\Phi'+\frac{k^2}{3}\Phi=0,
\quad y=\frac{k\eta}{\sqrt3},\quad
f(y)=3\frac{\sin y-y\cos y}{y^3},\quad f(0)=1.
\]

A second independent solution is `g_d(y)=(cos y+y sin y)/y³` (not the random field `g`). For `k>0`, a scalar retarded Green function is

\[
G_\Phi(\eta,s)=\mathbf1_{\eta\ge s}
\frac{f(ks/\sqrt3)g_d(k\eta/\sqrt3)
-g_d(ks/\sqrt3)f(k\eta/\sqrt3)}
{f(ks/\sqrt3)\,\partial_s g_d(ks/\sqrt3)
-\partial_s f(ks/\sqrt3)\,g_d(ks/\sqrt3)}.
\]

It has `G(s,s)=0` and `∂_ηG(s+,s)=1`. At `k=0`, use `1,η^−3` as the fundamental pair. The tensor homogeneous equation has operator `∂_η²+2∂_η/η+K²`, with solutions `sin(Kη)/(Kη)` and `cos(Kη)/(Kη)` for `K>0`, and `1,η^−1` at zero. Scalar, vector and tensor **sources and initial terms remain to be derived**. In the stated positive-log-volume ζ convention, the regular radiation scalar seed has `Φ_p=−2ζ_p/3`; verify the sign from the chosen constraint and slicing transformations rather than borrowing an incompatible curvature convention.

The eventual representation should have the explicit form

\[
X_A^{(r)}(\eta;\mathbf k_1,\ldots,\mathbf k_r)
=U_{AB}(\eta,\eta_i)I_B^{(r)}
+\int_{\eta_i}^{\eta}ds\,G_{AB}(\eta,s;\mathbf K)
S_B^{(r)}(s;\mathbf k_1,\ldots,\mathbf k_r),
\quad \mathbf K=\sum_j\mathbf k_j,
\]

with all entries and nesting specified, and any instantaneous constraint solution stated separately. Prove cancellation of dependence on the arbitrary finite matching time `η_i`. Never drop its boundary terms. One must not simply send `η_i→0` in an integral containing the singular decaying solution without proving the boundary limit exists.

## 4. Construct the observer and the ray before expanding the map

Solve the tagged fluid trajectory `dx_o^μ/dτ=u^μ`, including the perturbed location of its `T_i` event and its proper-time clock. For spatial tetrad vectors `e_(a)` impose

\[
g(e_{(a)},e_{(b)})=\delta_{ab},\quad u\cdot e_{(a)}=0,
\quad
\nabla_u e_{(a)}^\mu=(u^\mu a_\nu-a^\mu u_\nu)e_{(a)}^\nu.
\]

This is Fermi–Walker, not parallel transport along an accelerated worldline. Its initial orientation is part of the experiment.

Let `n^a n_a=1` point outward on the observed sky. Normalize the future-directed photon momentum at observation by

\[
p_o^\mu=u_o^\mu-n^a e_{(a)o}^\mu,\qquad E_o=1.
\]

Use affine parameter `λ_o=0` and integrate the future-directed momentum backward to negative `λ`:

\[
\frac{dx^\mu}{d\lambda}=p^\mu,\qquad
\frac{dp^\mu}{d\lambda}=-\Gamma^\mu_{\alpha\beta}p^\alpha p^\beta,
\qquad g_{\mu\nu}p^\mu p^\nu=0.
\]

The condition at emission is the scalar root

\[
F(\lambda,\epsilon)=\ln[T_{\rm phys}(x(\lambda,\epsilon),\epsilon)/T_\star]=0.
\]

All derivatives of `F` at fixed `λ` already include the perturbed ray. With `λ_e=λ_0+ελ_1+ε²λ_2/2+ε³λ_3/6`, and `F_r=∂_ε^rF|_0`, evaluated at `λ_0`, the exact implicit-function coefficients are

\[
\lambda_1=-F_1/F_0',
\]
\[
\lambda_2=-[F_2+2\lambda_1F_1'+\lambda_1^2F_0'']/F_0',
\]
\[
\lambda_3=-[F_3+3\lambda_1F_2'+3\lambda_2F_1'
+3\lambda_1^2F_1''+3\lambda_1\lambda_2F_0''
+\lambda_1^3F_0''']/F_0'.
\]

Use these rules also for shifted observer events with their appropriate defining scalars. Evaluate `E_e=−u_μ p^μ` at this **shifted** emission event. Temperature at emission is exactly `T_star`; there is no additional independent intrinsic-temperature fluctuation to add on this surface. Its effects enter through the shifted surface and emitter velocity.

Derive the ray/clock/tetrad coefficient equations by differentiating their exact ODEs through order three. Include evaluation shifts, time delay, transverse displacement, deflection–deflection and deflection–source products, nonlinear redshift, and both endpoints. This avoids a potentially incomplete hand-picked list of “Sachs–Wolfe”, “Doppler” and “lensing” terms. Every propagation integral must retain its actual finite limits; expansion of moving limits supplies boundary terms.

## 5. Normalization: an exactly checked part of the solution

Let `T_0=T_star a(η_e)/a(η_o)` be the background temperature, and define

\[
T_{\rm obs}/T_0=1+\epsilon t_1+\frac{\epsilon^2}{2}t_2+
\frac{\epsilon^3}{6}t_3+O(\epsilon^4),
\qquad b_r=\overline{t_r}:=\frac1{4\pi}\int d\Omega\,t_r.
\]

These `b_r` are random measured sky monopoles, **not ensemble means**. Direct series division gives

\[
\Theta^{(1)}=t_1-b_1,
\]
\[
\Theta^{(2)}=t_2-b_2-2b_1(t_1-b_1),
\]
\[
\Theta^{(3)}=t_3-b_3-3b_1(t_2-b_2)
+(6b_1^2-3b_2)(t_1-b_1).
\]

If redshift is accumulated as
`q=ln(T_obs/T_0)=εq_1+ε²q_2/2+ε³q_3/6`, then

\[
t_1=q_1,\qquad t_2=q_2+q_1^2,\qquad
t_3=q_3+3q_1q_2+q_1^3.
\]

Only after these operations apply

\[
\Theta_{\rm clean}(\hat n)=\Theta(\hat n)
-\sum_{m=-1}^{1}Y_{1m}(\hat n)\int d\Omega'\,
Y_{1m}^*(\hat n')\Theta(\hat n').
\]

No Lorentz deboosting is implied by dipole subtraction. The measured monopole is already zero; for the cleaned map `a_00=a_1m=0` at every order. All spectrum formulas below concern `ℓ≥2`. Products involving observer monopoles/dipoles must remain before this final projection; their products can contribute to higher multipoles.

## 6. Gaussian power: exact reduction once the observable kernels exist

Construct symmetric, generally complex harmonic kernels `K_lm^(r)` from the **complete normalized map**, including the observer:

\[
a_{\ell m}^{(r)}=\int_{\mathbf k_1}\cdots\int_{\mathbf k_r}
K_{\ell m}^{(r)}(\mathbf k_1,\ldots,\mathbf k_r)
\prod_{j=1}^r g(\mathbf k_j).
\]

These kernels already contain `∫dΩ Y_lm*`, observer-position phases and all convolution factors. If derived from an output-momentum kernel, integrate its `(2π)^3δ³(K−Σk_j)` with `d³K/(2π)^3`; do not leave an extra delta function in this expression.

Wick's theorem and connected subtraction give

\[
C_\ell^{(11)}=\frac1{2\ell+1}\sum_m\int_{\mathbf k}
P_g(k)|K_{\ell m}^{(1)}(\mathbf k)|^2,
\]
\[
C_\ell^{(22)}=\frac2{2\ell+1}\sum_m\int_{\mathbf k}\int_{\mathbf q}
P_g(k)P_g(q)|K_{\ell m}^{(2)}(\mathbf k,\mathbf q)|^2,
\]
\[
C_\ell^{(13)}=\frac3{2\ell+1}\sum_m\int_{\mathbf k}\int_{\mathbf q}
P_g(k)P_g(q)K_{\ell m}^{(1)}(\mathbf k)
K_{\ell m}^{(3)*}(\mathbf k,\mathbf q,-\mathbf q).
\]

In `22`, the one contraction internal to each quadratic factor is exactly the product of means and is subtracted. In `13`, all three pairings survive: the internal pair in the cubic factor is **not** disconnected as a covariance of the two composite observables. Normal-ordering the cubic map would incorrectly remove this physical term. Maintain internal momenta with exact total momentum zero.

Consequently

\[
C_\ell=\epsilon^2 C_\ell^{(11)}+
\epsilon^4\left[\frac14 C_\ell^{(22)}+
\frac13\operatorname{Re}C_\ell^{(13)}\right]+O(\epsilon^6).
\]

The `1/4` and `1/3` arise from Taylor factorials and the two cross terms, separately from Wick factors `2` and `3`. Odd orders vanish for the zero-mean Gaussian seed with fixed, unweighted observer labels. `O(ε^6)` also assumes continuation of the formal perturbation series: a bound on the map remainder alone would not prove this ensemble remainder. The fourth-order **correction** need not be positive because of `13` interference.

As an angular normalization check, define the linear transfer by
`K_lm^(1)(k)=4π i^l Δ_l(k)Y_lm*(khat)` at an origin observer. The addition theorem gives

\[
C_\ell^{(11)}=4\pi\int_0^\infty d\ln k\,
\mathcal P_g(k)|\Delta_\ell(k)|^2.
\]

The definition fixes the phase convention; the actual `Δ_l` still has to be derived with the stated clock and endpoint conditions.

## 7. Monochromatic spectrum and convergence

Only after constructing `K^(1..3)`, choose fixed primordial curvature variance `σ_g²`:

\[
\mathcal P_g(k)=\sigma_g^2\delta(\ln(k/k_0)),\quad
P_g(k)=\frac{2\pi^2\sigma_g^2}{k_0^2}\delta(k-k_0),
\qquad k_0>0.
\]

The physical seed variance is `ε²σ_g²`. Hold `σ_g²`, the background and the physical measurement labels fixed when varying `k_0`. Holding density variance instead defines a different experiment. For a sufficiently regular test function `f`,

\[
\int_{\mathbf k}P_g(k)f(\mathbf k)
=\frac{\sigma_g^2}{4\pi}\int d\Omega_{\hat k}\,f(k_0\hat k).
\]

Thus `11` is a two-angle integral and `22`/`13` are four-angle integrals before symmetry reductions. Do not replace the generated quadratic/cubic output momentum spectrum by another shell at `k_0`; sums of shell wavevectors are not confined to that shell. In particular the `13` contraction samples `K^(3)(k,q,−q)` exactly. An almost-everywhere-defined generic kernel does not suffice: this exact contraction restriction must itself be proved to exist, separately from integrability near zero partial sums.

For the general spectrum, first formulate exact results for explicitly declared smooth spectra supported on `0<k_min≤k≤k_max<∞`. This defines a controlled class of inputs, not an approximation to every spectrum. Even compact input support does **not** remove singularities at zero *sum* momentum. Show that the full gauge-invariant integrand has an integrable extension there, or retain an explicit regulator and demonstrate a regulator-independent limit. Prove ultraviolet/infrared convergence separately before extending the input class.

Use dominated convergence (or an appropriate proved distributional convergence result) for every interchange of primordial boundary limits, angular integrals, time integrals and shell/soft/high-frequency limits. Supply the majorant or other precise theorem hypotheses. Expanding individual gauge-dependent terms and then integrating is not sufficient. A delta-shell substitution also requires existence of the restricted full integrand on the shell; products of singular distributions must not be assumed defined.

Smooth-spectrum convergence also does not imply shell convergence. For illustration only, a full integrand behaving like `|k+q|^(−p)` has generic three-dimensional local threshold `p<3`, but the equal-radius antipodal shell has measure `θ dθ` and threshold `p<2`. Neither singularity is asserted for the actual observable. No large- or small-`k_0` scaling is established in this proposal.

## 8. Full analytical validation matrix

| Gate | Required evidence | Common error it detects |
|---|---|---|
| Background and linear seed | Exact radiation Friedmann/Euler solution; analytic regular scalar transfer, constraint-derived density and velocity; correct curvature sign | Wrong background, growing/decaying mixture |
| Initial data | Hamiltonian and all momentum constraints vanish through `ε³`; no forbidden independent integration constants | Setting all second/third-order data to zero |
| Dynamical closure | Substitute component kernels into every Einstein and fluid equation at orders 1–3; prove constraint propagation by Bianchi identities | Missing vector/tensor or source terms |
| Algebraic geometry | `u²=−1`, tetrad orthonormality and `u·e_a=0` through cubic | Missing lapse/velocity/observer products |
| Photon propagation | Null and all geodesic residuals vanish through cubic, including observed direction and endpoint initial conditions | Incomplete post-Born or time-delay terms |
| Physical endpoints | Emission temperature, initial clock temperature, and fixed elapsed proper time hold through cubic | Fixed-coordinate rather than physical event |
| Measurement composition | Exact log/exponential and sky-division identities; sky mean zero at each order; dipole removed last | Dropped monopole/dipole products |
| Stochastic normalization | Independently enumerate Wick pairings and derive Fourier/angular factors; retain `13` internal pair | Wrong factorials, erroneous normal ordering |
| Gauge invariance | Pull back fields, rays, tags, scalar surfaces and tetrad together by a general diffeomorphism, then differentiate exact scalar measurement through cubic; verify component kernels in a second gauge or explicit general-gauge transformation | Fixed-coordinate gauge comparison |
| Pure gauge tests | A coordinate-transformed homogeneous solution gives zero anisotropy with transformed physical labels; constant spatial dilation handled with the momentum convention transformed too | Spurious long-mode observer effects |
| Matching time and boundaries | `d/dη_i` of complete observable vanishes; all endpoint/initial terms accounted for | Lost initial-time or moving-limit terms |
| Rotations/reality | Real sky relation `a_l,-m=(−1)^m a_lm*`; isotropic unweighted ensemble covariance diagonal in `l,m` | Harmonic conjugation or observer-selection errors |
| Degenerate kinematics | Explicit treatment of zero total momentum, collinear and squeezed configurations, finite times | Illegitimate projector inverses |
| Integral legality | Convergence and limit-exchange proofs; exact lower limits and regular branch boundary prescriptions | Hidden approximation or divergent integral |

Do not demand that a physical finite-wavelength long mode vanish just because a coordinate-generated adiabatic mode does. Physical tidal effects are not a gauge failure. Isotropy checks assume the marked triad introduces no realization-dependent preferred external field.

Exact symbolic algebra is appropriate for coefficient extraction, tensor contractions, residuals, and elementary/special-function identities. Numerical quadrature, ODE solvers, fitting, stochastic sampling and floating-point agreement are not acceptance evidence under this statement. A CAS result must identify its assumptions, excluded singular surfaces, and boundary data; simplification to zero is not a convergence theorem.

## 9. Implementation and independent-review sequence

1. **Freeze conventions.** Write the gauge, seed, Fourier/factorial conventions, background zero-mode prescription and measurement labels. Review these before constructing kernels.
2. **Build an exact perturbation engine.** Use rational/symbolic coefficient arithmetic through order three, explicit spatial indices and momentum-conserving convolution symmetrization. Keep constraints separate from evolution. Begin with general component equations; scalar/vector/tensor projections are derived, not silently imposed on every order.
3. **Finish and validate initial data.** Produce order-one/two/three kernels at `η_i`, with regular primordial boundary proof and constraint residuals. This is a blocking milestone.
4. **Finish dynamical kernels.** Publish explicit sources and propagators, preserving finite time boundaries and general momenta. Verify every residual.
5. **Finish measurement kernels.** Solve observer, clock, triad, ray and surface recursions, assemble the scalar redshift, divide by the measured sky mean, then project the dipole.
6. **Perform a genuinely independent check.** A second derivation should use a different gauge or a general-gauge transformation and a separate implementation of endpoint/observer composition. Mere re-use of the same generated kernels is not independent validation.
7. **Contract and reduce.** Use the verified Wick formulas; preserve every interference term. Reduce angular/time integrals only with proofs. Then specialize the primordial spectrum.
8. **Publish a proof ledger.** For each kernel and integral list the defining expression, domain, boundary prescription, evaluated/unevaluated status, symbolic residual test, independent review result and remaining caveats. Native subagent attempts failed, but the user subsequently supplied an independent review of the original plan and algebra: **sound approach with stated gaps**. See `reviews/INDEPENDENT_REVIEW_01.md` for attribution and disposition. This review does not cover the newly derived quadratic construction.

Suggested initial tooling is Python with SymPy for exact identities and a literate mathematical derivation. If expression size warrants a tensor CAS, select and validate one before porting the generator. No large software stack or numerical simulator is necessary for the current planning and algebra-validation stage.

## 10. Delivered evidence versus unresolved work

`validation/check_algebra.py` checks normalization, log-redshift/endpoint composition, actual angular averages, Gaussian moments and pairings, scalar/tensor Green functions, and Fourier/shell measures. `validation/check_initial_data.py` adds checks of the exact fluid reduction, seed geometry, linear equations/constraints, quadratic source expansion, boundary coefficients, zero-total-momentum response, retarded matrix kernel and constraint-propagation energy identity. Neither constructs the observable temperature kernels `K^(2)` or `K^(3)`.

`derivations/QUADRATIC_INITIAL_DATA.md` gives the new modewise quadratic construction and its analytical constraint and convergence arguments. Its unevaluated finite-time integrals are explicit. Independent acceptance of this new work is pending.

Unresolved before a full solution can be claimed:

- independent review of the new quadratic boundary construction and its exact reduction;
- extension to cubic dynamical sources, constrained initial kernels and explicit finite-time representations;
- the complete cubic observer/clock/tetrad/ray/endpoint kernels;
- the independent component-level gauge check;
- explicit final finite-time kernels inside the spectrum integrals;
- convergence, shell restriction and any asymptotic-limit proofs.

The analytical route is defensible, but “full validation” of the requested spectrum remains conditional on passing **all** these gates. The current deliverable validates selected algebra and specifies how to validate the eventual solution; it must not be presented as the completed solution itself.
