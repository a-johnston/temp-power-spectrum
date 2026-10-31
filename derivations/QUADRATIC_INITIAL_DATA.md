# Quadratic regular initial data in a fixed comoving gauge

## Status and scope

This is the next analytical milestone after the user-supplied independent review. It gives a fixed gauge, exact reduced field equations, an explicit regular linear solution, an explicit general-momentum quadratic source, and a convergent retarded representation of the quadratic metric. It includes the exactly zero **total** momentum sector without scalar/vector/tensor projectors. It also derives the first nontrivial primordial time coefficient geometrically.

The representation is an **unevaluated analytical result**, not numerical evolution. The elementary time integrands and every nesting limit are specified below. The infinite series defining the retarded matrix kernel is absolutely convergent on each finite time interval. The source and integral construction are checked at selected algebraic seams by `validation/check_initial_data.py`; the full constraint argument uses the ADM propagation identities. This new derivation has not yet received an independent review.

This does **not** give cubic initial data, the cubic observed map, or a power spectrum. No finite-time answer is obtained by truncating the near-boundary time expansion.

## 1. Fix the gauge and the primordial class

Use the signature, background, Fourier measure and amplitude factorials of `APPROACH.md`. Write the seed as `ζ(x)=ε g(x)`; here `g` is the Gaussian field, never a metric component. Mark a Cartesian chart on the material-label space, with fixed origin, orientation and length convention. Declare the seed metric in this chart to be

\[
b_{ij}(x)=e^{2\zeta(x)}\delta_{ij}.
\]

Conformal flatness is a restriction on primordial tensor data, not a gauge argument for removing a physical tensor seed. Relabelling the material chart transforms the seed and its stochastic law together; it is not additional freedom to change the boundary data while leaving the experiment fixed.

A perfect barotropic radiation fluid with zero initial vorticity remains irrotational on its smooth evolution domain. Locally use its orthogonal hypersurfaces and transport the marked spatial labels with the fluid. Choose zero shift and write

\[
ds^2=a^2(\eta)\{-e^{2\chi(\eta,x)}d\eta^2
+q_{ij}(\eta,x)dx^idx^j\},\qquad
u^\mu=(a e^\chi)^{-1}(1,0,0,0).
\]

The spatial metric `q` is a general symmetric matrix; it is **not** assumed conformally flat at positive time, nor is its transverse or trace-free part discarded. The zero shift is a comoving-orthogonal gauge choice, not an assertion that scalar perturbations cannot generate vector/tensor metric contributions in other gauges.

Energy conservation and Euler's equation integrate to

\[
\rho\,(\sqrt{\det\gamma})^{4/3}=C(x),\qquad
\partial_i\ln(N\rho^{1/4})=0,
\quad \gamma=a^2q,\quad N=a e^\chi.
\]

Fix the homogeneous time normalization by `N ρ^(1/4)=a ρ̄^(1/4)` and put the common primordial boundary at `η=0`. The right-hand side is constant in the radiation background, so this removes time reparametrization freedom (the remaining additive constant is fixed by the boundary). The regular seed fixes `C(x)=ρ̄ a^4 e^(4ζ)`. Thus, **exactly**,

\[
\boxed{\chi=\frac16\ln\det q-\zeta,\quad
\rho=\bar\rho e^{-4\chi},\quad T_{\rm phys}=\bar T e^{-\chi}.}
\tag{1}
\]

There is no independently specifiable fluid density or velocity seed left in this gauge. These formulas satisfy the full fluid equations and `u²=−1`, not just their first-order truncations.

The regular perturbative boundary class is

\[
q(\eta,x)=b(x)+O(\eta^2),\qquad q'=O(\eta),\qquad
\chi=O(\eta^2),
\tag{2}
\]

with even analytic coefficient kernels near zero for each fixed finite set of input momenta. Exclude free singular modes, free seed anisotropies and independent higher-order corrections to `ζ`. This is a finite-order, modewise perturbative class; it is not a global existence assertion for every Gaussian realization. Spatial integral interchanges require separate bounds. Since `δρ/ρ̄=O(η²)`, transforming to uniform density changes the relative time and rescaled metric only at `O(η²)`. Hence the limit of the uniform-density logarithmic volume in the marked threading is the prescribed `ζ`.

This fixes the gauge for this milestone and permits subsequent measurement work directly in it. The previously suggested Poisson gauge was an option, not a requirement. A Poisson-gauge comparison, if made, must transform all physical and primordial labels.

## 2. Exact component equations before expansion

Let `H=1/η`, `M=q^(-1)q'/2`, `m=tr M=3χ'`, and use `D`, `R_ij`, `R` for the connection and curvatures of `q`. Spatial indices on these quantities are raised with `q`. The mixed extrinsic curvature is

\[
K^i{}_j=-\frac1{a e^\chi}(H\delta^i{}_j+M^i{}_j).
\]

The mixed ADM evolution equation for radiation is

\[
(K^i{}_j)'=-D_\gamma^iD_jN
+N\left({}^{(3)}R^i{}_j+K K^i{}_j-8\pi G\frac\rho3\delta^i{}_j\right).
\]

In terms of `q,χ` this becomes the fully explicit metric equation

\[
\begin{split}
q''={}&q'q^{-1}q'-2(H+\chi')q'\\
&+2[H^2(e^{-2\chi}-1)-2H\chi']q\\
&+2e^{2\chi}\left[\partial_i\partial_j\chi
-\Gamma^a{}_{ij}[q]\partial_a\chi
+\partial_i\chi\partial_j\chi-R_{ij}[q]\right].
\end{split}
\tag{3}
\]

The bracket on the last line is a covariant two-tensor; products on the other lines are ordinary matrix products. Explicitly,

\[
\Gamma^a{}_{ij}=\tfrac12q^{ab}(\partial_iq_{bj}+\partial_jq_{bi}-\partial_bq_{ij}),
\quad R_{ij}=\partial_a\Gamma^a{}_{ij}-\partial_j\Gamma^a{}_{ia}
+\Gamma^a{}_{ab}\Gamma^b{}_{ij}-\Gamma^a{}_{jb}\Gamma^b{}_{ia}.
\]

The exact constraints, in a form with background scale factors removed, are

\[
e^{2\chi}R+6H^2+4Hm+m^2-\operatorname{tr}M^2
=6H^2e^{-2\chi},
\tag{4}
\]
\[
D_jM^j{}_i-\partial_i m-(\partial_j\chi)M^j{}_i
+(2H+m)\partial_i\chi=0.
\tag{5}
\]

No inverse Laplacian or total-momentum projector occurs in these equations.

## 3. First nontrivial primordial coefficient: nonlinear in the seed

Set `q=b+η²c+O(η⁴)` and `C=b^(-1)c`, `t=tr C`. Equation (1) gives `χ=η²t/6+O(η⁴)`. Substitution into the leading part of (3), equivalently the mixed evolution equation, gives

\[
3C^i{}_j+t\delta^i{}_j=-R^i{}_j[b],\quad
 t=-R[b]/6,
\]
\[
\boxed{c_{ij}=-\frac13R_{ij}[b]+\frac1{18}R[b]b_{ij},\quad
\chi=-\frac{\eta^2}{36}R[b]+O(\eta^4).}
\tag{6}
\]

The Hamiltonian constraint at this order is `R[b]+6t=0`. The leading momentum constraint is
`D_j C^j_i−(2/3)∂_i t=0`, which follows from the contracted three-dimensional Bianchi identity `D_j R^j_i=(1/2)∂_i R` with the displayed coefficients. Thus all three momentum components, including the spatial connection terms, are retained.

For the conformal seed,

\[
R_{ij}[b]=-\zeta_{,ij}-\delta_{ij}\Delta\zeta
+\zeta_{,i}\zeta_{,j}-\delta_{ij}|\nabla\zeta|^2,
\quad R[b]=e^{-2\zeta}[-4\Delta\zeta-2|\nabla\zeta|^2],
\]
\[
c_{ij}=\frac13\zeta_{,ij}+\frac19\delta_{ij}\Delta\zeta
-\frac13\zeta_{,i}\zeta_{,j}+\frac29\delta_{ij}|\nabla\zeta|^2.
\tag{7}
\]

This coefficient is exact in `ζ`; only the time expansion is truncated here. In particular, with `q=I+εh+ε²j/2+…`,

\[
j_{ij}=4g^2\delta_{ij}
+\eta^2\left[-\frac23 g_{,i}g_{,j}
+\frac49\delta_{ij}|\nabla g|^2\right]+O(\eta^4).
\tag{8}
\]

Its symmetric quadratic Fourier kernel is

\[
J_{ij}(\eta;k,p)=4\delta_{ij}+\eta^2 B_{ij}(k,p)+O(\eta^4),
\quad B_{ij}=\frac13(k_ip_j+p_ik_j)-\frac49\delta_{ij}k\cdot p.
\tag{9}
\]

Also

\[
\chi^{(2)}(\eta;k,p)=\frac{\eta^2}{9}
[2(k^2+p^2)-k\cdot p]+O(\eta^4),
\quad (\rho/\bar\rho)^{(2)}=-4\chi^{(2)}+O(\eta^4).
\tag{10}
\]

At `p=−k`, `B_ij=−(2/3)k_i k_j+(4/9)δ_ij k²`, with trace `2k²/3` and trace-free part `−(2/3)(k_i k_j−δ_ij k²/3)`. The lapse kernel is `5η²k²/9`, and the density-contrast kernel is `−20η²k²/9`. These forced homogeneous responses are **not** removed or absorbed into a redefined background. A stochastic isotropic average of the trace-free part is a different operation from deleting the realization-dependent tensor.

Equations (6)–(10) alone are not finite-time initial data. The exact finite-time extension follows next.

## 4. Exact regular linear fields

For an input momentum `k`, put `y=|k|η/√3` and define the entire even functions

\[
s(y)=\sin y/y,\quad
 e(y)=[2s(y)-1-\cos y]/y^2,\quad s(0)=1,\ e(0)=1/6.
\]

The linear metric and lapse kernels are

\[
H_{ij}(\eta;k)=2s(y)\delta_{ij}-2\eta^2 e(y)k_i k_j,
\qquad c(\eta;k)=\cos y-s(y).
\tag{11}
\]

Here uppercase `H_ij` is a metric kernel, while the scalar `H=1/η` remains the conformal Hubble factor. In the source formulas below use `h(x)=∫_k H(k)g(k)e^(ik·x)` and `c(x)=χ^(1)(x)`.

For example write `h_ij=2U δ_ij+2∂_i∂_j V`, with
`U=s(y)g`, `V=η²e(y)g`. Then `U'=H c`, `c=U+ΔV/3−g`, and `U''+2H U'+(k²/3)U=0`. These identities verify the momentum constraint; the Hamiltonian constraint and all linear metric equations follow by substitution. The density and lapse have no singular zero-input-momentum limit. Time derivatives in subsequent formulas mean derivatives of (11), not independently specified initial velocities.

## 5. Explicit quadratic component source

Define the following operators on a symmetric spatial matrix `h`:

\[
\Gamma_1^a{}_{ij}[h]=\tfrac12(\partial_i h_{aj}+\partial_j h_{ai}-\partial_a h_{ij}),
\]
\[
R_{1,ij}[h]=\tfrac12(\partial_a\partial_i h_{aj}
+\partial_a\partial_j h_{ai}-\Delta h_{ij}-\partial_i\partial_j\operatorname{tr}h),
\]
\[
\Gamma_Q^a{}_{ij}[h,h]=-2h_{ab}\Gamma_1^b{}_{ij}[h],
\]
\[
R_{Q,ij}[h,h]=\partial_a\Gamma_Q^a{}_{ij}-\partial_j\Gamma_Q^a{}_{ia}
+2(\Gamma_1^a{}_{ab}\Gamma_1^b{}_{ij}
-\Gamma_1^a{}_{jb}\Gamma_1^b{}_{ia}).
\tag{12}
\]

`R[q]^(2)=R_1[j]+R_Q[h,h]`; these definitions fix the factorials. With

\[
v=\operatorname{tr}(h^2)/6,\qquad
\chi^{(2)}=\operatorname{tr}j/6-v,
\]

the order-two metric equation is

\[
\mathcal L j=\mathcal S[h,c],
\tag{13}
\]
\[
\begin{split}
\mathcal S_{ij}={}&2(h'h')_{ij}-4c'h'_{ij}
+\{4(H^2v+Hv')+8H^2c^2\}\delta_{ij}\\
&-8(H^2c+Hc')h_{ij}-2\partial_i\partial_j v
-4\Gamma_1^a{}_{ij}[h]\partial_a c
+4\partial_i c\partial_j c-2R_{Q,ij}[h,h]\\
&+8c(\partial_i\partial_j c-R_{1,ij}[h]).
\end{split}
\tag{14}
\]

This is a fully specified finite component expression, not an unnamed nonlinear source. Each repeated spatial index runs from 1 to 3, and all contractions in (12)–(14) use `δ`. Every `h,c` is the explicit linear field (11).

For total output momentum `K=k+p`, the linear operator is

\[
\mathcal L_K X=X''+\frac2\eta X'
+\frac{2}{3\eta}\delta\operatorname{tr}X'
+\frac{2}{3\eta^2}\delta\operatorname{tr}X+\mathcal D_K X,
\]
\[
(\mathcal D_K X)_{ij}=K^2X_{ij}-K_iK_aX_{aj}-K_jK_aX_{ai}
+\frac43K_iK_j\operatorname{tr}X.
\tag{15}
\]

For clarity, the symmetric Fourier kernel `S(η;k,p)` of (14) is defined by this **finite algebraic** prescription: replace `h` by `z H(η;k)+w H(η;p)` and `c` by `z c(η;k)+w c(η;p)`, use `∂_i z=ik_i z`, `∂_i w=ip_i w`, and take one half the coefficient of `zw`. Time derivatives act only on the displayed time functions. Products of two linear factors are therefore symmetrized with weight `1/2`. This fixes every derivative and convolution without an implicit PDE solve. `validation/initial_data.py` implements precisely this bilinear source with explicit index sums.

Explicitly, the quadratic metric field is

\[
j_{ij}(\eta,x)=\int_{\mathbf k}\int_{\mathbf p}
J_{ij}(\eta;k,p)g(\mathbf k)g(\mathbf p)e^{i(\mathbf k+\mathbf p)\cdot x},
\]

with the measure from `APPROACH.md` and no extra factorial in this convolution. The initial metric kernel is `J_0=4δ`. Write `J=J_0+Y`. The source for `Y` is

\[
F(\eta;k,p)=S(\eta;k,p)-\frac8{\eta^2}\delta
-4K^2\delta-8K_iK_j.
\tag{16}
\]

The apparent `η^(−2)` singularity cancels: `h→2gδ`, `v→2g²`, `c=O(η²)`, `h'=O(η)`. All remaining terms are even entire time functions with removable values at zero, for fixed finite input momenta. In particular

\[
F(0;k,p)=2(k_ip_j+p_ik_j)-4\delta_{ij}k\cdot p
=6B_{ij}+2\delta_{ij}\operatorname{tr}B.
\tag{17}
\]

Never integrate the two singular terms in (16) separately from zero.

## 6. Exact retarded representation, including K=0

Define trace and trace-free projectors on the six-dimensional space of symmetric matrices, **not** momentum-space longitudinal/transverse projectors:

\[
P_t X=\delta\operatorname{tr}X/3,\qquad P_f X=X-P_tX.
\]

For `0<s≤η`, define the elementary matrix kernel

\[
G_0(\eta,s)=s(1-s/\eta)P_f
+\frac{s^2}{\eta}(1-s/\eta)P_t.
\tag{18}
\]

Extend it by zero for `s>η`, and by its continuous value zero at `s=0`. It inverts `L_0=L_K−D_K` with zero regular seed: the trace-free scalar operator is `∂²+2∂/η`; the trace scalar operator is `∂²+4∂/η+2/η²`. Each Green function has zero coincidence value and unit first-derivative jump.

The full retarded kernel is the explicitly defined matrix-valued function

\[
\begin{split}
G_K(\eta,s)={}&G_0(\eta,s)\\
&+\sum_{n=1}^{\infty}(-1)^n
\int_s^\eta dt_1\int_s^{t_1}dt_2\cdots\int_s^{t_{n-1}}dt_n\,
G_0(\eta,t_1)\mathcal D_KG_0(t_1,t_2)\mathcal D_K\cdots
\mathcal D_KG_0(t_n,s).
\end{split}
\tag{19}
\]

All products are compositions of explicitly given six-dimensional linear maps, acting on the tensor to their right. This definition needs no unsupplied fundamental matrix. The exact quadratic metric is

\[
\boxed{J_{ij}(\eta;k,p)=4\delta_{ij}
+\left[\int_0^\eta ds\,G_K(\eta,s)F(s;k,p)\right]_{ij}.}
\tag{20}
\]

### Convergence and differentiation

In the Frobenius norm the projectors in (18) are orthogonal, and
`||G_0(η,s)||≤η−s`. Let `d=||D_K||`, which is finite and bounded by `7K²`. The norm of the `n`th term in (19) is at most

\[
d^n\frac{(\eta-s)^{2n+1}}{(2n+1)!}.
\]

If `M_T=sup_(0≤s≤T)||F(s;k,p)||`, the corresponding term in (20) is bounded by

\[
M_T d^n\frac{T^{2n+2}}{(2n+2)!}.
\]

Thus both series converge absolutely and uniformly on finite time triangles. The bound proves the nesting exchanges used here. For first derivatives, `|∂_η G_f|=s²/η²≤1` and `|∂_ηG_t|=|−s²/η²+2s³/η³|≤1`; the differentiated series is likewise uniformly convergent. The Volterra equation then gives the second derivative for positive `η`, verifies `L_KY=F`, and supplies the zero-time regular extension from the even source expansion. These bounds concern fixed momenta, not unbounded primordial integrations or the final observable.

Equivalently, if `F=Σ_m η^m F_m` and `Y=Σ_(n≥2)η^n Y_n`,

\[
Y_{m+2}=\left[\frac{P_f}{(m+2)(m+3)}
+\frac{P_t}{(m+3)(m+4)}\right]
(F_m-\mathcal D_KY_m),\quad Y_0=Y_1=0.
\tag{21}
\]

Odd coefficients vanish. Equation (19), rather than a truncated use of (21), supplies the finite-time answer. There are no positive-order resonances: the trace-free indicial factors are `n(n+1)`, and the trace factors are `(n+1)(n+2)`. Free singular modes and the free constant trace-free seed are excluded by (2). Uniqueness follows either from (21) in the analytic class or from the homogeneous Volterra equation with its convergent iteration bound.

At **exactly `K=0`**, `D_K=0`, so (20) reduces to

\[
J(\eta;k,-k)=4\delta+
\int_0^\eta ds\left[s(1-s/\eta)P_f
+\frac{s^2}{\eta}(1-s/\eta)P_t\right]F(s;k,-k).
\tag{22}
\]

The source is still (11)–(16) evaluated at `p=−k`, not discarded. This includes forced homogeneous trace and trace-free components. The only unevaluated operations here are explicitly defined time integrals (and the convergent series at nonzero `K`). No arbitrary direction at zero momentum is assigned. Further elementary or special-function reductions of these integrals have not been exhausted; no claim that (22) is an evaluated closed form is made.

## 7. Constraints and finite-time data

Equations (1) and (3) are the full fluid equations and ADM evolution in this gauge. To show that their regular perturbative solution is constrained at each positive time, let the unscaled ADM constraint residuals be

\[
\mathcal H_c={}^{(3)}R+K^2-K^i{}_jK^j{}_i-16\pi G\rho,
\quad \mathcal M_i=D_jK^j{}_i-D_iK.
\]

The Bianchi/ADM propagation identities at zero shift and conserved matter are

\[
\mathcal H_c'=-2N D_i\mathcal M^i+2NK\mathcal H_c
-4\mathcal M^iD_iN,
\quad
\mathcal M_i'=-\tfrac12N D_i\mathcal H_c+NK\mathcal M_i
-\mathcal H_cD_iN.
\tag{23}
\]

At perturbative order `r`, once all lower-order residuals vanish, put
`h_c=a² H_c^(r)` and `m_i=a M_i^(r)`. Their Fourier equations about radiation are

\[
h_c'=-4h_c/\eta-2iK_i m_i,\qquad
m_i'=-2m_i/\eta-\tfrac12iK_i h_c.
\tag{24}
\]

The regular boundary class gives bounded `h_c,m_i` as `η→0+`: cancellation of the universal background `η^(−2)` terms uses `χ=O(η²)` and `q'=O(η)`. Define the nonnegative Hermitian norm
`E=|h_c|²+4Σ_i|m_i|²`. Then, exactly,

\[
E'=-\frac8\eta|h_c|^2-\frac{16}\eta\sum_i|m_i|^2
\le-\frac4\eta E.
\]

Consequently `η⁴ E(η)≤s⁴ E(s)→0` as `s→0+`. All residuals vanish. Apply this first at order one, then at order two. This argument works also at `K=0` and does not assume an inverse Laplacian. The geometric substitution in §3 independently checks the leading nontrivial coefficient and signs of this construction. The test suite also builds `J` through `η⁴` from (21) and substitutes it into separately assembled full matrix-jet constraints: the Hamiltonian residual vanishes through `η²`, and every momentum component through `η³`, for symbolic general momentum pairs modulo rotation. This finite-coefficient check supplements rather than replaces the all-time argument.

At any specified positive `η_i`, the finite-time initial data are now explicitly

\[
q_i=I+\epsilon h_i+\frac{\epsilon^2}{2}j_i,
\quad q'_i=\epsilon h'_i+\frac{\epsilon^2}{2}j'_i,
\]

where `h` is (11), `j` is the double convolution of (20), and

\[
J'(\eta_i;k,p)=\int_0^{\eta_i}ds\,\partial_{\eta_i}G_K(\eta_i,s)F(s;k,p).
\]

There is no upper-end term because `G_K(η_i,η_i)=0`. Compute the lapse and density through quadratic from (1), and compute the extrinsic curvature from its exact formula in §2. In particular,

\[
\chi^{(2)}=\operatorname{tr}j/6-\operatorname{tr}(h^2)/6,
\quad (\rho/\bar\rho)^{(2)}=-4\chi^{(2)}+16(\chi^{(1)})^2,
\quad (N/a)^{(2)}=\chi^{(2)}+(\chi^{(1)})^2.
\]

The regular solution itself contains no arbitrary matching time. If it is split at `η_i`, the exact solution value and derivative above must be propagated together; uniqueness of the positive-time linear forced equation proves independence of that split. This is a statement about the quadratic dynamical data, not yet a check on an unconstructed cubic observable.

## 8. What this closes, and what it does not

This supplies an explicit modewise quadratic construction, rather than setting higher-order fields to zero. Its gauge, regular class, free seed convention, forced homogeneous sector, finite-time metric and velocity derivatives, and retarded boundary prescription are stated. The time integrals are unevaluated but fully specified. The constraint proof is analytical; exact symbolic tests supplement it.

Still required before accepting the **full** milestone independently: review of the exact reduction (3), the bilinear source and symmetrization (12)–(16), and the use of the propagation identities (23) for this boundary class. The new work has only parent-authored checks so far.

Still required for the original problem: extend the constrained construction to cubic order; construct the observer/clock/tetrad/ray/emission kernels; normalize the map; verify the explicit observable under gauge transformations; and prove the required spectral restrictions and convergence. The absence of total-momentum denominators in this quadratic gauge construction is helpful, but it does not prove convergence of the eventual full spectrum or legality of its cubic contraction restrictions.
