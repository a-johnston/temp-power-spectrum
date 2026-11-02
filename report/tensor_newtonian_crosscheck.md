# Independent second-order tensor cross-check

## Result and scope

**The monochromatic second-order temperature-frame tensor amplitude agrees exactly with the radiation-era Newtonian-gauge induced-GW kernel of Kohri–Terada, after the nonlinear gauge and variable transformations.** The agreement is at arbitrary finite time and general shell angle, not only at early or late times.

The comparison uses the report's field equations and variable definitions, but **no project source generator, cached kernel, numerical evolution, or quadrature**. The standalone check is:

```sh
.venv/bin/python -B solution/check_tensor_newtonian.py
```

It independently reconstructs the temperature-frame source from Christoffel/Ricci formulas, checks the gauge transformation, checks the literature's closed kernel, and proves equality of the finite-time differential equations and regular initial data. It does not import or test the existing expanded ancillary closed-form file. This is an external check of the analytical solution defining that file, not an independent audit of every stored coefficient. It also does not check the photon map, SH3, or the fourth-order temperature spectrum.

### Literature used

K. Kohri and T. Terada, *Semianalytic Calculation of Gravitational Wave Spectrum Nonlinearly Induced from Primordial Curvature Perturbations*, Phys. Rev. D **97**, 123532 (2018), [arXiv:1804.08577v2](https://arxiv.org/html/1804.08577).

Relevant equations: metric (1), tensor equation/source (7)–(9), retarded kernel (15), radiation Green function and scalar transfer (18)–(20), finite-time closed form (22). We use the **unaveraged finite-time amplitude**, not the late-time squared/oscillation-averaged kernel (25)–(26). The variable relativistic-degrees-of-freedom correction discussed around their (21) is not used: this problem has exactly `a proportional to conformal time` and constant radiation thermodynamics.

The independent check uses KT alone; it does not claim a second, separate comparison with Espinosa–Racco–Riotto.

## 1. Conventions: three otherwise easy factor-of-two mistakes

Write the Newtonian metric as in KT (1):

\[
ds^2=a^2(\eta)\left[-(1+2\Phi)d\eta^2+
\left((1-2\Phi)\delta_{ij}+\tfrac12h^{\rm KT}_{ij}\right)dy^i dy^j\right].
\]

Second-order scalars and vectors can be left implicit here: they do not contribute to the second-order spatial TT projection in the transformation below. `h^KT` is already second order; **the spatial metric contains `h^KT/2`, not `h^KT`**.

The report instead uses

\[
g=N^2(-dt^2+\gamma_{ij}dX^i dX^j),\qquad
\gamma=e^{2\zeta}e^H,\qquad
H=\epsilon H_1+\tfrac12\epsilon^2H_2+\cdots.
\]

Thus the desired tensor is in a **matrix logarithm after conformal rescaling**, not the spatial metric itself. For a seed `g=z_1 exp(ip.X)+z_2 exp(iq.X)`, the coefficient of `z_1 z_2` in the full `H` is exactly the report's symmetric kernel `A_2(p,q)`; there is no additional factorial.

KT uses Fourier measure `(2 pi)^(-3/2)`; the report uses `(2 pi)^(-3)`. Converting mode amplitudes by `g_KT=g_report/(2 pi)^(3/2)` converts the convolution to the report's measure without a residual numerical factor.

Take `K=p+q`, `|p|=|q|=k_0`, `mu=p.q/k_0^2`, and, for `-1<mu<1`,

\[
u=\hat K,\qquad v=\frac{p-q}{|p-q|},\qquad w=u\times v.
\]

The report's polarization is `v_i v_j-w_i w_j`, with squared norm 2. A unit-normalized KT plus polarization is this tensor divided by `sqrt(2)`. Working with the Cartesian TT projector avoids introducing a polarization-normalization error:

\[
\Lambda_{ij}{}^{ab}p_a p_b
=\frac{k_0^2(1-\mu)}4(v_i v_j-w_i w_j).
\]

## 2. First-order coordinate transformation

Use the embedding of temperature coordinates into Newtonian coordinates

\[
\eta=t+\alpha(t,X)+O(\zeta^2),\qquad
 y^i=X^i+\partial^i\beta(t,X)+O(\zeta^2).
\]

Here `alpha` and `beta` are first-order fields, not homogeneous background shifts. On the shell write `alpha=A(t) g(X)`, `beta=B(t) g(X)`. Define the **acoustic** variable and profile

\[
x=\frac{k_0t}{\sqrt3},\qquad s(x)=\frac{2\sin x}{x}-\cos x.
\]

Then

\[
B(t)=\frac3{k_0^2}[s(x)-1],\qquad
A(t)=B'(t)=\frac{\sqrt3}{k_0}s_x'(x),
\]

\[
\frac{\Phi}{g}=\phi(t)=-\frac{2}{x^2}\left(\frac{\sin x}{x}-\cos x\right),\qquad
\frac{\nu_1}{g}=n(t)=\cos x-\frac{\sin x}{x}.
\]

The minus sign of `Phi/g` follows from the report's seed `gamma -> exp(2 zeta) delta`; in particular `Phi/g -> -2/3`. KT's power normalization (13) supplies the same squared factor `4/9`, so a convention for the sign of primordial curvature does not affect the quadratic source.

These expressions satisfy exactly

\[
B'=A,\qquad n=\phi+A'+A/t,
\qquad A/t-\phi-n-k_0^2B/3=1,
\qquad A=-\frac t2(\phi+t\phi').
\]

Respectively these enforce zero shift, the lapse transformation, the primordial determinant, and following the radiation-fluid velocity. They also reproduce the report's first-order traceless metric:

\[
H_{1,ij}=2\left(\partial_i\partial_j\beta-
\tfrac13\delta_{ij}\Delta\beta\right).
\]

The integration constant in `B` is fixed by `B(0)=0`, so the primordial spatial labels agree. Since the induced tensor is second order, evaluating it at `eta=t+alpha` rather than `eta=t` changes it only at third order.

## 3. The nonlinear gauge map simplifies in logarithmic variables

In the following equations, a subscript `[2]` means the actual quadratic part of a field, **not** its second Taylor derivative.

Up to scalar multiples of `delta_ij` and symmetrized gradients (annihilated by TT projection), the spatial pullback divided by `a(t)^2` is

\[
\left(\frac{g_{ij}}{a(t)^2}\right)_{[2]}^{\rm TT}
=\frac12h^{\rm KT}_{ij}
+\Lambda_{ij}{}^{ab}\left[
 \beta_{,ac}\beta_{,bc}-\alpha_{,a}\alpha_{,b}
 +4(\alpha/t-\Phi)\beta_{,ab}\right].
\]

Dividing by `N^2/a(t)^2=exp(2 nu)` changes the last coefficient to `alpha/t-Phi-nu_1`. Extracting `H` from `gamma=exp(2 zeta) exp(H)` then subtracts

\[
\left(\tfrac12H_1^2+2gH_1\right)^{\rm TT}.
\]

Using `alpha/t-Phi-nu_1=g-Delta beta/3`, every mixed scalar-times-Hessian term cancels. The result is

\[
\boxed{
 H_{[2],ij}^{\rm TT}
 =\frac12h^{\rm KT}_{ij}
 -\Lambda_{ij}{}^{ab}\left[
  \beta_{,ac}\beta_{,bc}+\alpha_{,a}\alpha_{,b}
 \right].}
\]

Unknown second-order coordinate shifts enter the background spatial metric only as a trace or symmetrized gradient, so their TT projection is zero. Thus a full second-order scalar gauge solution is unnecessary for this comparison. In particular, **simply equating the report tensor with the Newtonian metric tensor would be wrong**.

## 4. Literature kernel and the explicit answer for `F_T`

Define the **tensor** phase and momentum ratio

\[
X=Kt,\qquad r=k_0/K=[2(1+\mu)]^{-1/2}.
\]

KT's radiation kernel on equal scalar momenta is

\[
I_{\rm RD}(r,r,X)=\int_0^X dY\,\frac YX\sin(X-Y)\,f_r(Y),
\]

\[
f_r(Y)=2\phi_r(Y)^2+
[\phi_r(Y)+Y\partial_Y\phi_r(Y)]^2,\qquad
\phi_r(Y)=-\frac2{z^2}\left(\frac{\sin z}{z}-\cos z\right),\quad z=rY/\sqrt3.
\]

All integration limits and the regular/no-free-tensor boundary prescription are explicit. In particular,

\[
I_{XX}+\frac2X I_X+I=f_r(X),\qquad
I=\frac29X^2+O(X^4).
\]

Equations (7)–(9) of KT imply, in the report's Fourier convention,

\[
h^{\rm KT}_{ij}(K,t)=4\int_p
\Lambda_{ij}{}^{ab}p_a p_b\,\frac{I_{\rm RD}(|p|/K,|K-p|/K,Kt)}{K^2}
 g_p g_{K-p}.
\]

Combining this with the gauge map and the report's Taylor convention gives the equal-shell kernel

\[
 A_{2,ij}^{\rm TT}(p,q;t)=
 \Lambda_{ij}{}^{ab}p_a p_b
 \left[\frac{4I_{\rm RD}(r,r,Kt)}{K^2}
       +2(p\cdot q)B(t)^2-2A(t)^2\right].
\]

Therefore the answer for the report's single amplitude is

\[
\boxed{
F_T(t,k_0,\mu)=k_0^2(1-\mu)
\left[\frac{I_{\rm RD}(r,r,Kt)}{K^2}
+\frac12\mu k_0^2 B(t)^2-\frac12 A(t)^2\right].}
\]

Equivalently, entirely dimensionlessly,

\[
\boxed{
F_T=
\frac{1-\mu}{2(1+\mu)}I_{\rm RD}(r,r,X)
+\frac32(1-\mu)\left[3\mu(s(x)-1)^2-s_x'(x)^2\right],
\quad X=\sqrt{6(1+\mu)}\,x.}
\]

This is a known, evaluated finite-time tensor response **plus an elementary gauge correction**. It is considerably shorter than the expanded temperature-frame tensor kernel.

### Closed form with no removable resonance logarithms

KT (22), specialized to `u=v=r`, can be written in the report's entire functions. Put

\[
d=2r^2-3,\quad a_\pm=1\pm2r/\sqrt3,\quad x=rX/\sqrt3,
\]

\[
\mathcal C=\operatorname{Cin}(a_-X)+\operatorname{Cin}(a_+X)-2\operatorname{Cin}(X),
\quad
\mathcal S=\operatorname{Si}(a_-X)+\operatorname{Si}(a_+X)-2\operatorname{Si}(X),
\]

\[
\mathcal E=-4r^2d\sin X+\frac{24r^2\cos^2x}{X}
-\frac{48\sqrt3\,r\cos x\sin x}{X^2}
+\frac{12(6+dX^2)\sin^2x}{X^3}.
\]

Then

\[
\boxed{
I_{\rm RD}(r,r,X)=\frac3{4r^6X}
\left\{\mathcal E+d^2[\sin X\,\mathcal C+\cos X\,\mathcal S]\right\}.}
\]

Here `Cin(z)=integral_0^z (1-cos u)/u du` and `Si(z)=integral_0^z sin(u)/u du` are entire, with their zero values defined by continuity. This follows from combining KT's Ci terms and explicit logarithm **before** taking a resonant limit. The standalone check differentiates this formula exactly and verifies both its wave equation and its `2 X^2/9` initial behavior. Apparent poles at `X=0` cancel in the assembled expression; it is not intended as a numerically stable small-X evaluation formula.

At the sound resonance `mu=-1/3`, `r=sqrt(3)/2`, so `a_-=0`. One simply sets `Cin(0)=Si(0)=0`. No frequency denominator, branch prescription, late-time limit, or resonant exclusion is needed.

## 5. Independent comparison with the temperature-frame equation

Rebuild the quadratic source using the report's exact mixed spatial equation,

\[
\operatorname{TF}\left[S'+2(t^{-1}+\nu')S+
 R^i{}_j[\gamma]-2D^iD_j\nu+2D^i\nu D_j\nu\right]=0.
\]

Choose `p=(b,0,c)`, `q=(-b,0,c)`, so `K` is along `z`,

\[
b^2=k_0^2(1-\mu)/2,\qquad c^2=k_0^2(1+\mu)/2.
\]

Let `h=-2B` be the first-order Hessian coefficient, `n=nu_1/g`, and let primes below mean `d/dt`. With the unknown `A_2` set to zero, expand `gamma=exp(2g) exp(H_1)` and compute Christoffels, Ricci curvature, raised lapse Hessian, and the lapse-gradient products directly. On the shell `H_1(t)` is proportional to a fixed spatial matrix, so its quadratic commutator shear vanishes. The TT equation is

\[
F_T''+\frac2tF_T'+K^2F_T=\mathcal S_T,
\]

\[
\boxed{
\mathcal S_T=\frac{k_0^2(1-\mu)}{36}
\left[k_0^4h^2(9\mu^2+12\mu+1)
-6k_0^2h(6\mu n+3\mu+2)
-36h'n'-72n^2+36\right].}
\]

Write `L_K=partial_t^2+2 partial_t/t+K^2`. Direct substitution of the radiation profiles proves the exact identity

\[
\boxed{
\mathcal S_T=k_0^2(1-\mu)f_r(Kt)
+L_K\left\{\frac{k_0^2(1-\mu)}2
[\mu k_0^2B^2-A^2]\right\}.}
\]

This proves that the literature-plus-gauge expression and the report's Green construction obey the **same complete finite-time equation**, including the velocity/lapse terms. No series truncation is used in this identity.

The initial data also agree:

\[
B=t^2/6+O(t^4),\quad A=t/3+O(t^3),\quad
I/K^2=2t^2/9+O(t^4),
\]

\[
\boxed{F_T=k_0^2(1-\mu)t^2/6+O(t^4).}
\]

This matches the TT projection of the report's nonlinear gradient initial data. Both solutions have zero constant primordial tensor and no decaying `1/t` tensor. Their difference solves the homogeneous tensor equation; regularity and zero constant term therefore make the difference identically zero. This uniqueness argument is what turns the source identity into an equality of the full solutions.

### Endpoints

- `mu=1`: the scalar pair is collinear and the tensor amplitude vanishes. The polarization basis itself degenerates; the vanishing Cartesian tensor is unambiguous.
- `mu=-1`: total momentum is zero and an SVT decomposition relative to `K` is undefined. Do **not** discard the separately specified forced homogeneous solution. A directional limit of the displayed TT projection is well defined through
  `I/K^2 -> integral_0^t ds (s/t)(t-s) f(s)`, but it is only that projection, not a replacement for the full homogeneous kernel.
- Finite resonance `mu=-1/3`: covered exactly by the entire closed expression above.

## 6. Verification ledger and suggested integration

The standalone symbolic script checks:

1. the temperature-frame TT source from elementary geometry, without the existing engine;
2. the quadratic pullback/conformal/matrix-log cancellation;
3. the first-order gauge, lapse, determinant, and comoving relations;
4. the equality of the radiation source to KT (20), including its normalization;
5. the exact finite-time source identity for symbolic angle and time;
6. matching regular initial data;
7. the entire version of KT (22), by exact differentiation and initial behavior;
8. negative controls: dropping the temporal gauge term, dropping the spatial/log correction, or confusing the metric tensor normalization fails.

No changes to the main report, run-all script, SH3 implementation, or shared symbolic engines are made here. The report author can integrate the boxed gauge map, the compact `F_T` formula, and the source/uniqueness proof; the literature comparison can then be marked completed **at the analytical solution level**, while keeping any separate audit of the stored expanded kernel distinct.
