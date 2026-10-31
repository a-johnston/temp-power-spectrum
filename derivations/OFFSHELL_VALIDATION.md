# Offshell ADM/Bianchi evidence retained after the quadratic review

## Scope

The second user-supplied independent review accepted the quadratic construction
at `ea6b053`. It recommended retaining reproducible checks **upstream** of the
source and constraint-energy tests. This note and `../validation/check_offshell.py`
record the reduction and offshell identities; the reviewer’s separate execution
is attributed in `../reviews/INDEPENDENT_REVIEW_02.md`.

These are local tensor identities. The coordinate existence assumptions and
regular boundary class remain those of `QUADRATIC_INITIAL_DATA.md` (Q).

## 1. Mixed ADM to the exact metric equation

Put `B=H I+M`, `M=q^(-1)q'/2`, `tr M=3χ'`, and `N=a exp(χ)`.
The mixed ADM equation uses `K^i_j=−B^i_j/N`. Multiplying it by `−N` gives

\[
B'-(H+\chi')B=e^{2\chi}q^{-1}
(D D\chi+d\chi\,d\chi-R)-(\operatorname{tr}B)B
+H^2e^{-2\chi}I.
\]

The last term uses the fixed radiation background and the exact fluid reduction:
`8πGN²ρ/3=H² exp(−2χ)`. The lapse conversion is

\[
N D_\gamma^iD_jN=e^{2\chi}q^{ia}
(D_aD_j\chi+\partial_a\chi\partial_j\chi).
\]

Using `H'=−H²` and `tr M=3χ'` gives

\[
M'=-2(H+\chi')M+
[H^2(e^{-2\chi}-1)-2H\chi']I
+e^{2\chi}q^{-1}(D D\chi+d\chi\,d\chi-R).
\]

Finally `q''=q'q^(-1)q'+2qM'` gives Q (3). Neither gravitational constraint has
been imposed. The determinant identity used here is a fluid/gauge relation,
not an on-constraint-surface substitution.

The regression evaluates this matrix identity at an arbitrary spatial
orthonormal point with fully symbolic symmetric `M` and spatial Hessian/curvature
combination. It checks all noncommuting products; it does not numerically test a
particular cosmological solution.

## 2. What spatial evolution imposes offshell

Let `E_μν=G_μν−8πGT_μν`. Gauss–Codazzi with the chosen extrinsic-curvature sign gives

\[
E_{nn}=\mathcal H_c/2,\qquad E_{ni}=-\mathcal M_i.
\]

The ADM evolution equations impose the **trace-reversed** spatial equations

\[
E_{ij}-\tfrac12\gamma_{ij}E^\mu{}_{\mu}=0.
\]

Writing `S_E=γ^ij E_ij`, their trace says
`S_E=(3/2)(−E_nn+S_E)`, hence `S_E=3E_nn` and

\[
E_{ij}=\gamma_{ij}E_{nn}=\tfrac12\gamma_{ij}\mathcal H_c.
\]

Setting `E_ij=0` before imposing the Hamiltonian constraint would be incorrect
and would change the propagation derivation. The tests keep `H_c` arbitrary.

## 3. Coordinate Bianchi calculation

At one spatial Riemann normal point take `γ_ij=δ_ij` and its spatial first
coordinate derivatives zero, but retain arbitrary symmetric `K_ij`, arbitrary
lapse gradient `N_i`, and arbitrary `N'`. With zero shift the needed spacetime
Christoffel symbols are

\[
\Gamma^0{}_{00}=N'/N,\quad \Gamma^0{}_{0i}=N_i/N,
\quad \Gamma^0{}_{ij}=-K_{ij}/N,
\]
\[
\Gamma^i{}_{00}=NN_i,\quad \Gamma^i{}_{0j}=-NK_{ij},
\quad \Gamma^i{}_{jk}=0.
\]

The mixed residual components are

\[
E^0{}_0=-\mathcal H_c/2,\quad E^0{}_i=\mathcal M_i/N,
\quad E^i{}_0=-N\mathcal M_i,
\quad E^i{}_j=\delta^i_j\mathcal H_c/2.
\]

Substitute these into

\[
\nabla_\mu E^\mu{}_\nu=\partial_\mu E^\mu{}_\nu
+\Gamma^\mu{}_{\mu a}E^a{}_\nu
-\Gamma^a{}_{\mu\nu}E^\mu{}_a.
\]

For `ν=0`, the result is minus one half of

\[
\mathcal H_c'+2ND_i\mathcal M^i-2NK\mathcal H_c
+4\mathcal M^iD_iN.
\]

For `ν=i`, it is `1/N` times

\[
\mathcal M_i'+\tfrac12ND_i\mathcal H_c-NK\mathcal M_i
+\mathcal H_cD_iN.
\]

The exact fluid equations ensure `∇_μ T^μ_ν=0`, so the Bianchi identity sets these
expressions to zero and recovers Q (23) offshell. The explicit `N'` cancels;
extrinsic-curvature components combine into the correct trace terms. The
symbolic test checks these cancellations with arbitrary independent residuals,
residual derivatives, lapse derivatives and extrinsic-curvature entries.

The pointwise normal-coordinate calculation loses no tensor content; restoring
spatial covariant derivatives gives the stated identities. It is not a claim
that spatial curvature or shear vanishes in the solution.

## 4. Additional retained review checks

`check_initial_data.py` now also compares the quadratic source to the separate
full metric expansion using **fully symbolic symmetric matrix functions**.
Their values and first time derivatives are independent, and all six input
momentum components remain symbolic. This strengthens the original fixed-matrix
regression without replacing the exact ADM derivation above.

A finite-time constant-dilation regression uses

\[
J(\eta;0,p)=(2-p\partial_p)H(\eta;p),
\qquad \chi^{(2)}(\eta;0,p)=-p\partial_p c(\eta;p).
\]

These follow by rescaling the spatial metric by `exp(2λ)` and the input wavevector
by `exp(−λ)`. This is a test of the dynamical kernels with the seed/chart convention
transformed consistently, not a proof of gauge invariance of a future observed map.
