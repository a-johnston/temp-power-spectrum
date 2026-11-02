# Exact redshift in longitudinal/Poisson variables

## Result and scope

The temperature-frame shear integral has an exact longitudinal-gauge reformulation with explicit lapse-potential endpoints, fluid Doppler factors, integrated scalar potentials, shift-vector terms, and tensor shear. This is the analytical basis for the report's proposed late-time recombination.

The formulation is **exact in the metric and photon trajectory**. Its derivation does not assume small angles, a Born ray, negligible emission time, or a vanishing induced shift. The explicit coordinate map is worked out at first order, sufficient to demonstrate why a spatial relabelling alone fails to remove the oscillating temperature-time slicing. Exact transformation formulas are supplied for the higher-order map, but its second- and third-order gauge generators are not solved here.

This is **not a proof of a late-time limit of the fourth-order spectrum**. Potentials, endpoint velocities, physical-surface displacements, and resonant field responses must still be combined before making such a claim.

Exact algebraic checks:

```sh
.venv/bin/python -B solution/check_longitudinal_redshift.py
```

## 1. Metric, shift sign, and measured energies

Use longitudinal/Poisson coordinates `(eta,y)` and write

\[
ds^2=a^2(\eta)\left[-\mathcal N^2d\eta^2+
 h_{ij}(dy^i+\beta^i d\eta)(dy^j+\beta^j d\eta)\right],
\]

\[
\mathcal N=e^\Phi,\qquad h=e^{-2\Psi}M,\qquad M=e^\chi,\quad\operatorname{tr}\chi=0.
\]

At linear scalar order this is the usual longitudinal metric. Beyond linear order a Poisson gauge fixes the scalar spatial shear and scalar shift and retains the induced transverse shift and tensors. The identities below are independent of the precise nonlinear completion of these coordinate conditions: they are valid for arbitrary lapse, shift, and positive spatial metric. In particular, no transverse-vector or tensor contribution is set to zero by the redshift derivation.

The ADM sign is important: with the displayed **plus** sign in `dy+beta d eta`, the conformal unit normal is

\[
\widetilde n=\mathcal N^{-1}(\partial_\eta-\beta^i\partial_i).
\]

Also `g_00/a^2=-N^2+h_ij beta^i beta^j`, not simply `-N^2` when the shift is nonzero.

Let the forward photon direction relative to this normal be `e^i`, with `h_ij e^i e^j=1`. If `p` is the physical affinely parametrized momentum, `p_tilde=a^2 p` is the conformal affinely parametrized momentum, and

\[
\widetilde p=\mathcal E(\widetilde n+e^i\partial_i),\qquad
\frac{dy^i}{d\eta}=\mathcal N e^i-\beta^i.
\]

Write the fluid velocity as

\[
u=\frac{\Gamma}{a}(\widetilde n+v^i\partial_i),\qquad
\Gamma=(1-h_{ij}v^iv^j)^{-1/2},
\]

and define the exact fluid Doppler factor

\[
\mathcal D=\Gamma(1-h_{ij}v^i e^j).
\]

Then the physical energy measured by the fluid is

\[
\boxed{E=-u\cdot p=\frac{\mathcal E}{a}\mathcal D.}
\]

The measured sky direction is fixed by the report's fluid-rest-frame tetrad, not by `e^i`. At the observer the correct initialization is

\[
p_o^\mu=E_o(u_o^\mu-\hat n^a E_{a,o}^\mu),
\]

followed by its decomposition relative to the ADM normal. Simply identifying `e_o=-hat n` would omit aberration when the fluid moves relative to the longitudinal slices. The specified Fermi–Walker tetrad and its initial orientation are transported as the same physical objects under the coordinate transformation.

## 2. Exact ADM energy equation

For the conformal metric, let `p_i` be the photon spatial covector and

\[
q=\sqrt{h^{ij}p_i p_j}=\mathcal E,\qquad
\mathcal H_\gamma=-p_\eta=\mathcal N q-\beta^i p_i.
\]

Hamilton's equations give

\[
y^{i\prime}=\mathcal N e^i-\beta^i,
\]

\[
p_i'=-q\partial_i\mathcal N+
\frac{\mathcal N q}{2}e^j e^k\partial_i h_{jk}
+p_j\partial_i\beta^j.
\]

Differentiating `q^2=h^ij p_i p_j` along the photon trajectory cancels the photon-advection derivatives of `h`, leaving

\[
\boxed{
\frac{d\ln\mathcal E}{d\eta}
=-e^i\partial_i\mathcal N
-\frac12e^i e^j\left(\partial_\eta h_{ij}-(\mathcal L_\beta h)_{ij}\right).
}
\]

The Lie derivative is

\[
(\mathcal L_\beta h)_{ij}=
\beta^k\partial_k h_{ij}+h_{kj}\partial_i\beta^k+h_{ik}\partial_j\beta^k.
\]

All partial derivatives are at fixed longitudinal spatial coordinates. All fields in the ray integral are then evaluated on the actual perturbed trajectory. The standalone script verifies the identity on independent symbolic spatial/time metric jets. Omitting the shift Lie derivative fails this check.

Define

\[
J=\int_{\eta_e}^{\eta_o}d\eta\left[
e^i\partial_i\mathcal N+
\frac12e^i e^j(\partial_\eta h_{ij}-\mathcal L_\beta h_{ij})\right].
\]

The exact measurement becomes

\[
\boxed{
T_{\rm obs}=T_\star\frac{a(\eta_e)}{a(\eta_o)}
\frac{\mathcal D_o}{\mathcal D_e}e^{-J}.
}
\]

Both limits are physical-event coordinates. In particular, `eta_e` is direction dependent and is not the background emission time.

## 3. Extracting the lapse-potential endpoints

With `h=e^(-2 Psi) M`, the energy integrand is

\[
e^i\partial_i\mathcal N
-(\partial_\eta-\beta\cdot\partial)\Psi
+\frac12e^{-2\Psi}e^i e^j
(\partial_\eta M_{ij}-(\mathcal L_\beta M)_{ij}).
\]

Along the photon,

\[
\mathcal N e^i\partial_i\Phi
=\frac{d\Phi}{d\eta}-\partial_\eta\Phi+\beta^i\partial_i\Phi.
\]

Consequently

\[
J=[\Phi]_e^o-\int_{\eta_e}^{\eta_o}\mathcal S\,d\eta,
\]

\[
\boxed{
\mathcal S=(\partial_\eta-\beta\cdot\partial)(\Phi+\Psi)
-\frac12e^{-2\Psi}e^i e^j
\left(\partial_\eta M_{ij}-(\mathcal L_\beta M)_{ij}\right).
}
\]

The desired exact longitudinal-gauge redshift formula is therefore

\[
\boxed{
T_{\rm obs}(\hat n)=T_\star
\frac{a_e}{a_o}
\frac{\mathcal D_o}{\mathcal D_e}
\exp\left[\Phi_e-\Phi_o+
\int_{\eta_e(\hat n)}^{\eta_o}\mathcal S(\eta,y(\eta))\,d\eta\right].
}
\]

This formula retains:

- scalar-potential endpoints, including their displaced evaluation;
- the time-dependent background scale factor at the actual endpoints;
- both fluid velocities and the nonlinear Doppler factors;
- scalar integrated redshift;
- shift advection of the scalar potentials;
- shift gradients in `L_beta M`;
- the full tensor matrix `M`, including products from its exponential.

In a purely scalar zero-shift metric (`M=1`, `beta=0`), it reduces **exactly**, even at nonlinear amplitude, to

\[
T_{\rm obs}=T_\star\frac{a_e}{a_o}
\frac{\mathcal D_o}{\mathcal D_e}
\exp\left[\Phi_e-\Phi_o+
\int_{\eta_e}^{\eta_o}(\Phi'+\Psi')\,d\eta\right].
\]

But imposing this purely scalar specialization on the actual second-/third-order radiation solution would discard required induced fields. For example, even with `M=1` the shift-gradient part of `S` contains `e^i e^j partial_i beta_j` at its first nonzero order.

## 4. Exact equivalence to the temperature-frame shear formula

Distinguish the temperature-frame conformal factor from the ADM lapse:

\[
N_T=\mathcal T/T_{\rm phys},\qquad
 g=N_T^2(-dt^2+\gamma_{ij}dX^i dX^j).
\]

The temperature-frame conformal photon energy is `N_T E`, whereas the longitudinal ADM energy is `a E/D`. Thus

\[
\frac{\widetilde E_T}{\mathcal E}=\frac{N_T}{a}\mathcal D.
\]

Since the original exact shear integral obeys `I=-ln[(N_T E)_o/(N_T E)_e]`,

\[
\boxed{
I=J-\left[\ln\frac{N_T}{a}+\ln\mathcal D\right]_e^o
=\left[\Phi-\ln\frac{N_T}{a}-\ln\mathcal D\right]_e^o
-\int_{\eta_e}^{\eta_o}\mathcal S\,d\eta.
}
\]

This is the exact recombination identity. It specifies which endpoint terms must accompany the replacement of the temperature-frame shear by longitudinal potentials and shift/tensor terms.

As a check, stay in temperature coordinates but choose the background conformal factor `a=a_1 t`. Then `N=e^nu`, `h=e^(2nu) gamma`, `beta=0`, `D=1`, and the ADM integral is

\[
J=\int(d\nu/dt+\sigma_{ij}e_T^i e_T^j)dt=[\nu]_e^o+I.
\]

The boxed identity reduces identically to the original redshift lemma.

## 5. The time shift and the induced shift are not optional

### Exact coordinate transformation formulas

Let `t=t(eta,y)` and `X=X(eta,y)` be the inverse map from longitudinal to temperature coordinates. With capital coordinate indices `A,B` define

\[
\mathcal B_{AB}=-\partial_A t\partial_B t+
\gamma_{ij}\partial_A X^i\partial_B X^j,\qquad \rho=N_T/a(\eta).
\]

All temperature-frame fields are evaluated at the transformed arguments. The transformed ADM variables are exactly

\[
h_{ij}=\rho^2\mathcal B_{ij},\qquad
\beta^i=(\mathcal B_{\rm sp}^{-1})^{ij}\mathcal B_{\eta j},
\]

\[
\mathcal N^2=\rho^2\left[-\mathcal B_{\eta\eta}+
\mathcal B_{\eta i}(\mathcal B_{\rm sp}^{-1})^{ij}\mathcal B_{\eta j}\right].
\]

These expressions keep the time-shift gradients and the shift vector at every order. A perturbative implementation can solve the chosen Poisson-gauge conditions for the coordinate generators and insert them here. It must not transform only the spatial metric and then declare the lapse and shift unchanged.

The direct pullback of the original ray integral is equally explicit. Define the derivative along the longitudinal photon by

\[
\mathscr D_\gamma=\partial_\eta+(\mathcal N e^i-\beta^i)\partial_i.
\]

Then

\[
\frac{dt}{d\eta}=\mathscr D_\gamma t,\qquad
e_T^i=\frac{\mathscr D_\gamma X^i}{\mathscr D_\gamma t},\qquad
I=\int_{\eta_e}^{\eta_o}\sigma_{ij}\,
\frac{(\mathscr D_\gamma X^i)(\mathscr D_\gamma X^j)}{\mathscr D_\gamma t}\,d\eta.
\]

In particular, the integration measure cannot be left as `dt=d eta` after a spatially dependent time shift. The velocity-potential identity `u_mu=-N_T partial_mu t` also gives

\[
\boxed{\mathcal D=\frac{N_T}{a\mathcal N}\mathscr D_\gamma t.}
\]

This relates the transformed time coordinate directly to the exact Doppler factor and provides a further consistency check on endpoint transformations.

### Explicit first-order radiation map

For a shell seed `g(X)`, use

\[
\eta=t+A(t)g(X),\qquad
y^i=X^i+B(t)\partial^i g(X),
\]

\[
x=k_0t/\sqrt3,\qquad
B=\frac3{k_0^2}\left(\frac{2\sin x}{x}-\cos x-1\right),\qquad A=B'.
\]

The spatial displacement is precisely the one in the report's Eulerian recombination lead: `Psi_displacement=grad(B g)=-Q grad(g)/2`. The additional `eta-t=A g` is the required time transformation.

Let `n(t)=cos x-sin x/x` be the temperature-frame lapse coefficient. The transformed fields are

\[
\frac{\Phi}{g}=n-A'-A/t,
\qquad
\frac{\Psi}{g}=A/t-n-1-k_0^2B/3,
\]

\[
\beta_i^{(1)}=\partial_i[(A-B')g]=0,
\]

\[
\boxed{
\frac\Phi g=\frac\Psi g=-\frac2{x^2}
\left(\frac{\sin x}{x}-\cos x\right).
}
\]

Both potentials therefore decay as `O(x^-2)` at first order. A purely spatial transformation (`eta=t`) instead retains the shift `beta=-grad(B' g)` and the unsuppressed temperature-frame lapse `n g`. Setting that shift to zero would not be a coordinate reformulation of the same experiment.

At second order and beyond a scalar time shift does not remove the physical transverse shift induced by scalar couplings. That shift must remain in the exact redshift expression. No explicit second-/third-order solution of these gauge generators is claimed in this note.

## 6. Physical emission, observer clock, and sky normalization

The emission event is the root

\[
N_T(\eta_e,y(\eta_e))=\mathcal T/T_\star,
\]

on the same branch as in the report. It is not a constant-longitudinal-time surface. At first order its time shift is

\[
\delta\eta_e=d_{e1}^{\rm temp}+\alpha_e
=-\bar t_\star\nu_1(\bar t_\star,\bar X_e)+\alpha_e.
\]

The observer remains the tagged fluid line, which moves in longitudinal coordinates. Its physical initial event is likewise the root `N_T=Tcal/T_I`. Along this worldline,

\[
\boxed{d\tau=\frac{a\mathcal N}{\Gamma}\,d\eta,\qquad
\int_{\eta_I}^{\eta_o}\frac{a\mathcal N}{\Gamma}\,d\eta=\tau_o.}
\]

The upper and lower coordinate endpoints both shift. At first order,

\[
\delta\eta_o=D_1^{\rm temp}+\alpha_o,\qquad
\delta\eta_I=d_{I1}^{\rm temp}+\alpha_I.
\]

This can also be derived directly from the longitudinal clock. Since `Phi=nu-alpha'-alpha/t`,

\[
\int_i^o t\Phi\,dt=\int_i^o t\nu\,dt-[t\alpha]_i^o,
\]

which turns the longitudinal proper-time equation into the temperature-frame clock equation plus the endpoint time transformation. The `alpha_I` terms cancel correctly only when the physical initial event is retained.

For comparison, the first-order log-temperature perturbation from the exact formula is

\[
\delta\ln T_{\rm obs}
=\mathcal H_e\delta\eta_e-\mathcal H_o\delta\eta_o
+\Phi_e-\Phi_o
+\hat n\cdot(v_o-v_e)
+\int_{\bar\eta_e}^{\bar\eta_o}(\Phi'+\Psi')\,d\eta,
\]

with forward background direction `e=-hat n`. This fixes the endpoint and Doppler signs. The nonlinear formula, not this linear expression, must be used to construct the higher-order map.

Finally use exactly the report's normalization:

\[
\Theta(\hat n)=\frac{T_{\rm obs}(\hat n)}{\langle T_{\rm obs}\rangle_\Omega}-1,
\]

and remove the measured dipole only afterwards. The factor `a_e` is direction dependent through the physical emission event. The observer's Doppler factor is direction dependent. Neither can be discarded as a monopole. Direction-independent observer factors may cancel from the completed exact ratio, but should not be selectively omitted from intermediate endpoint or clock expansions.

## 7. What this establishes for late-time recombination

The identity shows analytically how the large temperature-frame lapse and label-displacement contributions can be moved into:

1. decaying first-order longitudinal potentials;
2. physical endpoint temperature/time shifts;
3. observer and emitter Doppler/aberration factors;
4. the induced shift and tensor response.

It therefore supplies the missing derivation behind the proposed recombination. It does **not** imply that the complete redshift or its fourth-order spectrum decays like the scalar potentials: radiation velocities and temperature-surface shifts remain oscillatory, and the higher-order response may have nonuniform resonant limits.

To turn this into a late-time spectral result, remaining work is to construct the needed higher-order longitudinal variables or their contracted responses, expand the complete physical endpoints and Doppler factors, and control the resonant angular integrations. No such limit is interchanged with an integral in this derivation.
