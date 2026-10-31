# Unified derivation: observed temperature spectrum through fourth order

This document carries out [`PLAN.md`](../PLAN.md). Every displayed identity is labelled with the
script that checks it exactly. The status of each result is marked:

- **[evaluated]**: closed form;
- **[representation]**: an exact definite integral, with complete integrand, limits and nesting,
  not yet reduced;
- **[open]**.

Conventions are those of PLAN §1. A field is written `Y = Σ_r ε^r Y^(r)/r!`. Its symmetric kernel
`K_r` is the coefficient of `z₁⋯z_r` when `g = Σ_a z_a e^{ik_a·X}` is inserted and only
multilinear monomials are kept. This **multilinear rule** fixes every combinatorial factor: the
`r!` of the Taylor convention cancels the `r!` orderings of distinct labels. All products below
are evaluated by it (`solution/engine.py` implements it exactly).

---

## 1. Exact framework [evaluated]

For an irrotational radiation fluid (`P = ρ/3`, `ρ = αT⁴`), with `t` the velocity potential and
`X^i` comoving fluid labels:

```
g = N²(−dt² + γ_ij dX^i dX^j),   N = 𝒯/T_phys,   u = N⁻¹∂_t,   γ = e^{2ζ(X)} exp H,  tr H = 0.
```

The form is exact. It satisfies `g(u,u) = −1` and Euler's equation identically (lemma L1). The
continuity equation is equivalent to `∂_t det γ = 0`. So `ζ = ⅙ ln det γ` is time-independent at
every order, and the shear `σ_ij = ½∂_tγ_ij` is exactly traceless.

With `S = γ⁻¹σ = ½ e^{−H} ∂_t e^{H}` (mixed shear) and `N = a₁t e^{ν}`, Einstein's equations
are equivalent to:

| | Exact equation | Role |
|---|---|---|
| (N) | `3ν″ − Δ_γν − 2|Dν|² + tr S² + (3/t²)(e^{−2ν} − 1) = 0` | evolution (`R_00`) |
| (T) | `TF[ S′ + 2(1/t + ν′)S + R^k_j[γ] − 2D^kD_jν + 2D^kν D_jν ] = 0` | evolution (trace-free `R^k_j`) |
| (M) | `D_jS^j_i − 2∂_iν′ + (2/t)∂_iν + 2S^j_i∂_jν + 2ν′∂_iν = 0` | momentum constraint |
| (Ham) | `R − 4Δν − 2|Dν|² − tr S² + 6(1/t + ν′)² − 6e^{−2ν}/t² = 0` | Hamiltonian constraint |

(N), (M) and (Ham) are lemma L4 and its components (`verification/check_field_equations_generic.py`).
(T) is new. In mixed form `γ^{ki}σ′_ij = (S^k_j)′ + 2S^k_iS^i_j`, so the `σσ` terms of the spatial
equation cancel. The source `8πG T_ij ∝ γ_ij` is pure trace, so only the trace-free part evolves,
and the trace combines with (N) into (Ham).
**Check:** `solution/check_evolution_form.py` (generic 2-jets: six free functions in γ with
`t`-independent determinant, and arbitrary `N`).

**Exact redshift** (lemmas L2–L3). A photon emitted isotropically in the fluid frame on
`T_phys = T⋆` and received at observation event `o` has

```
T_obs(n̂) = T_fluid(o) · exp[−I(n̂)],    I(n̂) = ∫_{t_e(n̂)}^{t_o} σ_ij(t, X(t)) e^i e^j dt,
dX^i/dt = e^i,   γ_ij e^i e^j = 1,
de^i/dt = −2S^i_j e^j − Γ^i_jk e^j e^k + e^i σ_jk e^j e^k.
```

`T_fluid(o)` is direction-independent, so it cancels in `Θ` exactly. The emission event is the
root of `N(t_e, X(t_e)) = 𝒯/T⋆`, that is `ln t_e + ν(t_e, X(t_e)) = ln t̄⋆`.

---

## 2. Primordial state and initial data [evaluated]

**Definition.** The state is the solution of (N), (T), (M), (Ham) with:

- **(G)** `ζ = εg` exactly, with `g` zero-mean, homogeneous, isotropic and Gaussian, power `P(k)`;
- **(T)** the `t → 0` limit of `γ` is `e^{2ζ}δ_ij`. It is conformally flat in the marked material
  chart, so there is no tensor seed;
- **(V)** `d(Tu) = 0`;
- **(R)** `ν` and `H` are even and analytic in `t` at fixed momenta, with `ν(0) = 0` and
  `H(0) = 0`.

"No independent higher-order primordial modes" means `ζ^(r≥2) ≡ 0` together with (T), (V) and
(R) at every order. The residual relabellings of the chart are fixed as follows:

- the dilation is fixed by `⟨g⟩ = 0`;
- special conformal transformations are excluded by statistical homogeneity;
- translations and rotations are statistical symmetries.

**Identification.** At `t → 0`, `ν = O(t²)`, so constant-`t` and uniform-density slices coincide
in the limit. `ζ` is therefore the nonlinear uniform-density (δN) curvature perturbation. In both
branches' variables (`q = e^{2ν}γ`, `χ = ν`) the leading nonlinear data agree exactly:

```
ν = −(t²/36) R[e^{2ζ}δ] + O(t⁴) = (t²/18) e^{−2ζ}(2∇²ζ + |∇ζ|²) + O(t⁴),
H = −(t²/3) e^{−2ζ} TF R_ij[e^{2ζ}δ] + O(t⁴).
```

These are checked in `verification/check_initial_data_gradient.py` and orange's
`validation/check_initial_data.py`, and I confirmed by hand that they coincide.

**Initial data at finite `t_i`.** The data are the regular solution of §4, evaluated at `t_i`:

- the spatial metric `γ(t_i)`, the extrinsic curvature through `S(t_i)` and `N(t_i)`;
- the matter variables `ρ = 3a₁²/(8πG N⁴)` and `u = N⁻¹∂_t`.

They satisfy (M) and (Ham) at every order (§4.4). Higher-order data are **not** zero. For
example:

- `γ^(2)(t_i) = 4g²δ + O(t_i²)` at the `g²` level (the `ε²/2` coefficient of `e^{2εg}`);
- `ν^(2)` has kernel `(t_i²/9)[2(k₁²+k₂²) − k₁·k₂] + O(t_i⁴)`;
- at `k₂ = −k₁` this kernel is `5k²t_i²/9`, the forced homogeneous mode, which is kept.

Data posed at any `t_i` evolve back into the same regular solution, so `t_i` drops out of the
observable identically.

---

## 3. First order [evaluated]

For one mode (amplitude `g_k`, `x = kt/√3`):

```
ν₁ = −f(x),  f = sin x/x − cos x = x j₁(x);
H₁ = Q(t)(k_ik_j − δ_ij k²/3),  Q = −6(s(x) − 1)/k²,  s = 2 sin x/x − cos x;
ζ₁ = 1,   σ₁ = S₁ = ½H₁′.
```

These are verified in `verification/check_first_order.py` (all ten equations) and in
`solution/first_order_engine_check.py` (generic direction, through the unified equations).

```
Θ^(1)(n̂) = ∫_k g_k Σ_{ℓ≥1} (2ℓ+1) i^ℓ T_ℓ(k) P_ℓ(k̂·n̂),
T_ℓ = −3 ∫_{x⋆}^{x_o} dx s′(x) w_ℓ(√3(x_o − x)),   w_ℓ = j_ℓ/3 + j_ℓ″.
```

`T_ℓ` equals the Newtonian-gauge SW + Doppler + ISW transfer for every ℓ. Closed forms in
Si/Ci for ℓ = 2, 3 are in `verification/results/`.

```
C_ℓ^(11) = 4π ∫ d ln k 𝒫_g(k) T_ℓ(k)².
```

---

## 4. Field solution at every order

### 4.1 Linear operators (identical at every order) [evaluated]

At total momentum `K = Σk_a`, each order-`r` kernel solves `L_K[K_r] + c^(r) = 0`, where
`c^(r)` is the multilinear coefficient of the exact equations with `K_r` set to zero:

```
L_N = 3ν″ + K²ν − 6ν/t²
L_M = (i/2) K_j A′_ji − 2iK_i(ν′ − ν/t)
L_T = ½(A″ + 2A′/t) + TF R₁[A] + 2(K_iK_j − δ_ij K²/3)ν,
      R₁[A]_ij = ½(K²A_ij − K_aK_iA_aj − K_aK_jA_ai)
```

Here `A` is the order-`r` part of `H`. The operators are confirmed at second order in
`solution/second_order.py` and at third order (configuration `(k,q,−q)`) in
`solution/third_order.py`.

### 4.2 Triangular solution [representation]

Let `c = K/√3`, `x = ct`, `u₁(x) = x j₁(x)`, and `c_r` be the `t²` coefficient of
`(1/18) e^{−2ζ}(2∇²ζ + |∇ζ|²)` at order `r`:

- `c₂ = [2(k₁² + k₂²) − k₁·k₂]/9`;
- `c₃ = (2/9)[−2Σ_a k_a² + Σ_{a<b} k_a·k_b]`, which at `(k,q,−q)` is `−(2/9)(2k² + 5q²)`.

```
ν_r(t) = (9c_r/K²) u₁(Kt/√3) − (1/3) ∫₀^t G_N(t,t′) c_N^(r)(t′) dt′,
G_N(t,t′) = c⁻¹[(1 + 1/(xx′)) sin(x−x′) − ((x−x′)/(xx′)) cos(x−x′)].
```

At `K = 0` the first term is `c_r t²` and `G_N → (t³ − t′³)/(3tt′)`. The particular term is
`O(t⁴)`, because `c_N = O(t²)`. So the free regular mode is fixed *exactly* by the Hamiltonian
constraint at `O(t²)`; it is neither dropped nor set to zero.

For `K ≠ 0`, the trace-free metric is `H_r = (K̂K̂ − δ/3) h_S + (K ⊗ w + w ⊗ K) + h`, with

```
h_S(t) = ∫₀^t [6(ν_r′ − ν_r/t′) + 3i (K·c_M^(r))/K²] dt′            (M, longitudinal)
w(t)   = (2i/K²) ∫₀^t P_⊥ c_M^(r) dt′,   P_⊥ = δ − K̂K̂                (M, transverse)
h(t)   = −2 ∫₀^t (t′/t) sin(K(t−t′))/K · TT[c_T^(r)](t′) dt′         (T, TT part)
```

At `K = 0`, with no projector anywhere: `A(t) = −2∫₀^t (t′/t)(t − t′) c_T^(r)(t′) dt′`, and
`c_M^(r)` vanishes identically. For `r = 2` this is checked in `second_order_consistency.py`;
for `r = 3` it is not separately checked at `K = 0`, since only `K = k ≠ 0` occurs there.

The Green functions, `K → 0` limits and the `O(t⁴)` property are checked in
`solution/check_green.py`.

**Complete integrands.** `c^(2)` (`c_N`, `c_M`, `c_T`) are explicit elementary functions of
`t` and the momenta: products of `sin`, `cos` of `k_a t/√3` and rational functions. They are in
`solution/results/second_order_sources.txt`. `c^(3)` at `(k, q, −q)` is explicit in the
first-order closed forms and in the second-order kernels `K₂(k,q)`, `K₂(k,−q)` and `K₂(q,−q)`.
The last of these is the `K = 0` kernel, never a limit. All are in
`solution/results/third_order_sources.txt`. Second-order field kernels are therefore single
integrals of elementary functions. Third-order field kernels are double integrals.

**Boundedness near `K → 0` (needed in §7).** `c_M^(2)` is analytic in the momentum components
and vanishes at `K = 0`, so `c_M = O(K)`. Then `(K̂K̂ − δ/3)h_S` and `K ⊗ w` are bounded, with
a direction-dependent limit, and `ν₂` has no projector at all. The second-order kernels are
therefore bounded on a neighbourhood of `K = 0` and discontinuous only at `K = 0`. That is a set
of measure zero in `C^(22)`, and the exact `K = 0` sector is used where `K = 0` is exact (`C^(13)`).

### 4.3 Second order is consistent with all ten equations [checked]

The solution of §4.2 used only (N), the two parts of (M), and the TT part of (T). The scalar
and vector parts of (T) and the whole of (Ham) were not used. Inserting the solution as exact
Taylor series in `t`, with general symbolic momenta, every component of every equation vanishes
through `O(t⁴)`, for general `K` and for `K = 0` (`solution/second_order_consistency.py`; run with
`DEPTH=10` for `O(t⁸)`). §4.4 shows this check is a test of the implementation; the
consistency itself is proved.

### 4.4 Proof that the unused equations hold at every order

Let `𝔼_μν = G_μν − 8πG T_μν`. Matter conservation is exact by construction (L1), so
`∇^μ𝔼_μν = 0` identically. Suppose all residuals vanish below order `r`. Then the order-`r`
residual obeys the identity linearised about the background `a²(−dt² + dx²)`, with `ℋ = 1/t`.
Imposed at order `r`:

- `𝔼^t_i = 0`, from (M);
- `𝔼^t_t = ½𝔼`, from (N), because (N) is `R_tt = 8πG T_tt`, that is `𝔼_tt − ½g_tt𝔼 = 0`.
  Writing `τ` for the trace part (`𝔼^i_i = 3τ`), this gives `𝔼^t_t = 3τ`;
- for `K ≠ 0`, the TT part of TF `𝔼^i_j`, from (T).

Let the unknown scalar and vector trace-free parts be
`e_S (K̂K̂ − δ/3) + (K̂ ⊗ w + w ⊗ K̂)`, with `w ⊥ K`. Then:

- **Spatial component.** `∂_j𝔼^j_i = 0`, so `K_j𝔼^j_i = 0`. Its transverse part gives `K w = 0`,
  so `w = 0`. Its longitudinal part gives `⅔e_S + τ = 0`.
- **Time component.** `∂_t𝔼^t_t + 4ℋ𝔼^t_t − ℋ𝔼^t_t − ℋ𝔼^i_i = 0` gives `3τ′ + 6τ/t = 0`.
  So `τ ∝ t⁻²`.

The order-`r` residual is analytic at `t = 0`: the fields are even analytic with `ν = O(t²)`,
and the second-order series check confirms that negative powers cancel. Hence `τ = 0`, then
`e_S = 0`, and `(Ham) ∝ 𝔼^t_t = 0`.

At `K = 0`, (T) is imposed in full, so TF `𝔼^i_j = 0` and again `τ = 0`. The spatial component
gives `∂_t𝔼^t_i + 4ℋ𝔼^t_i = 0` (using `𝔼^i_t = −𝔼^t_i` on the background), so
`𝔼^t_i ∝ t⁻⁴` and `𝔼^t_i = 0` by regularity. This is the statement that the
momentum source vanishes at `K = 0`, checked explicitly at second order.

**So the triangular solution satisfies all ten equations at every order.** This includes
third order, given the exactly extracted sources. The Taylor-series check of §4.3 is an
independent test of the implementation, not a logical requirement.

### 4.5 Dilation Ward identity (residual relabelling) [checked]

A constant long mode is removed by `X → e^{ζ_L}X`. Hence `K₂(k, 0) = −k·∂_k K₁(k)` for `ν` and
for every component of `H`. Inserting this prediction into all ten equations gives exact
cancellation in closed form. It also lies in the regular class
(`ν₂ = 2k²t²/9 + …`, which equals `c₂` at `k₂ = 0`). Perturbed predictions fail
(`solution/check_ward_dilation.py`).

---

## 5. Observable through third order [evaluated at operator level]

### 5.1 Tetrad: symmetric root at the observation event

The Fermi–Walker transport of `g` along `u` is parallel transport of `g̃ = −dt² + γ dX²` along
`∂_t` (lemma L5): `dE_a^i/dt = −S^i_j E_a^j`. With initial orientation `E(t_I) = γ(t_I)^{−1/2}`,
at first order `E₁′ = −½γ₁′`, so `E₁(t) = −½γ₁(t)`. This is the first-order part of `γ(t)^{−1/2}`.

Hence `E_FW = γ_o^{−1/2} R` with `R = 1 + ε²Ω/2 + O(ε³)` and `Ω` antisymmetric. Then:

- `Θ_FW^(1)` and `Θ_FW^(2)` equal the symmetric-root maps;
- `Θ_FW^(3) = Θ_sym^(3) + 3𝒪_Ω Θ^(1)`, where `𝒪_Ω` is an infinitesimal rotation;
- the induced change of `a^(3)_ℓm` is `3D(Ω)a^(1)_ℓm`, with `D(Ω)` anti-Hermitian on each
  ℓ-multiplet, so `Re Σ_m a^(1)*_ℓm D(Ω)a^(1)_ℓm = 0` in every realisation.

So **`C^(22)` and `Re C^(13)` are individually unchanged.** All kernels are built with
`e(t_o) = −γ_o^{−1/2} n̂ = −e^{−ζ_o} exp(−H_o/2) n̂`. This direction satisfies the null
condition exactly.

### 5.2 Expansion about the background ray

Take `X̄(t) = n̂(t_o − t)` with `ē = −n̂`, and write `X = X̄ + εX₁ + ε²X₂/2`. All fields are
evaluated at `(t, X̄(t))`, and `σ(a,b) = σ_ij a^i b^j`.

**Integrand** `F = σ_ij e^i e^j` (`σ`, `S` and `Γ` vanish on the background):

```
F₁ = σ₁,nn
F₂ = σ₂,nn + 2X₁·∂σ₁,nn − 4σ₁(n, e₁)
F₃ = σ₃,nn + 3X₁·∂σ₂,nn − 6σ₂(n, e₁) + 3X₂·∂σ₁,nn + 3(X₁·∂)²σ₁,nn
     − 12(X₁·∂)σ₁(n, e₁) − 6σ₁(n, e₂) + 6σ₁(e₁, e₁)
```

`X₃` and `e₃` never enter.

**Ray**, `X_r″ = A_r`, with `X_r(t_o) = 0`, `X_r′(t_o) = e_r(t_o)`:

```
X_r(t) = e_r(t_o)(t − t_o) + ∫_t^{t_o} (t′ − t) A_r(t′) dt′,
A₁ = 2S₁n − Γ₁(n,n) − nσ₁,nn,
A₂ = [2S₂n − Γ₂(n,n) − nσ₂,nn] + 2(X₁·∂)[2S₁n − Γ₁(n,n) − nσ₁,nn]
     + 2[−2S₁e₁ + 2Γ₁(n, e₁) + e₁σ₁,nn + 2nσ₁(n, e₁)].
```

The initial directions are the Taylor coefficients of `−e^{−ζ} exp(−H/2) n̂` at `(t_o, 0)`.

**Emission root** (`t̄ = t̄⋆`; `φ(t) = ν(t, X(t))`, with `φ₁ = ν₁`, `φ₂ = ν₂ + 2X₁·∂ν₁`,
`φ₃ = ν₃ + 3X₁·∂ν₂ + 3X₂·∂ν₁ + 3(X₁·∂)²ν₁`; derivatives are along the ray parameter):

```
d₁ = −t̄φ₁;   u₁ = d₁/t̄
d₂ = t̄[u₁² − φ₂ − 2φ₁′d₁]
d₃ = t̄[3u₁u₂ − 2u₁³ − φ₃ − 3φ₂′d₁ − 3φ₁′d₂ − 3φ₁″d₁²]    (d₃ is not needed for I₃)
```

**Moving endpoint:**

```
I₁ = ∫F₁,   I₂ = ∫F₂ − 2F₁d₁,   I₃ = ∫F₃ − 3F₁d₂ − 3F₂d₁ − 3F₁′d₁²
```

The integrals run over `[t̄⋆, t_o]`, and the endpoint terms are evaluated at `t̄⋆`.

**Clock.** On the observer line `ω(t) = ν(t, 0)`; `W₁(T) = ∫₀^T tω₁`,
`W₂(T) = ∫₀^T t(ω₂ + ω₁²)`. The initial event shifts are `d_{I1} = −t̄_Iω₁(t̄_I)` and
`d_{I2} = t̄_I[u_{I1}² − ω₂ − 2ω₁′d_{I1}]`. Then:

```
D₁ = [t̄_I d_{I1} − W₁(t̄_o) + W₁(t̄_I)]/t̄_o
D₂ = [t̄_I d_{I2} + d_{I1}² − W₂(t̄_o) − 2t̄_oω₁(t̄_o)D₁ + W₂(t̄_I) + 2t̄_Iω₁(t̄_I)d_{I1} − D₁²]/t̄_o
```

`T_I → ∞` is the limit `t̄_I → 0`.

**Composition.** Let `Θ_c` be the normalised map for a coordinate observation event `(t_o, 0)`:

```
Θ^(1) = Θ_c^(1),   Θ^(2) = Θ_c^(2) + 2D₁∂Θ_c^(1),
Θ^(3) = Θ_c^(3) + 3D₂∂Θ_c^(1) + 3D₁∂Θ_c^(2) + 3D₁²∂²Θ_c^(1),     ∂ = ∂/∂t_o.
```

**Check:** `solution/check_observable_expansion.py`. Every line above is compared with the
direct ε-expansion of its exact definition. The fields used are arbitrary: general jets about
the ray, with arbitrary functions of `t` as coefficients.

### 5.3 Normalisation and dipole

With `T_obs ∝ e^{−I}` and `⟨·⟩` the sky mean (random, not an ensemble mean):

```
Θ^(1) = ⟨I₁⟩ − I₁
Θ^(2) = I₁² − I₂ − 2⟨I₁⟩I₁ + 2⟨I₁⟩² − ⟨I₁²⟩ + ⟨I₂⟩
Θ^(3) = −I₃ + 3I₁I₂ − I₁³ + 3⟨I₁⟩I₁² − 3⟨I₁⟩I₂ − 3⟨I₂⟩I₁ + 3⟨I₁²⟩I₁ − 6⟨I₁⟩²I₁
        + ⟨I₃⟩ − 3⟨I₁I₂⟩ + ⟨I₁³⟩ + 6⟨I₁⟩³ − 6⟨I₁⟩⟨I₁²⟩ + 6⟨I₁⟩⟨I₂⟩
```

These are checked in `verification/normalization_expansion.py`. They are equivalent to orange's
`t_r`/`b_r` form with `t₁ = −I₁`, `t₂ = I₁² − I₂`. Dipole subtraction is a final projection
removing `a_1m` only; `C_{ℓ≥2}` is unaffected, and no deboost is implied.

### 5.4 Fourier dictionary [representation]

Per mode `a` (momentum `k`, `μ_a = k̂·n̂`), on the ray the phase is `e^{iλ k·n̂}` with
`λ = t_o − t`:

```
σ_a,nn = ½Q_a′ k²(μ_a² − ⅓)
Γ_a(n,n) = i[2(k·n)n − k] + (i/2)[2(k·n)H_a n − k H_a,nn]
A_a = H_a′n − ½H_a′,nn n − Γ_a(n,n)
e_a(t_o) = (1 + ½H_a(t_o)) n̂
d_a = −t̄⋆ ν_a(t̄⋆) e^{iλ⋆k·n̂}
```

Second-order composite fields follow from the multilinear rule, for example

```
σ_ab = ½H_ab′ + ¼(H_aH_b + H_bH_a)′ + H_a′ + H_b′.
```

The second-order kernel of `I` is then

```
I_ab = ∫_{t̄⋆}^{t_o} dt { σ_ab,nn + i(k_b·X_a)σ_b,nn + i(k_a·X_b)σ_a,nn − 2σ_b(n, e_a) − 2σ_a(n, e_b) }
       − [F_a(t̄⋆) d_b + F_b(t̄⋆) d_a],
```

where `X_a` and `e_a = X_a′` are the single-mode ray solutions of §5.2 (one nested time
integral). `Θ_c,ab` follows from §5.3 by the multilinear rule:

```
Θ_ab = I_aI_b − I_ab − ⟨I_a⟩I_b − ⟨I_b⟩I_a + 2⟨I_a⟩⟨I_b⟩ − ⟨I_aI_b⟩ + ⟨I_ab⟩,
```

plus `D_a ∂Θ_b + D_b ∂Θ_a` from the clock. The third-order kernel at `(k, q, −q)` is the
multilinear coefficient of `Θ^(3)/6`. Every ingredient is defined above, and nothing is
truncated: post-Born displacements, time delay, the emitter and observer endpoints, the tetrad,
the clock, and the normalisation products, including the `ℓ' = 0` (sky-mean) couplings.

**Squeezed check of the observable.** At `q = 0` the long mode contributes `X_L = −ζ_L n̂λ`
(slower light), `e_L = ζ_L n̂`, `σ_{sL} = H_s′ + ½H_{sL}′`, and nothing else. With §4.5 this
gives `I_{sL} = −k·∂_k I_s` exactly, by the same steps written out by hand. So
`Θ[ζ_s + ζ_L] = Θ[ζ_s(e^{−ζ_L}·)]` at second order, including frame, ray speed and endpoints.

---

## 6. Gauge invariance at fixed physical labels

Every coordinate and label used is a physical scalar or frame:

- `t`: the fluid velocity potential, normalised at `T → ∞`;
- `X`: fluid labels in the marked conformally Cartesian primordial chart;
- the emission surface `T = T⋆`;
- the observer: the fluid line `X = 0`;
- the clock: proper time from `T = T_I`;
- the tetrad: the symmetric root of `γ` in the labels, then rotation-equivalent to Fermi–Walker.

A diffeomorphism acting on fields and labels together leaves all of these, and hence `Θ(n̂)`,
invariant. The temperature frame is just the representation in which they are coordinates. The
only residual freedom is the conformal group of the primordial chart, fixed in §2. Its dilation
generator is verified to act correctly (§4.5, §5.4). At first order the result equals the
Newtonian-gauge computation for all ℓ (purple). Orange's independent quadratic and cubic
constructions are in the same chart.

---

## 7. Spectrum

### 7.1 Gaussian reduction (general spectrum) [evaluated]

Let `𝒦_r(k₁..k_r; n̂)` be the kernel of `Θ^(r)` (observer at the origin, so there are no phases).
Define `A^(r)_ℓm = ∫dΩ Y*_ℓm 𝒦_r`. Then:

```
C^(11)_ℓ = (2ℓ+1)⁻¹ Σ_m ∫_k P(k)|A^(1)_ℓm(k)|²
C^(22)_ℓ = 2(2ℓ+1)⁻¹ Σ_m ∫_k∫_q P(k)P(q)|A^(2)_ℓm(k,q)|²
C^(13)_ℓ = 3(2ℓ+1)⁻¹ Σ_m ∫_k∫_q P(k)P(q) A^(1)_ℓm(k) A^(3)_ℓm(k,q,−q)*
C_ℓ = ε²C^(11) + ε⁴[¼C^(22) + ⅓ Re C^(13)] + O(ε⁶),   ℓ ≥ 2
```

- The disconnected part of 22 is `∝ δ_ℓ0` by isotropy.
- In 13 all three pairings survive, including the internal `(q, −q)` pair. It is a pairing of
  the fourth moment and is not removed by any normal ordering.
- Odd orders vanish.

The Wick factors and Fourier/shell measures are checked in orange's `validation/check_algebra.py`.

### 7.2 Monochromatic spectrum [representation]

Only after the complete observable is built, take

```
𝒫_g(k) = σ² δ(ln(k/k0)),   P(k) = (2π²σ²/k0²) δ(k − k0),   ⟨g²⟩ = σ² held fixed.
```

The background and the physical labels `T⋆`, `T_I`, `τ_o` are also held fixed. Then
`∫_k P f = (σ²/4π)∫dΩ f(k0k̂)`. Using joint rotation invariance
(`∫dΩ_k dΩ_q → 8π²∫_{−1}^{1}dμ`, with `k̂ = ẑ` and `q̂_μ = (√(1−μ²), 0, μ)`) and the addition
theorem:

```
C^(11)_ℓ = 4πσ² T_ℓ(k0)²                                                        [evaluated: ℓ=2,3 closed form]
C^(22)_ℓ = (σ⁴/4π) ∫_{−1}^{1}dμ ∫dΩ_n ∫dΩ_{n′} P_ℓ(n̂·n̂′) 𝒦₂(k0ẑ, k0q̂_μ; n̂) 𝒦₂(k0ẑ, k0q̂_μ; n̂′)*
C^(13)_ℓ = (3σ⁴/2) i^ℓ T_ℓ(k0) ∫_{−1}^{1}dμ [∫dΩ_n P_ℓ(ẑ·n̂) 𝒦₃(k0ẑ, k0q̂_μ, −k0q̂_μ; n̂)]*
```

The final integrals are:

- **`C^(22)`**: one angle `μ`, two sky directions, and in each `𝒦₂` the line-of-sight time
  plus one nested ray time and one field-kernel time;
- **`C^(13)`**: one angle, one sky direction, and up to five nested times
  (line of sight, `X₂`, `A₂`, and two field-kernel levels).

The domains are compact, and the integrands bounded:

- the field and ray integrands are continuous on compact time intervals;
- the `1/K²` projectors are bounded by §4.2;
- `K₂(q, −q)` is the exact `K = 0` kernel;
- the removable `ln|α|` resonance singularities appear only after Si/Ci reduction, never in
  these representations.

So both integrals exist at all finite `x⋆ < x_o`. The shell substitution is legitimate, because
the restricted integrand is defined everywhere except at `μ = ±1` and bounded. The result
depends on `k0` only through `x⋆ = k0t̄⋆/√3`, `x_o` and `x_I`.

---

## 8. Requirement-by-requirement status

| Problem-statement requirement | Status |
|---|---|
| Gaussian variable; no higher-order modes; constrained initial data at every order | §2, §4 [evaluated / representation]; constraints proved at 2nd order, operators at 3rd |
| Higher orders not set to zero | §2, §4.2 (forced `t²` data, `K = 0` modes kept) |
| Emission surface, collisionless propagation, observer on fluid line, proper time, FW tetrad | §1, §5 |
| Exact `T⋆E_o/E_e`; normalise; dipole last; keep observer/endpoint terms | §1, §5.3; `T_fluid(o)` cancels exactly; the remaining terms are all kept |
| Θ through third order incl. vector/tensor | §4.2 (`w`, `h`), §5.2–5.4 [operator level + Fourier dictionary] |
| `C_ℓ` with all interference and normalisations | §7 [representation] |
| Constraints, fluid normalisation, null/geodesic, gauge invariance verified | §1 (exact), §4.3–4.5, §5.2, §6 |
| Monochromatic only after constructing the observable; what is held fixed | §7.2 |
| Limits and their exchange justified; unresolved steps identified | §7.2 (finite times); §9 |

## 9. Ledger of unevaluated integrals and open steps

1. **Second-order field kernels.** Single integrals of elementary functions. Si/Ci reduction is
   mechanical with the antiderivative recursion of `verification/transfer_closed_form.py`. It is
   not generated. **[open: closed form]**
2. **Third-order field kernels at `(k,q,−q)`.** Double integrals. The source is explicit and the
   operators are confirmed. Consistency with all ten equations is *proved* by §4.4. An explicit
   implementation test, analogous to the second-order Taylor check, is **[open]**.
3. **Θ^(2), Θ^(3) kernels.** Defined above. A symbolic generator emitting them as nested
   `Integral` objects is **[open]**; the observable was verified at the operator level.
4. **Evaluation of `C^(22)`, `C^(13)` (monochromatic).** **[open]** beyond the representation.
   Resonance recombination (scalar `μ = ±1`, tensor/vector `μ = −1/3`) must be done before the
   μ-integral if closed forms are sought.
5. **Late-time and other asymptotic limits.** **[open]**. They are non-uniform at the
   resonances, and only the angle-first order of limits is admissible.
6. **Non-monochromatic spectra.** Infrared/ultraviolet convergence, and the separate `q → 0`
   sensitivity of 22 and 13. **[open]**. The dilation identity predicts the cancellation
   structure.
7. **Perturbative domain.** Single emission intersection and no caustics hold order by order
   only; there is no statement about convergence in `ε`.
