"""Exact, non-perturbative lemmas behind the "temperature-frame" formulation.

Coordinates (t, X^i): t = phi/Tref is the radiation-fluid velocity potential
(T u_mu = -d_mu phi, exact for an irrotational barotropic fluid) and X^i are
fluid labels.  Then, exactly,

    g = N^2 ( -dt^2 + gamma_ij dX^i dX^j ),   N = Tref / T_phys,

i.e. gtilde = g / N^2 is synchronous-comoving and u = N^{-1} d_t.

Checked here with ARBITRARY functions N(t,X), gamma_ij(t,X):
  L1  div T = 0 for the radiation fluid  <=>  d_t det(gamma) = 0
      (Euler is identically satisfied; continuity = "expansion-free in gtilde").
  L2  photon redshift law in gtilde:  d ln Etilde / dt = -(1/2) d_t gamma_ij e^i e^j,
      and the exact direction equation for e^i = p^i / p^0.
  L3  E = Omega * Etilde with Omega = 1/N, hence T_obs = T_fluid(obs) * Etilde_o/Etilde_e.
  L4  exact Hamiltonian constraint in these variables (det gamma = 1 family):
      N^2 (R[gamma] - 4 Lap ln N - 2 |D ln N|^2 - sigma_ij sigma^ij) + 6 (d_t N)^2 = 6 a1^2,
      where 8 pi G alpha Tref^4 = 3 a1^2 fixes the background a = a1 t.
"""
import sympy as sp
from gr import christoffel, ricci, assert_zero

t, x, y, z = X = sp.symbols('t x y z', real=True)
c = sp.Symbol('c', positive=True)  # rho = c / N^4  (c = alpha Tref^4)

# ---------------------------------------------------------------- L1
N = sp.Function('N')(*X)
gfun = {}
for i in range(3):
    for j in range(i, 3):
        gfun[i, j] = gfun[j, i] = sp.Function(f'g{i}{j}')(*X)
gam = sp.Matrix(3, 3, lambda i, j: gfun[i, j])
g = sp.zeros(4, 4)
g[0, 0] = -N**2
for i in range(3):
    for j in range(3):
        g[i + 1, j + 1] = N**2 * gam[i, j]
detg = gam.det()
ginv = sp.zeros(4, 4)
ginv[0, 0] = -1 / N**2
gaminv = gam.adjugate() / detg
for i in range(3):
    for j in range(3):
        ginv[i + 1, j + 1] = gaminv[i, j] / N**2
rho = c / N**4
p = rho / 3
u_up = [1 / N, 0, 0, 0]
u_dn = [-N, 0, 0, 0]
Tmix = sp.Matrix(4, 4, lambda m, n: (rho + p) * u_up[m] * u_dn[n] + (p if m == n else 0))
Tup = sp.Matrix(4, 4, lambda m, n: (rho + p) * u_up[m] * u_up[n] + p * ginv[m, n])
sqrtg = N**4 * sp.sqrt(detg)
# div of a symmetric tensor: nabla_mu T^mu_nu = (1/sqrt g) d_mu (sqrt g T^mu_nu) - 1/2 d_nu g_ab T^ab
div = []
for n_ in range(4):
    s = sum(sp.diff(sqrtg * Tmix[m, n_], X[m]) for m in range(4)) / sqrtg
    s -= sp.Rational(1, 2) * sum(sp.diff(g[a, b], X[n_]) * Tup[a, b] for a in range(4) for b in range(4))
    div.append(s)
for i in range(1, 4):
    assert_zero(div[i], f"L1 Euler (nu={i}) identically satisfied for arbitrary N, gamma_ij",
                simp=lambda e: sp.simplify(sp.together(e)))
expected0 = -sp.Rational(2, 3) * rho * sp.diff(detg, t) / detg
assert_zero(div[0] - expected0, "L1 continuity: nabla_mu T^mu_t = -(2/3) rho d_t ln det(gamma)",
            simp=lambda e: sp.simplify(sp.together(e)))

# ---------------------------------------------------------------- L2
# synchronous metric gtilde = -dt^2 + gamma_ij dX dX with arbitrary gamma
gt = sp.zeros(4, 4)
gt[0, 0] = -1
for i in range(3):
    for j in range(3):
        gt[i + 1, j + 1] = gam[i, j]
Gam, _ = christoffel(gt, X)
for i in range(3):
    for j in range(3):
        assert_zero(Gam[0][i + 1][j + 1] - sp.diff(gam[i, j], t) / 2, f"L2 Gamma^t_{{X{i+1}X{j+1}}} = (1/2) d_t gamma_ij")
assert all(sp.simplify(Gam[0][0][m]) == 0 for m in range(4))
print("[ok] L2 Gamma^t_{t mu} = 0, hence dp^t/dlambda = -(1/2) d_t gamma_ij p^i p^j (exact)")
# spatial: Gamma^i_{tj} = (1/2) gamma^{ik} d_t gamma_kj  (=> exact e^i equation)
for i in range(3):
    for j in range(3):
        rhs = sum(gaminv[i, k] * sp.diff(gam[k, j], t) for k in range(3)) / 2
        assert_zero(Gam[i + 1][0][j + 1] - rhs, f"L2 Gamma^X{i+1}_{{t X{j+1}}} = (1/2) gamma^ik d_t gamma_kj",
                    simp=lambda e: sp.simplify(sp.together(e)))
# spatial Christoffels of gtilde coincide with those of gamma (synchronous, g_tt = -1)
G3l, _ = christoffel(gam, list(X[1:]))
for i in range(3):
    for j in range(3):
        for k_ in range(3):
            assert sp.simplify(sp.together(Gam[i + 1][j + 1][k_ + 1] - G3l[i][j][k_])) == 0
print("[ok] L2 Gammat^i_jk = Gamma[gamma]^i_jk")
# exact ray equation from the geodesic equation, p = p^t (1, e^i):
#   de^i/dt = [ -Gamma^i_ab + e^i Gamma^t_ab ] (1,e)^a (1,e)^b
ev_ = [1] + list(sp.symbols('e1 e2 e3'))
sig_m = sp.Matrix(3, 3, lambda i, j: sp.diff(gam[i, j], t) / 2)
sig_up = gaminv * sig_m
for i in range(3):
    geo = (-sum(Gam[i + 1][a_][b_] * ev_[a_] * ev_[b_] for a_ in range(4) for b_ in range(4))
           + ev_[i + 1] * sum(Gam[0][a_][b_] * ev_[a_] * ev_[b_] for a_ in range(4) for b_ in range(4)))
    claim = (-2 * sum(sig_up[i, j] * ev_[j + 1] for j in range(3))
             - sum(G3l[i][j][k_] * ev_[j + 1] * ev_[k_ + 1] for j in range(3) for k_ in range(3))
             + ev_[i + 1] * sum(sig_m[j, k_] * ev_[j + 1] * ev_[k_ + 1] for j in range(3) for k_ in range(3)))
    assert sp.simplify(sp.together(sp.expand(geo - claim))) == 0, i
print("[ok] L2 exact ray equation de^i/dt = -2 sigma^i_j e^j - Gamma3^i_jk e^j e^k + e^i sigma_jk e^j e^k")
# redshift law: d ln p^t / dt = -Gamma^t_ab (1,e)(1,e) = -sigma_ij e^i e^j
red = -sum(Gam[0][a_][b_] * ev_[a_] * ev_[b_] for a_ in range(4) for b_ in range(4))
assert sp.simplify(red + sum(sig_m[j, k_] * ev_[j + 1] * ev_[k_ + 1] for j in range(3) for k_ in range(3))) == 0
print("[ok] L2 exact redshift law d ln Etilde/dt = -sigma_ij e^i e^j")

# ---------------------------------------------------------------- L3
# Conformal rescaling g = N^2 gt: null geodesics coincide as point sets; with affine
# parameters related by d lambda~ = N^{-2} d lambda one has p~^mu = N^2 p^mu, u~^mu = N u^mu.
# E = -g(u,p) = -N^2 gt(u~/N, p~/N^2) = N^{-1} * (-gt(u~,p~)) = Etilde / N  (Omega = 1/N).
# We check the affine-parameter claim: the geodesic acceleration of p~ = N^2 p vanishes in gt
# whenever p is a null geodesic tangent of g, via the conformal Christoffel identity.
Gg, _ = christoffel(g, X)
GamT = Gam
lnN = sp.log(N)
# Direct formula: Gamma[g]^a_bc = Gamma[gt]^a_bc + d^a_b w_c + d^a_c w_b - gt_bc gt^ad w_d, w = d ln N
gtinv = sp.zeros(4, 4)
gtinv[0, 0] = -1
for i in range(3):
    for j in range(3):
        gtinv[i + 1, j + 1] = gaminv[i, j]
w = [sp.diff(lnN, X[m]) for m in range(4)]
ok = True
for a in range(4):
    for b in range(4):
        for cc in range(b, 4):
            rhs = (GamT[a][b][cc] + (w[cc] if a == b else 0) + (w[b] if a == cc else 0)
                   - gt[b, cc] * sum(gtinv[a, d] * w[d] for d in range(4)))
            if sp.simplify(sp.together(Gg[a][b][cc] - rhs)) != 0:
                ok = False
assert ok
print("[ok] L3 conformal Christoffel identity Gamma[g] = Gamma[gt] + C(w)")
# (i) if p is an affinely parametrised null geodesic of g, then pt = N^2 p (dlambda~ = dlambda/N^2)
#     is an affinely parametrised geodesic of gt.  Pointwise identity at an arbitrary point with
#     arbitrary components p^a, using dp^a/dlambda = -Gamma[g]^a_bc p^b p^c and g(p,p) = 0.
pv = list(sp.symbols('p0 p1 p2 p3'))
null_p0 = sp.sqrt(sum(gt[i, j] * pv[i] * pv[j] for i in range(1, 4) for j in range(1, 4)))
psub = {pv[0]: null_p0}  # gt(p,p)=0 <=> g(p,p)=0 (future-directed)
for A_ in range(4):
    dp = -sum(Gg[A_][b_][c_] * pv[b_] * pv[c_] for b_ in range(4) for c_ in range(4))
    dlnN = sum(pv[b_] * w[b_] for b_ in range(4))
    # d pt^a / d lambda~ = N^2 d(N^2 p^a)/d lambda = N^4 (2 (p.w) p^a + dp^a/dlambda)
    acc = N**4 * (2 * dlnN * pv[A_] + dp) + sum(GamT[A_][b_][c_] * N**4 * pv[b_] * pv[c_] for b_ in range(4) for c_ in range(4))
    assert sp.simplify(sp.together(sp.expand(acc.subs(psub)))) == 0, A_
print("[ok] L3 null geodesics of g map to affinely parametrised geodesics of gt with pt = N^2 p")
# (ii) E = -g(u,p) with u = ut/N, p = pt/N^2  =>  E = Etilde / N
Et = -sum(gt[a_, b_] * [1, 0, 0, 0][a_] * pv[b_] for a_ in range(4) for b_ in range(4))  # Etilde for ut=d_t, pt=pv
E_ = -sum(g[a_, b_] * ([1, 0, 0, 0][a_] / N) * (pv[b_] / N**2) for a_ in range(4) for b_ in range(4))
assert sp.simplify(E_ - Et / N) == 0
print("[ok] L3 E = Etilde/N; with N_e = Tref/T_*, N_o = Tref/T_fluid(o):")
Ts, To_, Eo_t, Ee_t, Tref = sp.symbols('T_star T_o Etilde_o Etilde_e T_ref', positive=True)
Tobs = Ts * (Eo_t / (Tref / To_)) / (Ee_t / (Tref / Ts))
assert sp.simplify(Tobs - To_ * Eo_t / Ee_t) == 0
print("[ok] L3 T_obs = T_* E_o/E_e = T_fluid(o) Etilde_o/Etilde_e  (T_* and Tref cancel exactly)")

# ---------------------------------------------------------------- L4
a1 = sp.Symbol('a1', positive=True)
A = sp.Function('A')(*X)
B = sp.Function('B')(*X)
gam1 = sp.diag(sp.exp(2 * A), sp.exp(2 * B), sp.exp(-2 * A - 2 * B))  # det = 1: expansion-free family
gg = sp.zeros(4, 4)
gg[0, 0] = -N**2
for i in range(3):
    gg[i + 1, i + 1] = N**2 * gam1[i, i]
Gm, ggi = christoffel(gg, X)
Ric = ricci(Gm, X)
Rs = sum(ggi[m, n] * Ric[m, n] for m in range(4) for n in range(4))
G00_up = sum(ggi[0, m] * ggi[0, n] * (Ric[m, n] - gg[m, n] * Rs / 2) for m in range(4) for n in range(4))
# Hamiltonian constraint: 2 G_{mu nu} n^mu n^nu = 2 N^2 G^{tt} must equal 2 rho = 6 a1^2 / N^4 (8piG=1, c=3a1^2)
Hlhs = 2 * N**2 * G00_up * N**4
# right-hand formula built from 3D objects of gamma:
X3 = X[1:]
g3 = gam1
G3, g3i = christoffel(g3, X3)
R3 = ricci(G3, X3)
R3s = sum(g3i[i, j] * R3[i, j] for i in range(3) for j in range(3))
lnN3 = sp.log(N)
grad = [sp.diff(lnN3, X3[i]) for i in range(3)]
lap = sum(g3i[i, j] * (sp.diff(lnN3, X3[i], X3[j]) - sum(G3[k][i][j] * grad[k] for k in range(3)))
          for i in range(3) for j in range(3))
gradsq = sum(g3i[i, j] * grad[i] * grad[j] for i in range(3) for j in range(3))
sig = sp.Matrix(3, 3, lambda i, j: sp.diff(g3[i, j], t) / 2)
sigsq = sum(g3i[i, k] * g3i[j, l] * sig[i, j] * sig[k, l]
            for i in range(3) for j in range(3) for k in range(3) for l in range(3))
formula = N**2 * (R3s - 4 * lap - 2 * gradsq - sigsq) + 6 * sp.diff(N, t)**2
assert_zero(sp.expand(Hlhs - formula), "L4 exact Hamiltonian constraint: 2N^6 G^tt = N^2(R-4Lap lnN-2|D lnN|^2-sigma^2)+6 Ndot^2",
            simp=lambda e: sp.simplify(sp.expand(e)))
print("    (field equation: this equals 6 a1^2;  background N=a1 t, gamma=delta satisfies it)")
print("ALL EXACT LEMMAS VERIFIED")
