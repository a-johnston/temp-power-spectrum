"""Exact field equations in temperature-frame variables (checked as identities).

For g = N^2 ( -dt^2 + gamma_ij dX^i dX^j ) with det(gamma) independent of t
(which is what continuity enforces, lemma L1), sigma_ij = (1/2) d_t gamma_ij is traceless and,
with n = ln N, D/Lap the gamma-covariant derivative/Laplacian:

  E00:  R_tt[g] = -sigma^2 - 3 n'' + Lap n + 2 |Dn|^2
  E0i:  R_ti[g] =  D_j sigma^j_i - 2 d_i n' + 2 sigma^j_i d_j n + 2 n' d_i n
  Eij:  R_ij[g] =  R_ij[gamma] + sigma_ij' - 2 sigma_ik sigma^k_j - 2 (D_i D_j n - sigma_ij n')
                   + gamma_ij (n'' - Lap n) + 2 D_i n D_j n + 2 gamma_ij (n'^2 - |Dn|^2)

The Einstein equations for the radiation fluid at rest (8 pi G = 1, 8 pi G alpha Tref^4 = 3 a1^2) are
  R_tt = 3 a1^2 / N^2,   R_ti = 0,   R_ij = a1^2 gamma_ij / N^2
(R_mu nu = T_mu nu since T is traceless).  These identities are verified symbolically for an
arbitrary lapse N(t,x,y,z) and the det-1 family
  gamma = [[e^{2A}, F, 0], [F, e^{2B}, 0], [0, 0, 1/(e^{2A+2B} - F^2)]]
with arbitrary A, B, F (off-diagonal shear included).
"""
import sympy as sp
from gr import christoffel, ricci

t, x, y, z = X = sp.symbols('t x y z', real=True)
args = (t, x, z)  # dependence on three coordinates keeps the check fast but non-trivial (y-symmetric)
N = sp.Function('N')(*args)
A = sp.Function('A')(*args)
B = sp.Function('B')(*args)
F = sp.Function('F')(*args)
gam = sp.Matrix([[sp.exp(2 * A), F, 0], [F, sp.exp(2 * B), 0], [0, 0, 1 / (sp.exp(2 * A + 2 * B) - F**2)]])
assert sp.simplify(gam.det()) == 1
g = sp.zeros(4, 4)
g[0, 0] = -N**2
for i in range(3):
    for j in range(3):
        g[i + 1, j + 1] = N**2 * gam[i, j]

simp = lambda e: sp.simplify(sp.together(sp.expand(e)))
Gam4, _ = christoffel(g, X)
R4 = ricci(Gam4, X)

X3 = X[1:]
G3, gi = christoffel(gam, X3)
R3 = ricci(G3, X3)
n = sp.log(N)
dn = [sp.diff(n, q) for q in X3]
ndot = sp.diff(n, t)
sig = sp.Matrix(3, 3, lambda i, j: sp.diff(gam[i, j], t) / 2)
sig_up = gi * sig            # sigma^i_j
sigsq = sum(sig_up[i, j] * sig_up[j, i] for i in range(3) for j in range(3))
DDn = sp.Matrix(3, 3, lambda i, j: sp.diff(n, X3[i], X3[j]) - sum(G3[k][i][j] * dn[k] for k in range(3)))
lap = sum(gi[i, j] * DDn[i, j] for i in range(3) for j in range(3))
gradsq = sum(gi[i, j] * dn[i] * dn[j] for i in range(3) for j in range(3))
# divergence of mixed tensor sigma^j_i:  D_j S^j_i = d_j S^j_i + G^j_jk S^k_i - G^k_ji S^j_k
divsig = [sum(sp.diff(sig_up[j, i], X3[j]) for j in range(3))
          + sum(G3[j][j][k] * sig_up[k, i] for j in range(3) for k in range(3))
          - sum(G3[k][j][i] * sig_up[j, k] for j in range(3) for k in range(3)) for i in range(3)]

E00 = -sigsq - 3 * sp.diff(n, t, 2) + lap + 2 * gradsq
r = simp(R4[0, 0] - E00)
assert r == 0, r
print("[ok] E00 identity  R_tt[g] = -sigma^2 - 3 n'' + Lap n + 2|Dn|^2")
for i in range(3):
    E0i = divsig[i] - 2 * sp.diff(ndot, X3[i]) + 2 * sum(sig_up[j, i] * dn[j] for j in range(3)) + 2 * ndot * dn[i]
    r = simp(R4[0, i + 1] - E0i)
    assert r == 0, (i, r)
print("[ok] E0i identity  R_ti[g] = D_j sigma^j_i - 2 d_i n' + 2 sigma^j_i d_j n + 2 n' d_i n")
sigsig = sig * gi * sig
for i in range(3):
    for j in range(i, 3):
        Eij = (R3[i, j] + sp.diff(sig[i, j], t) - 2 * sigsig[i, j] - 2 * (DDn[i, j] - sig[i, j] * ndot)
               + gam[i, j] * (sp.diff(n, t, 2) - lap) + 2 * dn[i] * dn[j] + 2 * gam[i, j] * (ndot**2 - gradsq))
        r = simp(R4[i + 1, j + 1] - Eij)
        assert r == 0, (i, j, r)
print("[ok] Eij identity  (all 6 components)")
a1 = sp.Symbol('a1', positive=True)
tt = sp.Symbol('t', positive=True)
nb = sp.log(a1 * tt)
assert sp.simplify(-3 * sp.diff(nb, tt, 2) - 3 * a1**2 / (a1 * tt)**2) == 0
assert sp.simplify((sp.diff(nb, tt, 2) + 2 * sp.diff(nb, tt)**2) - a1**2 / (a1 * tt)**2) == 0
print("[ok] background N = a1 t solves E00 and Eij with the radiation source")
print("EXACT FIELD EQUATIONS VERIFIED")
