"""Lemma L5: Fermi-Walker transport in g  <=>  parallel transport in the temperature frame.

g = N^2 gt,  gt = -dt^2 + gamma_ij dX dX (arbitrary N, gamma), fluid u = N^{-1} d_t with
acceleration a_b = h_b^c d_c ln N (Euler, lemma L1).  For a spatial vector field
e = N^{-1} et (et orthogonal to u, i.e. et^t = 0), the Fermi-Walker operator

    FW[e]^a = u^b nabla_b e^a - u^a (a_b e^b)        (valid for e orthogonal to u)

satisfies, identically,

    FW[e]^a = N^{-2} ( ut^b nablat_b et^a ),   ut = d_t,

so FW transport of the observer tetrad in g is parallel transport along the geodesic ut in gt,
i.e. d E^i/dt = -sigma^i_j E^j (Gammat^i_{tj} = sigma^i_j).  Checked for arbitrary functions.
"""
import sympy as sp
from gr import christoffel

t, x, y, z = X = sp.symbols('t x y z', real=True)
N = sp.Function('N')(*X)
gf = {}
for i in range(3):
    for j in range(i, 3):
        gf[i, j] = gf[j, i] = sp.Function(f'g{i}{j}')(*X)
gam = sp.Matrix(3, 3, lambda i, j: gf[i, j])
et = [0] + [sp.Function(f'e{i}')(*X) for i in range(1, 4)]

gt = sp.zeros(4, 4); gt[0, 0] = -1
for i in range(3):
    for j in range(3):
        gt[i + 1, j + 1] = gam[i, j]
g = N**2 * gt
Gt, _ = christoffel(gt, X)
G, _ = christoffel(g, X)

u = [1 / N, 0, 0, 0]
e = [c / N for c in et]
# acceleration of u in g, lowered: a_b = g_bc (u^d nabla_d u^c)
acc_up = [sum(u[d] * sp.diff(u[c], X[d]) for d in range(4))
          + sum(G[c][d][f] * u[d] * u[f] for d in range(4) for f in range(4)) for c in range(4)]
acc_dn = [sum(g[b, c] * acc_up[c] for c in range(4)) for b in range(4)]
# Euler consistency: a_b = h_b^c d_c ln N  (spatial part d_i ln N, time part 0)
for b in range(4):
    expect = 0 if b == 0 else sp.diff(sp.log(N), X[b])
    assert sp.simplify(acc_dn[b] - expect) == 0, b
print("[ok] L5a acceleration of u = N^-1 d_t in g is a_b = h_b^c d_c ln N (matches Euler for radiation)")
a_e = sum(acc_dn[b] * e[b] for b in range(4))
for A in range(4):
    lhs = (sum(u[b] * sp.diff(e[A], X[b]) for b in range(4))
           + sum(G[A][b][c] * u[b] * e[c] for b in range(4) for c in range(4))
           - u[A] * a_e)
    rhs = (sp.diff(et[A], t) + sum(Gt[A][0][c] * et[c] for c in range(4))) / N**2
    assert sp.simplify(sp.together(lhs - rhs)) == 0, A
print("[ok] L5b FW_g[N^-1 et] = N^-2 (d_t et + Gammat^a_{t c} et^c)  for arbitrary N, gamma_ij, et")
print("     => Fermi-Walker in g  <=>  parallel transport along d_t in gt:  dE^i/dt = -sigma^i_j E^j")
