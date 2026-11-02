"""Nonlinear initial data (all orders in eps) at leading order in the small-t expansion.

Ansatz (Definition 1, conditions G, T, V, R):
    N = a1 t exp(nu),        nu  = t^2 nu2(X) + O(t^4)
    gamma = e^{2 zeta(X)} exp(H),  H = t^2 H2(X) + O(t^4),  tr H2 = 0
with zeta(X) an ARBITRARY function (no expansion in eps).  We solve the exact field
equations R_mu nu[g] = T_mu nu at the lowest non-trivial order in t for nu2, H2, and then
check that the momentum constraint at that order is satisfied identically (it is not used
to solve).  This is the explicit statement that the initial data are NOT obtained by
setting higher-order perturbations to zero: nu2, H2 are fixed nonlinear functionals of zeta.
"""
import sympy as sp
from gr import christoffel, ricci

t = sp.Symbol('t', positive=True)
x, y, z = Xs = sp.symbols('x y z', real=True)
X = [t, x, y, z]
a1 = sp.Symbol('a1', positive=True)
zeta = sp.Function('zeta')(*Xs)
nu2 = sp.Function('nu2')(*Xs)
h = {}
names = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2)]
for (i, j) in names:
    h[i, j] = h[j, i] = sp.Function(f'H{i}{j}')(*Xs)
h[2, 2] = -h[0, 0] - h[1, 1]
H2 = sp.Matrix(3, 3, lambda i, j: h[i, j])

lam = sp.Symbol('lam')  # bookkeeping: t^2 corrections carry lam, truncated at lam^1
N2 = a1**2 * t**2 * (1 + 2 * lam * t**2 * nu2)
gam = sp.exp(2 * zeta) * (sp.eye(3) + lam * t**2 * H2)
g = sp.zeros(4, 4)
g[0, 0] = -N2
for i in range(3):
    for j in range(3):
        g[i + 1, j + 1] = N2 * gam[i, j]

def trunc(e):
    e = sp.expand(e)
    return sp.expand(e.subs(lam, 0) + lam * sp.diff(e, lam).subs(lam, 0))

# inverse metric to O(lam)
ginv = sp.zeros(4, 4)
ginv[0, 0] = -1 / (a1**2 * t**2) * (1 - 2 * lam * t**2 * nu2)
gaminv = sp.exp(-2 * zeta) * (sp.eye(3) - lam * t**2 * H2)
for i in range(3):
    for j in range(3):
        ginv[i + 1, j + 1] = trunc(gaminv[i, j] / (a1**2 * t**2) * (1 - 2 * lam * t**2 * nu2))

n = 4
dg = [[[sp.diff(g[a, b], X[c]) for c in range(n)] for b in range(n)] for a in range(n)]
Gam = [[[0] * n for _ in range(n)] for _ in range(n)]
for a in range(n):
    for b in range(n):
        for c in range(b, n):
            s = sum(ginv[a, d] * (dg[d][b][c] + dg[d][c][b] - dg[b][c][d]) for d in range(n)) / 2
            Gam[a][b][c] = Gam[a][c][b] = trunc(s)
Ric = ricci(Gam, X, simp=trunc)
# sources: R_tt = 3 a1^2/N^2 * (N^2/N^2)... T_tt = rho N^2 = 3 a1^2 / N^2 ; T_ij = a1^2 gamma_ij / N^2 ; T_ti = 0
invN2 = (1 - 2 * lam * t**2 * nu2) / (a1**2 * t**2)
Ttt = trunc(3 * a1**2 * invN2)
E = {}
E['tt'] = trunc(Ric[0, 0] - Ttt)
for i in range(3):
    E[f't{i}'] = trunc(Ric[0, i + 1])
    for j in range(i, 3):
        E[f'{i}{j}'] = trunc(Ric[i + 1, j + 1] - a1**2 * gam[i, j] * invN2)

# background (lam^0) must vanish at the orders where it is O(t^-2); the gradient terms of the
# seed metric e^{2 zeta} delta appear at O(t^0) and must be balanced by the lam^1 t^2 terms.
eqs = []
for key, e in E.items():
    e0 = sp.expand(e.subs(lam, 0))
    e1 = sp.expand(sp.diff(e, lam).subs(lam, 0))
    tot = sp.expand(e0 + e1)
    # lowest power of t
    tot = sp.expand(tot * t**2)
    poly = sp.Poly(tot, t)
    lowest = min(m[0] for m in poly.monoms())
    coeffs = {m[0]: c for m, c in zip(poly.monoms(), poly.coeffs())}
    eqs.append((key, lowest, coeffs))
# tt and ij: lowest relevant is t^2 (i.e. O(t^0) before multiplying); t_i: t^3 (O(t^1))
unknowns = [nu2] + [h[i, j] for (i, j) in names]
lin = []
for key, lowest, coeffs in eqs:
    if key.startswith('t') and key != 'tt':
        continue
    c0 = coeffs.get(0, 0)
    assert sp.simplify(c0) == 0, (key, 'O(t^-2) background mismatch', c0)
    lin.append(coeffs.get(2, 0))
sol = sp.solve(lin, unknowns, dict=True)
assert len(sol) == 1, sol
sol = {k_: sp.simplify(v_) for k_, v_ in sol[0].items()}
print("[ok] unique solution at leading gradient order (Hamiltonian/evolution equations):")
for k_, v_ in sol.items():
    print(f"   {k_} =", v_)
# momentum constraint at its lowest order (t^1): must hold identically
for key, lowest, coeffs in eqs:
    if key.startswith('t') and key != 'tt':
        for pw, c in coeffs.items():
            r = sp.simplify(sp.expand(c.subs(sol).doit()))
            assert r == 0, (key, pw, r)
print("[ok] momentum constraint satisfied identically by that solution (all powers of t at this truncation)")
# compare with 3-curvature of the seed metric e^{2 zeta} delta
G3, g3i = christoffel(sp.exp(2 * zeta) * sp.eye(3), list(Xs))
R3 = ricci(G3, list(Xs), simp=sp.simplify)
R3s = sp.simplify(sum(g3i[i, j] * R3[i, j] for i in range(3) for j in range(3)))
cand_nu2 = sp.simplify(sol[nu2] / R3s)
print("   nu2 / R[e^{2zeta} delta] =", cand_nu2)
TF = sp.simplify(R3 - sp.exp(2 * zeta) * sp.eye(3) * R3s / 3)
ratio = sp.simplify(sol[h[0, 1]] / (sp.exp(-2 * zeta) * TF[0, 1]))
print("   H2_01 / (e^{-2zeta} traceless Ricci_01) =", ratio)
for (i, j) in names:
    assert sp.simplify(sol[h[i, j]] - ratio * sp.exp(-2 * zeta) * TF[i, j]) == 0
print("[ok] H2 is proportional to the traceless Ricci tensor of the seed metric, nu2 to its Ricci scalar")
