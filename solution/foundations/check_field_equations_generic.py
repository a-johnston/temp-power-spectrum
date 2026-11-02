"""Exact field equations E00/E0i/Eij and the Hamiltonian constraint (L4) for GENERIC gamma.

A strong test of Prop. 'Unified field equations' of the report.
gamma = w(X) L diag(a, b, 1/(ab)) L^T with L unit lower-triangular (P, Q, R): six free
functions; a, b, P, Q, R, N depend on (t,x,y,z), w on (x,y,z) only, so det gamma = w^3 is
t-independent but spatially varying -- exactly the class allowed by lemma L1
(det gamma = e^{6 zeta(X)}).  Every symmetric gamma with t-independent det can be written
this way locally (LDL^T decomposition).

All curvature quantities at a point depend only on the 2-jet of the metric functions at that
point.  Mode 'jet' (default) replaces every function value and derivative by an INDEPENDENT
SYMBOL, so a zero residual after exact cancellation is a proof of the identity for arbitrary
functions.  Mode 'random' substitutes random rational jets instead (fast, probabilistic,
for mutation testing only).
"""
import sys
import random
import sympy as sp
from sympy.core.function import AppliedUndef

MODE = sys.argv[1] if len(sys.argv) > 1 else 'jet'
t, x, y, z = X = sp.symbols('t x y z')
X3 = X[1:]
F = lambda nm, args: sp.Function(nm)(*args)
w = F('w', X3); a = F('a', X); b = F('b', X); N = F('N', X)
P_, Q_, R_ = F('P', X), F('Q', X), F('R', X)
L = sp.Matrix([[1, 0, 0], [P_, 1, 0], [Q_, R_, 1]])
Li = sp.Matrix([[1, 0, 0], [-P_, 1, 0], [P_ * R_ - Q_, -R_, 1]])
assert sp.simplify(L * Li - sp.eye(3)) == sp.zeros(3)
gam = w * L * sp.diag(a, b, 1 / (a * b)) * L.T
gi = (Li.T * sp.diag(1 / a, 1 / b, a * b) * Li) / w
assert sp.simplify(gam * gi - sp.eye(3)) == sp.zeros(3)
assert sp.simplify(gam.det() - w**3) == 0
g = sp.zeros(4, 4); g[0, 0] = -N**2
ginv = sp.zeros(4, 4); ginv[0, 0] = -1 / N**2
for i in range(3):
    for j in range(3):
        g[i + 1, j + 1] = N**2 * gam[i, j]
        ginv[i + 1, j + 1] = gi[i, j] / N**2


def make_ev(mode, seed=0):
    rnd = random.Random(seed)
    cache = {}

    def val(key, positive):
        if key not in cache:
            if mode == 'jet':
                nm = key[0] + ''.join(f'_{v}{c}' for v, c in key[1])
                cache[key] = sp.Symbol(nm, positive=True) if positive else sp.Symbol(nm)
            else:
                cache[key] = sp.Rational(rnd.randint(-30, 30), rnd.randint(7, 19)) + (2 if positive else 0)
        return cache[key]

    def ev(e):
        e = sp.sympify(e)
        rep = {}
        for d in e.atoms(sp.Derivative):
            vc = tuple(sorted((str(v), c) for v, c in d.variable_count))
            rep[d] = val((d.expr.func.__name__, vc), False)
        for f in e.atoms(AppliedUndef):
            rep[f] = val((f.func.__name__, ()), f.func.__name__ in ('w', 'a', 'b', 'N'))
        return sp.cancel(e.xreplace(rep))
    return ev


def chris(gm, gmi, xs):
    n = len(xs)
    dg = [[[sp.diff(gm[p, q], xs[r]) for r in range(n)] for q in range(n)] for p in range(n)]
    return [[[sum(gmi[A, d] * (dg[d][B][C] + dg[d][C][B] - dg[B][C][d]) for d in range(n)) / 2
              for C in range(n)] for B in range(n)] for A in range(n)]


G4 = chris(g, ginv, X)
G3 = chris(gam, gi, X3)


def ricci_at(G, xs, ev):
    n = len(xs)
    Gp = [[[ev(G[A][B][C]) for C in range(n)] for B in range(n)] for A in range(n)]
    Rm = sp.zeros(n, n)
    for B in range(n):
        for D in range(B, n):
            s = 0
            for A in range(n):
                s += ev(sp.diff(G[A][B][D], xs[A]) - sp.diff(G[A][B][A], xs[D]))
                for E in range(n):
                    s += Gp[A][A][E] * Gp[E][B][D] - Gp[A][D][E] * Gp[E][B][A]
            Rm[B, D] = Rm[D, B] = s
    return Rm


def residuals(ev):
    R4 = ricci_at(G4, X, ev)
    R3 = ricci_at(G3, X3, ev)
    n = sp.log(N)
    dn = [sp.diff(n, q) for q in X3]
    ndot = sp.diff(n, t)
    sig = sp.Matrix(3, 3, lambda i, j: sp.diff(gam[i, j], t) / 2)
    sig_up = gi * sig
    sigsq = sum(sig_up[i, j] * sig_up[j, i] for i in range(3) for j in range(3))
    DDn = sp.Matrix(3, 3, lambda i, j: sp.diff(n, X3[i], X3[j]) - sum(G3[k][i][j] * dn[k] for k in range(3)))
    lap = sum(gi[i, j] * DDn[i, j] for i in range(3) for j in range(3))
    gradsq = sum(gi[i, j] * dn[i] * dn[j] for i in range(3) for j in range(3))
    divsig = [sum(sp.diff(sig_up[j, i], X3[j]) for j in range(3))
              + sum(G3[j][j][k] * sig_up[k, i] for j in range(3) for k in range(3))
              - sum(G3[k][j][i] * sig_up[j, k] for j in range(3) for k in range(3)) for i in range(3)]
    res = {'E00': ev(R4[0, 0]) - ev(-sigsq - 3 * sp.diff(n, t, 2) + lap + 2 * gradsq)}
    for i in range(3):
        sj = sum(sig_up[j, i] * dn[j] for j in range(3))
        res[f'E0{i + 1}'] = ev(R4[0, i + 1]) - ev(divsig[i] - 2 * sp.diff(ndot, X3[i]) + 2 * sj + 2 * ndot * dn[i])
    sigsig = sig * gi * sig
    for i in range(3):
        for j in range(i, 3):
            Eij = (sp.diff(sig[i, j], t) - 2 * sigsig[i, j] - 2 * (DDn[i, j] - sig[i, j] * ndot)
                   + gam[i, j] * (sp.diff(n, t, 2) - lap) + 2 * dn[i] * dn[j] + 2 * gam[i, j] * (ndot**2 - gradsq))
            res[f'E{i + 1}{j + 1}'] = ev(R4[i + 1, j + 1]) - R3[i, j] - ev(Eij)
    gip = ginv.applyfunc(ev)
    gp = g.applyfunc(ev)
    Rs = sum(gip[m, q] * R4[m, q] for m in range(4) for q in range(4))
    Gtt = sum(gip[0, m] * gip[0, q] * (R4[m, q] - gp[m, q] * Rs / 2) for m in range(4) for q in range(4))
    R3s = sum(ev(gi[i, j]) * R3[i, j] for i in range(3) for j in range(3))
    res['Ham'] = 2 * ev(N)**6 * Gtt - ev(N**2) * (R3s - ev(4 * lap + 2 * gradsq + sigsq)) - ev(6 * sp.diff(N, t)**2)
    return res


if __name__ == '__main__':
    seeds = (0,) if MODE == 'jet' else (1, 2)
    for seed in seeds:
        res = residuals(make_ev(MODE, seed))
        for k_, v_ in res.items():
            r = sp.cancel(sp.together(v_))
            assert r == 0, (k_, MODE, seed)
            print(f"[ok] {k_} residual = 0  (mode={MODE}{'' if MODE == 'jet' else f', seed={seed}'})", flush=True)
    if MODE == 'jet':
        print("GENERIC FIELD EQUATIONS PROVED (independent symbolic 2-jets, t-independent spatially varying det gamma)")
