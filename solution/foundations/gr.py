"""Minimal exact (symbolic) tensor utilities used by every verification script.

Everything here is exact SymPy algebra: no numerical evaluation is used anywhere.
Index convention: coordinates x[0..3] = (t, x, y, z); metric g[mu, nu].
"""
import sympy as sp


def christoffel(g, x):
    """Gamma^a_{bc} for metric g in coordinates x."""
    n = len(x)
    ginv = g.inv()
    dg = [[[sp.diff(g[a, b], x[c]) for c in range(n)] for b in range(n)] for a in range(n)]
    Gam = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(b, n):
                s = 0
                for d in range(n):
                    if ginv[a, d] != 0:
                        s += ginv[a, d] * (dg[d][b][c] + dg[d][c][b] - dg[b][c][d])
                Gam[a][b][c] = Gam[a][c][b] = s / 2
    return Gam, ginv


def ricci(Gam, x, simp=lambda e: e):
    """R_{bd} = d_a Gam^a_{bd} - d_d Gam^a_{ba} + Gam^a_{ae} Gam^e_{bd} - Gam^a_{de} Gam^e_{ba}."""
    n = len(x)
    R = sp.zeros(n, n)
    for b in range(n):
        for d in range(b, n):
            s = 0
            for a in range(n):
                s += sp.diff(Gam[a][b][d], x[a]) - sp.diff(Gam[a][b][a], x[d])
                for e in range(n):
                    s += Gam[a][a][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][a]
            R[b, d] = R[d, b] = simp(s)
    return R


def series_eps(expr, eps, order):
    """Exact Taylor coefficients of expr in the bookkeeping parameter eps up to eps**order."""
    out = []
    e = expr
    fact = 1
    for n in range(order + 1):
        if n > 0:
            e = sp.diff(e, eps)
            fact *= n
        out.append(e.subs(eps, 0) / fact)
    return out


def assert_zero(expr, label, simp=sp.simplify):
    r = simp(expr)
    if r != 0:
        raise AssertionError(f"[FAIL] {label}: residual = {r}")
    print(f"[ok] {label}")
