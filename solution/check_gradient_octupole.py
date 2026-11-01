"""Exact checks of the leading-gradient octupole and the O(kappa^6) cancellation of M_l.

1. Conformal identity for the seed metric gamma = e^{2 zeta} delta, f = e^{-zeta}: in the orthonormal
   frame e^{-zeta} d_a,   (nabla_a R_bc)_STF = (f/2) [d_a d_b d_c (f^2)]_STF   (generic zeta).
2. Leading transfer coefficient for every l (small k):
   T_l = -k^l (t_o - t*)^{l-1} [t_o + (l-1) t*] / (3 (2l+1)!!) + O(k^{l+2}),   l = 2..7.
3. beta_3 = (1/3) int_{t*}^{t_o} t (t_o - t) dt = (t_o - t*)^2 (t_o + 2 t*)/18.
4. Gaussian average of the leading octupole, Q = e^{-3 eps g} STF[-eps g_abc + 6 eps^2 Sym(g_ab g_c)
   - 4 eps^3 g_a g_b g_c] (from (f/2) d^3 f^2), with exact isotropic Wick moments:
   <Q:Q> = (2/5) eps^2 e^{18 eps^2 m0} [m6 + 28 eps^2 m2 m4 + (560/9) eps^4 m2^3], and on the shell
   C^(22)/C^(11) = 148 sigma^2, C^(13)/C^(11) = 27 sigma^2 (fourth-order factor 1 + 46 eps^2 sigma^2).
5. M_l = O(kappa^8): with theta = -tau^2 A_nn and L2 = tau^2 (4 g A + 2 B)_mm at leading gradient,
   <theta(n) theta(n') L2(m)> = 0 for arbitrary directions (exact Wick).
"""
import itertools
import sympy as sp

X = sp.symbols('x1:4', real=True)
eps = sp.Symbol('epsilon')


# ------------------------------------------------------------------ 1. conformal identity
def check_conformal_identity():
    Z = sp.Function('Z')(*X)
    g = sp.exp(2 * Z) * sp.eye(3)
    gi = sp.exp(-2 * Z) * sp.eye(3)
    Gam = [[[sum(gi[i, l] * (sp.diff(g[l, j], X[k]) + sp.diff(g[l, k], X[j]) - sp.diff(g[j, k], X[l]))
                 for l in range(3)) / 2 for k in range(3)] for j in range(3)] for i in range(3)]
    Ric = sp.zeros(3)
    for i, j in itertools.product(range(3), repeat=2):
        Ric[i, j] = sum(sp.diff(Gam[k][i][j], X[k]) - sp.diff(Gam[k][i][k], X[j])
                        + sum(Gam[k][k][l] * Gam[l][i][j] - Gam[k][j][l] * Gam[l][i][k] for l in range(3))
                        for k in range(3))
    nab = {}
    for a, b, c in itertools.product(range(3), repeat=3):
        nab[a, b, c] = sp.diff(Ric[b, c], X[a]) - sum(Gam[d][a][b] * Ric[d, c] + Gam[d][a][c] * Ric[b, d]
                                                       for d in range(3))
    ortho = {key: sp.exp(-3 * Z) * v for key, v in nab.items()}          # three lowered frame indices
    f = sp.exp(-Z)
    rhs = {(a, b, c): f / 2 * sp.diff(f**2, X[a], X[b], X[c]) for a, b, c in itertools.product(range(3), repeat=3)}
    for a, b, c in [(0, 0, 0), (0, 0, 1), (0, 1, 2), (1, 1, 2)]:
        dlhs = STF(ortho, a, b, c) - STF(rhs, a, b, c)
        assert sp.simplify(sp.expand(dlhs)) == 0, (a, b, c)
    print('[ok] (nabla_<a R_bc>) = (f/2) [d_a d_b d_c f^2]_STF for the seed metric (generic zeta)')


def sym3(T, a, b, c):
    return sum(T[p] for p in set(itertools.permutations((a, b, c)))) / len(set(itertools.permutations((a, b, c)))) \
        if False else sum(T[p] for p in itertools.permutations((a, b, c))) / 6


def STF(T, a, b, c):
    S = lambda i, j, k: sym3(T, i, j, k)
    tr = lambda i: sum(S(i, d, d) for d in range(3))
    d = lambda i, j: 1 if i == j else 0
    return S(a, b, c) - sp.Rational(1, 5) * (d(a, b) * tr(c) + d(a, c) * tr(b) + d(b, c) * tr(a))


# ------------------------------------------------------------------ 2-3. transfer coefficient, beta_3
def check_transfer_and_beta():
    k, t, to, ts, x = sp.symbols('k t t_o t_star x', positive=True)
    for l in range(2, 8):
        y = sp.Symbol('y', positive=True)
        jl = sp.expand_func(sp.jn(l, y))
        w = sp.series(jl / 3 + sp.diff(jl, y, 2), y, 0, l + 1).removeO()
        sprime = sp.series(sp.diff(2 * sp.sin(x) / x - sp.cos(x), x), x, 0, 4).removeO()
        integrand = -3 * sprime.subs(x, k * t / sp.sqrt(3)) * w.subs(y, k * (to - t)) * k / sp.sqrt(3)
        T = sp.integrate(sp.expand(integrand), (t, ts, to))
        lead = sp.expand(T).coeff(k, l)
        target = -(to - ts)**(l - 1) * (to + (l - 1) * ts) / (3 * sp.factorial2(2 * l + 1))
        assert sp.simplify(lead - target) == 0, l
    beta3 = sp.integrate(t * (to - t), (t, ts, to)) / 3
    assert sp.simplify(beta3 - (to - ts)**2 * (to + 2 * ts) / 18) == 0
    print('[ok] leading T_l coefficient for l = 2..7; beta_3 time integral')


# ------------------------------------------------------------------ 4-5. exact isotropic Wick moments
m = {j: sp.Symbol(f'm{j}', positive=True) for j in (0, 2, 4, 6, 8)}


def sphere_avg(gamma):
    if any(e % 2 for e in gamma):
        return sp.S.Zero
    a, b, c = gamma
    num = sp.gamma(sp.Rational(a + 1, 2)) * sp.gamma(sp.Rational(b + 1, 2)) * sp.gamma(sp.Rational(c + 1, 2))
    return 2 * num / sp.gamma(sp.Rational(a + b + c + 3, 2)) / (4 * sp.pi)


def cov(al, be):
    """<d^al g(0) d^be g(0)> for an isotropic Gaussian field, al/be multi-indices (counts)."""
    n = sum(al) + sum(be)
    if n % 2:
        return sp.S.Zero
    gamma = tuple(x + y for x, y in zip(al, be))
    return sp.I**sum(al) * (-sp.I)**sum(be) * m[n] * sphere_avg(gamma)


DER = {}


def d(*idx):
    """symbol for the derivative of g with the given index list."""
    al = tuple(idx.count(i) for i in range(3))
    if al not in DER:
        DER[al] = sp.Symbol('g_' + ''.join(map(str, al)))
    return DER[al]


def wick(expr):
    expr = sp.expand(expr)
    inv = {v: k for k, v in DER.items()}
    syms = [s for s in expr.free_symbols if s in inv]
    if not syms:
        return expr
    P = sp.Poly(expr, *syms)
    out = 0
    for mon, coef in P.terms():
        items = [s for s, k in zip(syms, mon) for _ in range(k)]
        out += coef * pairings(items, inv)
    return sp.expand(out)


def pairings(items, inv):
    if not items:
        return sp.S.One
    if len(items) % 2:
        return sp.S.Zero
    first, rest = items[0], items[1:]
    tot = 0
    for i, other in enumerate(rest):
        c = cov(inv[first], inv[other])
        if c != 0:
            tot += c * pairings(rest[:i] + rest[i + 1:], inv)
    return tot


def check_octupole_and_M():
    g0 = d()
    gi = lambda i: d(i)
    gij = lambda i, j: d(i, j)
    gijk = lambda i, j, k: d(i, j, k)
    # independence claims at a point
    for i, j in itertools.product(range(3), repeat=2):
        A = gij(i, j) - (1 if i == j else 0) * sum(gij(a, a) for a in range(3)) / 3
        assert wick(g0 * A) == 0 and wick(gi(0) * A) == 0
    # rank-3 octupole tensor (eps-polynomial part, without the e^{-3 eps g} factor)
    T = {}
    for a, b, c in itertools.product(range(3), repeat=3):
        T[a, b, c] = (-eps * gijk(a, b, c) + 2 * eps**2 * (gij(a, b) * gi(c) + gij(a, c) * gi(b) + gij(b, c) * gi(a))
                      - 4 * eps**3 * gi(a) * gi(b) * gi(c))
    Q = {key: STF(T, *key) for key in T}
    QQ = sp.expand(sum(Q[key]**2 for key in Q))
    poly = sp.Poly(QQ, eps)
    coeffs = {k[0]: wick(c) for k, c in poly.terms()}
    assert sp.simplify(coeffs[2] - sp.Rational(2, 5) * m[6]) == 0
    assert sp.simplify(coeffs[4] - sp.Rational(2, 5) * 28 * m[2] * m[4]) == 0
    assert sp.simplify(coeffs[6] - sp.Rational(2, 5) * sp.Rational(560, 9) * m[2]**3) == 0
    assert coeffs.get(3, 0) == 0 and coeffs.get(5, 0) == 0
    # g is independent of every STF piece (pure-trace and odd covariances only)
    for key in [(0, 0, 0), (0, 1, 2)]:
        for e_ in sp.Poly(sp.expand(Q[key]), eps).coeffs():
            assert wick(g0 * e_) == 0
    print('[ok] <Q:Q> = (2/5) eps^2 e^{18 eps^2 m0} [m6 + 28 eps^2 m2 m4 + (560/9) eps^4 m2^3]')
    # fourth-order split on the shell: Theta^(r) Taylor coefficients of e^{-3 eps g} * poly
    full = {key: sp.expand(sp.series(sp.exp(-3 * eps * g0) * Q[key], eps, 0, 4).removeO()) for key in Q}
    Qr = {r: {key: sp.factorial(r) * full[key].coeff(eps, r) for key in Q} for r in (1, 2, 3)}
    shell = {m[j]: sp.Symbol('sigma2') for j in m}
    s2 = sp.Symbol('sigma2')
    c11 = wick(sum(Qr[1][k]**2 for k in Q)).subs(shell)
    c22 = wick(sum(Qr[2][k]**2 for k in Q)).subs(shell)
    c13 = wick(sum(Qr[1][k] * Qr[3][k] for k in Q)).subs(shell)
    assert sp.simplify(c22 / c11 - 148 * s2) == 0 and sp.simplify(c13 / c11 - 27 * s2) == 0
    assert sp.simplify(c22 / c11 / 4 + c13 / c11 / 3 - 46 * s2) == 0
    print('[ok] shell octupole: C22/C11 = 148 sigma^2, C13/C11 = 27 sigma^2, factor 1 + 46 eps^2 sigma^2')
    # M_l at O(kappa^6): <theta(n) theta(n') L2(m)> = 0
    n1, n2, n3 = [sp.symbols(f'{c}1:4', real=True) for c in 'uvw']
    A = lambda u: sum((gij(i, j) - (1 if i == j else 0) * sum(gij(a, a) for a in range(3)) / 3) * u[i] * u[j]
                      for i in range(3) for j in range(3))
    B = lambda u: sum((gi(i) * gi(j) - (1 if i == j else 0) * sum(gi(a)**2 for a in range(3)) / 3) * u[i] * u[j]
                      for i in range(3) for j in range(3))
    three = wick(A(n1) * A(n2) * (4 * g0 * A(n3) + 2 * B(n3)))
    assert sp.expand(three) == 0
    print('[ok] <theta theta L2> = 0 at leading gradient order for arbitrary directions => M_l = O(kappa^8)')


if __name__ == '__main__':
    check_transfer_and_beta()
    check_octupole_and_M()
    check_conformal_identity()
    print('GRADIENT OCTUPOLE CHECKS PASSED (exact symbolic algebra)')
