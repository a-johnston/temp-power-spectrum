"""Independent tests of formulas as printed in the report, sharing no code with the other solution modules.

Coefficients are retyped from the TeX. Field products use a small multilinear algebra: a function is a map from
sets of distinct plane-wave labels to coefficients, which is exactly the report's multilinear rule.

Exact parts: seed-curvature identities; the degree-six quadrupole from the printed three-tensor map; printed pair and
contracted cubic sources against exact-in-zeta t^2, t^4 data; acoustic amplitude; Green polynomial; V_2, V_3.
Numerical validation only (quadrature, not used in any derivation): Kohri-Terada amplitude, entire-function samples,
sphere moments, T_l table values and node positions.
"""

import itertools

import mpmath as mp
import sympy as sp

R = sp.Rational
rng = range(3)


# ---------------------------------------------------------------- multilinear algebra
class ML:
    __slots__ = ("c", "moms")

    def __init__(self, c, moms):
        self.c = {k: v for k, v in c.items() if v != 0}
        self.moms = moms

    def __add__(self, other):
        if not isinstance(other, ML):
            other = ML({0: other}, self.moms)
        d = dict(self.c)
        for k, v in other.c.items():
            d[k] = d.get(k, 0) + v
        return ML(d, self.moms)

    __radd__ = __add__

    def __neg__(self):
        return ML({k: -v for k, v in self.c.items()}, self.moms)

    def __sub__(self, other):
        return self + (-other if isinstance(other, ML) else -other)

    def __mul__(self, other):
        if not isinstance(other, ML):
            return ML({k: v * other for k, v in self.c.items()}, self.moms)
        d = {}
        for k, v in self.c.items():
            for l, w in other.c.items():
                if not k & l:
                    d[k | l] = d.get(k | l, 0) + v * w
        return ML({k: sp.expand(v) for k, v in d.items()}, self.moms)

    __rmul__ = __mul__

    def d(self, i):
        out = {}
        for k, v in self.c.items():
            total = sum(self.moms[j][i] for j in range(len(self.moms)) if k >> j & 1)
            out[k] = sp.expand(sp.I * total * v)
        return ML(out, self.moms)

    def coef(self, mask):
        return self.c.get(mask, sp.Integer(0))


def tf(m):
    tr = sum((m[i][i] for i in rng), 0 * m[0][0])
    return [[m[i][j] - (tr * R(1, 3) if i == j else 0) for j in rng] for i in rng]


def curvature_tensors(moms):
    """Seed tensors a, A, b, B of the printed gradient expansion, in the seed-orthonormal frame."""
    z = ML({1 << j: sp.Integer(1) for j in range(len(moms))}, moms)
    em2 = ML({0: sp.Integer(1)}, moms) + z * (-2) + z * z * 2 + z * z * z * R(-4, 3)  # at most three labels
    dz = [z.d(i) for i in rng]
    a = em2 * (sum((dz[i].d(i) for i in rng), 0 * z) * 2 + sum((dz[i] * dz[i] for i in rng), 0 * z)) * R(1, 18)
    raw = tf([[dz[i].d(j) - dz[i] * dz[j] for j in rng] for i in rng])
    A = [[em2 * raw[i][j] * R(1, 3) for j in rng] for i in rng]
    da = [a.d(i) for i in rng]
    dzda = sum((dz[i] * da[i] for i in rng), 0 * z)
    hess = [[em2 * (da[i].d(j) - dz[i] * da[j] - dz[j] * da[i] + (dzda if i == j else 0 * z)) for j in rng] for i in rng]
    A2 = [[sum((A[i][k] * A[k][j] for k in rng), 0 * z) for j in rng] for i in rng]
    b = (sum((hess[i][i] for i in rng), 0 * z) - a * a * 6 - sum((A2[i][i] for i in rng), 0 * z)) * R(1, 30)
    H, S = tf(hess), tf(A2)
    B = [[H[i][j] * R(3, 20) + S[i][j] * R(3, 20) - a * A[i][j] * R(7, 10) for j in rng] for i in rng]
    return {"a": a, "A": A, "H": H, "S": S, "b": b, "B": B}


def comps(m, mask):
    return [[m[i][j].coef(mask) for j in rng] for i in rng]


mu = sp.symbols("mu", real=True)
sn = sp.symbols("sn", positive=True)  # sqrt(1-mu^2)


def even_sn(x):
    """Replace sn^2 by 1-mu^2; odd powers remain and vanish only after an azimuthal or symmetry argument."""
    out = 0
    for (d,), c in sp.Poly(sp.expand(x), sn).terms():
        out += c * (1 - mu**2) ** (d // 2) * sn ** (d % 2)
    return sp.expand(out)


def avg(f):
    return sp.integrate(sp.expand(even_sn(f)), (mu, -1, 1)) / 2


# ---------------------------------------------------------------- curvature identities on a generic seed
def kd(i, j):
    return sp.Integer(1 if i == j else 0)


def check_seed_identities():
    x, y, z = X = sp.symbols("x y z")
    zeta = (
        R(1, 3) * x - R(2, 5) * y + R(1, 7) * z + R(1, 2) * x * y - R(3, 4) * z * z + R(2, 3) * x * x + R(1, 5) * x * y * z
        + R(1, 6) * y**3 - R(2, 7) * x * z * z + R(1, 9) * x**2 * y * z + R(1, 8) * z**4 - R(1, 4) * x * y**3
        + R(1, 11) * x**5 + R(1, 13) * y * z**4
    )  # fmt: skip
    g, gi = sp.exp(2 * zeta), sp.exp(-2 * zeta)
    pt = dict.fromkeys(X, 0)
    d = lambda f, i: sp.diff(f, X[i])
    gam = lambda c, a, b: kd(a, c) * d(zeta, b) + kd(b, c) * d(zeta, a) - kd(a, b) * d(zeta, c)

    def cov(T, rank):
        out = {}
        for idx in itertools.product(rng, repeat=rank):
            for c in rng:
                v = d(T[idx], c)
                for pos in range(rank):
                    for e in rng:
                        j = list(idx)
                        j[pos] = e
                        v -= gam(e, idx[pos], c) * T[tuple(j)]
                out[idx + (c,)] = v
        return out

    lap = sum(d(d(zeta, i), i) for i in rng)
    gr2 = sum(d(zeta, i) ** 2 for i in rng)
    a = gi * (2 * lap + gr2) / 18
    ric = {(i, j): -d(d(zeta, i), j) - kd(i, j) * lap + d(zeta, i) * d(zeta, j) - kd(i, j) * gr2 for i in rng for j in rng}
    rs = gi * sum(ric[(i, i)] for i in rng)
    assert sp.simplify((a + rs / 36).subs(pt)) == 0
    A = {k: -(v - kd(k[0], k[1]) * g * rs / 3) / 3 for k, v in ric.items()}
    dA, da = cov(A, 2), {(i,): d(a, i) for i in rng}
    ddA, ha = cov(dA, 3), cov(da, 1)

    def tfc(m):
        tr = gi * sum(m[(i, i)] for i in rng)
        return {(i, j): m[(i, j)] - kd(i, j) * g * tr / 3 for i in rng for j in rng}

    A2 = {(i, j): gi * sum(A[(i, k)] * A[(k, j)] for k in rng) for i in rng for j in rng}
    H, S = tfc(ha), tfc(A2)
    zero = lambda e: sp.simplify(e.subs(pt)) == 0
    assert all(zero(gi * sum(dA[(c, b, c)] for c in rng) - 2 * da[(b,)]) for b in rng), "Bianchi"
    for p, q, r in itertools.product(rng, repeat=3):
        assert zero(dA[(p, q, r)] + kd(p, q) * g * da[(r,)] - dA[(p, r, q)] - kd(p, r) * g * da[(q,)]), "Cotton"
    lapA = tfc({(p, q): gi * sum(ddA[(p, q, c, c)] for c in rng) for p in rng for q in rng})
    swapped = tfc({(p, q): gi * sum(ddA[(c, q, p, c)] for c in rng) for p in rng for q in rng})
    sym4 = lambda i, j, k, l: sum(ddA[p] for p in itertools.permutations((i, j, k, l))) / 24
    tr4 = tfc({(p, q): gi * sum(sym4(p, q, c, c) for c in rng) for p in rng for q in rng})
    for p, q in [(0, 0), (0, 1), (1, 2), (1, 0), (2, 2)]:
        common = -9 * S[(p, q)] - 18 * a * A[(p, q)]
        assert zero(lapA[(p, q)] - 3 * H[(p, q)] - common)
        assert zero(swapped[(p, q)] - 2 * H[(p, q)] - common)
        assert zero(tr4[(p, q)] - (R(11, 6) * H[(p, q)] - R(9, 2) * S[(p, q)] - 9 * a * A[(p, q)]))
    print("[ok] Bianchi, Cotton, Laplacian and symmetrised-trace identities on a generic conformally flat seed")


# ---------------------------------------------------------------- degree-six quadrupole
def check_degree_six():
    o, e, io, p, r = sp.symbols("o e iota p r", positive=True)
    t = sp.symbols("t")
    BH, BS, BA = sp.symbols("BH BS BA", real=True)
    D = o**2 - e**2
    L = sp.integrate(t * (o - t) ** 2, (t, e, o))
    E = o**4 - e**4
    bH = -R(11, 14) * L - R(3, 40) * E
    bS = R(27, 14) * L - R(3, 40) * E - D**2 / 28
    ba = R(27, 7) * L + R(7, 20) * E + (o**4 + 3 * io**4) / 4 - e**4
    assert sp.expand(bH + (o - e) * (59 * o**3 + 59 * o**2 * e - 106 * o * e**2 + 114 * e**3) / 420) == 0
    assert sp.expand(bS - (o - e) * (7 * o**3 + 7 * o**2 * e - 118 * o * e**2 + 62 * e**3) / 140) == 0
    assert sp.expand(ba - R(3, 140) * (43 * o**4 - 90 * o**2 * e**2 + 120 * o * e**3 - 108 * e**4 + 35 * io**4)) == 0
    q = (r * mu, r * sn, 0)
    c = curvature_tensors([(p, 0, 0), q, tuple(-w for w in q)])
    Q2 = [[c["A"][i][j] * (-D / 2) for j in rng] for i in rng]
    Q4 = [[c["H"][i][j] * BH + c["S"][i][j] * BS + c["a"] * c["A"][i][j] * BA for j in rng] for i in rng]
    pref = 8 * sp.pi / 75
    dot = lambda m, n: sp.expand(sum(m[i][j] * sp.conjugate(n[i][j]) for i in rng for j in rng))
    a2, a4 = comps(Q2, 1), comps(Q4, 1)
    p2, p4 = comps(Q2, 3), comps(Q4, 3)
    t2, t4 = comps(Q2, 7), comps(Q4, 7)
    q2, q4 = comps(Q2, 2), comps(Q4, 2)
    C11_4, C11_6 = pref * dot(a2, a2), 2 * pref * dot(a2, a4)
    C22_4, C13_4 = pref * avg(dot(p2, p2)) / 2, pref * avg(dot(a2, t2))
    C22_6, C13_6 = pref * avg(dot(p2, p4)), pref * avg(dot(a2, t4) + dot(a4, t2))

    assert sp.simplify(C22_4.subs(r, p) / C11_4 - R(17, 3)) == 0 and sp.simplify(C13_4.subs(r, p) / C11_4 - 4) == 0
    m0, m2, m4, m6 = sp.symbols("m0 m2 m4 m6")
    moments = {0: m0, 2: m2, 4: m4, 6: m6}
    to_m = lambda x: sp.expand(sum(co * moments[i] * moments[j] for (i, j), co in sp.Poly(sp.expand(x), p, r).terms()))
    pr22 = 32 * sp.pi / 18225 * D * (12 * BH * m0 * m6 + (18 * BH - BA) * m2 * m4)
    pr13 = 8 * sp.pi / 18225 * D * (60 * BH * m0 * m6 + (50 * BH - 9 * BA) * m2 * m4)
    assert sp.simplify(to_m(C22_6) - pr22) == 0 and sp.simplify(to_m(C13_6) - pr13) == 0
    assert sp.simplify(C11_6 - 16 * sp.pi / 6075 * D * BH * p**6) == 0
    print("[ok] printed general-moment coefficients; the (A^2)_TF tensor drops out")

    shell = lambda x: sp.expand(x.subs(r, p).subs(p, 1).subs({BH: bH, BS: bS, BA: ba}))
    pp11 = -4 * sp.pi / 637875 * (o - e) ** 2 * (o + e) * (59 * o**3 + 59 * o**2 * e - 106 * o * e**2 + 114 * e**3)
    pp22 = -8 * sp.pi / 637875 * D * (105 * io**4 + 719 * o**4 - 1920 * o**2 * e**2 + 2560 * o * e**3 - 1464 * e**4)
    pp13 = -2 * sp.pi / 1913625 * D * (2835 * io**4 + 9973 * o**4 - 25440 * o**2 * e**2 + 33920 * o * e**3 - 21288 * e**4)
    tot = -2 * sp.pi / 1913625 * D * (4095 * io**4 + 18601 * o**4 - 48480 * o**2 * e**2 + 64640 * o * e**3 - 38856 * e**4)
    assert sp.expand(shell(C11_6) - pp11) == 0 and sp.expand(shell(C22_6) - pp22) == 0 and sp.expand(shell(C13_6) - pp13) == 0
    assert sp.expand(pp22 + pp13 - tot) == 0
    print("[ok] printed explicit shell polynomials through k0^6")

    K3c = (4 * (e - o) * (12 * e**3 - 8 * e**2 * o + 7 * e * o**2 + 7 * o**3) + 5 * (o**4 - io**4)) / 270
    cross = 8 * sp.pi * D / 91125 * (204 * e**4 - 280 * e**3 * o + 210 * e**2 * o**2 - 89 * o**4 - 45 * io**4)
    assert sp.expand(18 * pp11 + R(28, 3) * pp11 + 4 * sp.pi / 25 * D / 9 * K3c + cross - tot) == 0
    J2, J4 = 4 * D / 9, (3984 * e**4 - 6256 * e**3 * o + 4692 * e**2 * o**2 - 1853 * o**4 - 567 * io**4) / 3402
    K12 = D / 9
    K14 = pp11 / (4 * sp.pi / 25 * 2 * K12)
    assert sp.simplify(4 * sp.pi / 25 * (K12 * J4 + K14 * J2) - pp13) == 0
    K2c4 = cross / (-12 * 4 * sp.pi / 25 * K12)
    assert sp.simplify(J4 + 8 * K2c4 - 16 * K14 - K3c) == 0
    print("[ok] four-piece dilation sum; stationary, centred and cross formulas are mutually consistent")

    kc2 = [[p2[i][j] + 2 * a2[i][j] + 2 * q2[i][j] for j in rng] for i in rng]
    kc4 = [[p4[i][j] + 4 * a4[i][j] + 4 * q4[i][j] for j in rng] for i in rng]
    assert sp.simplify(shell(pref * avg(dot(kc2, kc4))) / pp11 - R(28, 3)) == 0
    assert sp.simplify(pref * avg(dot(kc2, kc2)).subs(r, p) / 2 / C11_4 - R(5, 3)) == 0
    print("[ok] centred 22 ratios 5/3 and 28/3 from the map itself")


# ---------------------------------------------------------------- printed sources by power series
t = sp.symbols("t", positive=True)
NMAX = 14


def trunc(x, n=NMAX):
    x = sp.expand(x)
    return sum(x.coeff(t, i) * t**i for i in range(-4, n + 1))


def solve_lapse(K2, source, c0, nmax=NMAX):
    """3n'' + K2 n - 6n/t^2 = -source, regular with n = c0 t^2 + ..."""
    source = sp.expand(source)
    n = {2: c0}
    for m in range(4, nmax + 1, 2):
        n[m] = sp.expand((-source.coeff(t, m - 2) - K2 * n[m - 2]) / (3 * (m - 2) * (m + 1)))
    return sum(v * t**m for m, v in n.items())


def solve_wave(lam, rhs):
    """y'' + 2y'/t + lam y = rhs, regular with y(0) = 0."""
    rhs = sp.expand(rhs)
    y = {0: sp.Integer(0)}
    for m in range(2, NMAX + 1, 2):
        y[m] = sp.expand((rhs.coeff(t, m - 2) - lam * y[m - 2]) / (m * (m + 1)))
    return sum(v * t**m for m, v in y.items())


def check_sources():
    x = t / sp.sqrt(3)
    a = sp.expand(sp.series(sp.cos(x) - sp.sin(x) / x, t, 0, NMAX + 3).removeO())
    b = sp.expand(sp.series(-6 * (2 * sp.sin(x) / x - sp.cos(x) - 1), t, 0, NMAX + 3).removeO())
    da, db = sp.diff(a, t), sp.diff(b, t)
    assert trunc(db - 6 * (da - a / t), NMAX - 2) == 0 and trunc(3 * sp.diff(a, t, 2) + a - 6 * a / t**2, NMAX - 4) == 0
    K2 = 2 * (1 + mu)
    CN = trunc(
        (R(2, 3) * a * b - 4 * a + mu * (4 * a**2 - R(4, 3) * a * b + 2 * a) - 2 * mu**2 * a * b)
        + db**2 * (mu**2 - R(1, 3)) / 2 + 12 * a**2 / t**2
    )  # fmt: skip
    JS = trunc(
        2 * (3 + mu) * a**2 + (15 * mu**2 + 13 * mu - 8) * a * b / 3 + 16 * (1 + mu) * a
        + (5 - 5 * mu - 15 * mu**2 - 9 * mu**3) * b**2 / 36 + (-10 + 7 * mu + 15 * mu**2) * b / 6 + 5 + 7 * mu
        - (1 + 3 * mu) * da * db
    )  # fmt: skip
    JT = trunc(
        (1 - mu) / 36 * (b**2 * (9 * mu**2 + 12 * mu + 1) - 6 * b * (6 * mu * a + 3 * mu + 2) - 72 * a**2 + 36 - 36 * da * db)
    )
    cr = (4 - mu) / 9
    n2 = solve_lapse(K2, CN, cr)
    FS, FT = solve_wave(-K2 / 3, trunc(-4 * K2 * n2 + JS)), solve_wave(K2, JT)
    nh = solve_lapse(0, CN.subs(mu, -1), R(5, 9))

    mfrak = db * ((1 - 3 * mu) * (R(3, 2) + a) + (3 * mu**2 - 1) * b / 4) - 6 * a * da
    assert sp.simplify(trunc(sp.diff(FS, t) - 6 * (sp.diff(n2, t) - n2 / t) - mfrak, NMAX - 3)) == 0

    def Tj(j, A, d):
        return sum((-A) ** m * d ** (2 * m + j) / sp.factorial(2 * m + j) for m in range(9))

    aK, rr = K2 / 3, sp.symbols("r", positive=True)
    dd = t - rr
    G = 2 * Tj(1, aK, dd) + 2 * ((dd - rr) * Tj(2, aK, dd) - 2 * Tj(3, aK, dd)) / (t * rr)
    poly = sp.Poly(sp.expand(G * CN.subs(t, rr) * rr * t), rr)
    Gint = 0
    for (d,), co in poly.terms():
        assert d > 0 or sp.simplify(co) == 0
        if d > 0:
            Gint += co * t**d / d
    Gint = sp.expand(Gint / t)
    Iint = sum(sp.expand(a**2 / t).coeff(t, i) * t ** (i + 1) / (i + 1) for i in range(NMAX + 1))
    cin2 = sum((-1) ** (m + 1) * (2 * x) ** (2 * m) / (2 * m * sp.factorial(2 * m)) for m in range(1, 12))
    Iprint = sp.expand(
        sp.series(sp.sin(2 * x) / (2 * x) - sp.sin(x) ** 2 / (2 * x**2) - R(1, 2), t, 0, NMAX + 2).removeO() + cin2 / 2
    )
    assert trunc(Iint - Iprint, NMAX - 2) == 0
    acoustic = (
        18 * cr * (Tj(2, aK, t) - 2 * Tj(3, aK, t) / t) - Gint
        + (1 - 3 * mu) * (R(3, 2) * b + 3 * a**2 - 6 * Iint) + (3 * mu**2 - 1) * b**2 / 8 - 3 * a**2
    )  # fmt: skip
    assert sp.simplify(trunc(acoustic - FS, NMAX - 3)) == 0
    print("[ok] printed sources: momentum identity, acoustic amplitude and its integral I")

    def even(e):
        out = 0
        for (d,), co in sp.Poly(sp.expand(e), sn).terms():
            assert d % 2 == 0 or sp.simplify(co) == 0
            out += co * (1 - mu**2) ** (d // 2)
        return sp.expand(out)

    kv, qv = (1, 0, 0), (mu, sn, 0)
    c = curvature_tensors([kv, qv])
    assert even(c["a"].coef(3)) == sp.expand(n2.coeff(t, 2)) and sp.simplify(even(c["b"].coef(3)) - n2.coeff(t, 4)) == 0
    Kv, Vv = [kv[i] + qv[i] for i in rng], [kv[i] - qv[i] for i in rng]
    contr = lambda m, u: sum(m[i][j] * u[i] * u[j] for i in rng for j in rng)
    for key, power in (("A", 2), ("B", 4)):
        m = comps(c[key], 3)
        assert sp.simplify(even(R(3, 2) * contr(m, Kv)) - K2 * FS.coeff(t, power)) == 0
        assert sp.simplify(even(contr(m, Vv) - 2 * (1 - mu) * m[2][2]) - 4 * (1 - mu) * FT.coeff(t, power)) == 0
    assert sp.simplify(FT.coeff(t, 2) - (1 - mu) / 6) == 0
    c0 = curvature_tensors([kv, tuple(-w for w in kv)])
    assert c0["a"].coef(3) == nh.coeff(t, 2) and sp.simplify(c0["b"].coef(3) - nh.coeff(t, 4)) == 0
    print("[ok] printed pair sources reproduce exact-in-zeta t^2, t^4 data for lapse, F_S, F_T, homogeneous lapse")

    # contracted cubic sources
    mx = 8
    tr = lambda e: trunc(e, mx)
    mean = lambda e: sum(sp.integrate(sp.expand(tr(e)).coeff(t, i), (mu, -1, 1)) / 2 * t**i for i in range(-2, mx + 1))
    qS, qT = (1 + 3 * mu) / 6, (1 - mu) / 2
    dS, dT, dn, dnh = sp.diff(FS, t), sp.diff(FT, t), sp.diff(n2, t), sp.diff(nh, t)
    B0 = a * (16 + 8 * a + 2 * b**2 / 9) - 24 * a**3 / t**2
    D0 = 12 * a / t**2 - (1 + mu) * (5 + 4 * a + b * (mu - R(1, 3)))
    cN = mean(B0 + 12 * a / t**2 * nh + 2 * (D0 * n2 + a * ((3 + mu) / 6 * FS - qT * FT) + db / 2 * (qS * dS + qT * dT)))
    alS, alT = (R(3, 2) + a) * (3 + mu) / 6 - mu * b * qS / 4, -(R(3, 2) + a) * qT - mu * b * qT / 4
    beS, beT, ben = (1 + mu) * db * qS / 4, (1 + mu) * db * qT / 4, -(1 + mu) * (db * (mu - R(1, 3)) + 2 * da)
    mb = mean(-2 * a * dnh + 2 * (alS * dS + alT * dT + 2 * mu * a * dn + beS * FS + beT * FT + ben * n2))
    nu3 = solve_lapse(1, cN, R(-14, 9), mx)
    hp = sp.expand(6 * (sp.diff(nu3, t) - nu3 / t) + 3 * mb)
    h3 = sum(hp.coeff(t, i) * t ** (i + 1) / (i + 1) for i in range(mx))
    c3 = curvature_tensors([kv, qv, tuple(-w for w in qv)])
    Qk = [[R(2, 3), 0, 0], [0, R(-1, 3), 0], [0, 0, R(-1, 3)]]
    onto = lambda m: sum(m[i][j].coef(7) * Qk[i][j] for i in rng for j in rng)
    assert sp.simplify(avg(c3["a"].coef(7)) - nu3.coeff(t, 2)) == 0 and sp.simplify(avg(c3["b"].coef(7)) - nu3.coeff(t, 4)) == 0
    assert (
        sp.simplify(R(3, 2) * avg(onto(c3["A"])) - h3.coeff(t, 2)) == 0
        and sp.simplify(R(3, 2) * avg(onto(c3["B"])) - h3.coeff(t, 4)) == 0
    )
    Ups = mean(a * (qS * FS + qT * FT))
    BR = 8 * a - 8 * a**2 - R(8, 3) * a**2 * b + R(4, 9) * a * b**2 + 4 * a * da * db - 24 * a**3 / t**2
    DR = 12 * a / t**2 - (1 + mu) * (5 + 2 * a - 6 * a * mu + b * (mu - R(1, 3)))
    rN = BR + 12 * a / t**2 * nh + 2 * (DR * n2 + a / 6 * (2 - 3 * mu - 3 * mu**2) * FS + a / 2 * (2 + mu - 3 * mu**2) * FT)
    LU = 3 * sp.diff(Ups, t, 2) + Ups - 6 * Ups / t**2
    assert sp.simplify(trunc(mean(rN) - cN + LU, 6)) == 0
    assert sp.simplify(trunc(mean(qS * JS + qT * JT) - R(2, 27) * (12 * a * b + 72 * a - b**2 + 36 - 18 * da * db), 6)) == 0
    Pp = -2 * a * nh + 2 * (alS * FS + alT * FT + 2 * mu * a * n2)
    Qq = 2 * da * nh + 2 * (
        ((1 + 2 * mu) * db / 4 * qS - da * (3 + mu) / 6) * FS
        + ((1 + 2 * mu) * db / 4 * qT + da * qT) * FT
        + (-(1 + mu) * db * (mu - R(1, 3)) - 2 * da * (1 + 2 * mu)) * n2
    )  # fmt: skip
    assert sp.simplify(trunc(mean(sp.diff(Pp, t) + Qq) - mb, 6)) == 0
    print("[ok] contracted cubic sources: exact-in-zeta t^2, t^4 data, normal form, averaged forcing, P'+Q split")


# ---------------------------------------------------------------- operators and V_l, exact
def check_operators():
    K = sp.Matrix(sp.symbols("K1 K2 K3"))
    v = sp.symbols("a0:5")
    A0 = sp.Matrix([[v[0], v[1], v[2]], [v[1], v[3], v[4]], [v[2], v[4], -v[0] - v[3]]])
    z = K.dot(K)
    B = lambda A: (lambda M: M - sp.eye(3) * M.trace() / 3)(z * A - K * (A * K).T - (A * K) * K.T)
    assert sp.simplify(B(B(B(A0))) - R(2, 3) * z * B(B(A0)) - z**2 / 3 * B(A0)) == sp.zeros(3)
    d, zz = sp.symbols("d z")
    T = lambda j, a: sum((-a) ** m * d ** (2 * m + j) / sp.factorial(2 * m + j) for m in range(8))
    for lam in (zz, 0, -zz / 3):
        series = sum((-1) ** n * d ** (2 * n + 1) / sp.factorial(2 * n + 1) * lam**n for n in range(8))
        poly = d - (T(3, zz) + 3 * T(3, -zz / 3)) / 4 * lam + (3 * T(5, zz) + T(5, -zz / 3)) / 4 * lam**2
        diff = sp.expand(series - poly)
        assert all(sp.simplify(diff.coeff(d, i)) == 0 for i in range(12))
    y = sp.symbols("y", positive=True)
    printed = {
        2: -2 * sp.I * (y**4 + 6 * sp.I * y**3 - 24 * y**2 - 54 * sp.I * y + 54) / (3 * y**5),
        3: -2 * (y**5 + 9 * sp.I * y**4 - 54 * y**3 - 204 * sp.I * y**2 + 450 * y + 450 * sp.I) / (3 * y**6),
    }
    for l, pr in printed.items():
        h = sp.simplify(sp.expand_func(sp.jn(l, y)) + sp.I * sp.expand_func(sp.yn(l, y)))
        V = sp.exp(-sp.I * y) * (h / 3 + sp.diff(h, y, 2))
        assert sp.simplify(sp.simplify(V.rewrite(sp.exp)) - pr) == 0
    print("[ok] B_K cubic relation, spatial Green polynomial, printed V_2 and V_3")


# ---------------------------------------------------------------- numerical validation only
def validate_numerically():
    mp.mp.dps = 25
    sq3 = mp.sqrt(3)
    j = lambda l, y: mp.sqrt(mp.pi / (2 * y)) * mp.besselj(l + 0.5, y)
    w = lambda l, y: j(l, y) / 3 + mp.diff(lambda z: j(l, z), y, 2)
    dsx = lambda x: 2 * mp.cos(x) / x - 2 * mp.sin(x) / x**2 + mp.sin(x)

    def T(l, xo, xs):
        pts = mp.linspace(xs, xo, int(xo - xs) + 2)
        return -3 * mp.quad(lambda x: dsx(x) * w(l, sq3 * (xo - x)) if x < xo else 0, pts)

    for xo, tab in [(0.3, 3.414e-4), (1, 0.0258), (3, 0.116), (4, 0.958), (8, 0.00059), (16, 0.0533)]:
        val = 4 * mp.pi * T(2, mp.mpf(xo), mp.mpf(xo) / 3) ** 2
        assert abs(val / tab - 1) < 0.01, (xo, val, tab)
    for xo, tab in [(2, 0.0340), (8, 0.0014), (16, 0.00091)]:
        val = 4 * mp.pi * T(2, mp.mpf(xo), mp.mpf(1)) ** 2
        assert abs(val / tab - 1) < 0.01, (xo, val, tab)
    for l, ratio, nodes in [
        (2, True, [2.49867, 6.38035, 7.90609, 9.46677, 11.96208]),
        (2, False, [2.51338, 6.33088, 7.50682, 8.52178, 10.89081, 13.06476]),
        (3, True, [3.05037, 7.00073, 8.69752, 10.60041, 13.04624]),
        (3, False, [3.05208, 6.92329, 8.48546, 9.38590, 11.66731, 13.83506]),
    ]:
        for n in nodes:
            assert abs(T(l, mp.mpf(n), mp.mpf(n) / 3 if ratio else mp.mpf(1))) < 1e-5, (l, n)
    print("[ok] numerical: C_2^(11) table entries and all listed T_l node positions")

    a = lambda s: mp.cos(s / sq3) - mp.sin(s / sq3) / (s / sq3)
    b = lambda s: -6 * (2 * mp.sin(s / sq3) / (s / sq3) - mp.cos(s / sq3) - 1)

    def JT(s, m):
        A, B = a(s), b(s)
        return (
            (1 - m)
            / 36
            * (
                B**2 * (9 * m**2 + 12 * m + 1)
                - 6 * B * (6 * m * A + 3 * m + 2)
                - 72 * A**2
                + 36
                - 36 * mp.diff(a, s) * mp.diff(b, s)
            )
        )

    cin = lambda x: mp.euler + mp.log(abs(x)) - mp.ci(abs(x)) if x != 0 else mp.mpf(0)

    def ird(r, X):
        x = r * X / sq3
        d, ap, am = 2 * r**2 - 3, 1 + 2 * r / sq3, 1 - 2 * r / sq3
        C = cin(am * X) + cin(ap * X) - 2 * cin(X)
        S = mp.si(am * X) + mp.si(ap * X) - 2 * mp.si(X)
        E = (
            -4 * r**2 * d * mp.sin(X)
            + 24 * r**2 * mp.cos(x) ** 2 / X
            - 48 * sq3 * r * mp.cos(x) * mp.sin(x) / X**2
            + 12 * (6 + d * X**2) * mp.sin(x) ** 2 / X**3
        )
        return 3 / (4 * r**6 * X) * (E + d**2 * (mp.sin(X) * C + mp.cos(X) * S))

    for tv, mv in [(1.3, 0.4), (5, -0.2), (9, -mp.mpf(1) / 3), (3, 0.9), (7, -0.8)]:
        tv, mv = mp.mpf(tv), mp.mpf(mv)
        K = mp.sqrt(2 * (1 + mv))
        green = mp.quad(lambda s, tv=tv, mv=mv, K=K: (s / tv) * mp.sin(K * (tv - s)) / K * JT(s, mv), [0, tv])
        s_of = lambda y: 2 * mp.sin(y) / y - mp.cos(y)
        x = tv / sq3
        kt = (1 - mv) / (2 * (1 + mv)) * ird(1 / K, K * tv) + mp.mpf(3) / 2 * (1 - mv) * (
            3 * mv * (s_of(x) - 1) ** 2 - mp.diff(s_of, x) ** 2
        )
        assert abs(green - kt) < 1e-10, (tv, mv)
    print("[ok] numerical: printed Kohri-Terada form of F_T against the Green integral of J_T, including the resonance")

    al, ga = mp.mpf("1.7"), mp.mpf("0.6")
    half = mp.mpf(3) / 2
    B0 = sum(
        (-(al**2) / 4) ** r
        * (-(ga**2) / 4) ** h
        / (mp.rf(half, r) * mp.rf(half, h) * mp.factorial(r) * mp.factorial(h) * (r + h + 1))
        for r in range(50)
        for h in range(50)
    )
    assert abs(B0 - (cin(al + ga) - cin(al - ga)) / (al * ga)) < 1e-18
    assert abs(mp.hyp0f1(2.5, -(al**2) / 4) - 3 * (mp.sin(al) - al * mp.cos(al)) / al**3) < 1e-18
    print("[ok] numerical: entire-function samples")


if __name__ == "__main__":
    check_degree_six()
    check_sources()
    check_operators()
    check_seed_identities()
    validate_numerically()
    print("ALL PRINTED-FORMULA CHECKS PASSED")
