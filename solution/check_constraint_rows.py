"""Constraint rows remove the hyperbolic scalar Green family from physical pair responses.

For a nonzero pair momentum K the five trace-free spatial rows (T) can be replaced by the three
momentum rows (M) and the two TT rows of (T).  The scalar/vector spatial response then follows
from one cumulative integral,

    K.A(t) = V(t) = int_0^t [4K(nu' - nu/s) + 2i c_M] ds,     nonTT(A) = P_K[V],

instead of the hyperbolic (eigenvalue -K^2/3) Cartesian Green operator.  The replacement is
exact on every source family obeying the linearised Bianchi compatibility

    (C1)  c_Ham' - 4 c_N / t = 2i K.c_M,
    (C2)  c_M' + 2 c_M / t = i K.c_T - (i K / 6)(c_Ham + 4 c_N),

together with the regular lapse coefficient fixed by the initial Hamiltonian row.

Checks (exact symbolic algebra; the series test is an exact coefficient identity):
  1. (C1), (C2) for a physical second-order pair at general UNEQUAL momenta and angle, with the
     first-order fields reduced by their own field equations (no time expansion).
  2. The projector P_K, its transpose, and the projected weight used by the constraint-row adjoint.
  3. On the monochromatic shell: the explicit longitudinal source m(t, mu), the local scalar
     identity J_S = m' + 2m/t - 2C_N - c_Ham/2, the initial row, and the closed time integral of m.
  4. The evaluated acoustic kernel  int_r^t G_N^+(s,r)/s ds = (T_3 + r T_2)/(3 t r).
  5. The closed acoustic scalar amplitude S_2 agrees with the hyperbolic Cartesian solution
     coefficient by coefficient (symbolic k, mu), and K.A = V for a general unequal pair.
Run: .venv/bin/python -B solution/check_constraint_rows.py
"""

import sympy as s
from engine import Ctx, Field, I, equations, t

mu, y = s.symbols("mu y", real=True)
p, q, k = s.symbols("p q k", positive=True)


def reduce_y(expr):
    """Impose y^2 = 1 - mu^2 exactly after clearing denominators."""
    num, den = s.fraction(s.together(expr))
    out = 0
    for (power,), coeff in s.Poly(s.expand(num), y).terms():
        out += coeff * y ** (power % 2) * (1 - mu**2) ** (power // 2)
    return s.expand(out) / den


def pair_sources(k1, k2, profiles):
    """Exact second-order sources (c_N, c_M, c_T, c_Ham) of the pair (k1, k2).

    profiles = ((a1, b1, |k1|), (a2, b2, |k2|)): first fields zeta = 1, nu = a, H = b (k^k^ - 1/3).
    The unknown second-order fields are included and split off, so the returned sources are
    exactly the multilinear coefficients with the top-order fields set to zero.
    """
    Ctx.maxdeg, Ctx.trunc = 2, None
    Ctx.momenta = {1: tuple(k1), 2: tuple(k2)}
    nu, zeta = Field(), Field()
    H = [[Field() for _ in range(3)] for _ in range(3)]
    for lab, kv, (a, b, km) in ((1, k1, profiles[0]), (2, k2, profiles[1])):
        Q = kv * kv.T / km**2 - s.eye(3) / 3
        nu += Field({frozenset([lab]): a})
        zeta += Field({frozenset([lab]): s.Integer(1)})
        for i in range(3):
            for j in range(3):
                H[i][j] += Field({frozenset([lab]): b * Q[i, j]})
    n2 = s.Function("n2")(t)
    A = s.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            if (i, j) != (2, 2):
                A[i, j] = A[j, i] = s.Function(f"A{i}{j}")(t)
    A[2, 2] = -A[0, 0] - A[1, 1]
    nu += Field({frozenset([1, 2]): n2})
    H = [[H[i][j] + Field({frozenset([1, 2]): A[i, j]}) for j in range(3)] for i in range(3)]
    E = equations(nu, zeta, H)
    unknowns = [n2] + [A[i, j] for i in range(3) for j in range(i, 3) if (i, j) != (2, 2)]
    zero = {}
    for f in unknowns:
        zero.update({f.diff(t, 2): 0, f.diff(t): 0, f: 0})

    def src(e):
        return s.expand(e).subs(zero)

    K = k1 + k2
    K2 = reduce_y((K.T * K)[0])
    # The Hamiltonian row is normalised so that L_Ham = 4K^2 nu + 12 nu'/t + 12 nu/t^2 - K.A.K.
    lin_ham = s.expand(E["Ham"].coeff([1, 2]) - src(E["Ham"].coeff([1, 2])))
    expect = 4 * K2 * n2 + 12 * n2.diff(t) / t + 12 * n2 / t**2 - (K.T * A * K)[0]
    assert reduce_y(lin_ham - expect) == 0, "Hamiltonian-row normalisation"
    return (
        src(E["N"].coeff([1, 2])),
        [src(E["M"][i].coeff([1, 2])) for i in range(3)],
        s.Matrix(3, 3, lambda i, j: src(E["T"][i][j].coeff([1, 2]))),
        src(E["Ham"].coeff([1, 2])),
        K,
        K2,
    )


def on_shell(expr, profiles):
    """Reduce first-order profiles by their exact equations: L_k a = 0 and the first Hamiltonian
    row b = 6 + 6a + 18 (a'/t + a/t^2)/k^2 (which implies b' = 6(a' - a/t))."""
    for a, b, km in profiles:
        bexpr = 6 + 6 * a + 18 * (a.diff(t) / t + a / t**2) / km**2
        expr = expr.subs({b.diff(t, n): bexpr.diff(t, n) for n in (3, 2, 1)}).subs(b, bexpr)
        add = -(km**2) * a / 3 + 2 * a / t**2
        for _ in range(4):
            expr = s.expand(expr.subs({a.diff(t, 4): add.diff(t, 2), a.diff(t, 3): add.diff(t)}))
            expr = expr.subs(a.diff(t, 2), add)
    return reduce_y(s.expand(expr))


def check_compatibility():
    a1, b1, a2, b2 = (s.Function(n)(t) for n in ("a1", "b1", "a2", "b2"))
    prof = ((a1, b1, p), (a2, b2, q))
    k1, k2 = s.Matrix([p, 0, 0]), s.Matrix([q * mu, q * y, 0])
    cN, cM, cT, cH, K, _ = pair_sources(k1, k2, prof)
    c1 = cH.diff(t) - 4 * cN / t - 2 * I * sum(K[i] * cM[i] for i in range(3))
    assert reduce_y(c1) != 0, "control: (C1) should need the first-order field equations"
    assert on_shell(c1, prof) == 0, "(C1)"
    for i in range(3):
        KcT = sum(K[j] * cT[j, i] for j in range(3))
        c2 = cM[i].diff(t) + 2 * cM[i] / t - I * KcT + I * K[i] / 6 * (cH + 4 * cN)
        assert on_shell(c2, prof) == 0, ("(C2)", i)
    print("[ok] (C1), (C2) hold exactly for a physical pair at unequal |k1|, |k2| and arbitrary angle")


def projector(K):
    K2 = s.expand((K.T * K)[0])

    def P(V):
        out = (K * V.T + V * K.T) / K2
        out -= s.eye(3) * out.trace() / 3
        return out - (K.T * V)[0] * (K * K.T - s.eye(3) * K2 / 3) / (2 * K2**2)

    def PT(W):
        return 2 * W * K / K2 - (K.T * W * K)[0] * K / (2 * K2**2)

    return P, PT


def check_projector():
    K = s.Matrix(s.symbols("K1:4", real=True))
    V = s.Matrix(s.symbols("V1:4"))
    W = s.Matrix(3, 3, lambda i, j: s.Symbol(f"W{min(i, j)}{max(i, j)}"))
    W[2, 2] = -W[0, 0] - W[1, 1]
    P, PT = projector(K)
    A = P(V)
    assert s.simplify(A.trace()) == 0 and s.simplify(A - A.T) == s.zeros(3, 3)
    assert s.simplify(A * K - V) == s.zeros(3, 1), "K.P[V] = V"
    K2 = (K.T * K)[0]
    Pperp = s.eye(3) - K * K.T / K2
    tt = Pperp * A * Pperp - Pperp * (Pperp * A * Pperp).trace() / 2
    assert s.simplify(tt) == s.zeros(3, 3), "P[V] has no TT part"
    assert s.simplify(sum(W[i, j] * A[i, j] for i in range(3) for j in range(3)) - (PT(W).T * V)[0]) == 0
    assert s.simplify((K.T * PT(W))[0] - s.Rational(3, 2) * (K.T * W * K)[0] / K2) == 0
    print("[ok] projector P_K: K.P[V]=V, no TT part, transpose (2/K^2)WK - (KWK)K/(2K^4)")


def T(j, a, d, order=None):
    """Entire factorial tails T_j(a,d) = sum (-a)^m d^(2m+j)/(2m+j)!; closed form for a = c^2."""
    c = s.sqrt(a)
    closed = {
        0: s.cos(c * d),
        1: s.sin(c * d) / c,
        2: (1 - s.cos(c * d)) / a,
        3: (d - s.sin(c * d) / c) / a,
    }
    return closed[j]


def shell_m(a, b):
    ap, bp = a.diff(t), b.diff(t)
    return bp * ((1 - 3 * mu) * (s.Rational(3, 2) + a) + (3 * mu**2 - 1) * b / 4) - 6 * a * ap


def check_shell():
    a, b = s.Function("a")(t), s.Function("b")(t)
    prof = ((a, b, k), (a, b, k))
    k1, k2 = s.Matrix([k, 0, 0]), s.Matrix([k * mu, k * y, 0])
    cN, cM, cT, cH, K, K2 = pair_sources(k1, k2, prof)
    Kc = s.Matrix(cM)
    m = reduce_y(3 * I * (K.T * Kc)[0] / K2)
    assert s.simplify(m - shell_m(a, b)) == 0, "explicit longitudinal source"
    assert K.cross(Kc).applyfunc(lambda e: s.simplify(reduce_y(e))) == s.zeros(3, 1), "c_M parallel to K on the shell"
    Qh = K * K.T / K2 - s.eye(3) / 3
    JS = reduce_y(-3 * sum(Qh[i, j] * cT[j, i] for i in range(3) for j in range(3)))
    CN = reduce_y(cN)
    local = JS - (m.diff(t) + 2 * m / t - 2 * CN - reduce_y(cH) / 2)
    assert on_shell(local, prof) == 0, "local scalar identity"
    # Closed time integral of m, with I_a(t) = int_0^t a(s)^2/s ds.
    Ia = s.Function("Ia")(t)
    Mm = (1 - 3 * mu) * (s.Rational(3, 2) * b + 3 * a**2 - 6 * Ia) + (3 * mu**2 - 1) * b**2 / 8 - 3 * a**2
    dM = Mm.diff(t).subs(Ia.diff(t), a**2 / t)
    assert on_shell(dM - shell_m(a, b), prof) == 0, "d/dt of closed integral of m"
    # Initial Hamiltonian row on the physical profile, and M(0) = 0.
    x = s.Symbol("x", positive=True)
    aex, bex = s.cos(x) - s.sin(x) / x, -6 * (2 * s.sin(x) / x - s.cos(x) - 1)
    phys = {f: v.subs(x, k * t / s.sqrt(3)) for f, v in ((a, aex), (b, bex))}
    derivs = {f.diff(t): v.diff(t) for f, v in phys.items()}
    cH0 = s.limit(reduce_y(cH).subs(derivs).subs(phys), t, 0)
    assert s.simplify(-cH0 / 36 - k**2 * (4 - mu) / 9) == 0, "initial row c_r = k^2(4-mu)/9"
    assert all(s.limit(v, t, 0) == 0 for v in phys.values())
    # I_a in closed form: (1/2)Cin(2x) + sin(2x)/(2x) - sin(x)^2/(2x^2) - 1/2, with Cin' (z) = (1-cos z)/z.
    closed_ia = s.sin(2 * x) / (2 * x) - s.sin(x) ** 2 / (2 * x**2) - s.Rational(1, 2)
    d_ia = (1 - s.cos(2 * x)) / (2 * x) + closed_ia.diff(x)
    assert s.simplify(s.expand_trig(d_ia - aex**2 / x)) == 0 and s.limit(closed_ia, x, 0) == 0
    print("[ok] shell: c_M || K, explicit m(t,mu), J_S = m' + 2m/t - 2C_N - c_Ham/2, c_r = k^2(4-mu)/9")
    print("[ok] shell: int_0^t m = (1-3mu)[3b/2 + 3a^2 - 6 int a^2/s] + (3mu^2-1) b^2/8 - 3a^2")


def check_acoustic_kernel():
    cc, r = s.symbols("c r", positive=True)
    u1 = lambda X: s.sin(X) / X - s.cos(X)  # noqa: E731
    u2 = lambda X: -s.cos(X) / X - s.sin(X)  # noqa: E731
    G = lambda tt, ss: (u2(cc * tt) * u1(cc * ss) - u1(cc * tt) * u2(cc * ss)) / (3 * cc)  # noqa: E731
    d = t - r
    closed = (T(3, cc**2, d) + r * T(2, cc**2, d)) / (3 * t * r)
    assert s.simplify(s.diff(closed, t) - G(t, r) / t) == 0 and s.simplify(closed.subs(t, r)) == 0
    # The combined kernel of the scalar amplitude: 6 G - 6 int G/s.
    comb = 6 * G(t, r) - 6 * closed
    target = 2 * T(1, cc**2, d) + 2 * ((d - r) * T(2, cc**2, d) - 2 * T(3, cc**2, d)) / (t * r)
    assert s.simplify(s.expand_trig(comb - target)) == 0
    print("[ok] int_r^t G_N^+(s,r)/s ds = (T3 + r T2)/(3tr); scalar kernel 2T1 + 2[(d-r)T2 - 2T3]/(tr)")


def series(expr, n):
    return s.expand(s.series(expr, t, 0, n + 1).removeO())


def check_shell_series(order=16):
    x = k * t / s.sqrt(3)
    A = series(s.cos(x) - s.sin(x) / x, order)
    B = series(-6 * (2 * s.sin(x) / x - s.cos(x) - 1), order)
    cut = lambda e: sum(s.expand(e).coeff(t, n) * t**n for n in range(order + 1))  # noqa: E731
    Ap, Bp = A.diff(t), B.diff(t)
    third = s.Rational(1, 3)
    CN = cut(
        k**2 * (2 * third * A * B - 4 * A + mu * (4 * A**2 - 4 * third * A * B + 2 * A) - 2 * mu**2 * A * B)
        + Bp**2 * (mu**2 - third) / 2
        + 12 * s.expand(A**2 / t**2)
    )
    JS = cut(
        k**2
        * (
            2 * (3 + mu) * A**2
            + (15 * mu**2 + 13 * mu - 8) * A * B / 3
            + 16 * (1 + mu) * A
            + (5 - 5 * mu - 15 * mu**2 - 9 * mu**3) * B**2 / 36
            + (-10 + 7 * mu + 15 * mu**2) * B / 6
            + 5
            + 7 * mu
        )
        - (1 + 3 * mu) * Ap * Bp
    )
    K2 = 2 * k**2 * (1 + mu)
    n = {2: k**2 * (4 - mu) / 9}
    for m in range(4, order + 1, 2):
        n[m] = s.expand((-CN.coeff(t, m - 2) - K2 * n[m - 2]) / (3 * m * (m - 1) - 6))
    nser = sum(c * t**m for m, c in n.items())
    rhs = s.expand(-4 * K2 * nser + JS)
    S = {0: 0}
    for m in range(2, order + 1, 2):  # hyperbolic Cartesian scalar sector, regular branch
        S[m] = s.expand((rhs.coeff(t, m - 2) + K2 * S[m - 2] / 3) / (m * (m + 1)))
    cart = sum(c * t**m for m, c in S.items())
    Ia = s.integrate(s.expand(A**2 / t), (t, 0, t))
    In = s.integrate(s.expand(nser / t), (t, 0, t))
    mom = 6 * nser - 6 * In + (1 - 3 * mu) * (s.Rational(3, 2) * B + 3 * A**2 - 6 * Ia)
    mom += (3 * mu**2 - 1) * B**2 / 8 - 3 * A**2
    assert cut(mom - cart) == 0
    # Same amplitude from the closed acoustic kernel (series of the Green integral).
    aK = K2 / 3
    r = s.Symbol("r", positive=True)
    tail = lambda j, d: sum((-aK) ** i * d ** (2 * i + j) / s.factorial(2 * i + j) for i in range(order))  # noqa: E731
    d = t - r
    kern = 2 * tail(1, d) + 2 * s.expand(((d - r) * tail(2, d) - 2 * tail(3, d)) / (t * r))
    green = s.integrate(s.expand(kern * CN.subs(t, r)), (r, 0, t))
    init = 18 * n[2] * (tail(2, t) - 2 * s.expand(tail(3, t) / t))
    closed = cut(init - green + (1 - 3 * mu) * (s.Rational(3, 2) * B + 3 * A**2 - 6 * Ia))
    closed += cut((3 * mu**2 - 1) * B**2 / 8 - 3 * A**2)
    assert cut(closed - cart) == 0
    print(f"[ok] shell: acoustic closed S_2 equals the hyperbolic Cartesian solution through t^{order} (symbolic k, mu)")


def check_general_pair_series(order=10):
    """K.A_Cartesian = V for an unequal pair, coefficient by coefficient (exact rationals in t)."""
    pv, qv, mv = s.Integer(1), s.Rational(3, 2), s.Rational(1, 3)  # |k1|=1, |k2|=3/2, cos = 1/3
    yv = s.sqrt(1 - mv**2)
    a1, b1, a2, b2 = (s.Function(n)(t) for n in ("a1", "b1", "a2", "b2"))
    k1, k2 = s.Matrix([pv, 0, 0]), s.Matrix([qv * mv, qv * yv, 0])
    cN, cM, cT, cH, K, K2 = pair_sources(k1, k2, ((a1, b1, pv), (a2, b2, qv)))
    fields = {}
    for a, b, km in ((a1, b1, pv), (a2, b2, qv)):
        x = km * t / s.sqrt(3)
        fields[a] = series(s.cos(x) - s.sin(x) / x, order + 4)
        fields[b] = series(-6 * (2 * s.sin(x) / x - s.cos(x) - 1), order + 4)

    def ev(e):
        e = e.subs({f.diff(t, n): fields[f].diff(t, n) for f in fields for n in (2, 1)}).subs(fields)
        e = s.expand(e)
        return sum(e.coeff(t, n) * t**n for n in range(-2, order + 1))

    cN, cH = ev(cN), ev(cH)
    cM = [ev(c) for c in cM]
    cT = cT.applyfunc(ev)
    c_r = -cH.coeff(t, 0) / 36  # initial Hamiltonian row
    assert s.simplify(c_r - (2 * (pv**2 + qv**2) - (k1.T * k2)[0]) / 9) == 0, "report's c_2 kernel"

    def mismatch(c_reg):
        n = {2: c_reg}
        for m in range(4, order + 1, 2):
            n[m] = s.expand((-cN.coeff(t, m - 2) - K2 * n[m - 2]) / (3 * m * (m - 1) - 6))
        nser = sum(c * t**m for m, c in n.items())
        Bk = lambda X: X * K2 - K * (X * K).T - (X * K) * K.T + s.eye(3) * 2 * (K.T * X * K)[0] / 3  # noqa: E731
        QK = K * K.T - s.eye(3) * K2 / 3
        rhs = (-4 * QK * nser - 2 * cT).applyfunc(s.expand)
        Am = {0: s.zeros(3, 3)}
        for m in range(1, order + 1):
            prev = Am.get(m - 2, s.zeros(3, 3))
            Am[m] = ((rhs.applyfunc(lambda e, m=m: e.coeff(t, m - 2)) - Bk(prev)) / (m * (m + 1))).applyfunc(s.expand)
        Acart = sum((Am[m] * t**m for m in Am), s.zeros(3, 3))
        Vp = [4 * K[i] * (nser.diff(t) - s.expand(nser / t)) + 2 * I * cM[i] for i in range(3)]
        V = s.Matrix([s.integrate(s.expand(v), (t, 0, t)) for v in Vp])
        cut = lambda e: sum(s.expand(e).coeff(t, n) * t**n for n in range(order + 1))  # noqa: E731
        return (Acart * K - V).applyfunc(lambda e: s.simplify(cut(e)))

    assert mismatch(c_r) == s.zeros(3, 1)
    assert mismatch(c_r + s.Rational(1, 7)) != s.zeros(3, 1), "control: the initial row is required"
    print(f"[ok] unequal pair (|k2|/|k1|=3/2, cos=1/3): Cartesian K.A equals the momentum-row V through t^{order}")
    print("[ok] control: changing the initial Hamiltonian row breaks the equality")


if __name__ == "__main__":
    check_compatibility()
    check_projector()
    check_shell()
    check_acoustic_kernel()
    check_shell_series()
    check_general_pair_series()
    print("CONSTRAINT-ROW REPLACEMENT VERIFIED: no hyperbolic scalar kernel is needed for physical pair sources")
