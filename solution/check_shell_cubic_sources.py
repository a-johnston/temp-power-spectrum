"""Exact finite-time reductions of the shell-averaged cubic field sources.

Rebuilds N and the longitudinal momentum equation from matrix exponentials,
with arbitrary first-order time profiles and arbitrary second-order amplitudes.
No evolution, quadrature, truncated time series, or numerical evidence.
Run: .venv/bin/python -B solution/check_shell_cubic_sources.py
"""

import sys

import sympy as s
from engine import Ctx, Field, fexp, mexp, mmul, mscale, mtrace, t

mu, y = s.symbols("mu y", real=True)
k = s.Symbol("k", positive=True)
a, b = (s.Function(x)(t) for x in ("a", "b"))
identity = s.eye(3)


def reduce_shell(expr):
    """Clear denominators, then impose y^2=1-mu^2 exactly."""
    num, den = s.fraction(s.together(expr))
    out = 0
    for (power,), coeff in s.Poly(s.expand(num), y).terms():
        out += coeff * y ** (power % 2) * (1 - mu**2) ** (power // 2)
    return s.factor(out / den)


def assert_zero(expr):
    assert reduce_shell(expr) == 0


def check_sources():
    Ctx.maxdeg = 3
    p = s.Matrix([k, 0, 0])
    q = s.Matrix([k * mu, k * y, 0])
    vectors = {1: p, 2: q, 3: -q}
    Ctx.momenta = {j: tuple(v[i] for i in range(3)) for j, v in vectors.items()}
    Ctx.trunc = None
    zeta = Field({frozenset([j]): s.Integer(1) for j in vectors})
    nu = a * zeta
    H = [[Field() for _ in range(3)] for _ in range(3)]
    for lab, v in vectors.items():
        Q = v * v.T / k**2 - identity / 3
        for i in range(3):
            for j in range(3):
                H[i][j] += Field({frozenset([lab]): b * Q[i, j]})
    amps = {}
    for pair, sign in (((1, 2), 1), ((1, 3), -1)):
        n, fs, ft = (s.Function(f"{name}{'p' if sign == 1 else 'm'}")(t) for name in ("n", "S", "T"))
        total = p + sign * q
        diff = p - sign * q
        spatial = fs * (total * total.T / (2 * k**2 * (1 + sign * mu)) - identity / 3)
        spatial += ft * (diff * diff.T / (2 * k**2 * (1 - sign * mu)) - s.diag(0, 0, 1))
        amps[sign] = (n, fs, ft)
        nu += Field({frozenset(pair): n})
        for i in range(3):
            for j in range(3):
                H[i][j] += Field({frozenset(pair): spatial[i, j]})
    n0, z0 = (s.Function(x)(t) for x in ("n0", "Z0"))
    nu += Field({frozenset([2, 3]): n0})
    for i in range(3):
        for j in range(3):
            H[i][j] += Field({frozenset([2, 3]): z0 * (q[i] * q[j] / k**2 - identity[i, j] / 3)})

    M, Mi = mexp(H), mexp(mscale(H, -1))
    shear = mscale(mmul(Mi, [[v.dt() for v in row] for row in M]), s.Rational(1, 2))
    gi = mscale(Mi, fexp(-2 * zeta))
    density_inverse = mscale(Mi, fexp(zeta))
    lap = fexp(-3 * zeta) * sum((density_inverse[i][j] * nu.d_(j)).d_(i) for i in range(3) for j in range(3))
    grad2 = sum(gi[i][j] * nu.d_(i) * nu.d_(j) for i in range(3) for j in range(3))
    eq_n = nu.dt().dt() * 3 - lap - grad2 * 2 + mtrace(mmul(shear, shear)) + 3 / t**2 * (fexp(-2 * nu) - 1)
    cn = reduce_shell(eq_n.coeff([1, 2, 3]))

    # Exact divergence identity for det(M)=1 and tr(S)=0:
    # D_j S^j_i = partial_j S^j_i + 3 zeta_,j S^j_i
    #             - (1/2) tr(M^-1 partial_i M S).
    eq_m = []
    for i in range(3):
        Ci = mmul(Mi, [[v.d_(i) for v in row] for row in M])
        div = sum(shear[j][i].d_(j) + 3 * zeta.d_(j) * shear[j][i] for j in range(3))
        div -= mtrace(mmul(Ci, shear)) / 2
        eq_m.append(
            div
            - 2 * nu.dt().d_(i)
            + 2 / t * nu.d_(i)
            + 2 * sum(shear[j][i] * nu.d_(j) for j in range(3))
            + 2 * nu.dt() * nu.d_(i)
        )
    cm = reduce_shell(s.I * k * eq_m[0].coeff([1, 2, 3]))

    c111 = (
        a
        * k**2
        * (
            16
            + 8 * a
            + 4 * a * b * (mu**2 - s.Rational(1, 3))
            + b**2 * (mu**2 / 3 + s.Rational(1, 9))
            + 10 * b * (mu**2 - s.Rational(1, 3))
        )
        - 24 * a**3 / t**2
    )
    expect_n = c111 + 12 * a * n0 / t**2
    expect_n += (mu**2 - s.Rational(1, 3)) * (-a * k**2 * z0 + s.diff(b, t) * s.diff(z0, t) / 2)
    expect_m = -2 * a * k**2 * s.diff(n0, t)
    expect_m += k**2 * (mu**2 - s.Rational(1, 3)) * (-(s.Rational(3, 2) + a) * s.diff(z0, t) + b * s.diff(z0, t) / 4)
    for sign, (n, fs, ft) in amps.items():
        u = sign * mu
        qs, qt = (1 + 3 * u) / 6, (1 - u) / 2
        dn = 12 * a / t**2 - k**2 * (1 + u) * (5 + 4 * a + b * (u - s.Rational(1, 3)))
        expect_n += dn * n + a * k**2 * ((3 + u) * fs / 6 - (1 - u) * ft / 2)
        expect_n += s.diff(b, t) / 2 * (qs * s.diff(fs, t) + qt * s.diff(ft, t))
        expect_m += k**2 * (
            (s.Rational(3, 2) + a) * ((3 + u) * s.diff(fs, t) / 6 - (1 - u) * s.diff(ft, t) / 2)
            - u * b / 4 * (qs * s.diff(fs, t) + qt * s.diff(ft, t))
            + (1 + u) * s.diff(b, t) / 4 * (qs * fs + qt * ft)
            - (1 + u) * (s.diff(b, t) * (u - s.Rational(1, 3)) + 2 * s.diff(a, t)) * n
            + 2 * u * a * s.diff(n, t)
        )
    assert_zero(cn - expect_n)
    assert_zero(cm - expect_m)
    print("[ok] full cubic N and iK.cM sources: arbitrary time profiles and second-order amplitudes")
    if "--compare-engine" in sys.argv:
        from engine import equations

        original = equations(nu, zeta, H)
        assert_zero(cn - original["N"].coeff([1, 2, 3]))
        assert_zero(cm - s.I * k * original["M"][0].coeff([1, 2, 3]))
        print("[ok] independent comparison to original Christoffel-based field engine, arbitrary finite times")
    none = {f: 0 for values in amps.values() for f in values}
    none.update({n0: 0, z0: 0})
    assert_zero(cm.subs(none).doit())
    print("[ok] no pure-first-order cubic momentum source on the shell, pointwise in mu")
    c111_mean = s.integrate(c111, (mu, -1, 1)) / 2
    simple = a * k**2 * (16 + 8 * a + 2 * b**2 / 9) - 24 * a**3 / t**2
    assert s.expand(c111_mean - simple) == 0
    assert s.integrate(mu**2 - s.Rational(1, 3), (mu, -1, 1)) == 0
    assert s.expand(cn - expect_n + a * k**2 * amps[1][2]) != 0
    print("[ok] elementary pure-first-order mean; homogeneous tensor cancels, homogeneous lapse retained")
    print("[ok] omission control: deleting mixed tensor terms fails")


def check_normal_form():
    # Generic momenta, arbitrary second-order matrix jets. p is along x only
    # by rotation; the second-order momentum P is fully general.
    p, px, py, pz = s.symbols("p px py pz", real=True, nonzero=True)
    single = s.Matrix([p, 0, 0])
    pair = s.Matrix([px, py, pz])
    total = single + pair
    Q = single * single.T - identity * p**2 / 3
    mats = []
    for name in ("J", "V", "W"):
        xx, yy, xy, xz, yz = s.symbols(f"{name}xx {name}yy {name}xy {name}xz {name}yz")
        mats.append(s.Matrix([[xx, xy, xz], [xy, yy, yz], [xz, yz, -xx - yy]]))
    J, V, W = mats
    av, ad, add = s.symbols("av ad add")
    bdot = 6 / p**2 * (ad - av / t)
    dot = lambda X: s.trace(Q * X)
    phi = av * dot(J) / (2 * p**2)
    phidd = (add * dot(J) + 2 * ad * dot(V) + av * dot(W)) / (2 * p**2)
    ln_phi = 3 * phidd + (total.dot(total) - 6 / t**2) * phi
    ln_phi = ln_phi.subs(add, (-(p**2) * av + 6 * av / t**2) / 3)
    original = -av * (total.T * J * single)[0] + bdot * dot(V) / 2
    spatial_acceleration = W + 2 / t * V
    remainder = -av * (total.T * J * single)[0]
    remainder -= av / (2 * p**2) * ((total.dot(total) - p**2) * dot(J) + 3 * dot(spatial_acceleration))
    assert s.expand(original - ln_phi - remainder) == 0
    assert s.expand(original - ln_phi - remainder - av * dot(V) / t) != 0
    print("[ok] general-momentum cubic lapse normal form removes every H2' term")

    # Shell coefficients after substituting the complete scalar/tensor equations.
    n, fs, ft, js, jt = s.symbols("n S T JS JT")
    qs, qt = (1 + 3 * mu) / 6, (1 - mu) / 2
    kap2 = 2 * k**2 * (1 + mu)
    dn = 12 * a / t**2 - k**2 * (1 + mu) * (5 + 4 * a + b * (mu - s.Rational(1, 3)))
    reduced = dn * n + a * k**2 * ((3 + mu) * fs / 6 - (1 - mu) * ft / 2)
    reduced -= 3 * a / 2 * (qs * (kap2 * fs / 3 - 4 * kap2 * n + js) + qt * (-kap2 * ft + jt))
    target = (12 * a / t**2 - k**2 * (1 + mu) * (5 + 2 * a - 6 * a * mu + b * (mu - s.Rational(1, 3)))) * n
    target += a * k**2 * (2 - 3 * mu - 3 * mu**2) * fs / 6 + a * k**2 * (1 - mu) * (2 + 3 * mu) * ft / 2
    target -= 3 * a / 2 * (qs * js + qt * jt)
    assert s.expand(reduced - target) == 0
    for f in (n, fs, ft):
        assert s.Poly(s.expand(target).coeff(f), mu).degree() <= 2
    print("[ok] reduced shell lapse source uses only degree <= 2 moments, with no response derivatives")


def check_second_sources():
    from engine import equations

    Ctx.maxdeg = 2
    Ctx.momenta = {1: (k, 0, 0), 2: (k * mu, k * y, 0)}
    z = Field({frozenset([j]): s.Integer(1) for j in (1, 2)})
    H = [
        [
            Field(
                {
                    frozenset([j]): b * (s.sympify(Ctx.momenta[j][i]) * Ctx.momenta[j][l] / k**2 - s.Rational(int(i == l), 3))
                    for j in (1, 2)
                }
            )
            for l in range(3)
        ]
        for i in range(3)
    ]
    original = equations(a * z, z, H)
    total = s.Matrix([k * (1 + mu), k * y, 0])
    diff = s.Matrix([k * (1 - mu), -k * y, 0])
    U = total * total.T / (2 * k**2 * (1 + mu))
    V = diff * diff.T / (2 * k**2 * (1 - mu)) - s.diag(0, 0, 1)
    C = s.Matrix(3, 3, lambda i, j: original["T"][i][j].coeff([1, 2]))
    cn = k**2 * (2 * a * b / 3 - 4 * a + mu * (4 * a**2 - 4 * a * b / 3 + 2 * a) - 2 * mu**2 * a * b)
    cn += s.diff(b, t) ** 2 / 2 * (mu**2 - s.Rational(1, 3)) + 12 * a**2 / t**2
    js = k**2 * (
        2 * (3 + mu) * a**2
        + a * b * (15 * mu**2 + 13 * mu - 8) / 3
        + 16 * (1 + mu) * a
        + b**2 * (5 - 5 * mu - 15 * mu**2 - 9 * mu**3) / 36
        + b * (-10 + 7 * mu + 15 * mu**2) / 6
        + 5
        + 7 * mu
    )
    js -= (1 + 3 * mu) * s.diff(a, t) * s.diff(b, t)
    jt = (
        (1 - mu)
        / 36
        * (
            k**2 * (b**2 * (9 * mu**2 + 12 * mu + 1) - 6 * b * (6 * mu * a + 3 * mu + 2) - 72 * a**2 + 36)
            - 36 * s.diff(a, t) * s.diff(b, t)
        )
    )
    assert_zero(original["N"].coeff([1, 2]) - cn)
    assert_zero(-3 * s.trace(U * C) - js)
    assert_zero(-s.trace(V * C) - jt)
    print("[ok] complete elementary quadratic sources: N2 degree 2; JS and JT degree 3")
    forcing_mean = s.integrate(s.expand((1 + 3 * mu) * js / 6 + (1 - mu) * jt / 2), (mu, -1, 1)) / 2
    expected_mean = s.Rational(2, 27) * (k**2 * (12 * a * b + 72 * a - b**2 + 36) - 18 * s.diff(a, t) * s.diff(b, t))
    assert s.expand(forcing_mean - expected_mean) == 0
    reduced_local = k**2 * (8 * a - 8 * a**2 - s.Rational(8, 3) * a**2 * b + s.Rational(4, 9) * a * b**2)
    reduced_local += 4 * a * s.diff(a, t) * s.diff(b, t) - 24 * a**3 / t**2
    old_local = a * k**2 * (16 + 8 * a + 2 * b**2 / 9) - 24 * a**3 / t**2
    assert s.expand(reduced_local - old_local + 3 * a * forcing_mean) == 0
    print("[ok] forcing moments in reduced cubic lapse source evaluated completely")

    x = k * t / s.sqrt(3)
    actual_a = s.cos(x) - s.sin(x) / x
    actual_b = -6 * (2 * s.sin(x) / x - s.cos(x) - 1)
    assert s.simplify(3 * s.diff(actual_a, t, 2) + (k**2 - 6 / t**2) * actual_a) == 0
    assert s.simplify(s.diff(actual_b, t) - 6 * (s.diff(actual_a, t) - actual_a / t)) == 0
    print("[ok] exact first-order identities required by the normal form")


def check_momentum_adjoint():
    n, fs, ft, n0, w, phi, P = (s.Function(x)(t) for x in ("n", "S", "T", "n0", "w", "phi", "P"))
    qs, qt = (1 + 3 * mu) / 6, (1 - mu) / 2
    AS = (s.Rational(3, 2) + a) * (3 + mu) / 6 - mu * b * qs / 4
    AT = -(s.Rational(3, 2) + a) * qt - mu * b * qt / 4
    BS, BT = (1 + mu) * s.diff(b, t) * qs / 4, (1 + mu) * s.diff(b, t) * qt / 4
    BN = -(1 + mu) * (s.diff(b, t) * (mu - s.Rational(1, 3)) + 2 * s.diff(a, t))
    original = AS * s.diff(fs, t) + AT * s.diff(ft, t) + 2 * mu * a * s.diff(n, t) + BS * fs + BT * ft + BN * n
    local_p = AS * fs + AT * ft + 2 * mu * a * n
    local_q = ((1 + 2 * mu) * s.diff(b, t) * qs / 4 - s.diff(a, t) * (3 + mu) / 6) * fs
    local_q += ((1 + 2 * mu) * s.diff(b, t) * qt / 4 + s.diff(a, t) * qt) * ft
    local_q += (-(1 + mu) * s.diff(b, t) * (mu - s.Rational(1, 3)) - 2 * s.diff(a, t) * (1 + 2 * mu)) * n
    assert s.expand(original - s.diff(local_p, t) - local_q) == 0
    assert s.expand(-2 * a * s.diff(n0, t) - s.diff(-2 * a * n0, t) - 2 * s.diff(a, t) * n0) == 0
    Q = s.Function("Q")(t)
    original = -3 * w * (s.diff(phi, t) - phi / t) + s.Rational(3, 2) * w * (s.diff(P, t) + Q)
    Z = P - 2 * phi
    bulk = s.Rational(3, 2) * (w * Q - s.diff(w, t) * Z + 2 * w * phi / t)
    boundary_derivative = s.Rational(3, 2) * s.diff(w * Z, t)
    assert s.expand(original - bulk - boundary_derivative) == 0
    print("[ok] derivative-free momentum adjoint and combined normal-form boundary terms, arbitrary jets")


def check_green_functions():
    r, c = s.symbols("r c", positive=True)
    d = t - r
    G = -(s.sin(c * d) / c + (s.sin(c * d) - c * d * s.cos(c * d)) / (c**3 * t * r)) / 3
    assert s.simplify(3 * s.diff(G, t, 2) + (3 * c * c - 6 / t**2) * G) == 0
    assert s.simplify(G.subs(t, r)) == 0
    assert s.simplify(s.diff(G, t).subs(t, r)) == -s.Rational(1, 3)
    assert s.simplify(s.limit(G, c, 0) + (t**3 - r**3) / (9 * t * r)) == 0
    for sign in (-1, 1):
        f = s.sin(c * d) / c if sign == 1 else s.sinh(c * d) / c
        g = r / t * f
        assert s.simplify(s.diff(g, t, 2) + 2 / t * s.diff(g, t) + sign * c * c * g) == 0
        assert s.simplify(s.diff(g, t).subs(t, r)) == 1
    print("[ok] entire lapse, scalar and tensor Green kernels, including exact zero momentum")


def check_angular_functions():
    # Taylor coefficients are identities for every nonnegative integer order,
    # not a finite time-series test of the field solution.
    m, n, ell = s.symbols("m n ell", integer=True, nonnegative=True)
    beta, delta = s.symbols("beta delta", positive=True)
    coefficient = s.rf(m + 1, n) / (s.rf(m + 2, n) * (m + 1))
    assert s.simplify(s.combsimp(coefficient) - 1 / (m + n + 1)) == 0
    coefficient2 = s.rf(m + 1, n + ell) / (s.rf(m + 2, n + ell) * (m + 1))
    assert s.simplify(s.combsimp(coefficient2) - 1 / (m + n + ell + 1)) == 0
    print("[ok] exact all-order coefficient identities for 1F2 and two-variable entire angular kernels")

    alpha, gamma = s.symbols("alpha gamma", nonzero=True)

    def factors(order, freq):
        if order == 3:
            return [(freq, -1, 1 / (2 * s.I * freq)), (-freq, -1, -1 / (2 * s.I * freq))]
        return [
            (freq, -3, 3 / (2 * s.I * freq**3)),
            (-freq, -3, -3 / (2 * s.I * freq**3)),
            (freq, -2, -3 / (2 * freq**2)),
            (-freq, -2, -3 / (2 * freq**2)),
        ]

    maxpole = 0
    for moment in range(6):
        for left, right in ((3, None), (5, None), (3, 3), (3, 5)):
            terms = [(w, power + 2 * moment + 1, 2 * coeff) for w, power, coeff in factors(left, alpha)]
            if right is not None:
                terms = [(w + v, power + q, coeff * cf) for w, power, coeff in terms for v, q, cf in factors(right, gamma)]
            # Taylor subtractions in J remove exactly these principal parts.
            # Their sum must vanish, proving that no remainder is discarded.
            for power in range(min(p for _, p, _ in terms), 0):
                singular = sum(coeff * (s.I * w) ** (power - p) / s.factorial(power - p) for w, p, coeff in terms if p <= power)
                assert s.simplify(singular) == 0
            maxpole = max(maxpole, -min(p for _, p, _ in terms))
    assert maxpole == 3
    print("[ok] all required angular kernels reduce to finite existing Ephi/J sums, at most J3; no principal-part remainder")


if __name__ == "__main__":
    check_sources()
    check_normal_form()
    check_second_sources()
    check_momentum_adjoint()
    check_green_functions()
    check_angular_functions()
    print("SHELL CUBIC SOURCE REDUCTIONS VERIFIED (exact symbolic algebra)")
