"""Direct cubic measurement checks: exact off-shell fixtures, not spectra."""

from itertools import combinations

import sympy as s
from check_entire_sphere import reduce_radial
from check_quadrupole_adjoint_expansion import moment
from direct_cubic_measurement import (
    FourierFirst,
    RayFirst,
    compose_mixed_exponents,
    contract_direct_measurement,
    first_only_exponent_triple,
    measurement_phase_plan,
    ordered_radial_transport,
    polarized_bulk,
)
from entire_sphere import phase_average
from stationary_cubic import compose_first_only, gradient_linear_driver, linear_direction_driver, reduced_cubic_bulk

n = s.Matrix(s.symbols("n0:3", real=True))
time = s.Symbol("time", positive=True)


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, expr


def average_poly(expr):
    return s.Add(*(cf * moment(powers) for powers, cf in s.Poly(s.expand(expr), *n).terms()))


def canonical_homogeneous(seeds, C, E, O, I):
    """Conserved photon covector and scalar implicit roots; no rays/adjoints.

    Only the nilpotent label arithmetic is reused from engine. No field or
    geodesic equations, source kernels, or observable expansions are called.
    """
    import engine as eng

    saved = eng.Ctx.maxdeg, eng.Ctx.trunc
    eng.Ctx.maxdeg = 3
    eng.Ctx.trunc = s.expand
    try:

        def monomial(index):
            return eng.Field({frozenset([index + 1]): s.S.One})

        def power(F, k):
            out = eng.Field.const(1)
            for _ in range(k):
                out = out * F
            return out

        def substitute(F, T):
            result = eng.Field()
            for labels, expr in F.d.items():
                P = s.Poly(expr, time)
                for (degree,), coefficient in P.terms():
                    result += eng.Field({labels: coefficient}) * power(T, degree)
            return result

        a = sum((monomial(k) * s.sympify(seed.lapse) for k, seed in enumerate(seeds)), eng.Field())
        h = [
            [sum((monomial(k) * s.Matrix(seed.spatial)[i, j] for k, seed in enumerate(seeds)), eng.Field()) for j in range(3)]
            for i in range(3)
        ]

        def root(T0):
            T = eng.Field.const(T0)
            for _ in range(3):
                T = T0 * eng.fexp(-substitute(a, T))
            return T

        Te, Ti = root(E), root(I)
        proper = eng.fexp(a).apply(lambda f: s.integrate(time * f, time))
        To = eng.Field.const(O)
        for _ in range(3):
            To = To - (substitute(proper, To) - substitute(proper, Ti) - (O**2 - I**2) / 2) / O
        He = [[substitute(h[i][j], Te) for j in range(3)] for i in range(3)]
        Ho = [[substitute(h[i][j], To) for j in range(3)] for i in range(3)]
        hp = [[F.apply(lambda f: s.diff(f, time)) for F in row] for row in h]
        bracket = eng.madd(eng.mmul(h, hp), eng.mmul(hp, h), s=-1)
        rotation = [[F.apply(lambda f: s.integrate(f, (time, I, O)) / 8) for F in row] for row in bracket]
        R = eng.madd(eng.meye(), rotation)
        RT = eng.madd(eng.meye(), rotation, s=-1)
        M = eng.mmul(
            eng.mmul(
                eng.mmul(eng.mmul(RT, eng.mexp(eng.mscale(Ho, s.Rational(1, 2)))), eng.mexp(eng.mscale(He, -1))),
                eng.mexp(eng.mscale(Ho, s.Rational(1, 2))),
            ),
            R,
        )
        delta = sum((n[i] * n[j] * (M[i][j] - (1 if i == j else 0)) for i in range(3) for j in range(3)), eng.Field())
        temperature = 1 - delta / 2 + 3 * power(delta, 2) / 8 - 5 * power(delta, 3) / 16
        mean = temperature.apply(average_poly) - 1
        inverse = 1 - mean + power(mean, 2) - power(mean, 3)
        theta = temperature * inverse - 1
        return average_poly(s.Rational(15, 2) * n.dot(C * n) * theta.coeff([1, 2, 3]))
    finally:
        eng.Ctx.maxdeg, eng.Ctx.trunc = saved


def check_radial_ordering():
    u, v, w, t = s.symbols("u v w t", real=True)
    # All three directions equal: the 6 ordered labels reproduce -3 R v^2.
    original = -3 * s.integrate((s.exp(t) - 1) * ((s.exp(2 * t) - s.exp(2)) / 2) ** 2, (t, 0, 1))
    ordered = ordered_radial_transport(lambda a, b, c, U, V, W: s.exp(U + 2 * V + 2 * W), u, v, w, 0, 1).doit()
    zero(original - ordered)
    assert ordered != 0
    print("[ok] common fourth ray time evaluated exactly; six ordered sectors and sign/factor verified by exponential fixture")


def check_fourier_polarization():
    axis = s.Matrix([0, 0, 1])
    labels = s.symbols("label0:3")
    states = []
    for a in range(3):
        p = s.Matrix(s.symbols(f"p{a}_0:3"))
        G = s.Matrix([[a + 2, 1, 2], [1, a - 1, 1], [2, 1, 3 - a]]) / s.Integer(a + 3)
        hp = s.Matrix([[1, a, 1], [a, -2, 2], [1, 2, 1]]) / s.Integer(a + 5)
        X = s.Matrix(s.symbols(f"X{a}_0:3"))
        v = s.Matrix([*s.symbols(f"v{a}_0:2"), 0])
        R = s.Matrix(s.symbols(f"R{a}_0:3"))
        K = s.Matrix(s.symbols(f"K{a}_0:3"))
        driver = linear_direction_driver(axis, hp, [s.I * p[i] * G for i in range(3)])
        states.append(RayFirst(p, 1, G, hp, driver, X, v, R, K))

    def total(attr):
        return sum((labels[a] * getattr(A, attr) for a, A in enumerate(states)), s.zeros(*getattr(states[0], attr).shape))

    G, hp, X, v, R, K = (total(attr) for attr in ("G", "hp", "X", "v", "R", "K"))
    gradG = [sum((labels[a] * s.I * A.p[i] * A.G for a, A in enumerate(states)), s.zeros(3)) for i in range(3)]
    gradHp = [sum((labels[a] * s.I * A.p[i] * A.hp for a, A in enumerate(states)), s.zeros(3)) for i in range(3)]
    grad2G = [
        [sum((-labels[a] * A.p[i] * A.p[j] * A.G for a, A in enumerate(states)), s.zeros(3)) for j in range(3)] for i in range(3)
    ]
    grad2Hp = [
        [sum((-labels[a] * A.p[i] * A.p[j] * A.hp for a, A in enumerate(states)), s.zeros(3)) for j in range(3)] for i in range(3)
    ]
    old = reduced_cubic_bulk(axis, G, hp, gradG, gradient_linear_driver(axis, gradHp, grad2G), gradHp, grad2Hp, X, v, R, K)
    mixed = s.diff(old, *labels).subs(dict.fromkeys(labels, 0)) / 6
    radial = -sum(axis.dot(states[a].R) * states[b].v.dot(states[c].v) for a, (b, c) in enumerate(((1, 2), (0, 2), (0, 1))))
    zero(polarized_bulk(states, axis) + radial - mixed)
    # The cubic identity checks the R and K contractions of both quadratic
    # forces at arbitrary probe vectors and nine symbolic momentum components.
    print(
        "[ok] exact Fourier-specialized bilinear forcing and all cubic polarization factors versus the prior generic-jet formula"
    )


def check_composition_polarization():
    o = s.Symbol("observer", positive=True)
    labels = s.symbols("L0:3")

    def first(name):
        return tuple(s.Function(name + str(i))(o) for i in range(3))

    def pairs(name):
        return {pair: s.Function(name + "".join(map(str, pair)))(o) for pair in combinations(range(3), 2)}

    U1, F1, Dt, d1, D1 = (first(name) for name in ("U", "F", "T", "d", "D"))
    U2, F2, d2, D2 = (pairs(name) for name in ("V", "G", "e", "E"))
    U3 = s.Symbol("U3")
    second, third = compose_mixed_exponents(U1, U2, U3, F1, F2, Dt, d1, d2, D1, D2, o)

    def linear(values):
        return sum(labels[i] * values[i] for i in range(3))

    def quadratic(values):
        return 2 * sum(labels[i] * labels[j] * value for (i, j), value in values.items())

    jets = compose_first_only(
        linear(U1),
        quadratic(U2),
        6 * s.prod(labels) * U3,
        linear(F1),
        quadratic(F2),
        linear(Dt),
        linear(d1),
        quadratic(d2),
        linear(D1),
        quadratic(D2),
        linear(U1).diff(o),
        linear(U1).diff(o, 2),
        quadratic(U2).diff(o),
        linear(F1).diff(o),
        linear(d1).diff(o),
    )
    zero(s.diff(jets.third, *labels).subs(dict.fromkeys(labels, 0)) / 6 - third)
    for pair, value in second.items():
        zero(s.diff(jets.second, *(labels[i] for i in pair)).subs(dict.fromkeys(labels, 0)) / 2 - value)
    print("[ok] exact arbitrary-function endpoint/clock polarization, including both mixed observer-emission derivatives")


def check_noncoplanar_plan():
    from quadrupole_adjoint_expansion import Q

    ps = (
        s.Matrix([1, 0, 0]),
        s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0]),
        s.Matrix([s.Rational(1, 3), s.Rational(2, 3), s.Rational(2, 3)]),
    )
    # Exact off-shell polynomial profiles, not truncated radiation transfers.
    seeds = tuple(FourierFirst(p, 1, 0, time**2 * Q(p), time) for p in ps)
    data = first_only_exponent_triple(seeds, n, 1, 2, s.Rational(1, 2))
    plan = measurement_phase_plan(data)
    assert plan and max(len(term.limits) for term in plan) <= 3
    assert any(len(term.limits) == 3 for term in plan)
    for term in plan:
        assert not any(k.has(*n) for k in term.phase)
        assert len({limit[0] for limit in term.limits}) == len(term.limits)
    print(
        f"[ok] complete exact finite-phase noncoplanar triple plan: {len(plan)} sectors, at most three ray times; dense Fr expansion intentionally deferred"
    )


def check_double_dilation():
    p = s.Matrix(s.symbols("p0:3", real=True))
    A = s.Matrix([[1, 1, 0], [1, -2, 1], [0, 1, 1]]) / 7
    C = s.diag(2, -1, -1) / 3
    E, O, I = s.symbols("emission observer initial", positive=True)
    profile = s.Function("f")
    seeds = (
        FourierFirst(p, 1, s.Function("a")(time), profile(time) * A, time),
        FourierFirst(s.zeros(3, 1), 1, 0, s.zeros(3), time),
        FourierFirst(s.zeros(3, 1), 1, 0, s.zeros(3), time),
    )
    data = first_only_exponent_triple(seeds, n, E, O, I)
    plan = measurement_phase_plan(data)
    assert len(plan) == 1 and len(plan[0].limits) == 1
    u = plan[0].limits[0][0]
    hpnn = profile(u).diff(u) * n.dot(A * n)
    x = p.dot(n) * (O - u)
    zero(plan[0].polynomial - hpnn * (x * x - s.I * x) / 2)
    result = contract_direct_measurement(data, C)
    assert len(result.sectors) == 1 and not result.expression.has(*n)
    limit = result.sectors[0].limits[0]
    u = limit[0]
    alpha = s.Symbol("alpha", real=True)
    linear = -s.Rational(1, 2) * phase_average(
        profile(u).diff(u) * n.dot(A * n), n, alpha * p * (O - u), s.Rational(15, 2) * n.dot(C * n)
    )
    expected = (linear.diff(alpha, 2) + linear.diff(alpha)).subs(alpha, 1)
    got = result.expression.function
    for_value = got - expected
    zero(reduce_radial(for_value, (p * (O - u)).dot(p * (O - u))))
    print(
        "[ok] exact finite-wavelength double-dilation measurement identity at arbitrary momentum; two-time term collapsed to one"
    )


def check_homogeneous():
    E, O, I = s.S.One, s.Integer(2), s.Rational(1, 2)
    C = s.Matrix([[2, 1, -1], [1, 0, 2], [-1, 2, -2]]) / 3
    matrices = [
        s.Matrix([[1, 1, 0], [1, 0, 1], [0, 1, -1]]) / 7,
        s.Matrix([[2, 0, 1], [0, -1, 1], [1, 1, -1]]) / 11,
        s.Matrix([[1, -1, 1], [-1, -2, 0], [1, 0, 1]]) / 13,
    ]
    seeds = tuple(
        FourierFirst(s.zeros(3, 1), s.Rational(i + 1, 7), time ** (i + 2) / s.Integer(11 + i), time ** (i + 2) * A, time)
        for i, A in enumerate(matrices)
    )
    data = first_only_exponent_triple(seeds, n, E, O, I)
    # At zero phase the explicit remaining time integrals are polynomial and
    # can be evaluated before contraction for this independent exact fixture.
    got = contract_direct_measurement(
        type(data)(
            tuple(x.doit() for x in data.first),
            {k: x.doit() for k, x in data.second.items()},
            data.third.doit(),
            data.radial_transport.doit(),
            data.n,
        ),
        C,
    ).expression.doit()
    expected = canonical_homogeneous(seeds, C, E, O, I)
    zero(got - expected)
    assert expected != 0
    assert data.radial_transport == 0 or data.radial_transport.doit() == 0
    print(
        "[ok] full arbitrary-triple normalized measurement versus conserved photon energy, physical roots/initial clock and noncommuting FW frame"
    )


if __name__ == "__main__":
    check_radial_ordering()
    check_fourier_polarization()
    check_composition_polarization()
    check_double_dilation()
    check_noncoplanar_plan()
    check_homogeneous()
    print("DIRECT FIRST-ONLY CUBIC MEASUREMENT CHECKS PASSED; full Einstein cubic map and spectrum not evaluated")
