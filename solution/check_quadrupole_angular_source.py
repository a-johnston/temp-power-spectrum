"""Exact source/Volterra identities and gradient regression; finite-phase
quadrature is separately labeled numerical validation, not a physics solution.
"""

import time
from functools import lru_cache

import numpy as np
import sympy as s
from check_entire_sphere import numerical, quadrature, reduce_radial, zero
from entire_sphere import RadialBessel, phase_moment, polynomial_rank, sphere_reduce
from quadrupole_adjoint_expansion import Q, QuadrupoleFirstAdjoint
from quadrupole_angular_source import QuadrupoleAngularSource, ray_polynomials
from stationary_cubic import direct_field_residuals, linear_direction_driver
from unified_adjoint_source import (
    FieldVariationCoefficients,
    FirstFieldAdjointSource,
    frame_terminal_covector,
    inner,
    measurement_source,
)

n = s.Matrix(s.symbols("n0:3", real=True))
t, u = s.symbols("t u", positive=True)
e, o, ti = s.symbols("e o iota", positive=True)


def zero_sphere(expr):
    if isinstance(expr, s.MatrixBase):
        for z in expr:
            zero_sphere(z)
    else:
        zero(sphere_reduce(expr, n))


def sky_polynomial(expr):
    from check_quadrupole_adjoint_expansion import moment

    return s.expand(sum(c * moment(powers) for powers, c in s.Poly(s.expand(expr), *n).terms()))


def sym_matrix(prefix):
    a = s.symbols(prefix + "0:5")
    return s.Matrix([[a[0], a[2], a[3]], [a[2], a[1], a[4]], [a[3], a[4], -a[0] - a[1]]])


def polynomial_identities():
    p = s.Matrix(s.symbols("p0:3"))
    q = s.Matrix(s.symbols("q0:3"))
    C = s.diag(2, -1, -1) / 3
    h, hp = sym_matrix("h"), sym_matrix("d")
    z = s.Symbol("z")
    obj = QuadrupoleAngularSource(p, q, C, e, o, ti, n)
    atom = ray_polynomials(p, n, z, h, hp, t, u)
    G = 2 * z * s.eye(3) + h
    P = s.eye(3) - n * n.T
    zero_sphere(atom.v_upper - linear_direction_driver(n, hp, [s.I * p[i] * G for i in range(3)]))
    zero_sphere(atom.X_upper.diff(t) + atom.v_upper)
    zero_sphere(atom.X_upper.subs(u, t) + G * n / 2)
    zero_sphere(atom.K_lower.diff(t) + atom.R_lower)
    zero_sphere(atom.K_lower.subs(u, t) + P * hp * n)
    B0, Bup, Blo = obj.polynomial_blocks(z, h, hp, t, u)
    phase = s.exp(s.I * q.dot(n) * (o - t))

    def old(X, v, R, K):
        return measurement_source(q, n, 1, X, v, R, K, 0, 0, 0, 0, 0, e, o, ti, t)

    upper = old(atom.X_upper, -atom.v_upper, s.zeros(3, 1), s.zeros(3, 1))
    lower = old(s.zeros(3, 1), s.zeros(3, 1), atom.R_lower, atom.K_lower)
    zero_sphere(upper.H_prime / phase - Bup)
    zero_sphere(lower.H_prime / phase - Blo)
    zero_sphere(lower.H / phase - B0)
    ranks = [polynomial_rank(obj.weight * b, n) for b in (B0, Bup, Blo)]
    assert max(ranks) == 7 and all(rank <= 7 for rank in ranks), ranks
    print("[ok] arbitrary-momentum Volterra signs/boundaries and exact old/new ray-source polynomial identities")
    print("[ok] unit-sphere-reduced weighted angular ranks:", ranks, "(F0,...,F7 suffice)")


def frame_ward_identity():
    p = s.Matrix(s.symbols("p0:3"))
    C = sym_matrix("C")
    hp = sym_matrix("H")
    obj = QuadrupoleAngularSource(p, s.zeros(3, 1), C, 1, 2, s.Rational(1, 2), n)
    P = s.eye(3) - n * n.T
    grad = -P * (hp * n + s.I * p * n.dot(hp * n) / 2)
    raw = frame_terminal_covector(obj.weight, grad, n)
    observer_kernel = ray_polynomials(p, n, 0, s.zeros(3), hp, s.Integer(2), s.S.One).K_lower
    zero(raw - frame_terminal_covector(obj.weight, observer_kernel, n))

    def average_unreduced(poly):
        return s.Add(*(cf * phase_moment(powers, p) for powers, cf in s.Poly(s.expand(poly), *n).terms()))

    frame = raw.applyfunc(average_unreduced)
    S = (-s.Rational(15, 4) * n.dot(hp * n) * obj.N).applyfunc(average_unreduced)
    zero(frame - (S * C - C * S))
    # Compare the unit-sphere-reduced production representation as well.
    difference = obj.frame_density(hp, s.S.One) - frame
    for value in difference:
        zero(reduce_radial(value, p.dot(p)))
    print("[ok] exact finite-wavelength rotation Ward identity K_a=[linear quadrupole,C] at generic momentum and STF matrices")


def finite_phase_validation():
    C = s.Matrix([[2, 1, -1], [1, 0, 2], [-1, 2, -2]]) / 3
    h = s.Matrix([[1, 2, 0], [2, -2, 1], [0, 1, 1]]) / 5
    hp = s.Matrix([[2, -1, 1], [-1, 1, 2], [1, 2, -3]]) / 7
    p = s.Matrix([1, s.Rational(-1, 2), s.Rational(1, 3)])
    ordinary = s.Matrix([s.Rational(2, 3), 1, s.Rational(-3, 4)])
    cases = [
        (ordinary, s.Rational(5, 4), s.Rational(7, 4)),
        (-p, s.Rational(5, 4), s.Rational(5, 4)),
        (-p + s.Matrix([s.Rational(1, 10**12), 0, 0]), s.Rational(5, 4), s.Rational(5, 4)),
        (s.zeros(3, 1), s.Rational(5, 4), s.Rational(7, 4)),
        (ordinary, s.Rational(7, 4), s.Rational(5, 4)),
    ]
    for case, (q, T, U) in enumerate(cases):
        E, O, I = s.S.One, s.Rational(5, 2), s.Rational(1, 2)
        obj = QuadrupoleAngularSource(p, q, C, E, O, I, n)
        atom = ray_polynomials(p, n, 1, h, hp, T, U)
        blocks = obj.blocks(1, h, hp, T, U)
        Kp, Kq, _ = obj.phases(T, U)
        phi = s.exp(s.I * Kp.dot(n))

        def old(X, v, R, K, lam=0, q=q, obj=obj, E=E, O=O, I=I, T=T):
            return measurement_source(q, n, obj.weight / (4 * s.pi), X, v, R, K, lam, 0, 0, 0, 0, E, O, I, T)

        upper = old(atom.X_upper * phi, -atom.v_upper * phi, s.zeros(3, 1), s.zeros(3, 1))
        lower = old(s.zeros(3, 1), s.zeros(3, 1), atom.R_lower * phi, atom.K_lower * phi)
        tests = [
            (blocks.Hprime_upper, 4 * s.pi * upper.H_prime),
            (blocks.Hprime_lower, 4 * s.pi * lower.H_prime),
            (blocks.H_lower, 4 * s.pi * lower.H),
        ]
        Ua = n.dot(hp * n) * phi / 2
        mean = quadrature(Ua, n, s.zeros(3, 1), 56)
        Ja = -quadrature(obj.weight * Ua, n, s.zeros(3, 1), 56)
        lam = (obj.weight * (-Ua + s.sympify(mean)) - s.sympify(Ja)) / (4 * s.pi)
        norm = old(s.zeros(3, 1), s.zeros(3, 1), s.zeros(3, 1), s.zeros(3, 1), lam)
        tests.append((blocks.Hprime_normalization, 4 * s.pi * norm.H_prime))
        for closed, raw in tests:
            got = numerical(closed)
            direct = quadrature(raw, n, s.zeros(3, 1), 40)
            refined = quadrature(raw, n, s.zeros(3, 1), 56)
            assert np.max(np.abs(direct - refined)) < 2e-12
            assert np.max(np.abs(got - refined)) < 2e-12, (case, np.max(np.abs(got - refined)))
        if case == 0:
            # Both separate-mean terms are nonzero at finite wavelength.
            hpnn = n.dot(hp * n)
            first = obj.plain(hpnn, Kp) * obj.weighted(obj.N, Kq) / 4
            second = obj.weighted(hpnn, Kp) * obj.plain(obj.N, Kq) / 4
            assert np.max(np.abs(numerical(first))) > 1e-6
            assert np.max(np.abs(numerical(second))) > 1e-6
    print("[ok] complete source kernels versus independent finite-phase sphere quadrature, including K=0, near K=0 and q=0")
    print("      numerical validation only; both random-mean normalization products are nonzero")


def normalization_identity():
    p = s.Matrix(s.symbols("p0:3"))
    q = s.Matrix(s.symbols("q0:3"))
    C = s.diag(2, -1, -1) / 3
    hp = s.Matrix([[2, -1, 1], [-1, 1, 2], [1, 2, -3]]) / 7
    obj = QuadrupoleAngularSource(p, q, C, 1, 3, s.Rational(1, 2), n)
    T, U = s.S.One, s.S.One
    Kp, Kq, K = obj.phases(T, U)
    label = s.Symbol("first_phase")
    Ua = n.dot(hp * n) / 2
    mean = obj.plain(Ua, Kp)
    Ja = -obj.weighted(Ua, Kp)
    old = measurement_source(
        q,
        n,
        obj.weight / (4 * s.pi),
        s.zeros(3, 1),
        s.zeros(3, 1),
        s.zeros(3, 1),
        s.zeros(3, 1),
        (obj.weight * (-Ua * label + mean) - Ja) / (4 * s.pi),
        0,
        0,
        0,
        0,
        1,
        3,
        s.Rational(1, 2),
        T,
    )
    raw = 4 * s.pi * old.H_prime / s.exp(s.I * Kq.dot(n))
    coupled = raw.applyfunc(lambda z: s.expand(z).coeff(label, 1))
    separated = raw.applyfunc(lambda z: s.expand(z).coeff(label, 0))
    closed = obj.plain(coupled, K) + obj.plain(separated, Kq)
    zero(closed - obj.normalization(hp, T, U))
    print("[ok] exact finite-phase random-mean normalization identity at six independent momentum components")


def gradient_regression():
    kap = s.Symbol("kappa")
    p = s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0])
    q = s.Matrix([s.Rational(2, 5), -s.Rational(4, 5), 0])
    C = s.Matrix([[2, 1, -1], [1, 0, 2], [-1, 2, -2]]) / 3
    obj = QuadrupoleAngularSource(kap * p, kap * q, C, e, o, ti, n)
    h = -(kap**2) * u**2 * Q(p) / 3
    a = -(kap**2) * p.dot(p) * u**2 / 9
    source = obj.assemble(1, h, a, t, u)
    adj = QuadrupoleFirstAdjoint(p, q, C, e, o, ti)

    @lru_cache(None)
    def jet(expr):
        # Coefficient extraction commutes with the specified finite integrals.
        # Truncate at each algebra operation instead of differentiating giant
        # nested Integral expressions or expanding unneeded high-rank powers.
        if not expr.has(kap):
            return (expr, s.S.Zero, s.S.Zero)
        if expr == kap:
            return (s.S.Zero, s.S.One, s.S.Zero)
        if isinstance(expr, RadialBessel):
            r, z = expr.args
            zjet = jet(z)
            assert zjet[0] == zjet[1] == 0
            return (RadialBessel(r, 0), s.S.Zero, -RadialBessel(r + 1, 0) * zjet[2] / 2)
        if isinstance(expr, s.Integral):
            assert not any(limit.has(kap) for limit in expr.limits)
            return tuple(s.Integral(c, *expr.limits) if c != 0 else s.S.Zero for c in jet(expr.function))
        if expr.is_Add:
            terms = [jet(a) for a in expr.args]
            return tuple(s.Add(*(a[k] for a in terms)) for k in range(3))

        def multiply(a, b):
            return tuple(s.Add(*(a[j] * b[k - j] for j in range(k + 1))) for k in range(3))

        if expr.is_Mul:
            result = (s.S.One, s.S.Zero, s.S.Zero)
            for a in expr.args:
                result = multiply(result, jet(a))
            return result
        if expr.is_Pow and expr.exp.is_Integer and expr.exp >= 0:
            result = (s.S.One, s.S.Zero, s.S.Zero)
            for _ in range(int(expr.exp)):
                result = multiply(result, jet(expr.base))
            return result
        raise AssertionError(("unsupported gradient dependency", expr))

    coeff = lambda expr: s.expand(jet(expr)[2].doit())
    M = source.coefficients
    zero(M.H.applyfunc(coeff))
    zero(M.H_prime.applyfunc(coeff) - adj.measurement_tensor(t))
    Cp = inner(C, Q(p))
    p2 = p.dot(p)
    zero(coeff(M.nu_emission) + e**2 * Cp / 3)
    zero(M.H_prime_emission.applyfunc(coeff) + p2 * e**3 * C / 18)
    zero(M.H_prime_observer.applyfunc(coeff) - p2 * (o**4 + 3 * ti**4) * C / (72 * o))
    zero(coeff(M.nu_clock_bulk) - t * Cp / 3)
    zero(coeff(M.nu_clock_initial) - ti**2 * Cp / 3)
    assert not n.free_symbols & set().union(
        *(
            z.free_symbols
            for z in (M.H, M.H_prime, M.nu_emission, M.H_prime_emission, M.H_prime_observer, M.nu_clock_bulk, M.nu_clock_initial)
        )
    )
    print("[ok] exact finite-phase assembled source recovers every previous leading-gradient coefficient at arbitrary times")


def constant_seed_dilation():
    q = s.Matrix(s.symbols("q0:3"))
    C = s.diag(2, -1, -1) / 3
    obj = QuadrupoleAngularSource(s.zeros(3, 1), q, C, e, o, ti, n)
    result = obj.assemble(1, s.zeros(3), 0, t, u).coefficients
    B = obj.background_tensor(q, t)
    derivative = sum((q[i] * B.diff(q[i]) for i in range(3)), s.zeros(3))
    difference = result.H_prime.applyfunc(lambda z: z.doit()) + derivative / 2
    for value in difference:
        zero(reduce_radial(value, (q * (o - t)).dot(q * (o - t))))
    zero(result.H)
    assert all(z == 0 for z in (result.nu_emission, result.nu_clock_bulk, result.nu_clock_initial))
    zero(result.H_prime_emission)
    zero(result.H_prime_observer)
    print("[ok] exact finite-wavelength p=0 measurement dilation identity, including the observer ray boundary")


def homogeneous_pair_covariance():
    import engine as eng

    alpha = s.Matrix(s.symbols("alpha0:3", real=True))
    radius = s.Symbol("radius", positive=True)
    a = s.Function("a")(eng.t)
    b = s.Function("b")(eng.t)
    saved = eng.Ctx.maxdeg, eng.Ctx.momenta, eng.Ctx.trunc
    try:
        eng.Ctx.maxdeg = 3
        eng.Ctx.momenta = {1: [radius * alpha[i] for i in range(3)], 2: [-radius * alpha[i] for i in range(3)]}
        eng.Ctx.trunc = s.expand
        first = lambda value: eng.Field({frozenset([1]): value, frozenset([2]): value})
        H = [[first(b * Q(alpha)[i, j]) for j in range(3)] for i in range(3)]
        sourceN, sourceH, initial = direct_field_residuals(first(a), first(s.S.One), H)
        fixed = {alpha[0]: 1, alpha[1]: 0, alpha[2]: 0}
        tensor = s.Matrix(3, 3, lambda i, j: sourceH[i][j].coeff([1, 2]))
        radial = s.expand(3 * tensor[0, 0].subs(fixed) / 2)
        for value in tensor - radial * Q(alpha):
            zero(sphere_reduce(value, alpha))
        scalar = sourceN.coeff([1, 2])
        zero(sphere_reduce(scalar - scalar.subs(fixed), alpha))
        assert s.expand(scalar.subs(fixed)) != 0
        zero(sphere_reduce(initial.coeff([1, 2]) - s.Rational(5, 9) * radius**2, alpha))
    finally:
        eng.Ctx.maxdeg, eng.Ctx.momenta, eng.Ctx.trunc = saved
    print("[ok] exact homogeneous Einstein pair source is scalar times STF(alpha alpha); lapse/initial sources remain nonzero")


def homogeneous_assembly():
    E, O, I = s.S.One, s.Integer(2), s.Rational(1, 2)
    C = s.Matrix([[2, 1, -1], [1, 0, 2], [-1, 2, -2]]) / 3
    A = s.Matrix([[1, 1, 0], [1, 0, 1], [0, 1, -1]]) / 7
    B = s.Matrix([[2, 0, 1], [0, -1, 1], [1, 1, -1]]) / 11
    h = u**2 * A + u**3 * B
    a = u**2 / 7 - u**3 / 11
    obj = QuadrupoleAngularSource(s.zeros(3, 1), s.zeros(3, 1), C, E, O, I, n)
    result = obj.assemble(s.Rational(1, 3), h, a, t, u)
    Z = FieldVariationCoefficients(0, 0, s.zeros(3), s.zeros(3))
    combined = result.with_field(Z)
    arbitrary = s.Symbol("already_contracted")
    assert combined.sky(arbitrary) == arbitrary
    P = s.eye(3) - n * n.T
    ht = h.subs(u, t)
    v = P * (ht - h.subs(u, O)) * n / 2
    X = ((u - t) * P * h.diff(u) * n / 2 - (s.Rational(2, 3) * s.eye(3) + h) * n / 2).applyfunc(
        lambda z: s.integrate(z, (u, t, O))
    )
    K = -P * (ht - h.subs(u, E)) * n
    theta = -n.dot((h.subs(u, O) - h.subs(u, E)) * n) / 2
    Ja = -inner(C, h.subs(u, O) - h.subs(u, E)) / 2
    D = -(I**2 * a.subs(u, I) + s.integrate(u * a, (u, I, O))) / O
    Ae = inner(C, h.diff(u).subs(u, O)) / 2
    old = measurement_source(
        s.zeros(3, 1),
        n,
        obj.weight / (4 * s.pi),
        X,
        v,
        s.zeros(3, 1),
        K,
        (obj.weight * theta - Ja) / (4 * s.pi),
        D,
        a.subs(u, E),
        n.dot(h.diff(u).subs(u, E) * n) / 2,
        Ae,
        E,
        O,
        I,
        t,
    )
    old_source = FirstFieldAdjointSource(Z, old, t, E, O, I, lambda expr: 4 * s.pi * sky_polynomial(expr))
    probe = s.Matrix([[2, 1, 0], [1, -1, 1], [0, 1, -1]]) / 5
    got = combined.action(t**2 + t**3 / 3, t**2 * probe).doit()
    expected = old_source.action(t**2 + t**3 / 3, t**2 * probe).doit()
    zero(got - expected)
    assert got != 0 and not got.has(s.Integral)
    nonzero_field = FieldVariationCoefficients(t / 13, t**2 / 17, t * A, (t + 1) * B)
    full = result.with_field(nonzero_field).action(t**2 + t**3 / 3, t**2 * probe).doit()
    field_action = s.integrate(nonzero_field.action(t**2 + t**3 / 3, 2 * t + t**2, t**2 * probe, 2 * t * probe), (t, 0, O))
    zero(full - got + field_action)
    Delta = h.subs(u, O) - h.subs(u, E)
    zero(obj.frame_covector(h, u).doit() + (Delta * C - C * Delta) / 2)
    print("[ok] complete clock/emission/observer action and no-double-sky adapter in exact homogeneous off-shell fixture")


if __name__ == "__main__":
    for check in (
        polynomial_identities,
        frame_ward_identity,
        normalization_identity,
        finite_phase_validation,
        gradient_regression,
        constant_seed_dilation,
        homogeneous_pair_covariance,
        homogeneous_assembly,
    ):
        started = time.monotonic()
        check()
        print(f"      {check.__name__}: {time.monotonic() - started:.2f}s", flush=True)
    print("EXACT QUADRUPOLE ANGULAR SOURCE CLOSURE VERIFIED; Green/time integrals and full spectrum remain unevaluated")
