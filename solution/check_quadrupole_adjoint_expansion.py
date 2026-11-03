"""Exact gradient-coefficient checks; NOT a finite-wavelength identity.

The first-adjoint coefficients are tested against the full source implementation,
then against causal Green probes, boundary traces and the singular origin.
The first-adjoint cubic sector is only one part of the fourth-gradient map.
"""

from functools import lru_cache

import sympy as s
from quadrupole_adjoint_expansion import (
    Q,
    QuadrupoleFirstAdjoint,
    angular_A,
    angular_B,
    averaged_shell_first_correction,
    conformal_ricci_dual,
    first_cubic_correction,
    first_pair_correction,
    homogeneous_lapse_correction,
    initial_pair,
    leading_cubic_source,
    physical_pair_source,
)
from unified_adjoint_source import field_coefficients, inner, measurement_source, ricci_mixed, stf, stf_dual, sym_outer

n = s.Matrix(s.symbols("nx ny nz", real=True))
t, x, kap = s.symbols("t x kappa", positive=True)
e, o, ti = s.symbols("e o iota", positive=True)


@lru_cache(None)
def moment(powers):
    if any(a % 2 for a in powers):
        return s.S.Zero
    return s.prod(s.factorial2(a - 1) for a in powers) / s.factorial2(sum(powers) + 1)


def sky(expr):
    return s.expand(sum(c * moment(powers) for powers, c in s.Poly(s.expand(expr), *n).terms()))


def sky_matrix(A):
    return A.applyfunc(sky)


def zero(expr):
    if isinstance(expr, s.MatrixBase):
        for c in expr:
            zero(c)
    else:
        assert s.simplify(s.expand_log(s.expand(expr), force=True)) == 0, expr


def coefficient2(expr):
    if isinstance(expr, s.MatrixBase):
        return expr.applyfunc(coefficient2)
    return s.expand(s.diff(expr, kap, 2).subs(kap, 0) / 2)


def generic_tensor_checks():
    u = s.Matrix(s.symbols("u0:3"))
    v = s.Matrix(s.symbols("v0:3"))
    c = s.symbols("c0:5")
    C = s.Matrix([[c[0], c[2], c[3]], [c[2], c[1], c[4]], [c[3], c[4], -c[0] - c[1]]])
    w = s.Rational(15, 2) * n.dot(C * n)
    zero(sky_matrix(w * stf(n * n.T)) - C)
    zero(sky_matrix(w * u.dot(n) * v.dot(n) * stf(n * n.T)) - angular_A(C, u, v))
    zero(sky_matrix(w * v.dot(n) * stf(sym_outer(u, n))) - angular_B(C, u, v))
    actual = stf_dual(lambda j: inner(C, ricci_mixed(u, v, 2 * s.eye(3), j)))
    zero(actual - conformal_ricci_dual(u, v, C))
    print("[ok] all quadrupole moment and conformal curvature-dual formulas at arbitrary vectors and STF covector")


def source_checks(p, q, C):
    adj = QuadrupoleFirstAdjoint(p, q, C, e, o, ti)
    p2 = p.dot(p)
    P = s.eye(3) - n * n.T
    delta = t - o
    X = delta * n + s.I * kap * delta**2 * (p - 2 * p.dot(n) * n) / 2
    v = s.I * kap * delta * P * p + kap**2 * p.dot(n) * delta * (t - 2 * o) * P * p / 3
    K = kap**2 * (t * t - e * e) * P * Q(p) * n / 3
    w = s.Rational(15, 2) * n.dot(C * n)
    theta = kap**2 * adj.tau2 * n.dot(Q(p) * n)
    Ja = kap**2 * adj.tau2 * inner(C, Q(p))
    D = kap**2 * p2 * (o**4 + 3 * ti**4) / (36 * o)
    meas = measurement_source(
        kap * q,
        n,
        w / (4 * s.pi),
        X,
        v,
        s.zeros(3, 1),
        K,
        (w * theta - Ja) / (4 * s.pi),
        D,
        -(kap**2) * p2 * e**2 / 9,
        -(kap**2) * e * n.dot(Q(p) * n) / 3,
        -(kap**2) * o * inner(C, Q(p)) / 3,
        e,
        o,
        ti,
        t,
    )
    fullsky = lambda z: 4 * s.pi * sky(z)
    zero(coefficient2(meas.H))
    got = coefficient2(meas.H_prime).applyfunc(fullsky)
    zero(got - adj.measurement_tensor(t))
    Cp = inner(C, Q(p))
    zero(fullsky(coefficient2(meas.nu_emission)) + e**2 * Cp / 3)
    zero(coefficient2(meas.H_prime_emission).applyfunc(fullsky) + p2 * e**3 * C / 18)
    zero(coefficient2(meas.H_prime_observer).applyfunc(fullsky) - p2 * (o**4 + 3 * ti**4) * C / (72 * o))
    zero(coefficient2(meas.nu_clock_bulk) - t * Cp / 3)
    zero(coefficient2(meas.nu_clock_initial) - ti**2 * Cp / 3)
    W = s.Symbol("W")
    fld = field_coefficients(
        kap * p,
        kap * q,
        1,
        -(kap**2) * p2 * t**2 / 9,
        -2 * kap**2 * p2 * t / 9,
        -(kap**2) * t**2 * Q(p) / 3,
        -2 * kap**2 * t * Q(p) / 3,
        0,
        W * C,
        t,
    )
    zero(coefficient2(fld.H) - 2 * W * adj.rho)
    zero(coefficient2(fld.H_prime) + 4 * p2 * t * W * C / 9)
    zero(coefficient2(fld.nu) + 8 * W * (p + q).dot(C * q))
    zero(coefficient2(fld.nu_prime) + 4 * t * W * Cp / 3)
    # The normalization part of M is generally nonzero and cannot be deleted.
    norm = adj.tau2 * (angular_A(C, p, p) - p2 * C / 3) / 2
    assert any(s.expand(z) != 0 for z in norm)
    print("[ok] closed weak source matches complete source API for p=" + str(list(p)) + ", q=" + str(list(q)))


def reconstruction_checks(p, q, C):
    E, O, I = s.S.One, s.Integer(3), s.Rational(1, 2)
    adj = QuadrupoleFirstAdjoint(p, q, C, E, O, I)
    W = lambda z, early: z * z * (1 / E - 1 / O) / 2 if early else z * (O - z) / (2 * O)
    integ = lambda f, lo, hi: s.integrate(s.expand(f), (t, lo, hi))

    def field_integral(f, at):
        if at < E:
            return integ(f * W(t, True), at, E) + integ(f * W(t, False), E, O)
        return integ(f * W(t, False), at, O)

    Cp = inner(C, Q(p))
    D = (p + q).dot(C * q)
    p2 = p.dot(p)
    for at in (s.Rational(1, 4), s.Rational(3, 4), s.Integer(2)):
        early, before = at < E, at < I
        GH = at * (t - at) / t
        GHp = at**2 / t**2
        L = max(at, E)
        Hprobe = adj.measurement_tensor(t).applyfunc(lambda z, GHp=GHp, L=L: integ(z * GHp, L, O))
        Hprobe += -2 * field_integral(GH, at) * adj.rho + 4 * p2 * field_integral(t * GHp, at) * C / 9
        if early:
            Hprobe -= p2 * E**3 * (at**2 / E**2) * C / 18
        Hprobe += p2 * (O**4 + 3 * I**4) * (at**2 / O**2) * C / (72 * O)
        zero(Hprobe - adj.spatial(at, early))
        GN = (t**3 - at**3) / (9 * t * at)
        GNp = s.diff(GN, t)
        Nprobe = 8 * D * field_integral(GN, at) + 4 * Cp * field_integral(t * GNp, at) / 3
        Nprobe += Cp * integ(t * GN, max(at, I), O) / 3
        if early:
            Nprobe -= E**2 * Cp * GN.subs(t, E) / 3
        if before:
            Nprobe += I**2 * Cp * GN.subs(t, I) / 3
        zero(Nprobe - adj.lapse(at, early, before))
    symbolic = QuadrupoleFirstAdjoint(p, q, C, e, o, ti)
    zero(symbolic.spatial(o, False) - p2 * (o**4 + 3 * ti**4) * C / (72 * o))
    zero(symbolic.spatial(e, True) - symbolic.spatial(e, False) + p2 * e**3 * C / 18)
    zero(s.limit(9 * x * symbolic.lapse(x, True, True), x, 0, dir="+") - symbolic.origin())
    # Independent exact quadrature of the explicit logarithmic adjoint.
    integrated = adj.spatial(t, True).applyfunc(lambda z: integ(z, 0, E))
    integrated += adj.spatial(t, False).applyfunc(lambda z: integ(z, E, O))
    zero(integrated - adj.integrated_spatial())
    assert not any(z.has(s.log) for z in adj.integrated_spatial())
    print("[ok] both causal Green reconstructions, observer trace, emission jump, singular origin and exact log cancellation")


def symbolic_time_integrals():
    """Prove log cancellation for arbitrary times, not just chosen values."""
    adj = QuadrupoleFirstAdjoint(s.Matrix([1, 0, 0]), s.zeros(3, 1), s.diag(2, -1, -1) / 3, e, o, ti)
    F = lambda z: o**2 / z - z - 2 * o * s.log(o / z)
    A = o**4 / 12 - o**2 * e**2 / 2 + 2 * o * e**3 / 3 - e**4 / 4
    zero(s.integrate(x * x * F(e), (x, 0, e)) + s.integrate(x * x * F(x), (x, e, o)) - A / 3)
    zero(s.integrate(x * x * (1 / e - 1 / o), (x, 0, e)) + s.integrate(x * x * (1 / x - 1 / o), (x, e, o)) - (o**2 - e**2) / 6)
    J2 = (o**4 - e**4) / 40
    for kind in (0, 1):
        expr = lambda early, kind=kind: (
            x * (adj.moment(0, x, early) - x * adj.moment(-1, x, early)) if kind == 0 else x * x * adj.moment(-1, x, early)
        )
        val = s.integrate(expr(True), (x, 0, e)) + s.integrate(expr(False), (x, e, o))
        zero(val - J2 / (6 if kind == 0 else 3))
    print("[ok] complete logarithm cancellation and field-moment integrals at arbitrary symbolic e,o")


def physical_leading_comparison():
    a = s.Matrix(s.symbols("a0:3"))
    b = s.Matrix(s.symbols("b0:3"))
    c = s.Matrix(s.symbols("c0:3"))
    z = s.symbols("z0:3")
    vectors = (a, b, c)
    Z = sum(z)
    grad = sum((s.I * zj * k for zj, k in zip(z, vectors, strict=True)), s.zeros(3, 1))
    Hess = -sum((zj * k * k.T for zj, k in zip(z, vectors, strict=True)), s.zeros(3))
    # Independent local-curvature observable; no adjoint or field residual.
    curvature = s.exp(-2 * Z) * (Hess - grad * grad.T)
    independent = -((o**2 - e**2) / 6) * stf(curvature.applyfunc(lambda f: s.diff(f, *z).subs(dict.fromkeys(z, 0))))
    assembled = ((o**2 - e**2) / 12) * leading_cubic_source(a, b, c)
    zero(independent - assembled)
    print("[ok] full arbitrary-triple physical cubic quadrupole at TWO gradients agrees with independent curvature map")


def shell_checks():
    C = Q(s.Matrix([1, 0, 0]))
    k = s.Matrix([1, 0, 0])
    q = s.Matrix(s.symbols("qx qy qz", real=True))

    # Full isotropic q average, with |q|=1. Polynomial moments enforce this
    # without any square-root/angular parametrization or momentum division.
    def qaverage(z):
        return s.expand(sum(cf * moment(powers) for powers, cf in s.Poly(s.expand(z), *q).terms()))

    hom = QuadrupoleFirstAdjoint(k, s.zeros(3, 1), C, e, o, ti)
    zero(qaverage(inner(hom.integrated_spatial(), physical_pair_source(q, -q))))
    homlapse = qaverage(-hom.origin() * initial_pair(q, -q))
    zero(homlapse + (3 * o**4 + 5 * ti**4 - 8 * e**4) * inner(C, Q(k)) / 36)
    assert s.expand(homlapse) != 0
    full = qaverage(first_cubic_correction(k, q, -q, C, e, o, ti))
    zero(full - averaged_shell_first_correction(C, k, e, o, ti))
    zero(homlapse - homogeneous_lapse_correction(C, k, e, o, ti))
    print("[result] shell-averaged first-adjoint kappa^4 block:", s.factor(full))
    print("[result] homogeneous-lapse partition:", s.factor(homlapse))
    print("[ok] homogeneous tensor averages to zero; homogeneous lapse and physical initial clock remain")
    # No extra factor 3 multiplying an individual partition.
    direct = sum(first_pair_correction(p, u, v, C, e, o, ti) for p, u, v in ((k, q, -q), (q, k, -q), (-q, k, q)))
    zero(direct - first_cubic_correction(k, q, -q, C, e, o, ti))


if __name__ == "__main__":
    generic_tensor_checks()
    C = s.Matrix([[2, 1, -1], [1, 0, 2], [-1, 2, -2]]) / 3
    cases = [
        (s.Matrix([1, 0, 0]), s.zeros(3, 1)),
        (s.Matrix([0, 1, 0]), s.Matrix([1, -1, 0])),
        (s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0]), s.Matrix([s.Rational(2, 5), -s.Rational(4, 5), 0])),
    ]
    for p, q in cases:
        source_checks(p, q, C)
    reconstruction_checks(*cases[-1], C)
    symbolic_time_integrals()
    physical_leading_comparison()
    shell_checks()
    print("QUADRUPOLE ADJOINT GRADIENT COEFFICIENTS VERIFIED; full finite-wavelength spectrum NOT evaluated")
