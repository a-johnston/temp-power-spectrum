"""Exact finite-wavelength Green certificates, plus independent quadrature.

These are operator identities, not an evaluated nonlinear spectrum. Rational
momentum-zero duality fixture is exact; finite-phase quadrature is numerical.
"""

import sympy as s
from quadrupole_angular_source import ContractedMeasurement
from regular_field_green import (
    RegularFieldGreen,
    TimeTail,
    WeakGreenAdjoint,
    clock_origin_primitive,
    clock_origin_weight,
    initial_lapse,
    physical_first_clock,
    physical_first_fields,
    scalar_convolutions,
    spatial_operator,
)
from unified_adjoint_source import FieldVariationCoefficients, MeasurementSource, inner


def reduce_tails(expr):
    """Prove at nonzero parameters using a*T_(r+2)=d^r/r!-T_r.

    All asserted kernels are entire, so polynomial/analytic continuation also
    proves the zero parameter. Zero is additionally tested without this reducer.
    """
    if isinstance(expr, s.MatrixBase):
        return expr.applyfunc(reduce_tails)
    while True:
        atoms = [T for T in expr.atoms(TimeTail) if T.args[0] >= 2]
        if not atoms:
            return s.cancel(s.expand(expr))
        expr = expr.xreplace(
            {
                T: (T.args[2] ** (T.args[0] - 2) / s.factorial(T.args[0] - 2) - TimeTail(T.args[0] - 2, T.args[1], T.args[2]))
                / T.args[1]
                for T in atoms
            }
        )


def zero(expr):
    if isinstance(expr, s.MatrixBase):
        for value in expr:
            zero(value)
    else:
        assert reduce_tails(expr) == 0, expr


def check_entire_and_operators():
    t, r, a, z = s.symbols("t r a z", positive=True)
    for rank in range(6):
        zero(TimeTail(rank, a, t) + a * TimeTail(rank + 2, a, t) - t**rank / s.factorial(rank))
        assert TimeTail(rank, 0, t) == t**rank / s.factorial(rank)
        assert TimeTail(rank, a, 0) == (1 if rank == 0 else 0)
    G = RegularFieldGreen(s.Matrix([s.sqrt(3 * a), 0, 0]), t, r)
    LN = lambda f: 3 * f.diff(t, 2) + (3 * a - 6 / t**2) * f
    LS = lambda f: f.diff(t, 2) + 2 * f.diff(t) / t - a * f
    zero(LN(G.lapse))
    zero(G.lapse.subs(t, r))
    zero(G.lapse.diff(t).subs(t, r) - s.Rational(1, 3))
    zero(LS(G.coupled_scalar) - G.lapse)
    for derivative in range(3):
        zero(G.coupled_scalar.diff(t, derivative).subs(t, r))
    m, mp, j, jp = G.origin_probe()
    zero(LN(m))
    zero(j.diff(t, 2) + 2 * jp / t + spatial_operator(G.p, j) + 4 * G.Q * m)
    zero(s.cancel(9 * r * G.lapse).subs(r, 0) - m)
    zero(-4 * G.Q * s.cancel(9 * r * G.coupled_scalar).subs(r, 0) - j)
    zero(G.clock_primitive(t).diff(t) - t * G.lapse)
    zero(clock_origin_primitive(a, t).diff(t) - t * initial_lapse(a, t))
    zero(s.cancel(9 * r * G.clock_primitive(t)).subs(r, 0) - clock_origin_primitive(a, t))
    # Cartesian operator polynomial, independently verified at generic momentum.
    p = s.Matrix(s.symbols("p0:3", real=True))
    b = s.symbols("b0:5")
    B = s.Matrix([[b[0], b[2], b[3]], [b[2], b[1], b[4]], [b[3], b[4], -b[0] - b[1]]])
    BB = spatial_operator(p, B)
    BBB = spatial_operator(p, BB)
    zero(spatial_operator(p, BBB) - 2 * p.dot(p) * BBB / 3 - p.dot(p) ** 2 * BB / 3)
    # Its two nonzero eigenvalue evaluations are EXACT, not finite series.
    H = RegularFieldGreen(s.Matrix([s.sqrt(z), 0, 0]), t, r)
    g0, g1, g2 = H.spatial_coefficients
    zero(g0 + z * g1 + z * z * g2 - r * TimeTail(1, z, t - r) / t)
    zero(g0 - z * g1 / 3 + z * z * g2 / 9 - r * TimeTail(1, -z / 3, t - r) / t)
    zero(g0.subs(t, r))
    zero(g0.diff(t).subs(t, r) - 1)
    for g in (g1, g2):
        zero(g.subs(t, r))
        zero(g.diff(t).subs(t, r))
    # Physical source+initial matching: the hyperbolic family cancels exactly.
    nu, h = physical_first_fields(G.p, t)
    D = scalar_convolutions(a, t)[2]
    direct = -2 * G.Q * TimeTail(3, -a, t) / t
    zero(direct + (-a / 3) * (-12 * G.Q * D / t) - h)
    zero(LN(nu))
    zero(h.diff(t, 2) + 2 * h.diff(t) / t + spatial_operator(G.p, h) + 4 * G.Q * nu + 2 * G.Q)
    zero(nu - (TimeTail(0, a, t) - TimeTail(1, a, t) / t))
    zero(physical_first_clock(G.p, r, t) - a * clock_origin_weight(a, r, t) / (3 * t))
    # Exact zero momentum: retain the singular 1/r and the clock origin row.
    G0 = RegularFieldGreen(s.zeros(3, 1), t, r)
    zero(G0.lapse - (t**3 - r**3) / (9 * t * r))
    zero(G0.spatial(B) - r * (t - r) * B / t)
    zero(G0.coupled_scalar - (r * (t - r) ** 3 / (18 * t) + (t - r) ** 4 / (36 * t) + (t - r) ** 5 / (180 * r * t)))
    zero(clock_origin_weight(0, r, t) - (t**4 + 3 * r**4) / 4)
    print(
        "[ok] entire retarded lapse/Cartesian spatial/coupled kernels: exact equations, jumps, zero momenta and regular-origin limits"
    )
    print(
        "[ok] lapse-spatial and proper-clock Green convolutions EVALUATED; physical first scalar hyperbolic terms cancel with initial matching"
    )


def check_weak_transpose():
    t, r = s.symbols("t r", positive=True)
    o, e, i = s.Integer(4), s.Integer(2), s.Integer(1)
    N0, N1 = t, 1 + t
    H0, H1 = 1 + t, 2 - t
    bnu, bhe, bho, A = s.symbols("bnu bhe bho A")
    B = s.diag(1, -1, 0)
    norm = inner(B, B)
    F = FieldVariationCoefficients(-N0, -N1, -H0 * B / norm, -H1 * B / norm)
    M = MeasurementSource(s.zeros(3), s.zeros(3), bnu, bhe * B / norm, bho * B / norm, -A * t / o, -A * i**2 / o)
    measure = ContractedMeasurement(M, t, e, o, i)
    adj = WeakGreenAdjoint(s.zeros(3, 1), F, measure)
    cov = adj.at(r)
    nu = t**2 + 2 * t**4
    h = 3 * t**2 - t**4
    fnu = 3 * nu.diff(t, 2) - 6 * nu / t**2
    fh = h.diff(t, 2) + 2 * h.diff(t) / t
    # All integrals here are elementary. This is an exact weak duality test.
    dual = s.integrate(s.simplify(cov.nu.doit()) * fnu.subs(t, r), (r, 0, o)) + cov.origin.doit()
    direct = s.integrate(N0 * nu + N1 * nu.diff(t), (t, 0, o)) + bnu * nu.subs(t, e)
    direct -= A * (s.integrate(t * nu, (t, i, o)) + i**2 * nu.subs(t, i)) / o
    zero(dual - direct)
    lh = inner(cov.H, B).doit()
    dual = s.integrate(lh * fh.subs(t, r), (r, 0, o))
    direct = s.integrate(H0 * h + H1 * h.diff(t), (t, 0, o)) + bhe * h.diff(t).subs(t, e) + bho * h.diff(t).subs(t, o)
    zero(dual - direct)
    zero(inner(adj.spatial(o), B) - bho)
    # SymPy limit(Piecewise(..., r<e), r,e,dir='-') can incorrectly
    # select the value AT the event. Evaluate the two explicit branches.
    branches = s.piecewise_fold(lh).args
    zero(branches[0].expr.subs(r, e) - branches[1].expr.subs(r, e) - bhe)
    nu_low = s.piecewise_fold(cov.nu.doit().subs(s.Max(i, r), i)).args[0]
    assert nu_low.cond.subs(r, 0) == s.true
    zero(s.limit(9 * r * nu_low.expr, r, 0, dir="+") - cov.origin.doit())
    assert not adj.clock_lapse(r).has(s.Integral)
    # Differentiating the source coefficient would violate the transpose.
    G = RegularFieldGreen(s.Matrix([1, 0, 0]), t, r)
    zero(G.spatial(t * B, 1) - t * G.spatial(B, 1))
    assert G.spatial(t * B, 1) != G.spatial(t * B).diff(t)
    print(
        "[ok] complete weak transpose, nonzero observer trace, emission jump, singular-origin row and kernel-only time derivatives"
    )


def check_quadrature():
    import mpmath as mp

    mp.mp.dps = 35

    def TT(r, a, d):
        return d**r * mp.hyper([1], [(r + 1) / 2, (r + 2) / 2], -a * d * d / 4) / mp.factorial(r)

    def gn(t, r, a):
        d = t - r
        return (TT(1, a, d) + (d * TT(2, a, d) - TT(3, a, d)) / (t * r)) / 3

    t, r = s.symbols("t r", positive=True)
    q = s.Symbol("q", nonnegative=True)
    G = RegularFieldGreen(s.Matrix([q, 0, 0]), t, r)
    fun = s.lambdify((q, t, r), G.coupled_scalar, modules=[{"TimeTail": TT}, "mpmath"])
    clock = s.lambdify((q, t, r), G.clock_primitive(t), modules=[{"TimeTail": TT}, "mpmath"])
    z = s.Symbol("z", nonnegative=True)
    parameter_lapse = RegularFieldGreen(s.Matrix([s.sqrt(z), 0, 0]), t, r).lapse
    lapse_resolvent = s.lambdify((z, t, r), -parameter_lapse.diff(z), modules=[{"TimeTail": TT}, "mpmath"])
    spectral = s.Symbol("spectral", real=True)
    spatial_resolvent = s.lambdify(
        (spectral, t, r), -(r * TimeTail(1, spectral, t - r) / t).diff(spectral), modules=[{"TimeTail": TT}, "mpmath"]
    )
    physical = s.lambdify((q, r, t), physical_first_clock(s.Matrix([q, 0, 0]), r, t), modules=[{"TimeTail": TT}, "mpmath"])
    for kval in [mp.mpf(0), mp.mpf("1e-9"), mp.mpf("0.7"), mp.mpf("4.1")]:
        end, start = mp.mpf("2.3"), mp.mpf("0.4")
        a = kval * kval / 3
        direct = mp.quad(lambda u, a=a, end=end, start=start: u * TT(1, -a, end - u) * gn(u, start, a) / end, [start, end])
        assert abs(fun(kval, end, start) - direct) < mp.mpf("1e-29")
        direct = mp.quad(lambda u, a=a, start=start: u * gn(u, start, a), [start, end])
        assert abs(clock(kval, end, start) - direct) < mp.mpf("1e-29")
        lapse = lambda u, a=a: TT(0, a, u) - TT(1, a, u) / u
        direct = -(start**2 * lapse(start) + mp.quad(lambda u, lapse=lapse: u * lapse(u), [start, end])) / end
        assert abs(physical(kval, start, end) - direct) < mp.mpf("1e-29")
        direct = mp.quad(lambda u, a=a, start=start, end=end: gn(end, u, a) * gn(u, start, a), [start, end])
        assert abs(lapse_resolvent(kval * kval, end, start) - direct) < mp.mpf("1e-29")
        direct = mp.quad(
            lambda u, a=a, start=start, end=end: start * TT(1, -a, end - u) * TT(1, -a, u - start) / end, [start, end]
        )
        assert abs(spatial_resolvent(-a, end, start) - direct) < mp.mpf("1e-29")
    print(
        "[ok] independent finite-frequency quadrature validates coupled/clock integrals, exact and near-zero momenta (numerical validation)"
    )


if __name__ == "__main__":
    check_entire_and_operators()
    check_weak_transpose()
    check_quadrature()
    print("REGULAR FIELD GREEN TRANSPOSE VERIFIED; remaining generic source/time pairings and spectrum unevaluated")
