"""Exact weak-source insertion and conservative time-arity certificates.

Generic pairings are CONSTRUCTED, not evaluated here. The homogeneous lapse
origin is additionally checked against the independent leading gradient result.
"""

from functools import lru_cache
from unittest.mock import patch

import sympy as s
from quadrupole_angular_source import QuadrupoleAngularSource
from regular_field_green import TimeTail, physical_quadrupole_adjoints
from unified_adjoint_source import inner, stf


@lru_cache(None)
def _arity(expr):
    if isinstance(expr, s.Integral):
        return len(expr.limits) + _arity(expr.function)
    if isinstance(expr, s.Add):
        return max(map(_arity, expr.args), default=0)
    if isinstance(expr, s.Mul):
        return sum(map(_arity, expr.args))
    if isinstance(expr, s.Pow):
        base, power = expr.args
        if power.is_Integer and power >= 0:
            return int(power) * _arity(base)
        if _arity(base):
            raise ValueError("Nonpolynomial function of a definite integral")
    if isinstance(expr, s.Piecewise):
        return max((_arity(case.expr) for case in expr.args), default=0)
    value = max((_arity(arg) for arg in expr.args), default=0)
    if value:
        raise ValueError(("Hidden nonlinear integral", expr.func))
    return 0


def arity(expr):
    if isinstance(expr, s.MatrixBase):
        return max(map(_arity, expr), default=0)
    return _arity(s.sympify(expr))


def zero(expr):
    if isinstance(expr, s.MatrixBase):
        for x in expr:
            zero(x)
    else:
        assert s.expand(expr) == 0, expr


def lazy_build(operation):
    # sympy.trace(Matrix(...)) used to call Integral.doit behind our backs.
    with patch.object(s.Integral, "doit", side_effect=AssertionError("Unexpected automatic integration")):
        return operation()


def physical_pair_source(time, kb, kc):
    import engine as eng
    from regular_field_green import physical_first_fields
    from stationary_cubic import direct_field_residuals

    saved = eng.Ctx.maxdeg, eng.Ctx.momenta, eng.Ctx.trunc
    try:
        eng.Ctx.maxdeg = 3
        eng.Ctx.momenta = {1: list(kb), 2: list(kc)}
        eng.Ctx.trunc = s.expand
        ab, hb = physical_first_fields(kb, eng.t)
        ac, hc = physical_first_fields(kc, eng.t)
        field = lambda b, c: eng.Field({frozenset([1]): b, frozenset([2]): c})
        N, H, c = direct_field_residuals(
            field(ab, ac), field(1, 1), [[field(hb[i, j], hc[i, j]) for j in range(3)] for i in range(3)]
        )
        N = N.coeff([1, 2]).subs(eng.t, time)
        H = s.Matrix(3, 3, lambda i, j: H[i][j].coeff([1, 2]).subs(eng.t, time))
        c = c.coeff([1, 2])
        zero(c - (2 * (kb.dot(kb) + kc.dot(kc)) - kb.dot(kc)) / 9)
        return N, H, c
    finally:
        eng.Ctx.maxdeg, eng.Ctx.momenta, eng.Ctx.trunc = saved


def check_source_insertion():
    t, u, r = s.symbols("t u r", positive=True)
    p = s.Matrix([1, 0, 0])
    q = s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0])
    C = s.Matrix(
        [
            [s.Rational(2, 3), s.Rational(1, 7), s.Rational(1, 5)],
            [s.Rational(1, 7), -s.Rational(1, 3), s.Rational(1, 11)],
            [s.Rational(1, 5), s.Rational(1, 11), -s.Rational(1, 3)],
        ]
    )
    background, first = lazy_build(lambda: physical_quadrupole_adjoints(p, q, C, t, u, 1, 2, s.Rational(1, 2)))
    for label, adjoint, maximum in [("background", background, 1), ("first", first, 2)]:
        cov = lazy_build(lambda adjoint=adjoint: adjoint.at(r))
        counts = [arity(value) for value in (cov.nu, cov.H, cov.origin)]
        assert max(counts) <= maximum, (label, counts)
        assert all(count > 0 for count in counts)
        assert not (cov.nu.free_symbols | cov.H.free_symbols | cov.origin.free_symbols) & set(s.symbols("n0:3", real=True))
        print(f"[ok] exact {label} adjoint inserted at noncollinear momenta: time arities (N,H,origin)={counts}")
    # The inferred source already carries its one linear sky/observer pairing;
    # closing the clock Green integral introduces NO new definite time integral.
    A = first.clock_amplitude
    assert first.clock_lapse(r).atoms(s.Integral) <= A.atoms(s.Integral)
    assert not s.Matrix(first.measurement.coefficients.H_prime_observer).has(s.Integral)
    kb = s.Matrix([0, 0, 1])
    kc = q - kb
    assert s.Matrix.hstack(p, kb, kc).det() != 0
    CN, CH, c_pair = physical_pair_source(r, kb, kc)
    assert arity(CN) == arity(CH) == 0 and c_pair != 0
    pair = lazy_build(lambda: first.stationary_pair(CN, CH, c_pair, r))
    assert arity(pair) <= 3
    # Actual dense kernel values are not claimed to have been integrated.
    assert pair.has(s.Integral) and pair.has(TimeTail)
    print(
        "[ok] actual Einstein single-pair functional at a noncoplanar unequal-frequency triple has time arity at most three, with its nonzero initial row"
    )
    # Supplying an independently evaluated clock does not change the generic API.
    angular = QuadrupoleAngularSource(p, q, C, 1, 2, s.Rational(1, 2))
    h = u * u * stf(p * p.T)
    a = u * u / 7
    generic = angular.assemble(1, h, a, t, u)
    D = -(s.Rational(1, 2) ** 4 / 7 + s.integrate(u * a, (u, s.Rational(1, 2), 2))) / 2
    closed = angular.assemble(1, h, a, t, u, clock_shift=D)
    for name in vars(generic.coefficients):
        x = getattr(generic.coefficients, name)
        y = getattr(closed.coefficients, name)
        zero((x - y).doit())
    print(
        "[ok] optional evaluated first clock agrees with the full off-shell assembly, including both physical event coefficients"
    )


def check_complete_three_time_bound():
    from direct_cubic_measurement import FourierFirst, first_only_exponent_triple, normalized_cubic_raw
    from regular_field_green import physical_first_fields

    t = s.Symbol("field_time", positive=True)
    ps = (
        s.Matrix([1, 0, 0]),
        s.Matrix([s.Rational(3, 5), s.Rational(4, 5), 0]),
        s.Matrix([s.Rational(2, 3), s.Rational(4, 3), s.Rational(4, 3)]),
    )
    seeds = tuple(FourierFirst(p, 1, *physical_first_fields(p, t), t) for p in ps)
    n = s.Matrix(s.symbols("n0:3", real=True))
    data = first_only_exponent_triple(seeds, n, 1, 2, s.Rational(1, 2))
    raw = normalized_cubic_raw(data)
    count = arity(raw)
    assert count <= 3 and count == 3
    print(
        "[ok] complete normalized direct measurement at noncoplanar unequal-frequency physical first inputs has time arity three (frame and clock retained)"
    )


def check_homogeneous_probe_origin():
    t, u = s.symbols("t u", positive=True)
    k = s.Symbol("k", real=True)
    direction = s.Matrix([1, 0, 0])
    C = stf(direction * direction.T)
    background, first = physical_quadrupole_adjoints(k * direction, s.zeros(3, 1), C, t, u, 1, 2, s.Rational(1, 2))
    M = first.measurement.coefficients
    zero(first.spatial(2) - s.Matrix(M.H_prime_observer))
    assert first.origin != 0 and arity(first.origin) <= 2
    # The coefficient is exact; it does NOT assert a finite-wavelength field
    # comparison beyond the already proved four-gradient sector.
    leading = first.origin.diff(k, 2).subs(k, 0).doit() / 2
    expected = inner(C, stf(direction * direction.T)) * (3 * 2**4 + 5 * s.Rational(1, 2) ** 4 - 8) / 20
    zero(s.simplify(leading - expected))
    assert leading != 0
    assert arity(background.lapse(t)) == 1  # query time may equal source API symbol
    print(
        "[ok] zero probe momentum: nonzero origin covector reproduces the physical gradient reference; observer trace and fresh-time binding preserved"
    )


if __name__ == "__main__":
    check_source_insertion()
    check_homogeneous_probe_origin()
    check_complete_three_time_bound()
    print("EXACT QUADRUPOLE GREEN SOURCE INSERTION VERIFIED; generic time integrals and fourth-order spectrum remain unevaluated")
