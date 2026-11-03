"""Exact phase/Fubini/binding fixtures, not a radiation-spectrum evaluation."""

import sympy as s
from entire_sphere import RadialBessel, phase_average
from phase_integral_contraction import contract_phase_integrals, phase_sectors

n = s.Matrix(s.symbols("n0:3", real=True))
u, t = s.symbols("u t", real=True)
e, o = s.symbols("e o", positive=True)


def zero(expr):
    assert s.expand(expr) == 0, expr


def canonical(expression):
    return expression.as_dummy()


def checks():
    K = s.Matrix(s.symbols("K0:3", real=True))
    w = s.Rational(15, 2) * (n[0] ** 2 - n[1] ** 2)
    poly = n[0] ** 2 * n[2] + 3 * n[1]
    direct = contract_phase_integrals(poly * s.exp(s.I * K.dot(n)), n, w)
    zero(direct.expression - phase_average(poly, n, K, w))
    assert direct.maximum_rank == 5
    wave = s.exp(s.I * u * n[0])
    primitive = s.Integral(wave, (u, e, o))
    out = contract_phase_integrals(primitive**2, n)
    assert len(out.sectors) == 1
    a, b = (limit[0] for limit in out.sectors[0].limits)
    assert a != b
    zero(out.sectors[0].phase[0] - a - b)
    zero(out.expression - s.Integral(RadialBessel(0, (a + b) ** 2), (a, e, o), (b, e, o)))
    # Causal nesting and an independent copy sharing the SAME source name.
    nested = s.Integral(s.Integral(s.exp(s.I * (u + t) * n[0]), (u, t, o)), (t, e, o))
    out = contract_phase_integrals(nested * primitive, n)
    assert len(out.sectors) == 1
    limits = out.sectors[0].limits
    assert len({limit[0] for limit in limits}) == 3
    assert limits[0][1] == limits[1][0] or limits[1][1] == limits[2][0]
    # A product of means must not become a mean of a product.
    mean = contract_phase_integrals(primitive, n).expression
    separated = contract_phase_integrals(mean * primitive, n).expression
    assert separated.has(RadialBessel(0, u * u)) is False  # bound names were freshened
    assert len(separated.atoms(s.Integral)) >= 2
    assert canonical(separated) != canonical(contract_phase_integrals(primitive**2, n).expression)
    # Opposite phases: exact zero argument, not a small-argument expansion.
    zero(contract_phase_integrals((n.dot(n) - 1) * s.exp(s.I * K.dot(n)), n).expression)
    zero(contract_phase_integrals(s.exp(s.I * K.dot(n)) * s.exp(-s.I * K.dot(n)), n).expression - 1)
    # Integral scopes: the inner bound t shadows the outer t; its upper limit
    # is independent. The phase has two independent times, not twice one time.
    shadow = s.Integral(s.exp(s.I * t * n[0]) * s.Integral(s.exp(s.I * t * n[0]), (t, e, o)), (t, e, o))
    sectors = phase_sectors(shadow, n)
    assert len(sectors) == 1 and len(sectors[0].limits) == 2
    a, b = (limit[0] for limit in sectors[0].limits)
    zero(sectors[0].phase[0] - a - b)
    # A previously contracted mean depending on an outer time must not be
    # captured by the new canonical dummy when the outer integral is lifted.
    mean_t = contract_phase_integrals(s.Integral(s.exp(s.I * u * t * n[0]), (u, e, o)), n).expression
    outer = contract_phase_integrals(s.Integral(mean_t * s.exp(s.I * t * n[1]), (t, e, o)), n)
    assert len(outer.sectors) == 1
    z = outer.sectors[0].limits[0][0]
    embedded = next(iter(outer.sectors[0].polynomial.atoms(s.Integral)))
    v = embedded.limits[0][0]
    assert v != z and z in embedded.free_symbols
    # An arbitrary inner profile is retained while the phase-free polynomial
    # outer time is integrated exactly. A trigonometric outer time is NOT.
    f = s.Function("f")
    triangle = s.Integral(t**2 * f(u) * wave, (u, e, t), (t, e, o))
    sectors = phase_sectors(triangle, n)
    assert len(sectors) == 1 and len(sectors[0].limits) == 1
    v = sectors[0].limits[0][0]
    zero(sectors[0].polynomial - (o**3 - v**3) * f(v) / 3)
    zero(sectors[0].phase[0] - v)
    nonpolynomial = phase_sectors(s.Integral(s.sin(t) * f(u) * wave, (u, e, t), (t, e, o)), n)
    assert len(nonpolynomial) == 1 and len(nonpolynomial[0].limits) == 2
    print("[ok] polynomial triangular convolution evaluated with arbitrary inner profile; oscillatory outer integral retained")
    # Exact polynomial-in-time fixture allows all remaining integrals to be
    # evaluated, independently of the Fubini lift.
    fixture = s.Integral((1 + n[0] * u) * (n[1] + t), (u, t, o), (t, e, o)) ** 2
    zero(contract_phase_integrals(fixture, n, w).expression.doit() - phase_average(fixture.doit(), n, s.zeros(3, 1), w))
    for bad in (s.exp(n[0] ** 2), 1 / (1 + n[0]), s.Integral(wave, (u, e, n[1]))):
        try:
            contract_phase_integrals(bad, n)
        except ValueError:
            pass
        else:
            raise AssertionError(("invalid angular dependence accepted", bad))
    print("[ok] exact polynomial/plane-wave moments; independent repeated-integral variables and nested causal bounds")
    print("[ok] products of means, shadowed scopes and outer-time-dependent contracted scalars")
    print("[ok] exact zero phases, evaluated polynomial-time fixture and invalid-input guards")


if __name__ == "__main__":
    checks()
