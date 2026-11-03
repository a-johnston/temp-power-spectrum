"""Exact sky contraction of polynomial plane waves and finite Volterra integrals.

Independent copies of an integral are alpha-renamed before multiplying. Already
sky-contracted factors remain scalar coefficients, so products of means are NOT
replaced by a single combined-phase average. Only phase-free polynomial times
and their elementary triangular convolutions are evaluated; oscillatory/Green
integrals are left explicit.
"""

from dataclasses import dataclass
from itertools import product

import sympy as s
from entire_sphere import angular_polynomial, contract_angular_polynomial
from sympy.polys.polyerrors import PolynomialError


@dataclass(frozen=True)
class PhaseSector:
    polynomial: object
    phase: tuple
    limits: tuple = ()

    def integral(self, expression):
        return s.Integral(expression, *self.limits) if self.limits else expression


@dataclass(frozen=True)
class ContractedPhases:
    expression: object
    sectors: tuple
    maximum_rank: int


def collapse_triangle(sector):
    """Evaluate a polynomial outer time with no phase dependence.

    int_L^O dt int_L^t du P(t,u) exp(i K(u).n)
      = int_L^O du [int_u^O P(t,u) dt] exp(i K(u).n).
    No oscillatory or Green integral is evaluated by this operation.
    """
    if len(sector.limits) != 2:
        return sector
    (u, L, upper), (t, lower, O) = sector.limits
    if upper != t or lower != L or L.has(u, t) or O.has(u, t) or any(k.has(t) for k in sector.phase):
        return sector
    try:
        poly = s.Poly(sector.polynomial, t)
    except PolynomialError:
        return sector
    result = s.Add(*(coefficient * (O ** (power + 1) - u ** (power + 1)) / (power + 1) for (power,), coefficient in poly.terms()))
    return PhaseSector(result, sector.phase, ((u, L, O),))


def polynomial_time_reduction(sector):
    """Evaluate only polynomial times absent from the oscillatory phase."""
    while True:
        changed = False
        for i, (var, lo, hi) in enumerate(sector.limits):
            if any(k.has(var) for k in sector.phase):
                continue
            if any(bound.has(var) for j, limit in enumerate(sector.limits) if i != j for bound in limit[1:]):
                continue
            try:
                poly = s.Poly(sector.polynomial, var)
            except PolynomialError:
                continue
            coefficient = s.Add(*(cf * (hi ** (power + 1) - lo ** (power + 1)) / (power + 1) for (power,), cf in poly.terms()))
            sector = PhaseSector(coefficient, sector.phase, sector.limits[:i] + sector.limits[i + 1 :])
            changed = True
            break
        if changed:
            continue
        reduced = collapse_triangle(sector)
        if reduced != sector:
            sector = reduced
            continue
        return sector


def phase_sectors(expression, n):
    """Fubini-expand finite products, preserving nested causal limits.

    The inner-to-outer order of Integral.limits is retained. A repeated integral
    in a power is parsed anew for each multiplicative copy, never identified
    with the same time variable. Bounds must be independent of sky direction.
    """
    n = tuple(n)
    expression = s.sympify(expression)

    def visit(expr):
        if expr == 0:
            return []
        if expr.is_polynomial(*n):
            return [PhaseSector(expr, (s.S.Zero,) * 3)]
        if isinstance(expr, s.Integral):
            if any(bound.has(*n) for limit in expr.limits for bound in limit[1:]):
                raise ValueError("Sky-dependent integration bounds are not supported")
            if any(len(limit) != 3 for limit in expr.limits):
                raise ValueError("Only specified definite integrals are supported")
            var, lo, hi = expr.limits[-1]
            fresh = s.Dummy(str(var), real=True)
            inside = s.Integral(expr.function, *expr.limits[:-1]) if len(expr.limits) > 1 else expr.function
            # subs respects scopes of nested integrals; xreplace would not.
            inside = inside.subs(var, fresh)
            return [PhaseSector(term.polynomial, term.phase, term.limits + ((fresh, lo, hi),)) for term in visit(inside)]
        if expr.func == s.exp:
            arg = expr.args[0]
            offset = arg.subs(dict.fromkeys(n, 0))
            phase = tuple(s.expand(s.diff(arg, x) / s.I) for x in n)
            if (
                any(k.has(*n) for k in phase)
                or s.expand(arg - offset - s.I * sum(k * x for k, x in zip(phase, n, strict=True))) != 0
            ):
                raise ValueError("The angular exponential is not a Cartesian plane wave")
            return [PhaseSector(s.exp(offset), phase)]
        if expr.is_Add:
            return [term for arg in expr.args for term in visit(arg)]

        def multiply(left, right):
            return [
                PhaseSector(
                    a.polynomial * b.polynomial, tuple(x + y for x, y in zip(a.phase, b.phase, strict=True)), a.limits + b.limits
                )
                for a, b in product(left, right)
            ]

        if expr.is_Mul:
            out = [PhaseSector(s.S.One, (s.S.Zero,) * 3)]
            for arg in expr.args:
                out = multiply(out, visit(arg))
            return out
        if expr.is_Pow and expr.exp.is_Integer and expr.exp >= 0:
            out = [PhaseSector(s.S.One, (s.S.Zero,) * 3)]
            for _ in range(int(expr.exp)):
                out = multiply(out, visit(expr.base))
            return out
        raise ValueError(f"Non-polynomial, non-plane-wave sky dependence: {expr.func}")

    raw = [polynomial_time_reduction(term) for term in visit(expression)]
    # Fresh per invocation: avoid capturing bound variables in scalar means
    # created by an earlier contraction, even when they depend on an outer time.
    canonical = tuple(s.Dummy(f"tau{i}", real=True) for i in range(max((len(a.limits) for a in raw), default=0)))
    grouped = {}
    for term in raw:
        renaming = {limit[0]: canonical[i] for i, limit in enumerate(term.limits)}
        limits = tuple(tuple(s.sympify(x).xreplace(renaming) for x in limit) for limit in term.limits)
        phase = tuple(s.expand(x.xreplace(renaming)) for x in term.phase)
        key = (phase, limits)
        grouped.setdefault(key, []).append(term.polynomial.xreplace(renaming))
    return tuple(
        PhaseSector(s.Add(*coefficients), phase, limits)
        for (phase, limits), coefficients in grouped.items()
        if s.Add(*coefficients) != 0
    )


def contract_phase_integrals(expression, n, weight=1, maximum_rank=None):
    """Evaluate sky integrals; reduce phase-free polynomial times only.

    The input contains no sky integral operator: it is a sky function whose
    normalized average, optionally weighted, is requested. Mean factors already
    contracted by an earlier call are n-independent and remain coefficients.
    """
    if isinstance(expression, s.MatrixBase):
        return expression.applyfunc(lambda z: contract_phase_integrals(z, n, weight, maximum_rank).expression)
    sectors = phase_sectors(expression, n)
    rank = 0
    terms = []
    for sector in sectors:
        poly, restore = angular_polynomial(weight * sector.polynomial, n)
        rank = max(rank, 0 if poly.is_zero else int(poly.total_degree()))
        terms.append(sector.integral(contract_angular_polynomial(poly, restore, sector.phase, maximum_rank)))
    return ContractedPhases(s.Add(*terms), sectors, rank)
