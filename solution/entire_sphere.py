"""Exact Cartesian spherical phase moments, with no division by phase length.

F_r(z)=j_r(sqrt(z))/sqrt(z)^r is defined ENTIRELY in z by 0F1, including
F_r(0)=1/(2r+1)!!. Polynomial sphere integrals are finite sums of these
functions. No angular or small-phase approximation is used.
"""

from functools import lru_cache
from itertools import product
from math import factorial

import sympy as s


class RadialBessel(s.Function):
    """Entire F_r(z); the first argument is a nonnegative integer index."""

    nargs = 2

    @classmethod
    def eval(cls, r, z):
        if not isinstance(r, s.Integer) or r < 0:
            raise ValueError("RadialBessel index must be a nonnegative integer")
        if z.is_zero:
            return 1 / s.factorial2(2 * r + 1)

    def fdiff(self, argindex=2):
        if argindex != 2:
            raise ValueError("The RadialBessel index is discrete")
        r, z = self.args
        return -RadialBessel(r + 1, z) / 2

    def _eval_rewrite_as_hyper(self, r, z, **hints):
        return s.hyper((), (r + s.Rational(3, 2),), -z / 4) / s.factorial2(2 * r + 1)

    def _eval_evalf(self, prec):
        return self.rewrite(s.hyper)._eval_evalf(prec)


@lru_cache(None)
def moment_terms(powers):
    """(radial index, powers of K, coefficient) for a Cartesian monomial."""
    assert len(powers) == 3 and all(isinstance(m, int) and m >= 0 for m in powers)
    degree = sum(powers)
    out = []
    for pairs in product(*(range(m // 2 + 1) for m in powers)):
        remaining = tuple(m - 2 * r for m, r in zip(powers, pairs, strict=True))
        count = s.prod(
            s.Rational(factorial(m), 2**r * factorial(r) * factorial(m - 2 * r)) for m, r in zip(powers, pairs, strict=True)
        )
        out.append((degree - sum(pairs), remaining, s.I**degree * (-1) ** sum(pairs) * count))
    return tuple(out)


def phase_moment(powers, K):
    """Normalized <n_x^a n_y^b n_z^c exp(i K.n)>, exact for all K."""
    K = tuple(map(s.sympify, K))
    assert len(K) == 3
    z = sum(k * k for k in K)
    return s.Add(
        *(
            coef * RadialBessel(r, z) * s.prod(k**m for k, m in zip(K, remaining, strict=True))
            for r, remaining, coef in moment_terms(tuple(powers))
        )
    )


def sphere_reduce(expr, n):
    """Canonical polynomial representative modulo n.n=1 (eliminate n_z^2)."""
    return s.expand(
        sum(
            coef * n[2] ** (m % 2) * (1 - n[0] ** 2 - n[1] ** 2) ** (m // 2)
            for (m,), coef in s.Poly(s.expand(expr), n[2]).terms()
        )
    )


def polynomial_rank(expr, n):
    if isinstance(expr, s.MatrixBase):
        return max((polynomial_rank(z, n) for z in expr), default=0)
    poly = s.Poly(sphere_reduce(expr, n), *n)
    return 0 if poly.is_zero else int(poly.total_degree())


def phase_average(polynomial, n, K, weight=1, maximum_rank=None):
    """Normalized <weight(n) polynomial(n) exp(i K.n)>; scalar or matrix.

    Both input factors must be polynomials in n. The result contains at most
    F_0,...,F_rank at K.K; the coefficients are finite polynomials in K.
    Sphere reduction avoids redundant P_n^2 terms without changing the average.
    """
    if isinstance(polynomial, s.MatrixBase):
        return polynomial.applyfunc(lambda z: phase_average(z, n, K, weight, maximum_rank))
    poly = s.Poly(sphere_reduce(weight * polynomial, n), *n)
    if poly.is_zero:
        return s.S.Zero
    rank = int(poly.total_degree())
    if maximum_rank is not None and rank > maximum_rank:
        raise ValueError(f"Angular rank {rank} exceeds the specified bound {maximum_rank}")
    K = tuple(map(s.sympify, K))
    assert len(K) == 3
    coefficients = {}
    for powers, scalar in poly.terms():
        for r, remaining, count in moment_terms(tuple(map(int, powers))):
            coefficients.setdefault(r, []).append(scalar * count * s.prod(k**m for k, m in zip(K, remaining, strict=True)))
    z = sum(k * k for k in K)
    return s.Add(*(s.Add(*terms) * RadialBessel(r, z) for r, terms in sorted(coefficients.items())))
