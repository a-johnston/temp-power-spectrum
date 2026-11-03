"""Exact entire moment identities, plus separate finite-phase quadrature tests."""

from functools import lru_cache
from itertools import product

import numpy as np
import sympy as s
from entire_sphere import RadialBessel, phase_average, phase_moment, sphere_reduce

n = s.Matrix(s.symbols("nx ny nz", real=True))
K = s.Matrix(s.symbols("Kx Ky Kz", real=True))
z = K.dot(K)


def zero(expr):
    if isinstance(expr, s.MatrixBase):
        for value in expr:
            zero(value)
    else:
        assert s.simplify(s.expand(expr)) == 0, expr


def powers_through(rank):
    return (p for p in product(range(rank + 1), repeat=3) if sum(p) <= rank)


@lru_cache(None)
def characteristic_derivative(powers):
    if sum(powers) == 0:
        return RadialBessel(0, z)
    i = next(i for i, m in enumerate(powers) if m)
    previous = list(powers)
    previous[i] -= 1
    return s.diff(characteristic_derivative(tuple(previous)), K[i])


def reduce_radial(expr, argument=z):
    """Use F_r=(2r+3)F_(r+1)-z F_(r+2) in a linear radial expression."""
    atoms = expr.atoms(RadialBessel)
    if not atoms:
        return s.expand(expr)
    high = max(int(f.args[0]) for f in atoms)
    if high < 2:
        return s.expand(expr)

    @lru_cache(None)
    def rec(r):
        if r >= high - 1:
            return RadialBessel(r, argument)
        return (2 * r + 3) * rec(r + 1) - argument * rec(r + 2)

    return s.expand(expr.xreplace({f: rec(int(f.args[0])) for f in atoms}))


def quadrature(poly, direction, phase, order=40):
    """Independent Gauss-Legendre/azimuth quadrature, normalized to <1>=1."""
    if isinstance(poly, s.MatrixBase):
        return np.array([[quadrature(poly[i, j], direction, phase, order) for j in range(poly.cols)] for i in range(poly.rows)])
    zz, ww = np.polynomial.legendre.leggauss(order)
    az = np.arange(2 * order) * np.pi / order
    X = np.sqrt(1 - zz[:, None] ** 2) * np.cos(az)[None, :]
    Y = np.sqrt(1 - zz[:, None] ** 2) * np.sin(az)[None, :]
    Z = np.broadcast_to(zz[:, None], X.shape)
    points = (X, Y, Z)
    val = s.lambdify(tuple(direction), poly, "numpy")(*points)
    phi = np.exp(1j * sum(complex(a) * b for a, b in zip(phase, points, strict=True)))
    return np.sum(val * phi * ww[:, None]) / (4 * order)


def numerical(expr):
    expr = expr.xreplace({f: f.evalf(50) for f in expr.atoms(RadialBessel)})
    if isinstance(expr, s.MatrixBase):
        return np.array(expr.applyfunc(lambda z: complex(z.evalf(40))).tolist(), dtype=complex)
    return complex(s.sympify(expr).evalf(40))


def exact_checks():
    count = 0
    for p in powers_through(7):
        zero(phase_moment(p, K) - characteristic_derivative(p) / s.I ** sum(p))
        got = phase_moment(p, s.zeros(3, 1))
        expected = 0 if any(m % 2 for m in p) else s.prod(s.factorial2(m - 1) for m in p) / s.factorial2(sum(p) + 1)
        zero(got - expected)
        count += 1
    for p in powers_through(4):
        trace = sum(phase_moment(tuple(m + 2 * (j == i) for j, m in enumerate(p)), K) for i in range(3)) - phase_moment(p, K)
        zero(reduce_radial(trace))
    x = s.Symbol("x", positive=True)
    zero(s.hyperexpand(RadialBessel(0, x * x).rewrite(s.hyper)) - s.sin(x) / x)
    argument = s.Symbol("argument")
    for r in range(8):
        zero(s.diff(RadialBessel(r, argument).rewrite(s.hyper), argument) + RadialBessel(r + 1, argument).rewrite(s.hyper) / 2)
    # Sphere reduction is algebraic, not a direction approximation.
    poly = (n.dot(n) - 1) * (n[0] + 2 * n[1] * n[2]) ** 3
    zero(sphere_reduce(poly, n))
    zero(phase_average(poly, n, K))
    print(f"[ok] {count} exact Cartesian derivative/zero-phase moments through rank seven and raw unit-sphere traces")
    print("[ok] entire hypergeometric base equals j0; no phase-length division or piecewise zero prescription")


def finite_phase_checks():
    polynomials = [s.S.One, n[0], n[0] ** 2 * n[1] ** 2, n[0] ** 3 * n[1] * n[2] ** 3]
    phases = [
        s.zeros(3, 1),
        s.Matrix([s.Rational(1, 10**10), 0, 0]),
        s.Matrix([1, s.Rational(-2, 3), s.Rational(5, 4)]),
        s.Matrix([3, -2, 4]),
    ]
    for phase in phases:
        for poly in polynomials:
            got = numerical(phase_average(poly, n, phase))
            check = quadrature(poly, n, phase, 40)
            refined = quadrature(poly, n, phase, 56)
            assert abs(check - refined) < 3e-13, (phase, poly, check - refined)
            assert abs(got - refined) < 3e-13, (phase, poly, got - refined)
    print("[ok] independent finite-phase quadrature, including exact and near-zero phases (numerical validation, not the proof)")


if __name__ == "__main__":
    exact_checks()
    finite_phase_checks()
    print("ENTIRE SPHERICAL PHASE MOMENTS VERIFIED")
