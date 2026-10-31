"""Exact component helpers for derivations/QUADRATIC_INITIAL_DATA.md.

This module constructs symbolic quadratic sources. It does not evaluate the
retarded time integrals, evolve numerically, or compute the temperature map.
Indices are Cartesian 0..2; all amplitude kernels use the factorial convention.
"""

import sympy as s

I3 = s.eye(3)


def trace_part(x):
    return I3 * s.trace(x) / 3


def spatial_operator(momentum, x):
    k = s.Matrix(momentum)
    return ((k.dot(k)) * x - k * (k.T * x) - (x * k) * k.T
            + s.Rational(4, 3) * (k * k.T) * s.trace(x))


def time_operator(eta, x):
    return (x.diff(eta, 2) + 2 * x.diff(eta) / eta
            + 2 * I3 * s.trace(x.diff(eta)) / (3 * eta)
            + 2 * I3 * s.trace(x) / (3 * eta**2))


def green_zero(eta, time, x):
    """L_0 Green action for 0 <= time <= eta, eta > 0."""
    trace = trace_part(x)
    return time * (1 - time / eta) * (x - trace) + (
        time**2 / eta * (1 - time / eta) * trace
    )


def indicial_inverse(n, x):
    """n>=2, inverse of eta^n coefficient of eta^2 L_0."""
    trace = trace_part(x)
    return (x - trace) / (n * (n + 1)) + trace / ((n + 1) * (n + 2))


def linear_mode(eta, momentum):
    """Exact metric/lapse kernels (11), with the zero input handled explicitly."""
    k = s.Matrix(momentum)
    k2 = s.simplify(k.dot(k))
    if k2 == 0:
        return 2 * I3, s.S.Zero
    y = s.sqrt(k2 / 3) * eta
    sinc = s.sin(y) / y
    e = (2 * sinc - 1 - s.cos(y)) / y**2
    metric = 2 * sinc * I3 - 2 * eta**2 * e * (k * k.T)
    return metric, s.cos(y) - sinc


def linear_boundary_jet(eta, momentum):
    """Only for checking the leading time coefficient, not finite-time data."""
    k = s.Matrix(momentum)
    h2 = -(k * k.T) / 3 - I3 * k.dot(k) / 9
    return 2 * I3 + eta**2 * h2, -eta**2 * k.dot(k) / 9


def gamma_linear(momentum, h):
    k = s.Matrix(momentum)
    return [[[s.I * (k[i] * h[a, j] + k[j] * h[a, i] - k[a] * h[i, j]) / 2
              for j in range(3)] for i in range(3)] for a in range(3)]


def ricci_linear(momentum, h):
    k = s.Matrix(momentum)
    return (-k * (k.T * h) - (h * k) * k.T + k.dot(k) * h
            + (k * k.T) * s.trace(h)) / 2


def ricci_quadratic(k, p, h, f):
    """Symmetric Fourier kernel of R_Q in (12), not coefficient of zw."""
    k, p = s.Matrix(k), s.Matrix(p)
    total = k + p
    gh, gf = gamma_linear(k, h), gamma_linear(p, f)
    gq = [[[-sum(h[a, b] * gf[b][i][j] + f[a, b] * gh[b][i][j]
                  for b in range(3)) for j in range(3)]
            for i in range(3)] for a in range(3)]
    out = s.zeros(3)
    for i in range(3):
        for j in range(3):
            out[i, j] = sum(s.I * total[a] * gq[a][i][j]
                            - s.I * total[j] * gq[a][i][a] for a in range(3))
            out[i, j] += sum(
                gh[a][a][b] * gf[b][i][j] + gf[a][a][b] * gh[b][i][j]
                - gh[a][j][b] * gf[b][i][a] - gf[a][j][b] * gh[b][i][a]
                for a in range(3) for b in range(3)
            )
    return out


def quadratic_source(eta, k, p, h, f, c, d):
    """S in (14), polarized with 1/2 weight; h,c at k and f,d at p."""
    k, p = s.Matrix(k), s.Matrix(p)
    total = k + p
    hp, fp = h.diff(eta), f.diff(eta)
    cp, dp = s.diff(c, eta), s.diff(d, eta)
    v = s.trace(h * f) / 6
    gh, gf = gamma_linear(k, h), gamma_linear(p, f)
    gamma_grad = s.Matrix(3, 3, lambda i, j: sum(
        gh[a][i][j] * s.I * p[a] * d + gf[a][i][j] * s.I * k[a] * c
        for a in range(3)
    ))
    return (
        hp * fp + fp * hp - 2 * (cp * fp + dp * hp)
        + (4 * (v / eta**2 + s.diff(v, eta) / eta) + 8 * c * d / eta**2) * I3
        - 4 * ((c / eta**2 + cp / eta) * f + (d / eta**2 + dp / eta) * h)
        + 2 * (total * total.T) * v - 2 * gamma_grad
        - 2 * (k * p.T + p * k.T) * c * d
        - 2 * ricci_quadratic(k, p, h, f)
        + 4 * (c * (-(p * p.T) * d - ricci_linear(p, f))
               + d * (-(k * k.T) * c - ricci_linear(k, h)))
    )


def regular_source(eta, k, p, h, f, c, d):
    return (quadratic_source(eta, k, p, h, f, c, d)
            - 8 * I3 / eta**2 - spatial_operator(s.Matrix(k) + s.Matrix(p), 4 * I3))


def quadratic_boundary_coefficient(k, p):
    k, p = s.Matrix(k), s.Matrix(p)
    return (k * p.T + p * k.T) / 3 - s.Rational(4, 9) * I3 * k.dot(p)
