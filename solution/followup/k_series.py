"""Exact t-series of the second-order pair (q1, q2) in the plane-symmetric system, at K = q1 + q2 != 0.

Needed for forward_K.py: at small K, (N) has the homogeneous solution nu ~ t^2 (the constraint-violating curvature mode),
so truncation error in the K-channel initial data grows like t^2 and enters f' through (M). Starting from this
series removes that error to the order retained.

nu_2 from (N), with the t^2 coefficient fixed by the exact gradient data; f_2' from (M) at K != 0.
"""

import sympy as sp

t, y, e, z1, z2 = sp.symbols("t y e z1 z2", positive=True)


def first_order(q, order):
    """Series of n1 = cos x - sin x/x and h = 6(1 + cos x - 2 sin x/x), x = q t / sqrt3."""
    u = q**2 * t**2 / 3  # x^2
    n1 = 0
    h = 0
    for m in range(1, order // 2 + 2):
        c = (-1) ** m * u**m / sp.factorial(2 * m)
        j = (-1) ** m * u**m / sp.factorial(2 * m + 1)
        n1 += c - j
        h += 6 * (c - 2 * j)
    return n1, h


def trunc(expr, order):
    p = sp.Poly(sp.expand(expr * t**2), t)
    return sum(c * t ** (m[0] - 2) for m, c in p.terms() if m[0] - 2 <= order)


def pair_series(q1, q2, order=8):
    """Return (coefficient lists) nu2 = {m: c}, f2 = {m: d} for even m in 2..order, as sympy rationals."""
    q1, q2 = sp.Integer(q1), sp.Integer(q2)
    K = q1 + q2
    E1, E2 = sp.exp(sp.I * q1 * y), sp.exp(sp.I * q2 * y)
    n1a, ha = first_order(q1, order)
    n1b, hb = first_order(q2, order)
    g = z1 * E1 + z2 * E2
    Z = e * g
    nu1 = e * (n1a * z1 * E1 + n1b * z2 * E2)
    f1 = e * (ha * z1 * E1 + hb * z2 * E2)
    D = sp.diff

    def pick(expr):
        return sp.simplify(sp.expand(D(expr, e, 2).subs(e, 0) / 2).coeff(z1, 1).coeff(z2, 1) / (E1 * E2))

    # (N): 3 nu2'' + K^2 nu2 - 6 nu2 / t^2 + src = 0, src from first-order products.
    fy, ny, zy = D(f1, y), D(nu1, y), D(Z, y)
    ex = sp.exp(-2 * f1 / 3 - 2 * Z)
    lap = ex * (D(ny, y) - (zy + fy / 3) * ny + 2 * (zy - fy / 6) * ny)
    N_src = 3 * D(nu1, t, 2) - lap - 2 * ex * ny**2 + D(f1, t) ** 2 / 6 + 3 / t**2 * (sp.exp(-2 * nu1) - 1)
    src = trunc(pick(N_src), order)

    # Exact gradient value of the t^2 coefficient: nu = (t^2/9) e^{-2 zeta}(zeta_yy + zeta_y^2/2).
    grad = sp.exp(-2 * Z) * (D(Z, y, 2) + D(Z, y) ** 2 / 2) / 9
    c = {2: sp.nsimplify(pick(grad))}
    nu2 = c[2] * t**2
    for m in range(4, order + 1, 2):
        residual = sp.expand(3 * D(nu2, t, 2) + K**2 * nu2 - 6 * nu2 / t**2 + src)
        coeff = residual.coeff(t, m - 2)
        c[m] = -coeff / (3 * m * (m - 1) - 6)
        nu2 += c[m] * t**m

    # (M): (E w)_y = 3 E (2 nu'_y - 2 nu_y / t - 2 nu' nu_y), E = exp(3 zeta - f/2 + 2 nu), w = f'.
    nu_tot = nu1 + e**2 * nu2 * z1 * z2 * E1 * E2
    w1 = D(f1, t)
    Ee = sp.exp(3 * Z - f1 / 2 + 2 * nu1)
    nt = D(nu_tot, t)
    nyt = D(nu_tot, y)
    rhs = 3 * Ee * (2 * D(nt, y) - 2 * nyt / t - 2 * nt * nyt)
    lhs_known = D(Ee * w1, y)
    w2 = trunc(pick(rhs - lhs_known) / (sp.I * K), order)
    w2 = sp.expand(w2)
    d = {m: sp.nsimplify(w2.coeff(t, m - 1) / m) for m in range(2, order + 1, 2)}
    return c, d


if __name__ == "__main__":
    # Control 1: K -> 0 values. Generic near-cancelling pair approaches the closed K=0 series with q = 1 by scaling.
    # In x = q t/sqrt3: nu2 = 5/3 x^2 - 13/45 x^4, a = 2 x^2 - 11/30 x^4.
    q = 400
    c, d = pair_series(q, -(q - 1), order=4)
    s = lambda v, m: float(v) * (3 / (q - 0.5) ** 2) ** (m / 2)
    print("nu2 x^2, x^4:", s(c[2], 2), s(c[4], 4), " expect 5/3 =", 5 / 3, " -13/45 =", -13 / 45)
    print("f2  x^2, x^4:", s(d[2], 2), s(d[4], 4), " expect 2, -11/30 =", -11 / 30)
    # Control 2: the t^2 coefficient of f2 against the gradient formula f = (t^2/3) e^{-2 zeta}(zeta_yy - zeta_y^2).
    q1, q2 = 7, -5
    c, d = pair_series(q1, q2, order=4)
    exact = sp.Rational(1, 3) * (2 * (q1**2 + q2**2) + 2 * q1 * q2)
    print("f2 t^2 at (7,-5):", d[2], " gradient formula:", exact)
    assert d[2] == exact
    print("PASS")
