"""Exact K=0 pair fields (q,-q), closed forms, and their checks.

Conventions: q = 1 (restore with t -> q t), x = q t / sqrt3. zeta = g = z1 E + z2 / E, E = exp(i q y).
All quantities are FULL-SERIES coefficients of z1 z2 (report convention K_2, so F_2 = K_2 / 2).
"""

import sympy as sp
from plane_reduction import f, nu, reduced_equations, t, y, zeta

X = sp.symbols("x", positive=True)
s3 = sp.sqrt(3)


def Cin2(x):
    """Cin(2x) = int_0^{2x} (1-cos u)/u du, entire."""
    return sp.Integer(0) if x == 0 else sp.Function("Cin2")(x)


u1 = sp.sin(X) / X - sp.cos(X)
n1 = -u1
h1 = -6 * (2 * sp.sin(X) / X - sp.cos(X) - 1)  # H_1 = h1 (qhat qhat - 1/3) g
n2 = -sp.cos(2 * X) - 6 * sp.cos(X) + 6 * sp.sin(X) / X + 3 * sp.sin(2 * X) / X - 5 * sp.sin(X) ** 2 / X**2
Cin = sp.Function("C2")(X)  # stands for Cin(2x); C2'(x) = (1-cos 2x)/x
a2 = (
    -78 - 36 * sp.cos(X) - 6 * sp.cos(2 * X) + 72 * sp.sin(X) / X + 36 * sp.sin(2 * X) / X - 24 * sp.sin(X) ** 2 / X**2 + 24 * Cin
)


def dx(e):
    return sp.diff(e, X).subs(sp.Derivative(Cin, X), (1 - sp.cos(2 * X)) / X)


def to_t(e):
    return e.subs(X, t / s3)


def check_first_order():
    eq = reduced_equations()
    eps = sp.Symbol("eps")
    E = sp.exp(sp.I * y)
    res = {}
    for key in ("N", "T", "Ham", "M"):
        e = eq[key].subs({f: eps * to_t(h1) * E, nu: eps * to_t(n1) * E, zeta: eps * E}).doit()
        res[key] = sp.simplify(sp.diff(e, eps).subs(eps, 0))
    return res


def check_pair():
    """Insert the closed forms into all reduced equations and take the z1 z2 coefficient."""
    eq = reduced_equations()
    z1, z2, eps = sp.symbols("z1 z2 eps")
    E = sp.exp(sp.I * y)
    g = z1 * E + z2 / E
    A2 = sp.Function("A2")(t)
    N2 = sp.Function("N2")(t)
    res = {}
    for key in ("T", "N", "Ham", "M"):
        e = (
            eq[key]
            .subs({f: eps * to_t(h1) * g + eps**2 * A2 * z1 * z2, nu: eps * to_t(n1) * g + eps**2 * N2 * z1 * z2, zeta: eps * g})
            .doit()
        )
        e = sp.expand(sp.diff(e, eps, 2).subs(eps, 0) / 2)
        c = e.coeff(z1, 1).coeff(z2, 1)
        c = c.subs(sp.Derivative(A2, (t, 2)), dx(dx(a2)) / 3).subs(sp.Derivative(A2, t), dx(a2) / s3)
        c = c.subs(sp.Derivative(N2, (t, 2)), dx(dx(n2)) / 3).subs(sp.Derivative(N2, t), dx(n2) / s3)
        c = c.subs({A2: a2, N2: n2}).subs(t, s3 * X)
        res[key] = sp.simplify(sp.expand(sp.expand_trig(sp.expand(c))))
    return res


def series_checks():
    """Small-x series of a2 and n2 against the exact gradient data, and the late-time forms."""
    A = a2.subs(Cin, sum((-1) ** (m + 1) * (2 * X) ** (2 * m) / (2 * m * sp.factorial(2 * m)) for m in range(1, 8)))
    a_series = sp.series(A, X, 0, 6).removeO()
    n_series = sp.series(n2, X, 0, 6).removeO()
    return a_series, n_series


def KR_parts():
    """Isotropic and quadrupolar parts of 2 F_2 (K->0): 2F2 = I(x) - (mu^2-1/3)/4 * W(x)."""
    psi1 = 1 + n1
    Iso = n2 + h1 * psi1 - h1**2 / 12
    W = a2 + h1**2 / 2
    return Iso, W


if __name__ == "__main__":
    print("first-order residuals:", check_first_order())
    print("pair residuals (T,N,Ham,M):", check_pair())
    a_s, n_s = series_checks()
    print("a2 series:", a_s)
    print("n2 series:", n_s)
    Iso, W = KR_parts()
    print("Iso =", sp.simplify(sp.expand(sp.expand_trig(Iso))))
    print("W   =", sp.simplify(sp.expand(sp.expand_trig(W))))
    print("R1 = psi1 - h1/6 =", sp.simplify(1 + n1 - h1 / 6))
