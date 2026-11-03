"""Exact cubic orthonormal integrand, endpoint and clock identities.

Generic symmetric non-commuting matrices; no shell, gradient, or time-series
approximation. Matrix jets are pulled back to the ray before coefficient
extraction. Negative controls detect the off-shell commutator and boundary terms.
"""

import sympy as s
from check_orthonormal_map import sym


def check_driver():
    zero = s.zeros(3)
    eye = s.eye(3)
    H1, H2, H3 = [sym(x) for x in ("H", "J", "K")]
    P1, P2, P3 = [sym(x) for x in ("P", "Q", "R")]
    H, P = [zero, H1, H2 / 2, H3 / 6], [zero, P1, P2 / 2, P3 / 6]

    def mul(A, B):
        return [sum((A[j] * B[i - j] for j in range(i + 1)), zero) for i in range(4)]

    Hsq = mul(H, H)
    Hcu = mul(Hsq, H)
    # e^H and its partial time derivative; zeta is time independent at a
    # fixed point. In the contracted driver its conformal factors cancel.
    dM = [
        P[i] + (mul(P, H)[i] + mul(H, P)[i]) / 2 + (mul(P, Hsq)[i] + mul(mul(H, P), H)[i] + mul(Hsq, P)[i]) / 6 for i in range(4)
    ]
    R = [(eye if i == 0 else zero) - H[i] / 2 + Hsq[i] / 8 - Hcu[i] / 48 for i in range(4)]
    B = mul(mul(R, dM), R)
    comm = H1 * (H1 * P1 - P1 * H1) - (H1 * P1 - P1 * H1) * H1
    assert (B[3] / 2 - P3 / 12 - comm / 48).applyfunc(s.expand) == zero
    # Fix n along z without loss of generality; matrices and direction jets
    # remain arbitrary. Unit-vector constraint fixes the second radial jet.
    n = s.Matrix([0, 0, 1])
    x, y, u, w = s.symbols("x y u w")
    v1 = s.Matrix([x, y, 0])
    v2 = s.Matrix([u, w, x * x + y * y])
    vv = [-n, v1, v2 / 2, s.zeros(3, 1)]
    old = 3 * sum((vv[i].T * B[j] * vv[3 - i - j])[0] for i in range(4) for j in range(4 - i))
    bil = lambda A, a, b: (a.T * A * b)[0]
    new = bil(P3 / 2 + comm / 8, n, n) - 3 * bil(P2, n, v1) - 3 * bil(P1, n, v2) + 3 * bil(P1, v1, v1)
    assert s.expand(old - new) == 0
    assert s.expand(bil(comm, n, n)) != 0  # dropping it off-shell must fail
    print("[ok] cubic driver with generic non-commuting matrices and unit direction, including commutator")
    # Pulled-back partial-time jets supply every position product. Use a
    # generic scalar component and independent directional derivative jets.
    X1, X2, h1x, h1xx, h2x = s.symbols("X1 X2 h1x h1xx h2x")
    eps, X = s.symbols("eps X")
    h1, h2, h3 = s.symbols("h1 h2 h3")
    field = eps * (h1 + X * h1x + X**2 * h1xx / 2) + eps**2 * (h2 + X * h2x) / 2 + eps**3 * h3 / 6
    third = s.expand(field.subs(X, eps * X1 + eps**2 * X2 / 2)).coeff(eps, 3) * 6
    assert s.expand(third - h3 - 3 * X1 * h2x - 3 * X2 * h1x - 3 * X1**2 * h1xx) == 0
    print("[ok] pulled-back cubic jets supply all ray-position products")


def check_boundaries():
    eps, x, z = s.symbols("eps x z")
    d1, d2, D1, D2 = s.symbols("d1 d2 D1 D2")
    f1, f1t, f2, f1o, f2o, f1oo, d1o = s.symbols("f1 f1t f2 f1o f2o f1oo d1o")
    U1, U1o, U1oo, U2, U2o, U3 = s.symbols("U1 U1o U1oo U2 U2o U3")
    F = eps * (f1 + x * f1t + z * f1o) + eps**2 * (f2 + z * f2o) / 2
    lower = eps * (d1 + z * d1o) + eps**2 * d2 / 2
    boundary = s.integrate(F, (x, 0, lower))
    integral = eps * (U1 + z * U1o + z * z * U1oo / 2) + eps**2 * (U2 + z * U2o) / 2 + eps**3 * U3 / 6
    composed = s.expand((integral - boundary).subs(z, eps * D1 + eps**2 * D2 / 2)).coeff(eps, 3) * 6
    target = (
        U3
        - 3 * f1 * d2
        - 3 * f2 * d1
        - 3 * f1t * d1 * d1
        + 3 * D1 * (U2o - 2 * (f1o * d1 + f1 * d1o))
        + 3 * D2 * U1o
        + 3 * D1 * D1 * U1oo
    )
    assert s.expand(composed - target) == 0
    for omitted in (
        -3 * f1 * d2,
        -3 * f2 * d1,
        -3 * f1t * d1 * d1,
        -6 * D1 * (f1o * d1 + f1 * d1o),
        3 * D2 * U1o,
        3 * D1 * D1 * U1oo,
    ):
        assert s.expand(composed - (target - omitted)) != 0
    print("[ok] moving emission and physical clock cubic identity; all six omission controls detected")


if __name__ == "__main__":
    check_driver()
    check_boundaries()
    print("CUBIC PRODUCT IDENTITIES PASSED")
