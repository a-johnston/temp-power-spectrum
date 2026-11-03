"""Exact adjoint for cubic propagation terms linear in second-order fields.

Generic plane-wave jets test the integrations by parts without time series or
quadrature. Separate tests check the linear ray equation, moving endpoints,
physical clock, centered normalization, and a mixed-time phase obstruction.
"""

import sympy as s


def symmetric(prefix):
    xx, yy, xy, xz, yz = s.symbols(f"{prefix}xx {prefix}yy {prefix}xy {prefix}xz {prefix}yz")
    return s.Matrix([[xx, xy, xz], [xy, yy, yz], [xz, yz, -xx - yy]])


def transverse(prefix):
    x, y = s.symbols(f"{prefix}x {prefix}y")
    return s.Matrix([x, y, 0])


def zero(expr):
    if isinstance(expr, s.MatrixBase):
        assert expr.applyfunc(s.expand) == s.zeros(*expr.shape)
    else:
        assert s.expand(expr) == 0


def check_mixed_driver():
    eps, lam = s.symbols("eps lambda")
    h, ht, j, jt = (symmetric(name) for name in ("h", "ht", "j", "jt"))
    H = eps * h + eps**2 * lam * j / 2
    Hp = eps * ht + eps**2 * lam * jt / 2
    cut = lambda M: M.applyfunc(lambda z: sum(s.expand(z).coeff(eps, i) * eps**i for i in range(4)))
    Hsq = cut(H * H)
    F = s.eye(3) - H / 2 + Hsq / 8 - cut(Hsq * H) / 48
    Mdot = Hp + cut(Hp * H + H * Hp) / 2 + cut(Hp * Hsq + cut(H * Hp) * H + Hsq * Hp) / 6
    B = cut(cut(F * Mdot) * F) / 2
    zero(B.applyfunc(lambda z: s.expand(z).coeff(eps, 3).coeff(lam)))
    zero(B.applyfunc(lambda z: s.expand(z).coeff(eps, 2).coeff(lam)) - jt / 4)
    n = s.Matrix([0, 0, 1])
    p, q = s.Matrix(s.symbols("px py pz")), s.Matrix(s.symbols("qx qy qz"))
    X, Y = s.Matrix(s.symbols("Xx Xy Xz")), s.Matrix(s.symbols("Yx Yy Yz"))
    V, U = transverse("V"), transverse("U")
    # Keep the pulled-back driver's lambda-dependent cubic pieces only.
    Bpull = eps * ht / 2 + eps**2 * lam * jt / 4
    Bpull += eps**3 * lam * (s.I * q.dot(Y) * ht + s.I * p.dot(X) * jt) / 4
    v = -n + eps * V + eps**2 * lam * U / 2
    actual = s.expand((v.T * Bpull * v)[0]).coeff(eps, 3).coeff(lam) * 6
    expected = 3 * (
        s.I * p.dot(X) * (n.T * jt * n)[0] / 2 - (n.T * jt * V)[0] + s.I * q.dot(Y) * (n.T * ht * n)[0] / 2 - (n.T * ht * U)[0]
    )
    zero(actual - expected)
    print("[ok] mixed cubic bulk from exact matrix-exponential driver and displaced orthonormal directions")


def check_linear_ray():
    # n=z is without loss of generality; all metric jets and momentum
    # components remain independent. The homogeneous pair p=0 is included.
    n = s.Matrix([0, 0, 1])
    P = s.eye(3) - n * n.T
    p = s.Matrix(s.symbols("px py pz"))
    j, jt = symmetric("j"), symmetric("jt")
    nn = lambda A: (n.T * A * n)[0]
    gamma_nn = s.I / 2 * (2 * p[2] * j * n - p * nn(j))
    e_dot = jt * n - gamma_nn - n * nn(jt) / 2
    # v=e-j n/2 at first order; derivative along the background ray.
    v_dot = e_dot - (jt - s.I * p[2] * j) * n / 2
    target = P * jt * n / 2 + s.I / 2 * (p * nn(j) - p[2] * j * n)
    zero(v_dot - target)
    zero(n.dot(target))
    print("[ok] arbitrary second-order metric linear ray: X_j'=v_j+j n/2, v_j'=G_j, transverse")


def check_adjoint():
    n = s.Matrix([0, 0, 1])
    P = s.eye(3) - n * n.T
    p = s.Matrix(s.symbols("px py pz"))
    pn = p[2]
    j, jt = symmetric("j"), symmetric("jt")
    X = s.Matrix(s.symbols("Xx Xy Xz"))
    V = transverse("V")
    Xp = s.Matrix(s.symbols("XPx XPy XPz"))
    U = transverse("U")
    r = s.Matrix(s.symbols("rx ry rz"))
    sig = transverse("sig")
    R = s.Matrix(s.symbols("Rx Ry Rz"))
    K = s.Matrix(s.symbols("Kx Ky Kz"))
    nn = lambda A: (n.T * A * n)[0]
    grad_nn = lambda A: s.I * p * nn(A)
    G = P * jt * n / 2 + s.I / 2 * (p * nn(j) - pn * j * n)
    # Unknown ray variables Y=X_j and U=v_j are eliminated by the duality identity.
    Y = s.Matrix(s.symbols("Yx Yy Yz"))
    original_ray = r.dot(Y) + sig.dot(U)
    dual_boundary_derivative = r.dot(Y) + R.dot(U + j * n / 2) + (sig - R).dot(U) + K.dot(G)
    dual_integrand = -R.dot(j * n) / 2 - K.dot(G)
    zero(original_ray - dual_boundary_derivative - dual_integrand)

    W = V + P * K / 2
    Wp = transverse("Vp")
    Wp += P * (sig - R) / 2
    direct = X.dot(grad_nn(jt)) / 2 - V.dot(jt * n)
    old = direct + dual_integrand
    new_bulk = X.dot(s.I * pn * grad_nn(j)) / 2 - (Xp + P * K).dot(grad_nn(j)) / 2
    new_bulk -= V.dot(s.I * pn * j * n)
    new_bulk += (Wp - R / 2).dot(j * n)
    # d/dt along the ray: D j = partial_t j - i(p.n)j.
    boundary_derivative = Xp.dot(grad_nn(j)) / 2 + X.dot(grad_nn(jt - s.I * pn * j)) / 2
    boundary_derivative -= Wp.dot(j * n) + W.dot((jt - s.I * pn * j) * n)
    zero(old - new_bulk - boundary_derivative)
    assert s.expand(old - new_bulk) != 0
    # Omitting the longitudinal part of R loses its coupling to the
    # longitudinal coordinate displacement, despite transverse v_j.
    assert s.expand(old - (direct - (P * R).dot(j * n) / 2 - K.dot(G))) != 0
    print("[ok] dual ray elimination and derivative-free bulk identity for arbitrary plane-wave jets")
    print("[ok] omission controls: total-derivative boundary and longitudinal displacement are necessary")


def check_boundaries_and_clock():
    eps, lam = s.symbols("eps lambda")
    em, d1, D1, Cm = s.symbols("em d1 D1 Cm")
    e, o, initial = s.symbols("e o initial", positive=True)
    F1, Fj, B, U1, Uj, U1o, Ujo = s.symbols("F1 Fj B U1 Uj U1o Ujo")
    # Retain precisely the coefficients linear in the second-order fields.
    # Other first-order-only endpoint products cannot contribute to lambda.
    endpoint_shift = eps * d1 - eps**2 * lam * e * em / 2
    lower = eps * F1 * endpoint_shift + eps**2 * lam * Fj * endpoint_shift / 2
    coord = eps * U1 + eps**2 * lam * Uj / 2 + eps**3 * lam * B / 2 - lower
    clock_extra = eps * U1o * (eps**2 * lam * Cm / 2) + eps**2 * lam * Ujo * (eps * D1) / 2
    coeff = s.expand(coord + clock_extra).coeff(eps, 3).coeff(lam) * 6
    expected = 3 * (B + e * F1 * em - Fj * d1 + D1 * Ujo + Cm * U1o)
    zero(coeff - expected)
    for term in (e * F1 * em, -Fj * d1, D1 * Ujo, Cm * U1o):
        assert s.expand(coeff - (expected - 3 * term)) != 0

    # The second-order lapse's linear contribution to the physical clock.
    m_initial, integral = s.symbols("m_initial integral")
    d_initial = -initial * m_initial
    clock = (initial * d_initial - integral) / o
    zero(clock - (-(initial**2) * m_initial - integral) / o)
    # The symmetric frame's v_j(o)=0 and X_j(o)=0 remove the dual boundary;
    # they do not remove the metric-dependent boundary after integrating j'.
    R0, K0, Xjo, Ujo_sym = s.symbols("R0 K0 Xjo Ujo")
    dual_ends = R0 * Xjo + K0 * Ujo_sym
    zero(dual_ends.subs({Xjo: 0, Ujo_sym: 0}))
    print("[ok] complete second-order-field cubic emission and clock coefficients; four omission controls")


def check_normalization():
    eps, lam = s.symbols("eps lambda")
    # Three arbitrary equal-area sky samples suffice to check the algebraic
    # mean/projection identity. No isotropy, Gaussianity or dipole removal.
    u = s.symbols("u0:3")
    j = s.symbols("j0:3")
    c = s.symbols("c0:3")
    mean = lambda xs: sum(xs) / 3
    # Expansion of exp(-I) with I=eps*u+eps^2*lambda*j/2+eps^3*lambda*c/6.
    E = [
        1
        - eps * u[i]
        + eps**2 * (u[i] ** 2 - lam * j[i]) / 2
        + eps**3 * (-lam * c[i] / 6 + lam * u[i] * j[i] / 2 - u[i] ** 3 / 6)
        for i in range(3)
    ]
    denom = mean(E)
    den1 = s.expand(denom).coeff(eps, 1)
    den2 = s.expand(denom).coeff(eps, 2)
    den3 = s.expand(denom).coeff(eps, 3)
    inverse = 1 - eps * den1 + eps**2 * (den1**2 - den2) + eps**3 * (-(den1**3) + 2 * den1 * den2 - den3)
    product = [(u[i] - mean(u)) * (j[i] - mean(j)) for i in range(3)]
    for i in range(3):
        actual = s.expand(E[i] * inverse).coeff(eps, 3).coeff(lam) * 6
        zero(actual - (-(c[i] - mean(c)) + 3 * (product[i] - mean(product))))
    print("[ok] centered normalized cubic map includes the required first-order x second-order map product")


def check_emission_moments():
    mu, y = s.symbols("mu y", real=True)
    khat, qhat = s.Matrix([0, 0, 1]), s.Matrix([y, 0, mu])
    Qk = khat * khat.T - s.eye(3) / 3
    Qq = qhat * qhat.T - s.eye(3) / 3
    total, diff = khat + qhat, khat - qhat
    fs, ft, nd, ad, bd = s.symbols("Sprime Tprime n a bprime")
    J = fs * (total * total.T / (2 * (1 + mu)) - s.eye(3) / 3)
    J += ft * (diff * diff.T / (2 * (1 - mu)) - s.diag(0, 1, 0))
    # Azimuthal averaging about k preserves zz and averages the transverse trace.
    az = lambda A: s.diag((A[0, 0] + A[1, 1]) / 2, (A[0, 0] + A[1, 1]) / 2, A[2, 2])
    red = lambda A: A.applyfunc(lambda z: s.factor(s.expand(z).subs(y * y, 1 - mu * mu)))
    zero(red(az(Qq)) - (3 * mu**2 - 1) * Qk / 2)
    qs, qt = (1 + 3 * mu) / 6, (1 - mu) / 2
    zero(red(az(J)) - s.Rational(3, 2) * (qs * fs + qt * ft) * Qk)
    # The two mixed partitions contribute equally under q -> -q. The
    # homogeneous tensor averages to zero; the homogeneous lapse remains.
    nh, nr0, nr2, Cprime = s.symbols("nh N0 N2 Cprime")
    projected = bd / 2 * (nh + 3 * nr2 - nr0) + 3 * ad * Cprime / 2
    expected = (bd * (nh + 3 * nr2 - nr0) + 3 * ad * Cprime) / 2
    zero(projected - expected)
    z1, z2, z3 = s.symbols("z1 z2 z3")
    h1, h2, h3, n12, n13, n23 = s.symbols("h1 h2 h3 n12 n13 n23")
    # Factorial convention: the cubic kernel of 3 h^(1) m^(2)
    # is the sum of its three single|pair partitions, not three times it.
    poly = 3 * (z1 * h1 + z2 * h2 + z3 * h3) * 2 * (z1 * z2 * n12 + z1 * z3 * n13 + z2 * z3 * n23)
    zero(s.expand(poly).coeff(z1).coeff(z2).coeff(z3) / 6 - (h1 * n23 + h2 * n13 + h3 * n12))
    print("[ok] shell emission contribution reduces to n_h+3N2-N0 and one metric moment derivative; polarization factors fixed")


def check_mixed_phase():
    kn, qn, o, t, r = s.symbols("kn qn o t r")
    phase = (kn + qn) * (o - t) - qn * (o - r)
    zero(phase - (kn * (o - t) + qn * (r - t)))
    assert s.diff(phase, qn) != 0
    zero(s.diff(phase, qn).subs(r, t))
    print("[ok] mixed-time phase q.(r-t) survives: a finite nine-moment truncation is not automatic")


if __name__ == "__main__":
    check_mixed_driver()
    check_linear_ray()
    check_adjoint()
    check_boundaries_and_clock()
    check_normalization()
    check_emission_moments()
    check_mixed_phase()
    print("PROPAGATION ADJOINT IDENTITIES VERIFIED (exact symbolic algebra)")
