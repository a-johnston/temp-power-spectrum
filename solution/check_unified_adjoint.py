"""Exact checks for a unified stationary third-response construction.

Checks a polarized nonlinear implicit system, the exact normalization adjoint,
and a boundary/clock/redshift/normalization fixture. The fixture is an algebraic
verification model, NOT an Einstein-radiation solution or a physical approximation.
No quadrature, numerical evolution, floating point, or time-series evidence.
"""

from itertools import combinations

import sympy as s


def zero(expr):
    if isinstance(expr, s.MatrixBase):
        assert expr.applyfunc(s.simplify) == s.zeros(*expr.shape)
    else:
        assert s.simplify(expr) == 0


def check_polarized_identity():
    y = s.Matrix(s.symbols("y1 y2"))
    x = s.Matrix(s.symbols("x1 x2 x3"))
    u, v = y
    a, b, c = x
    F = s.Matrix(
        [
            2 * u + v + a + 2 * b - c + u**2 + 2 * u * v + a * v - b * u + b * c + u**3 + a * u * v + c**2 * v + a * b * c,
            -u + 3 * v + 2 * a - b + c + 2 * v**2 - u * v + c * u + a * b + v**3 - b * u**2 + a * c * v + 2 * a * b * c,
        ]
    )
    J = u - 2 * v + a * u + b * v + u * v + b * c + u**3 / 3 + u * v * c + a * b * c
    origin = dict.fromkeys([*y, *x], s.Integer(0))
    at0 = lambda z: z.subs(origin)
    A = at0(F.jacobian(y))
    assert A != A.T
    tangent = [A.inv() * (-at0(F.diff(z))) for z in x]
    lam0 = -A.T.inv() * at0(s.Matrix([J]).jacobian(y).T)
    L0 = J + (lam0.T * F)[0]

    def D(expr, index):
        result = s.diff(expr, x[index])
        for i in range(2):
            result += tangent[index][i] * s.diff(expr, y[i])
        return result

    gradL = s.Matrix([L0]).jacobian(y).T
    lam1 = [-A.T.inv() * at0(D(gradL, i)) for i in range(3)]
    cubic = at0(D(D(D(L0, 0), 1), 2))
    terms = []
    for i in range(3):
        j, k = [q for q in range(3) if q != i]
        terms.append((lam1[i].T * at0(D(D(F, j), k)))[0])
    cubic += sum(terms)

    # Independent recursive substitution in the actual implicit equations.
    # Distinct-label coefficients are polarized derivatives, with no factorial.
    state = s.zeros(2, 1)
    for degree in (1, 2, 3):
        for labels in combinations(range(3), degree):
            substituted = F.subs(dict(zip(y.flat(), state.flat(), strict=True)), simultaneous=True)
            monomial = s.prod(x[i] for i in labels)
            coeff = []
            for expr in substituted:
                value = s.expand(expr)
                for i in labels:
                    value = value.coeff(x[i], 1)
                coeff.append(value.subs(dict.fromkeys(x.flat(), s.Integer(0))))
            state += -A.inv() * s.Matrix(coeff) * monomial
    direct = s.expand(J.subs(dict(zip(y.flat(), state.flat(), strict=True)), simultaneous=True))
    direct = direct.coeff(a, 1).coeff(b, 1).coeff(c, 1)
    zero(cubic - direct)
    for term in terms:
        assert term != 0
        assert s.simplify(cubic - term - direct) != 0
    # At second order no first adjoint response is needed.
    second = at0(D(D(L0, 0), 1))
    direct2 = (
        s.expand(J.subs(dict(zip(y.flat(), state.flat(), strict=True)), simultaneous=True)).coeff(a, 1).coeff(b, 1).subs(c, 0)
    )
    zero(second - direct2)
    print(
        "[ok] polarized stationary cubic identity: nonsymmetric implicit system, nonlinear objective, explicit input dependence"
    )
    print("[ok] each first-adjoint term is necessary; second and third state responses eliminated")


def check_normalization_adjoint():
    # Finite exact weighted measure: the same identity holds on the sphere
    # with q_i replaced by dOmega/(4 pi). No angular quadrature is involved.
    q = s.Matrix([s.Rational(1, 6), s.Rational(1, 3), s.Rational(1, 2)])
    w = s.Matrix([1, -2, 1])  # objective weights, sum w_i=0
    Z = s.Matrix(s.symbols("Z0:3", positive=True))
    M = (q.T * Z)[0]
    theta = Z / M - s.ones(3, 1)
    J = (w.T * theta)[0]
    lam = (Z / M).multiply_elementwise(w - J * q)
    derivative = s.Matrix([s.diff(J, z) for z in Z])
    zero(derivative.multiply_elementwise(-Z) + lam)
    zero(sum(lam))

    eps = s.Symbol("eps")
    I1 = s.Matrix(s.symbols("I10:13"))
    Z1 = -I1
    M1 = (q.T * Z1)[0]
    theta1 = Z1 - M1 * s.ones(3, 1)
    J1 = (w.T * theta1)[0]
    lam_series = lam.subs(dict(zip(Z.flat(), (s.ones(3, 1) + eps * Z1).flat(), strict=True)), simultaneous=True)
    lam_first = lam_series.diff(eps).subs(eps, 0)
    zero(lam_first - (w.multiply_elementwise(theta1) - J1 * q))
    print("[ok] exact normalized-map adjoint is mean-free; its first variation has no uncentered random monopole")


def check_ray_adjoint_block():
    w = s.Symbol("w")
    n = s.Matrix([0, 0, 1])
    P = s.eye(3) - n * n.T
    r = s.Matrix(s.symbols("rx ry rz"))
    sig = s.Matrix([*s.symbols("sx sy"), 0])
    R = s.Matrix(s.symbols("Rx Ry Rz"))
    K = s.Matrix(s.symbols("Kx Ky Kz"))
    lamX = -w * R
    lamV = -w * P * K
    lamXp = -w * r
    lamVp = -w * P * (sig - R)
    zero(-w * r - lamXp)
    zero(-w * sig - P * lamX - lamVp)
    zero(n.dot(lamV))
    print("[ok] unified first ray-adjoint block reproduces the existing R,K propagation adjoint")


def check_field_transpose():
    t, r, T, c, k = s.symbols("t r T c k", positive=True)
    p = s.Matrix(s.symbols("p0:3"))
    matrices = []
    for name in ("U", "V"):
        xx, yy, xy, xz, yz = s.symbols(f"{name}xx {name}yy {name}xy {name}xz {name}yz")
        matrices.append(s.Matrix([[xx, xy, xz], [xy, yy, yz], [xz, yz, -xx - yy]]))

    def cartesian(M):
        value = p.dot(p) * M - p * (M * p).T - (M * p) * p.T
        return value - s.eye(3) * s.trace(value) / 3

    zero(s.trace(matrices[0] * cartesian(matrices[1])) - s.trace(cartesian(matrices[0]) * matrices[1]))
    u1 = lambda x: s.sin(x) / x - s.cos(x)
    u2 = lambda x: -s.cos(x) / x - s.sin(x)
    GN = (u2(c * t) * u1(c * r) - u1(c * t) * u2(c * r)) / (3 * c)
    zero(3 * s.diff(GN, t, 2) + (3 * c**2 - 6 / t**2) * GN)
    zero(GN.subs(t, r))
    zero(s.diff(GN, t).subs(t, r) - s.Rational(1, 3))
    vk = 3 * u1(c * t) / c**2
    zero(9 * s.limit(r * GN, r, 0) - vk)
    gh = r / t * s.sin(c * (t - r)) / c
    zero(s.diff(gh, r, 2) - 2 / r * s.diff(gh, r) + (2 / r**2 + c**2) * gh)
    zero(gh.subs(r, t))
    zero(s.diff(gh, r).subs(r, t) + 1)

    # A coupled scalar-sector example of the Cartesian field operator.
    # Lambda_N has the physical 1/t adjoint behavior at the regular origin.
    b, u, h, v = s.symbols("b u h v")
    nu = b * t**2 + u * t**4
    H = h * t**2 + v * t**4
    lN = (T - t) ** 2 / t
    lH = t * (T - t) ** 2
    LN = lambda f: 3 * s.diff(f, t, 2) + (k**2 - 6 / t**2) * f
    LH = lambda f: s.diff(f, t, 2) + 2 / t * s.diff(f, t) - k**2 * f / 3
    LHstar = lambda f: s.diff(f, t, 2) - 2 / t * s.diff(f, t) + (2 / t**2 - k**2 / 3) * f
    fN = LN(nu)
    fH = LH(H) + 4 * k**2 * nu
    pN = LN(lN) + 4 * k**2 * lH
    pH = LHstar(lH)
    difference = s.integrate(s.expand(pN * nu + pH * H - lN * fN - lH * fH), (t, 0, T))
    boundary = 9 * b * s.limit(t * lN, t, 0)
    zero(difference - boundary)
    assert s.simplify(difference) != 0
    x = k * t / s.sqrt(3)
    physical_nu = s.cos(x) - s.sin(x) / x
    physical_H = -6 * (2 * s.sin(x) / x - s.cos(x) - 1)
    zero(LN(physical_nu))
    zero(LH(physical_H) + 4 * k**2 * physical_nu + 2 * k**2)
    zero(s.limit(physical_nu / t**2, t, 0) + k**2 / 9)
    print("[ok] regular field Green transpose, spatial adjoint and diagonal jumps")
    print("[ok] exact physical linear field branch agrees with Cartesian coupling and initial-source normalization")
    print("[ok] regular-origin boundary covector Lambda_c = 9 lim(t Lambda_N); omission changes the coupled response")


def check_curvature_limit():
    # The established physical leading-gradient map, exact in zeta.
    # A and B are its trace-free Hessian and gradient-square sky contractions.
    eps, f = s.symbols("eps f")
    g, A, B, tau2 = s.symbols("g A B tau2")
    F = f * s.exp(eps * g) - 1
    J = -tau2 * f**2 * (eps * A - eps**2 * B)
    y1 = -g
    L0 = J  # lambda_0=0 since J_f vanishes at the background.
    D = lambda expr: s.diff(expr, eps) + y1 * s.diff(expr, f)
    at0 = lambda expr: s.simplify(expr.subs({eps: 0, f: 1}))
    lam1 = -at0(D(s.diff(L0, f)))
    predicted = at0(D(D(D(L0)))) + 3 * lam1 * at0(D(D(F)))
    direct = s.diff(J.subs(f, s.exp(-eps * g)), eps, 3).subs(eps, 0)
    expected = -12 * tau2 * (g**2 * A + g * B)
    zero(predicted - direct)
    zero(predicted - expected)
    assert s.simplify(predicted - at0(D(D(D(L0))))) != 0
    print("[ok] physical leading-gradient curvature cubic map recovered without f2 or f3; nonlinear seed term is essential")


def coefficient(expr, eps, degree):
    return s.expand(expr).coeff(eps, degree)


def truncate(expr, eps):
    return sum(coefficient(expr, eps, j) * eps**j for j in range(4))


def check_boundary_fixture():
    eps = s.Symbol("eps")
    initial, emission, observation = s.Integer(1), s.Integer(2), s.Integer(3)
    proper = (observation**2 - initial**2) / 2
    co = s.Rational(2, 5)
    ce = [s.Rational(1, 3), -s.Rational(2, 7), s.Rational(3, 4)]
    a = [s.Rational(2, 3), -s.Rational(3, 5), s.Rational(4, 7)]
    b = [s.Rational(1, 5), s.Rational(2, 9), -s.Rational(1, 4)]
    w = s.Matrix([1, -2, 1])
    ti, to = s.symbols("ti to")
    te = s.Matrix(s.symbols("te0:3"))
    I = s.Matrix(s.symbols("I0:3"))
    Z = s.Matrix(s.symbols("Z0:3"))
    M = s.Symbol("M")
    theta = s.Matrix(s.symbols("theta0:3"))
    state = s.Matrix([ti, to, *te, *I, *Z, M, *theta])
    E = [ti + eps * co * ti**2 - initial, (to**2 - ti**2) / 2 + eps * co * (to**3 - ti**3) / 3 - proper]
    E += [te[i] + eps * ce[i] * te[i] ** 2 - emission for i in range(3)]
    E += [I[i] - eps * a[i] * (to - te[i]) - eps**2 * b[i] * (to**2 - te[i] ** 2) / 2 for i in range(3)]
    E += [Z[i] - s.exp(-I[i]) for i in range(3)]
    E += [M - sum(Z.flat()) / 3]
    E += [M * (1 + theta[i]) - Z[i] for i in range(3)]
    F = s.Matrix(E)
    state0 = s.Matrix([initial, observation, *([emission] * 3), *([0] * 3), *([1] * 3), 1, *([0] * 3)])
    origin = dict(zip(state.flat(), state0.flat(), strict=True)) | {eps: s.Integer(0)}
    zero(F.subs(origin))
    A = F.jacobian(state).subs(origin)
    Jy = s.Matrix([0] * (len(state) - 3) + w.flat())
    y1 = -A.inv() * F.diff(eps).subs(origin)
    lam0 = -A.T.inv() * Jy
    linear_path = dict(zip(state.flat(), (state0 + eps * y1).flat(), strict=True))
    A1 = F.jacobian(state).subs(linear_path, simultaneous=True).diff(eps).subs(eps, 0)
    lam1 = -A.T.inv() * A1.T * lam0
    path = F.subs(linear_path, simultaneous=True)
    F2 = path.diff(eps, 2).subs(eps, 0)
    F3 = path.diff(eps, 3).subs(eps, 0)
    cubic = (lam0.T * F3)[0] + 3 * (lam1.T * F2)[0]

    # Independent exact amplitude expansion of the physical-form roots,
    # clock integral, redshift exponential and random-mean normalization.
    root = lambda value, c: value - c * value**2 * eps + 2 * c**2 * value**3 * eps**2 - 5 * c**3 * value**4 * eps**3
    ti_series = root(initial, co)
    to_series = observation
    for degree in (1, 2, 3):
        unknown = s.Symbol(f"o{degree}")
        trial = to_series + unknown * eps**degree
        residual = (trial**2 - ti_series**2) / 2 + eps * co * (trial**3 - ti_series**3) / 3 - proper
        answer = s.solve(coefficient(residual, eps, degree), unknown)[0]
        to_series += answer * eps**degree
    te_series = [root(emission, c) for c in ce]
    I_series = [
        truncate(eps * a[i] * (to_series - te_series[i]) + eps**2 * b[i] * (to_series**2 - te_series[i] ** 2) / 2, eps)
        for i in range(3)
    ]
    Z_series = [truncate(1 - z + z**2 / 2 - z**3 / 6, eps) for z in I_series]
    M_series = sum(Z_series) / 3
    d = M_series - 1
    inverse = truncate(1 - d + d**2 - d**3, eps)
    theta_series = s.Matrix([truncate(z * inverse - 1, eps) for z in Z_series])
    J_series = (w.T * theta_series)[0]
    zero(cubic - 6 * coefficient(J_series, eps, 3))
    zero((lam0.T * F2)[0] - 2 * coefficient(J_series, eps, 2))
    zero((Jy.T * y1)[0] - coefficient(J_series, eps, 1))

    # The first normalization adjoint matches the full block solve.
    theta1 = y1[-3:, :]
    J1 = (w.T * theta1)[0]
    zero(lam0[5:8, :] - w)  # rows of I - integral = 0
    zero(lam1[5:8, :] - (w.multiply_elementwise(theta1) - s.ones(3, 1) * J1 / 3))
    # Proper-clock and emission adjoint signs, with the fixture's unlogged roots.
    zero(lam1[1] - sum(w[i] * a[i] for i in range(3)) / observation)
    zero(lam1[0] - initial * lam1[1])
    for i in range(3):
        zero(lam1[2 + i] + w[i] * a[i])
    assert (lam1.T * F2)[0] != 0
    assert s.simplify((lam0.T * F3)[0] - 6 * coefficient(J_series, eps, 3)) != 0
    # Omitting physical boundary or random-mean blocks changes the answer.
    for indices in ([0, 1], list(range(2, 5)), list(range(11, 15))):
        omitted = sum(lam0[i] * F3[i] + 3 * lam1[i] * F2[i] for i in indices)
        assert s.simplify(omitted) != 0
    print("[ok] boundary fixture: temperature roots, proper clock, integrated redshift and random-mean normalization")
    print("[ok] full stationary cubic response equals independent exact root/clock/map amplitude expansion")
    print("[ok] clock, emission and normalization omission controls fail as required")


if __name__ == "__main__":
    check_polarized_identity()
    check_normalization_adjoint()
    check_ray_adjoint_block()
    check_field_transpose()
    check_curvature_limit()
    check_boundary_fixture()
    print("UNIFIED ADJOINT IDENTITIES VERIFIED (exact algebra; fixture is not a radiation solution)")
