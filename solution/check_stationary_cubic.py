"""Exact stationary cubic assembly checks, including a full null-energy cross-check.

Generic local metric jets prove ray identities at finite wavelength. A separate
homogeneous anisotropic OFF-SHELL metric has conserved canonical photon momentum,
so its full cubic normalized observable is computed independently from endpoint
energies. Polynomial time profiles there are exact functions, not a time-series
approximation and not an Einstein-radiation solution.
"""

import sympy as s
from stationary_cubic import (
    assemble_polarized_cubic,
    clock_first_only,
    comm,
    compose_first_only,
    cubic_field_integrand,
    direct_field_kinetic_cubic,
    direct_field_residuals,
    emission_first_only,
    first_only_driver_cubic,
    frame_exponent_boundary,
    initial_cubic_pairing,
    initial_diagonal_sources,
    linear_direction_driver,
    normalized_map_cubic,
    polarize_diagonal_cubic,
    reduced_cubic_bulk,
    second_ray_forcing,
)
from unified_adjoint_source import (
    FieldVariationCoefficients,
    FirstFieldAdjointSource,
    inner,
    measurement_source,
)


def zero(expr):
    if isinstance(expr, s.MatrixBase):
        assert expr.applyfunc(s.expand) == s.zeros(*expr.shape)
    else:
        assert s.expand(expr) == 0


def symmetric(name, tracefree=False):
    xx, yy, zz, xy, xz, yz = s.symbols(f"{name}xx {name}yy {name}zz {name}xy {name}xz {name}yz")
    return s.Matrix([[xx, xy, xz], [xy, yy, yz], [xz, yz, -xx - yy if tracefree else zz]])


def check_second_ray_jets():
    eps = s.Symbol("eps")
    n = s.Matrix([0, 0, 1])
    P = s.eye(3) - n * n.T
    G, hp = symmetric("G"), symmetric("hp", True)
    dg = [symmetric("d" + str(i)) for i in range(3)]
    ddriver = [s.Matrix(s.symbols(f"driver{i}_0:3")) for i in range(3)]
    X = s.Matrix(s.symbols("X0:3"))
    v = s.Matrix([*s.symbols("vx vy"), 0])
    cut = lambda M: M.applyfunc(lambda x: sum(s.expand(x).coeff(eps, k) * eps**k for k in range(3)))
    M = s.eye(3) + eps * G / 2 + eps**2 * G * G / 8
    Mi = s.eye(3) - eps * G / 2 + eps**2 * G * G / 8
    gi = s.eye(3) - eps * G + eps**2 * G * G / 2
    dgamma = [eps * d + eps**2 * (d * G + G * d) / 2 for d in dg]
    dM = [eps * d / 2 + eps**2 * (d * G + G * d) / 8 for d in dg]
    Mt = eps * hp / 2 + eps**2 * (hp * G + G * hp) / 8
    gt = eps * hp + eps**2 * (hp * G + G * hp) / 2
    S = cut(gi * gt / 2)
    Gam = [
        s.Matrix(
            3, 3, lambda i, j, l=l: sum(gi[l, r] * (dgamma[i][r, j] + dgamma[j][r, i] - dgamma[r][i, j]) / 2 for r in range(3))
        )
        for l in range(3)
    ]
    Gam = [cut(g) for g in Gam]
    vv = -n + eps * v + eps**2 * n * v.dot(v) / 2
    ee = cut(Mi * vv)
    F = cut(ee.T * gt * ee / 2)[0]
    coordinate = -2 * cut(S * ee) - s.Matrix([cut(ee.T * g * ee)[0] for g in Gam]) + cut(ee * F)
    spatial = sum((cut(ee[i] * dM[i] * ee) for i in range(3)), s.zeros(3, 1))
    exact = cut(Mt * ee) + spatial + cut(M * coordinate)
    raw2 = exact.applyfunc(lambda x: 2 * s.expand(x).coeff(eps, 2))
    pull = sum((2 * X[i] * ddriver[i] for i in range(3)), s.zeros(3, 1))
    QX, QV = second_ray_forcing(n, G, hp, dg, ddriver, X, v)
    zero(QV - P * (raw2 + pull))
    zero(n.dot(raw2) - 2 * v.dot(linear_direction_driver(n, hp, dg)))
    # Coordinate velocity contains the radial unit-direction contribution.
    xexact = ee.applyfunc(lambda x: 2 * s.expand(x).coeff(eps, 2))
    xexact += sum((X[i] * dg[i] * n for i in range(3)), s.zeros(3, 1))
    zero(QX - xexact)
    print("[ok] exact coordinate geodesic/metric jets give the complete transverse second-ray forcing")
    print("[ok] unit-direction radial jet and noncommuting frame-connection term retained")


def check_cubic_ray_duality():
    eps = s.Symbol("eps")
    n = s.Matrix([0, 0, 1])
    P = s.eye(3) - n * n.T
    G, hp = symmetric("G"), symmetric("hp", True)
    dh = [symmetric("dh" + str(i), True) for i in range(3)]
    ddh = [[symmetric("ddh" + str(min(i, j)) + str(max(i, j)), True) for j in range(3)] for i in range(3)]
    X, Y, R, K = (s.Matrix(s.symbols(name + "0:3")) for name in ("X", "Y", "R", "K"))
    v = s.Matrix([*s.symbols("vx vy"), 0])
    u = s.Matrix([*s.symbols("ux uy"), 0])
    r = s.Matrix([inner(n * n.T, d) / 2 for d in dh])
    sig = -P * hp * n
    pulled = eps * hp / 2 + eps**2 * sum((X[i] * dh[i] for i in range(3)), s.zeros(3)) / 2
    pulled += (
        eps**3
        * (
            sum((Y[i] * dh[i] for i in range(3)), s.zeros(3))
            + sum((X[i] * X[j] * ddh[i][j] for i in range(3) for j in range(3)), s.zeros(3))
        )
        / 4
    )
    pulled += eps**3 * comm(G, comm(G, hp)) / 48
    vv = -n + eps * v + eps**2 * (u + n * v.dot(v)) / 2
    raw = 6 * s.expand((vv.T * pulled * vv)[0]).coeff(eps, 3)
    known = first_only_driver_cubic(n, G, hp, dh, ddh, X, v)
    zero(raw - known - 3 * (r.dot(Y) + sig.dot(u)))
    QX, QV = s.Matrix(s.symbols("QX0:3")), s.Matrix([*s.symbols("QVx QVy"), 0])
    boundary_derivative = r.dot(Y) + R.dot(u + QX) + (sig - R).dot(u) + K.dot(QV)
    reduced = known - 3 * R.dot(QX) - 3 * K.dot(QV)
    zero(raw - reduced - 3 * boundary_derivative)
    assert s.expand(3 * v.dot(v) * inner(n * n.T, hp)) != 0
    assert s.expand(inner(n * n.T, comm(G, comm(G, hp)))) != 0
    print("[ok] complete cubic redshift driver and second-ray adjoint elimination for arbitrary spatial jets")


def check_field_kinetic_and_polarization():
    eps, t = s.symbols("eps t")
    H, D, E = s.symbols("H D E", commutative=False)
    c = lambda A, B: A * B - B * A
    S = eps * D / 2 + eps**2 * c(D, H) / 4 + eps**3 * c(c(D, H), H) / 12
    Sp = eps * E / 2 + eps**2 * c(E, H) / 4 + eps**3 * (c(c(E, H), H) + c(c(D, H), D)) / 12
    cut = lambda x: sum(s.expand(x).coeff(eps, k) * eps**k for k in range(4))
    plus = sum((eps * H / 2) ** k / s.factorial(k) for k in range(4))
    minus = sum((-eps * H / 2) ** k / s.factorial(k) for k in range(4))
    kinetic = cut(2 * plus * cut(Sp + 2 * S / t) * minus)
    actual = 6 * s.expand(kinetic).coeff(eps, 3)
    predicted = c(H, c(H, E)) / 4 + c(D, c(H, D)) + c(H, c(H, D)) / (2 * t)
    zero(actual - predicted)
    h, hp, hpp = symmetric("h", True), symmetric("hp", True), symmetric("hpp", True)
    assert any(s.expand(x) != 0 for x in direct_field_kinetic_cubic(h, hp, hpp, t))
    zero(s.trace(hp * comm(hp, h)))
    f, fp, fpp = s.symbols("f fp fpp")
    zero(direct_field_kinetic_cubic(f * h, fp * h, fpp * h, t))
    a, b, d = s.symbols("a b d")
    x0, x1, x2, l0, l1, l2 = s.symbols("x0 x1 x2 l0 l1 l2")
    linear = a * l0 + b * l1 + d * l2
    quadratic = (a * x0 + b * x1 + d * x2) ** 2
    result = polarize_diagonal_cubic(3 * linear * quadratic, [a, b, d])
    zero(result - (l0 * x1 * x2 + l1 * x0 * x2 + l2 * x0 * x1))
    Jy, Jyy, Jyyy, c1, c2, c3 = s.symbols("Jy Jyy Jyyy c1 c2 c3")
    assembled = Jyyy * c1**3 + initial_cubic_pairing(-Jy, c3, -Jyy * c1, c2)
    zero(assembled - (Jy * c3 + 3 * Jyy * c1 * c2 + Jyyy * c1**3))
    print("[ok] direct cubic field kinetic commutators; lapse tr(S^2) has no direct cubic term")
    print("[ok] diagonal polarization divided by 3! and nonlinear initial-data signs/factors")
    print("[ok] monochromatic first field h=f(t) B(X) kills every direct cubic temporal field commutator")
    labels = [a, b, d]
    values = [x0, x1, x2]
    g = sum(x * y for x, y in zip(labels, values, strict=True))
    Ns = s.symbols("N0:3")
    Cs = s.symbols("C0:3")
    Ls = [symmetric("La" + str(i), True) for i in range(3)]
    H0 = symmetric("H0", True)
    Q = symmetric("Q", True)
    R = symmetric("R", True)
    N0, C0, meas = s.symbols("Nback Cback meas")
    Nfirst = sum(x * y for x, y in zip(labels, Ns, strict=True))
    Cfirst = sum(x * y for x, y in zip(labels, Cs, strict=True))
    Hfirst = sum((labels[i] * Ls[i] for i in range(3)), s.zeros(3))
    diagonal = meas * g**3 + cubic_field_integrand(N0, H0, g**3, g**3 * R, Nfirst, Hfirst, g**2, g**2 * Q)
    diagonal += initial_cubic_pairing(C0, g**3, Cfirst, g**2)
    partitions = []
    initial_parts = []
    for i in range(3):
        j, k = [j for j in range(3) if j != i]
        pair = values[j] * values[k]
        partitions.append((Ns[i], Ls[i], pair, pair * Q))
        initial_parts.append((Cs[i], pair))
    xyz = x0 * x1 * x2
    assembled = assemble_polarized_cubic(meas * xyz, N0, H0, xyz, xyz * R, partitions, C0, xyz, initial_parts, t, 1).doit()
    zero(assembled - polarize_diagonal_cubic(diagonal, labels))
    print("[ok] executable arbitrary-triple assembly agrees with diagonal polarization, including initial-data pairings")


def check_endpoint_composition():
    eps, x, z = s.symbols("eps x z")
    e, i, o = s.symbols("e i o", positive=True)
    ae, ate, Xae = s.symbols("ae ate Xae")
    d1, d2 = emission_first_only(e, ae, ate, Xae)
    shift = eps * d1 + eps**2 * d2 / 2
    root = shift / e - shift**2 / (2 * e**2) + eps * (ae + shift * ate + eps * Xae)
    zero(s.expand(root).coeff(eps, 1))
    zero(s.expand(root).coeff(eps, 2))
    ai, api, ao, integral1, integral2 = s.symbols("ai api ao integral1 integral2")
    Di, Dii = emission_first_only(i, ai, api, 0)
    D1, D2 = clock_first_only(i, o, ai, api, ao, integral1, integral2)
    zero(s.expand(o * D1 - i * Di + integral1))
    zero(s.expand(o * D2 + D1**2 - i * Dii - Di**2 + integral2 + 2 * o * ao * D1 - 2 * i * ai * Di))
    d1, d2, D1, D2 = s.symbols("d1 d2 D1 D2")
    f1, ft, f2, fo, do = s.symbols("f1 ft f2 fo do")
    U1, U2, U3, Uo, Uoo, U2o = s.symbols("U1 U2 U3 Uo Uoo U2o")
    F = eps * (f1 + x * ft + z * fo) + eps**2 * f2 / 2
    lower = eps * (d1 + z * do) + eps**2 * d2 / 2
    raw = eps * (U1 + z * Uo + z * z * Uoo / 2) + eps**2 * (U2 + z * U2o) / 2 + eps**3 * U3 / 6
    raw -= s.integrate(F, (x, 0, lower))
    raw = s.expand(raw.subs(z, eps * D1 + eps**2 * D2 / 2))
    jets = compose_first_only(U1, U2, U3, f1, f2, ft, d1, d2, D1, D2, Uo, Uoo, U2o, fo, do)
    zero(2 * raw.coeff(eps, 2) - jets.second)
    zero(6 * raw.coeff(eps, 3) - jets.third)
    assert s.expand(6 * D1 * (fo * d1 + f1 * do)) != 0
    print(
        "[ok] emission/initial-temperature roots, simplified second proper-clock shift and complete mixed clock-emission composition"
    )


def check_direct_spatial_sources():
    import engine as eng

    t = eng.t
    vectors = [s.Matrix(s.symbols(name + "0:3")) for name in ("p", "q", "r")]
    eng.Ctx.maxdeg = 3
    eng.Ctx.momenta = {i + 1: v.flat() for i, v in enumerate(vectors)}
    eng.Ctx.trunc = s.expand
    z = eng.Field({frozenset([i]): s.Integer(1) for i in (1, 2, 3)})
    a = eng.Field({frozenset([i]): s.Function("nu" + str(i))(t) for i in (1, 2, 3)})
    N, H, c = direct_field_residuals(a, z, eng.mzero())
    lap = lambda f: sum((f.d_(i).d_(i) for i in range(3)), eng.Field())
    gz = [z.d_(i) for i in range(3)]
    ga = [a.d_(i) for i in range(3)]
    za = sum((gz[i] * ga[i] for i in range(3)), eng.Field())
    aa = sum((x * x for x in ga), eng.Field())
    expectedN3 = (-12 * z * z * lap(a) + 12 * z * za + 24 * z * aa - 24 * a * a * a / t**2) / 6
    expectedN2 = (4 * z * lap(a) - 2 * za - 4 * aa + 12 * a * a / t**2) / 2
    zero(N.coeff([1, 2, 3]) - expectedN3.coeff([1, 2, 3]))
    zero(N.coeff([1, 2]) - expectedN2.coeff([1, 2]))
    Hess = [[z.d_(i).d_(j) + 2 * a.d_(i).d_(j) for j in range(3)] for i in range(3)]
    quadratic = [[gz[i] * gz[j] + 2 * (gz[i] * ga[j] + gz[j] * ga[i]) + 2 * ga[i] * ga[j] for j in range(3)] for i in range(3)]
    expectedH3 = eng.tf(eng.madd(eng.mscale(Hess, -4 * z * z), eng.mscale(quadratic, -4 * z)))
    expectedH2 = eng.tf(eng.madd(eng.mscale(Hess, 4 * z), eng.mscale(quadratic, 2)))
    for i in range(3):
        for j in range(3):
            zero(H[i][j].coeff([1, 2, 3]) - expectedH3[i][j].coeff([1, 2, 3]))
            zero(H[i][j].coeff([1, 2]) - expectedH2[i][j].coeff([1, 2]))
    p, q, r = vectors
    c2 = (2 * (p.dot(p) + q.dot(q)) - p.dot(q)) / 9
    c3 = 2 * (-2 * (p.dot(p) + q.dot(q) + r.dot(r)) + p.dot(q) + p.dot(r) + q.dot(r)) / 9
    zero(c.coeff([1, 2]) - c2)
    zero(c.coeff([1, 2, 3]) - c3)
    diag2, diag3 = initial_diagonal_sources(z, gz, lap(z))
    zero(diag2.coeff([1, 2]) / 2 - c2)
    zero(diag3.coeff([1, 2, 3]) / 6 - c3)
    eng.Ctx.trunc = None
    print("[ok] physical-input second/third spatial residuals in conformal sector at three arbitrary momenta")
    print("[ok] nonlinear initial Hamiltonian coefficients and both distinct polarization conventions")


def homogeneous_measurement_check():
    t, O, eps = s.symbols("t O eps", positive=True)
    n = s.Matrix(s.symbols("nx ny nz", real=True))
    nx, ny, nz = n
    e, initial = s.Integer(2), s.Integer(1)
    sphere_relation = s.Poly(nz**2 + nx**2 + ny**2 - 1, nz)

    def unit(expr):
        return s.Poly(s.expand(expr), nz).rem(sphere_relation).as_expr()

    def avg(expr):
        poly = s.Poly(s.expand(expr), nx, ny, nz)
        value = 0
        for powers, coefficient in poly.terms():
            if any(k % 2 for k in powers):
                continue
            moment = s.prod(s.factorial2(k - 1) for k in powers) / s.factorial2(sum(powers) + 1)
            value += coefficient * moment
        return s.expand(value)

    def integral(expr, lo, hi):
        primitive = s.Poly(s.expand(expr), t).integrate().as_expr()
        return s.expand(primitive.subs(t, hi) - primitive.subs(t, lo))

    def mint(M, lo, hi):
        return M.applyfunc(lambda x: integral(x, lo, hi))

    A = s.diag(1, -1, 0)
    B = s.Matrix([[0, 1, 1], [1, 0, 0], [1, 0, 0]])
    C = s.Matrix([[0, 0, 1], [0, 1, 1], [1, 1, -1]])
    D = s.Matrix([[1, 0, 1], [0, -1, 0], [1, 0, 0]])
    h = t**2 * A / 3 + t**4 * B / 17
    hp = h.diff(t)
    a = t**2 / 7
    p = s.eye(3) - n * n.T
    v = (p * (h - h.subs(t, O)) * n / 2).applyfunc(unit)
    K = (-p * (h - h.subs(t, e)) * n).applyfunc(unit)
    grad = [s.zeros(3) for _ in range(3)]
    grad_driver = [s.zeros(3, 1) for _ in range(3)]
    grad2 = [[s.zeros(3) for _ in range(3)] for _ in range(3)]
    bulk = unit(reduced_cubic_bulk(n, h, hp, grad, grad_driver, grad, grad2, s.zeros(3, 1), v, s.zeros(3, 1), K))
    frame2 = mint(comm(h, hp) / 4, initial, O)
    frame = unit(frame_exponent_boundary(K.subs(t, O), frame2, n))
    U1 = inner(n * n.T, h.subs(t, O) - h.subs(t, e)) / 2
    F2 = unit(-2 * v.dot(hp * n))
    U2 = integral(F2, e, O)
    U3 = integral(bulk, e, O) + frame
    F1e = inner(n * n.T, hp.subs(t, e)) / 2
    d1, d2 = emission_first_only(e, a.subs(t, e), s.diff(a, t).subs(t, e), 0)
    D1, D2 = clock_first_only(
        initial,
        O,
        a.subs(t, initial),
        s.diff(a, t).subs(t, initial),
        a.subs(t, O),
        integral(t * a, initial, O),
        integral(t * a * a, initial, O),
    )
    jets = compose_first_only(
        U1,
        U2,
        U3,
        F1e,
        F2.subs(t, e),
        inner(n * n.T, hp.diff(t).subs(t, e)) / 2,
        d1,
        d2,
        D1,
        D2,
        s.diff(U1, O),
        s.diff(U1, O, 2),
        s.diff(U2, O),
        0,
        0,
    )
    w = nx * nz
    first_only = avg(w * unit(normalized_map_cubic(jets, avg))).subs(O, 3)
    # Every polynomial in t above is an exact off-shell profile, not a short-time
    # solution. Fixed-time endpoint energies give an independent observable.
    zero_matrix = s.zeros(3)

    def multiply(L, R):
        return [sum((L[j] * R[k - j] for j in range(k + 1)), zero_matrix) for k in range(4)]

    def exponential(H):
        out = [s.eye(3), zero_matrix, zero_matrix, zero_matrix]
        term = out
        for k in range(1, 4):
            term = [x / k for x in multiply(term, H)]
            out = [out[j] + term[j] for j in range(4)]
        return out

    def endpoint_answer(h2, h3, m2):
        # Temperature roots/proper clock expanded in amplitude only.
        de2 = d2 - e * m2.subs(t, e)
        do2 = D2 - (initial**2 * m2.subs(t, initial) + integral(t * m2, initial, O)) / O

        def at(T, shift1, shift2):
            return [
                zero_matrix,
                h.subs(t, T),
                (h2.subs(t, T) + 2 * shift1 * hp.subs(t, T)) / 2,
                (
                    h3.subs(t, T)
                    + 3 * shift1 * h2.diff(t).subs(t, T)
                    + 3 * shift2 * hp.subs(t, T)
                    + 3 * shift1**2 * hp.diff(t).subs(t, T)
                )
                / 6,
            ]

        Ho, He = at(O, D1, do2), at(e, d1, de2)
        R = [s.eye(3), zero_matrix, frame2 / 2, zero_matrix]
        Eo = exponential([x / 2 for x in Ho])
        Ee = exponential([-x for x in He])
        Q = multiply(multiply(multiply(multiply([x.T for x in R], Eo), Ee), Eo), R)
        q1, q2, q3 = [unit(inner(n * n.T, Q[j])) for j in (1, 2, 3)]
        z1 = -q1 / 2
        z2 = unit(-q2 / 2 + 3 * q1**2 / 8)
        z3 = unit(-q3 / 2 + 3 * q1 * q2 / 4 - 5 * q1**3 / 16)
        answer = 6 * avg(w * unit(z3 - avg(z1) * z2 + (avg(z1) ** 2 - avg(z2)) * z1))
        return s.simplify(answer.subs(O, 3))

    direct_first = endpoint_answer(s.zeros(3), s.zeros(3), s.Integer(0))
    zero(s.simplify(first_only - direct_first))
    assert avg(w * frame).subs(O, 3) != 0
    print("[ok] complete first-only cubic measurement equals independent conserved-photon-momentum endpoint energy")
    # Nonzero second/third metric and second lapse test the new measurement
    # Hessian's compatibility with the full independent cubic observable.
    h2 = t**2 * C / 5
    h3 = t**4 * D / 11
    m2 = -(t**2) / 13
    theta1 = -U1 + avg(U1)
    J1 = avg(w * theta1)
    lambdaI = (w * theta1 - J1) / (4 * s.pi)
    sourceM = measurement_source(
        s.zeros(3, 1),
        n,
        w / (4 * s.pi),
        s.zeros(3, 1),
        v,
        s.zeros(3, 1),
        K,
        lambdaI,
        D1,
        a.subs(t, e),
        F1e,
        avg(w * s.diff(U1, O)),
        e,
        O,
        initial,
        t,
    )
    source = FirstFieldAdjointSource(
        FieldVariationCoefficients(0, 0, s.zeros(3), s.zeros(3)), sourceM, t, e, O, initial, lambda x: 4 * s.pi * avg(unit(x))
    )
    second_coupling = source.action(m2, h2).doit().subs(O, 3)
    third_linear = -avg(w * inner(n * n.T, h3.subs(t, O) - h3.subs(t, e)) / 2).subs(O, 3)
    assembled = first_only - 3 * second_coupling + third_linear
    direct_full = endpoint_answer(h2, h3, m2)
    zero(s.simplify(assembled - direct_full))
    assert second_coupling != 0
    print("[ok] full cubic measurement including higher fields agrees with first-only + first-source + linear-third assembly")
    print("[ok] noncommuting FW rotation, both physical endpoints, second clock and exact sphere normalization included")
    print("      fixture is OFF SHELL; this is not an Einstein-radiation spectrum computation")


if __name__ == "__main__":
    check_second_ray_jets()
    check_cubic_ray_duality()
    check_field_kinetic_and_polarization()
    check_endpoint_composition()
    check_direct_spatial_sources()
    homogeneous_measurement_check()
    print("STATIONARY CUBIC ASSEMBLY CHECKS PASSED (exact algebra; remaining spectrum integrals unevaluated)")
