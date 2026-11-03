"""Exact checks for the complete first field-adjoint weak source.

Independent Christoffel-based field equations test arbitrary functions/momenta
OFF SHELL. Raw coordinate driver jets, moving endpoints, the proper clock and
normalized-map algebra test the measurement Hessian separately. A regular
zero-momentum Green example tests derivative-valued boundary sources and the
origin covector. No quadrature, time-series truncation or physical toy model.
"""

import engine as eng
import sympy as s
from unified_adjoint_source import (
    FieldVariationCoefficients,
    FirstFieldAdjointSource,
    MeasurementSource,
    field_coefficients,
    field_variation,
    frame_driver_ab,
    frame_terminal_covector,
    inner,
    measurement_source,
    stf,
)


def zero(expr):
    if isinstance(expr, s.MatrixBase):
        assert expr.applyfunc(s.expand) == s.zeros(*expr.shape)
    else:
        assert s.expand(expr) == 0


def symmetric(name, function=False):
    values = [
        s.Function(name + suffix)(eng.t) if function else s.Symbol(name + suffix) for suffix in ("xx", "yy", "xy", "xz", "yz")
    ]
    xx, yy, xy, xz, yz = values
    return s.Matrix([[xx, xy, xz], [xy, yy, yz], [xz, yz, -xx - yy]])


def check_field_engine():
    t = eng.t
    p = s.Matrix(s.symbols("p0:3", real=True))
    q = s.Matrix(s.symbols("q0:3", real=True))
    z = s.Symbol("z")
    a, m = s.Function("a")(t), s.Function("m")(t)
    h, j = symmetric("h", True), symmetric("j", True)
    hp, jp = h.diff(t), j.diff(t)
    eng.Ctx.maxdeg = 2
    eng.Ctx.momenta = {1: p.flat(), 2: q.flat()}
    eng.Ctx.trunc = s.expand
    nf = lambda x, y: eng.Field({frozenset([1]): x, frozenset([2]): y})
    H = [[nf(h[i, k], j[i, k]) for k in range(3)] for i in range(3)]
    eq = eng.equations(nf(a, m), nf(z, 0), H)
    half = eng.mexp(eng.mscale(H, s.Rational(1, 2)))
    minus = eng.mexp(eng.mscale(H, -s.Rational(1, 2)))
    transformed = eng.mscale(eng.mmul(eng.mmul(half, eq["T"]), minus), 2)
    exactN = eq["N"].coeff([1, 2])
    exactH = s.Matrix(3, 3, lambda i, k: transformed[i][k].coeff([1, 2]))
    CN, CH = field_variation(p, q, z, a, a.diff(t), h, hp, m, m.diff(t), j, jp, t)
    zero(CN - exactN)
    zero(CH - exactH)
    rawH = s.Matrix(3, 3, lambda i, k: 2 * eq["T"][i][k].coeff([1, 2]))
    assert any(s.expand(v) != 0 for v in (rawH - exactH))
    lh = symmetric("L")
    ln = s.Symbol("Lnu")
    coeff = field_coefficients(p, q, z, a, a.diff(t), h, hp, ln, lh, t)
    zero(coeff.action(m, m.diff(t), j, jp) - ln * exactN - inner(lh, exactH))
    print("[ok] arbitrary-time off-shell mixed field operator equals independent Christoffel engine at general p,q")
    print("[ok] all field transpose coefficients match the exact mixed residual; orthonormal conversion retained")

    # The frame connection, derived independently from the exact metric/frame ODE.
    S = eq["S"]
    halfprime = [[half[i][k].dt() for k in range(3)] for i in range(3)]
    Omega = eng.madd(eng.mmul(halfprime, minus), eng.mmul(eng.mmul(half, S), minus), s=-1)
    for label in (1, 2):
        zero(s.Matrix(3, 3, lambda i, k, label=label: Omega[i][k].coeff([label])))
    actual = s.Matrix(3, 3, lambda i, k: Omega[i][k].coeff([1, 2]))
    zero(actual - frame_driver_ab(h, hp, j, jp))
    zero(actual + actual.T)
    print("[ok] frame rotation driver has no linear term and its exact mixed term is ([h,j']+[j,h'])/8")
    eng.Ctx.trunc = None


def check_zero_and_soft_modes():
    t = s.Symbol("t", positive=True)
    p = s.Matrix(s.symbols("p0:3"))
    q = s.Matrix(s.symbols("q0:3"))
    h, hp, j, jp = (symmetric(n) for n in ("h", "hp", "j", "jp"))
    z, a, ap, m, mp = s.symbols("z a ap m mp")
    CN, CH = field_variation(p, s.zeros(3, 1), z, a, ap, h, hp, m, mp, s.zeros(3), s.zeros(3), t)
    zero(CN - 12 * a * m / t**2)
    zero(CH - 2 * mp * hp)
    # Constant primordial dilation: the spatial operator rescales by exp(-2z).
    CN, CH = field_variation(s.zeros(3, 1), q, z, 0, 0, s.zeros(3), s.zeros(3), m, mp, j, jp, t)
    B = stf(q.dot(q) * j - q * (j * q).T - (j * q) * q.T)
    Q = stf(q * q.T)
    zero(CN + 2 * z * q.dot(q) * m)
    zero(CH + 2 * z * B + 8 * z * Q * m)
    print("[ok] homogeneous STATE lapse probe retained; no momentum-projector singularity")
    print("[ok] finite-wavelength constant-seed dilation of the off-shell spatial operator")


def check_measurement_ray_and_endpoints():
    t, e, o, initial = s.symbols("t e o initial", positive=True)
    n = s.Matrix([0, 0, 1])
    P = s.eye(3) - n * n.T
    p = s.Matrix(s.symbols("p0:3"))
    q = s.Matrix(s.symbols("q0:3"))
    h, hp, j, jp = (symmetric(v) for v in ("h", "hp", "j", "jp"))
    X, Y, R, K = (s.Matrix(s.symbols(v + "0:3")) for v in ("X", "Y", "R", "K"))
    v = s.Matrix([*s.symbols("vx vy"), 0])
    u = s.Matrix([*s.symbols("ux uy"), 0])
    w, ae, Fe, Da, A, lia = s.symbols("w ae Fe Da A lia")
    src = measurement_source(q, n, w, X, v, R, K, lia, Da, ae, Fe, A, e, o, initial, t)
    phi = s.exp(s.I * q.dot(n) * (o - t))
    nn = lambda M: inner(n * n.T, M)
    r = s.I * p * nn(hp) / 2
    sig = -P * hp * n
    Gj = P * jp * n / 2 + s.I * (q * nn(j) - q.dot(n) * j * n) / 2
    direct = phi * (s.I * q.dot(X) * nn(jp) / 2 - v.dot(jp * n))
    ray_other = r.dot(Y) + sig.dot(u)
    total_derivative = r.dot(Y) + R.dot(u + phi * j * n / 2) + (sig - R).dot(u) + K.dot(phi * Gj)
    # The raw coordinate driver has the same bilinear local products, before
    # transforming either direction into its symmetric orthonormal chart.
    z = s.Symbol("z")
    sig_ab = (hp * j + h * jp + jp * h + j * hp) / 4 + z * jp
    ea = v + z * n + h * n / 2
    eb = u + j * n / 2
    rawlocal = nn(sig_ab) - ea.dot(jp * n) - eb.dot(hp * n)
    zero(rawlocal + v.dot(jp * n) + u.dot(hp * n))
    base = measurement_source(q, n, 1, X, v, R, K, 0, 0, 0, 0, 0, e, o, initial, t)
    zero(direct + ray_other - base.ray_bulk(j, jp) - total_derivative)
    norm_clock_bulk = phi * (lia / 2 + s.I * w * Da * q.dot(n) / 2) * nn(jp)
    zero(src.ray_bulk(j, jp) - w * base.ray_bulk(j, jp) - norm_clock_bulk)
    # The dual boundary [R.Y+K.u]_e^o vanishes by initial adjoint / observer ray
    # data, not by dropping the metric emission/observer boundary terms below.
    jpe, jpo = symmetric("jpe"), symmetric("jpo")
    me, mI, intm = s.symbols("me mI intm")
    phase_e = s.exp(s.I * q.dot(n) * (o - e))
    expected = w * (e * Fe * me * phase_e + e * ae * nn(jpe) * phase_e / 2 + Da * nn(jpo) / 2)
    zero(src.sky_endpoints(me, jpe, jpo) - expected)
    Cm = -(initial**2 * mI + intm) / o
    zero(src.nu_clock_initial * mI + src.nu_clock_bulk / t * intm - A * Cm)
    for term in (src.nu_emission * me, inner(src.H_prime_emission, jpe), inner(src.H_prime_observer, jpo)):
        assert s.simplify(term) != 0
    print("[ok] complete measurement ray source equals raw coordinate driver modulo the exact dual total derivative")
    print("[ok] normalization, observer phase derivative, both emission terms and clock initial event retained")


def check_normalized_hessian_and_clock():
    alpha, beta = s.symbols("alpha beta")
    U = s.symbols("U0:3")
    V = s.symbols("V0:3")
    B = s.symbols("B0:3")
    ah = s.symbols("ah0:3")
    mb = s.symbols("mb0:3")
    Fh = s.symbols("Fh0:3")
    Fj = s.symbols("Fj0:3")
    Uo = s.symbols("Uo0:3")
    Vo = s.symbols("Vo0:3")
    Da, Db, e = s.symbols("Da Db e")
    w = [1, -2, 1]
    Iab = [B[i] + e * Fh[i] * mb[i] + e * Fj[i] * ah[i] + Da * Vo[i] + Db * Uo[i] for i in range(3)]
    # Independently expand moving lower limit / coordinate-observer composition.
    for i in range(3):
        delta_e = e * (1 - alpha * ah[i] - beta * mb[i]) - e
        raw = alpha * (U[i] + beta * Db * Uo[i]) + beta * (V[i] + alpha * Da * Vo[i]) + alpha * beta * B[i]
        raw -= (alpha * Fh[i] + beta * Fj[i]) * delta_e
        zero(s.expand(raw).coeff(alpha, 1).coeff(beta, 1) - Iab[i])
    I = [alpha * U[i] + beta * V[i] + alpha * beta * Iab[i] for i in range(3)]
    Z = [1 - x + x * x / 2 for x in I]
    d = sum(Z) / 3 - 1
    inverse = 1 - d + d * d
    J = sum(w[i] * (Z[i] * inverse - 1) for i in range(3))
    raw_hessian = s.expand(J).coeff(alpha, 1).coeff(beta, 1)
    theta_a = [-u + sum(U) / 3 for u in U]
    Ja = sum(w[i] * theta_a[i] for i in range(3))
    lia = [w[i] * theta_a[i] - Ja / 3 for i in range(3)]
    source = sum(w[i] * Iab[i] + lia[i] * V[i] for i in range(3))
    zero(source + raw_hessian)
    assert s.expand(source - sum(w[i] * Iab[i] for i in range(3))) != 0
    initial, o, aI, mI, inta, intm = s.symbols("initial o aI mI inta intm")
    dIa, dIb = -initial * aI, -initial * mI
    ca = (initial * dIa - inta) / o
    cb = (initial * dIb - intm) / o
    zero(ca - (-(initial**2) * aI - inta) / o)
    zero(cb - (-(initial**2) * mI - intm) / o)
    print("[ok] assembled measurement source is minus the independently expanded normalized physical Hessian")
    print("[ok] physical-temperature roots and proper-time clock fix both clock shifts including T_I")


def check_frame_covector():
    n = s.Matrix([0, 0, 1])
    gradient = s.Matrix([*s.symbols("gx gy"), 0])
    w = s.Symbol("w")
    x, y, z = s.symbols("rx ry rz")
    rotation = s.Matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    K = frame_terminal_covector(w, gradient, n)
    zero(inner(K, rotation) - w * gradient.dot(rotation * n))
    zero(K + K.T)
    # Polarized cubic composition: one first sky times each complementary
    # second rotation, not an extra factor three for each partition.
    alpha, beta, gamma = s.symbols("alpha beta gamma")
    grads = [gradient, s.Matrix([*s.symbols("hx hy"), 0]), s.Matrix([*s.symbols("jx jy"), 0])]
    h, hp, j, jp, k, kp = (symmetric(v) for v in ("h", "hp", "j", "jp", "k", "kp"))
    rotations = [frame_driver_ab(j, jp, k, kp), frame_driver_ab(h, hp, k, kp), frame_driver_ab(h, hp, j, jp)]
    first = alpha * grads[0] + beta * grads[1] + gamma * grads[2]
    second = beta * gamma * rotations[0] + alpha * gamma * rotations[1] + alpha * beta * rotations[2]
    direct = s.expand(w * first.dot(second * n)).coeff(alpha, 1).coeff(beta, 1).coeff(gamma, 1)
    predicted = sum(inner(frame_terminal_covector(w, g, n), r) for g, r in zip(grads, rotations, strict=True))
    zero(direct - predicted)
    assert s.expand(predicted) != 0
    print("[ok] first frame adjoint and polarized cubic rotation composition, including single-pair factors")


def check_boundary_green_solution():
    # Exact homogeneous-momentum field Green map. This is an operator test,
    # not a time-series approximation or a separate physical-universe model.
    t, r = s.symbols("t r", positive=True)
    o, e, initial = s.Integer(4), s.Integer(2), s.Integer(1)
    gN = (t**3 - r**3) / (9 * t * r)
    gH = r * (t - r) / t
    N0, N1 = t, 1 + t
    H0, H1 = 1 + t, 2 - t
    bnu, bhe, bho, A = s.symbols("bnu bhe bho A")
    # Split the source domain at its actual event points.
    ln_bulk = s.integrate(N0 * gN + N1 * s.diff(gN, t), (t, r, o))
    lh_bulk = s.integrate(H0 * gH + H1 * s.diff(gH, t), (t, r, o))
    clock_low = s.integrate(-A * t * gN / o, (t, initial, o)) - A * initial**2 * gN.subs(t, initial) / o
    clock_high = s.integrate(-A * t * gN / o, (t, r, o))
    ln_low = ln_bulk + bnu * gN.subs(t, e) + clock_low
    ln_mid = ln_bulk + bnu * gN.subs(t, e) + clock_high
    ln_high = ln_bulk + clock_high
    lh_low = lh_bulk + bhe * s.diff(gH, t).subs(t, e) + bho * s.diff(gH, t).subs(t, o)
    lh_high = lh_bulk + bho * s.diff(gH, t).subs(t, o)
    c = s.integrate(N0 * t**2 + N1 * 2 * t, (t, 0, o)) + bnu * e**2
    c += s.integrate(-A * t * t**2 / o, (t, initial, o)) - A * initial**4 / o
    zero(s.simplify(9 * s.limit(r * ln_low, r, 0) - c))
    zero(s.simplify((lh_low - lh_high).subs(r, e) - bhe))
    zero(s.simplify(lh_high.subs(r, o) - bho))
    # Pair with arbitrary regular test states through their exact residuals.
    nu = t**2 + 2 * t**4
    h = 3 * t**2 - t**4
    fN = 3 * s.diff(nu, t, 2) - 6 * nu / t**2
    fH = s.diff(h, t, 2) + 2 / t * s.diff(h, t)
    dualN = (
        sum(
            s.integrate(s.simplify(expr * fN.subs(t, r)), (r, lo, hi))
            for expr, lo, hi in ((ln_low, 0, initial), (ln_mid, initial, e), (ln_high, e, o))
        )
        + c
    )
    dualH = s.integrate(s.simplify(lh_low * fH.subs(t, r)), (r, 0, e)) + s.integrate(
        s.simplify(lh_high * fH.subs(t, r)), (r, e, o)
    )
    directN = s.integrate(N0 * nu + N1 * s.diff(nu, t), (t, 0, o)) + bnu * nu.subs(t, e)
    directN += s.integrate(-A * t * nu / o, (t, initial, o)) - A * initial**2 * nu.subs(t, initial) / o
    directH = s.integrate(H0 * h + H1 * s.diff(h, t), (t, 0, o)) + bhe * s.diff(h, t).subs(t, e) + bho * s.diff(h, t).subs(t, o)
    zero(s.simplify(dualN - directN))
    zero(s.simplify(dualH - directH))
    # Exercise the actual complete-source assembly API. A nonunit formal sky
    # measure catches accidental reintegration of the already contracted clock.
    D = s.diag(1, -1, 0)
    norm = inner(D, D)
    field = FieldVariationCoefficients(-N0, -N1, -H0 * D / norm, -H1 * D / norm)
    measure = MeasurementSource(
        s.zeros(3), s.zeros(3), bnu / 3, bhe * D / (3 * norm), bho * D / (3 * norm), -A * t / o, -A * initial**2 / o
    )
    source = FirstFieldAdjointSource(field, measure, t, e, o, initial, lambda expr: 3 * expr)
    zero(s.simplify(source.action(nu, h * D).doit() - directN - directH))
    zero(s.simplify(source.action(t**2, s.zeros(3)).doit() - c))
    for point, ln, lh in (
        (s.Rational(1, 2), ln_low, lh_low),
        (s.Rational(3, 2), ln_mid, lh_low),
        (s.Integer(3), ln_high, lh_high),
    ):
        causalN = s.Piecewise((gN.subs(r, point), t >= point), (0, True))
        causalH = s.Piecewise((gH.subs(r, point), t >= point), (0, True))
        zero(s.simplify(source.action(causalN, s.zeros(3)).doit() - ln.subs(r, point)))
        zero(s.simplify(source.action(0, causalH * D).doit() - lh.subs(r, point)))
    print("[ok] complete source assembly API and causal Green action; clock is not integrated over sky twice")
    print("[ok] weak-source Green solve including derivative endpoint evaluations and physical clock point source")
    print(
        "[ok] first spatial adjoint has a nonzero observer terminal value and an emission jump; origin covector still equals 9 lim(r Lambda_N)"
    )


if __name__ == "__main__":
    check_field_engine()
    check_zero_and_soft_modes()
    check_measurement_ray_and_endpoints()
    check_normalized_hessian_and_clock()
    check_frame_covector()
    check_boundary_green_solution()
    print("FIRST FIELD-ADJOINT SOURCE VERIFIED (exact algebra; full cubic spectrum not evaluated)")
