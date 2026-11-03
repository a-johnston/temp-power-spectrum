"""Complete cubic quadrupole through four spatial gradients.

Default: shell contraction and two noncoplanar triples, all five components.
--symbolic: prove the uncontracted identity at six independent momentum
components, covering arbitrary triples by rotation and polynomial continuity.

This is an exact gradient-coefficient comparison with an independent covariant
curvature calculation. It is NOT an arbitrary-triple finite-wavelength result.
The observable_fast module is used only for its rational polynomial arithmetic,
not for its geometry, geodesic, root, clock, or observable implementations.
"""

import sys
import time
from functools import lru_cache

import check_cubic_curvature_hand as hand
import engine as eng
import observable_fast as alg
import sympy as s
from quadrupole_adjoint_expansion import Q, angular_A, first_cubic_correction, leading_cubic_source
from stationary_cubic import direct_field_residuals
from unified_adjoint_source import inner, ricci_linear, stf

start = time.monotonic()


def log(message):
    print(f"[{time.monotonic() - start:.1f}s] {message}", flush=True)


def qaverage(expr, mu, v):
    expr = s.Poly(s.expand(expr), v).rem(s.Poly(v * v + mu * mu - 1, v)).as_expr().expand()
    assert s.diff(expr, v) == 0, "axisymmetric scalar must not retain odd transverse momentum"
    return s.expand(sum(cf / (m + 1) for (m,), cf in s.Poly(expr, mu).terms() if m % 2 == 0))


@lru_cache(None)
def first_only_U3(wavevalues=None):
    """Leading-gradient specialization of reduced_cubic_bulk, not a ray solve."""
    A = alg
    A.KMAX = 4
    A.Ctx.maxdeg = 3
    zero, one = A.ZERO, A.ONE
    k = [one, zero, zero]
    q = [A.B1, A.B2, zero]
    # The old ring's names ip,u,v are independent momentum components here,
    # NOT inverse momenta. Do not use observable_fast.to_sympy to convert them.
    wave = [k, q, [-z for z in q]] if wavevalues is None else [[A.R.from_expr(s.sympify(z)) for z in p] for p in wavevalues]
    A.Ctx.momenta = {i + 1: tuple(A.KAP * z for z in p) for i, p in enumerate(wave)}
    F = A.Fld
    n = A.NV
    T, O, E = A.T, A.TO, A.TS
    P = [[((one if i == j else zero) - n[i] * n[j]) for j in range(3)] for i in range(3)]
    dot = lambda u, v: sum((u[i] * v[i] for i in range(3)), F())
    matvec = lambda M, v: [sum((M[i][j] * v[j] for j in range(3)), F()) for i in range(3)]
    integral = lambda f: f.map(lambda z: A.integ(z, A._IT, E, T))
    amp = F({frozenset([i]): one for i in (1, 2, 3)})
    g = A.on_ray(amp)
    hp = [[F() for _ in range(3)] for _ in range(3)]
    v = [F() for _ in range(3)]
    for lab, p in enumerate(wave, 1):
        p2 = sum((c * c for c in p), zero)
        for i in range(3):
            v[i] += F({frozenset([lab]): A.IU * A.KAP * (T - O) * sum((P[i][j] * p[j] for j in range(3)), zero)})
            for j in range(3):
                Qij = p[i] * p[j] - (p2 * A.rat(1, 3) if i == j else zero)
                hp[i][j] += F({frozenset([lab]): -2 * A.KAP**2 * T * Qij * A.rat(1, 3)})
    hp = [[A.on_ray(z) for z in row] for row in hp]
    X = [amp * ((T - O) * ni) for ni in n]
    hpn = matvec(hp, [F.const(ni) for ni in n])
    hpnn = dot([F.const(ni) for ni in n], hpn)
    gradg = [g.d_(i) for i in range(3)]
    R = [integral(hpnn.d_(i) * A.rat(1, 2)) for i in range(3)]
    sray = [-z for z in matvec(P, hpn)]
    K = [integral(sray[i] - R[i]) for i in range(3)]
    QX = [-g * g * n[i] - 2 * g * v[i] + 2 * dot(X, gradg) * n[i] for i in range(3)]
    HessX = [sum((g.d_(i).d_(j) * X[j] for j in range(3)), F()) for i in range(3)]
    PHessX = matvec(P, HessX)
    Pgrad = matvec(P, gradg)
    ngrad = dot([F.const(ni) for ni in n], gradg)
    QV = [2 * PHessX[i] + 2 * ngrad * v[i] - 2 * g * Pgrad[i] for i in range(3)]
    local = sum((X[i] * X[j] * hpnn.d_(i).d_(j) * A.rat(3, 2) for i in range(3) for j in range(3)), F())
    local -= 6 * sum((v[i] * X[j] * hpn[i].d_(j) for i in range(3) for j in range(3)), F())
    local += 3 * dot(v, matvec(hp, v)) - 3 * dot(v, v) * hpnn
    reduced = local - 3 * dot(R, QX) - 3 * dot(K, QV)
    # reduced is a diagonal third derivative: divide distinct-label extraction
    # by 6, unlike direct residual coefficients below.
    U3 = A.integ(reduced.coeff([1, 2, 3]), A._IT, E, O)
    return U3


def direct_measurement(wavevalues=None, C=None):
    A = alg
    U3 = first_only_U3(wavevalues)
    if C is None:
        C = Q(s.Matrix([1, 0, 0]))
    weight = sum((A.as_ring(C[i, j]) * A.NV[i] * A.NV[j] for i in range(3) for j in range(3)), A.ZERO) * A.rat(15, 2)
    result = A.sky_mean(-U3 * weight * A.rat(1, 6)).as_expr().expand()
    out = result.coeff(s.Symbol("kap"), 4)
    if wavevalues is None:
        out = qaverage(out, s.Symbol("b1"), s.Symbol("b2"))
    return s.expand(out.subs({s.Symbol("ts"): hand.e, s.Symbol("to"): hand.o}))


@lru_cache(None)
def field_coefficients(wavevalues=None):
    kap = s.Symbol("kappa")
    mu, v = hand.mu, hand.v
    t = eng.t
    k = s.Matrix([1, 0, 0])
    q = s.Matrix([mu, v, 0])
    wave = [k, q, -q] if wavevalues is None else [s.Matrix(p) for p in wavevalues]
    eng.Ctx.maxdeg = 3
    eng.Ctx.momenta = {i + 1: list(kap * p) for i, p in enumerate(wave)}

    @lru_cache(None)
    def trunc(expr):
        return s.Add(*(cf * kap**power for (power,), cf in s.Poly(s.expand(expr), kap).terms() if power <= 4))

    eng.Ctx.trunc = trunc
    nu = eng.Field({frozenset([i + 1]): -(kap**2) * p.dot(p) * t**2 / 9 for i, p in enumerate(wave)})
    z = eng.Field({frozenset([i]): s.S.One for i in (1, 2, 3)})
    H = [
        [eng.Field({frozenset([lab + 1]): -(kap**2) * t**2 * Q(p)[i, j] / 3 for lab, p in enumerate(wave)}) for j in range(3)]
        for i in range(3)
    ]
    _, source, _ = direct_field_residuals(nu, z, H)
    full = s.Matrix(3, 3, lambda i, j: s.expand(source[i][j].coeff([1, 2, 3])))
    source2 = full.applyfunc(lambda z: s.expand(z).coeff(kap, 2))
    assert (source2 - leading_cubic_source(*wave)).applyfunc(s.expand) == s.zeros(3)
    coeff = full.applyfunc(lambda z: s.expand(z.coeff(kap, 4) / t**2))
    assert not any(z.has(t) for z in coeff)
    eng.Ctx.trunc = None
    return coeff


def direct_field():
    C = Q(s.Matrix([1, 0, 0]))
    return qaverage(inner(C, field_coefficients()), hand.mu, hand.v)


def main():
    e, o, ti = hand.e, hand.o, hand.ti
    k = s.Matrix([1, 0, 0])
    C = Q(k)
    q = s.Matrix([hand.mu, hand.v, 0])
    A = o**4 / 12 - o**2 * e**2 / 2 + 2 * o * e**3 / 3 - e**4 / 4
    diff = o**4 - e**4
    delta = o**2 - e**2
    first = qaverage(first_cubic_correction(k, q, -q, C, e, o, ti), hand.mu, hand.v)
    # Background adjoint paired with C_H,3^[2]=8(Q_k+Q_q), followed by the
    # complete q angular average. Its next Green/phase coefficient is explicit.
    source2 = 8 * Q(k)
    phase = -A * inner(angular_A(C, k, k), source2) / 12
    Bsource = 2 * stf(ricci_linear(k, source2))
    green = -diff * inner(C, Bsource) / 240
    c3 = s.Rational(-14, 9)
    initial = diff * inner(C, Q(k)) * c3 / 10
    direct4 = diff * direct_field() / 40
    measurement = direct_measurement()
    total = s.expand(first + phase + green + initial + direct4 + measurement)
    log("all stationary cubic blocks assembled")
    # Independent full Einstein-radiation local-curvature map, not an adjoint
    # or field-engine formula. Apply the objective BEFORE averaging q.
    data = hand.geometry()
    aA = data[5]
    HessaTF = data[7]
    A2TF = data[8]
    compact = (-11 * A / 14 - 3 * diff / 40) * hand.nn(HessaTF)
    compact += (27 * A / 14 - 3 * diff / 40 - delta**2 / 28) * hand.nn(A2TF)
    compact += (27 * A / 7 + 7 * diff / 20 + (o**4 + 3 * ti**4) / 4 - e**4) * hand.nn(aA)
    weight = s.Rational(15, 2) * (hand.nx**2 - s.Rational(1, 3))
    reference = hand.qaverage(hand.sky(weight * compact.coeff(7)))
    difference = s.factor(total - reference)
    for name, value in [
        ("first adjoint", first),
        ("direct field", direct4),
        ("background phase", phase),
        ("background Green", green),
        ("background initial", initial),
        ("direct measurement", measurement),
        ("total", total),
        ("independent", reference),
    ]:
        print(name + ": " + str(s.factor(value)), flush=True)
    print("difference:", difference, flush=True)
    assert difference == 0
    # Detect dropped homogeneous lapse/initial-clock and direct contributions.
    assert s.diff(total, ti) != 0
    assert all(s.expand(value) != 0 for value in (first, direct4, measurement, initial, phase, green))
    J12 = delta / 9
    J32 = 4 * delta / 9
    J14 = -A * inner(angular_A(C, k, k), Q(k)) / 6 - diff * inner(C, Q(k)) / 120
    # w=5 P2, hence the 4pi/25 conversion. Keep BOTH 1[2]*3[4]
    # and 1[4]*3[2] orderings in the fourth-order spectrum coefficient.
    spectrum6 = s.expand(4 * s.pi * (J12 * total + J14 * J32) / 25)
    expected6 = -2 * s.pi * delta * (2835 * ti**4 + 9973 * o**4 - 25440 * o**2 * e**2 + 33920 * o * e**3 - 21288 * e**4) / 1913625
    assert s.expand(spectrum6 - expected6) == 0
    assert s.expand(J14 * J32) != 0
    print("[ok] known weighted C13/3 coefficient through kappa^6, including both linear/cubic gradient orderings")
    print("FULL SHELL-AVERAGED CUBIC QUADRUPOLE AGREES THROUGH FOUR GRADIENTS; not a finite-wavelength identity")


def arbitrary_triple_checks(symbolic=False):
    # These triples are noncoplanar and need not have equal magnitudes. All
    # five independent quadrupole covectors are checked BEFORE Wick/sky-q
    # contraction. --symbolic uses a rotation-complete polynomial chart for
    # arbitrary triples, not a sample. Wavelengths are still gradient-expanded.
    cases = [((1, 0, 0), (0, 1, 0), (0, 0, 1)), ((1, 2, 0), (0, 1, 1), (1, -1, 2))]
    if symbolic:
        p, b1, b2, ip, u, v = s.symbols("p b1 b2 ip u v")
        cases = [((p, 0, 0), (b1, b2, 0), (ip, u, v))]
    basis = [s.diag(1, 0, -1), s.diag(0, 1, -1)]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        B = s.zeros(3)
        B[i, j] = B[j, i] = 1
        basis.append(B)
    saved = hand.SHELL_REDUCTION, hand.momenta, hand.K
    try:
        for values in cases:
            wave = [s.Matrix(p) for p in values]
            P = sum(wave, s.zeros(3, 1))
            e, o, ti = hand.e, hand.o, hand.ti
            diff = o**4 - e**4
            delta = o**2 - e**2
            A = o**4 / 12 - o**2 * e**2 / 2 + 2 * o * e**3 / 3 - e**4 / 4
            M2 = leading_cubic_source(*wave)
            M4 = field_coefficients(values)
            c3 = 2 * (-2 * sum(p.dot(p) for p in wave) + sum(wave[i].dot(wave[j]) for i, j in ((0, 1), (0, 2), (1, 2)))) / 9
            hand.SHELL_REDUCTION = False
            hand.momenta = values
            hand.K = [tuple(sum(values[a][j] for a in range(3) if mask & (1 << a)) for j in range(3)) for mask in range(8)]
            data = hand.geometry()
            compact = (-11 * A / 14 - 3 * diff / 40) * hand.nn(data[7])
            compact += (27 * A / 14 - 3 * diff / 40 - delta**2 / 28) * hand.nn(data[8])
            compact += (27 * A / 7 + 7 * diff / 20 + (o**4 + 3 * ti**4) / 4 - e**4) * hand.nn(data[5])
            for C in basis:
                got = first_cubic_correction(*wave, C, e, o, ti) + diff * inner(C, M4) / 40
                got -= A * inner(angular_A(C, P, P), M2) / 12
                got -= diff * inner(C, 2 * stf(ricci_linear(P, M2))) / 240
                got += diff * inner(C, Q(P)) * c3 / 10 + direct_measurement(values, C)
                weight = s.Rational(15, 2) * sum(C[i, j] * hand.n[i] * hand.n[j] for i in range(3) for j in range(3))
                expected = hand.sky(weight * compact.coeff(7))
                assert s.expand(got - expected) == 0, (values, C, s.factor(got - expected))
            log("full uncontracted cubic quadrupole agrees at FOUR gradients for " + str(values) + "; all five components")
    finally:
        hand.SHELL_REDUCTION, hand.momenta, hand.K = saved


if __name__ == "__main__":
    main()
    arbitrary_triple_checks("--symbolic" in sys.argv)
