"""Exact first-field-only cubic measurement, with explicit mixed-time phases.

This is J_111 in stationary_cubic, NOT the full Einstein cubic map. Fourier
spatial jets are contracted as momenta before assembly. The sole four-time ray
term has its common time evaluated exactly. Remaining time integrals are not
evaluated; no gradient expansion or early dipole projection is used.
"""

from dataclasses import dataclass
from itertools import combinations, permutations

import sympy as s
from phase_integral_contraction import contract_phase_integrals, phase_sectors
from quadrupole_angular_source import ray_polynomials
from stationary_cubic import comm, linear_direction_driver


@dataclass(frozen=True)
class FourierFirst:
    momentum: s.MatrixBase
    zeta: s.Expr | int
    lapse: s.Expr | int
    spatial: s.MatrixBase
    time: s.Symbol

    def at(self, time):
        return s.sympify(self.lapse).subs(self.time, time), s.Matrix(self.spatial).subs(self.time, time)

    def prime(self, time):
        return s.Matrix(self.spatial).diff(self.time).subs(self.time, time)


@dataclass(frozen=True)
class RayFirst:
    p: s.MatrixBase
    phase: s.Expr | int
    G: s.MatrixBase
    hp: s.MatrixBase
    driver: s.MatrixBase
    X: s.MatrixBase
    v: s.MatrixBase
    R: s.MatrixBase
    K: s.MatrixBase


@dataclass(frozen=True)
class DirectExponentTriple:
    first: tuple
    second: dict
    third: s.Expr
    radial_transport: s.Expr
    n: s.MatrixBase


def integral(expr, var, lo, hi):
    if isinstance(expr, s.MatrixBase):
        return expr.applyfunc(lambda z: integral(z, var, lo, hi))
    return s.S.Zero if expr == 0 or lo == hi else s.Integral(expr, (var, lo, hi))


def ordered_radial_transport(terms, u, v, w, e, o):
    """Mixed I3 from -3 integral (n.R) |v|^2 dt, with t integrated out.

    terms(a,b,c,u,v,w) supplies (n.r_a(u)) V_b(v).V_c(w), including phases.
    On e<=u<=v<=w<=o the eliminated time has length v-u. Six permutations
    include both V-time orderings and the three choices of the R direction.
    """
    return -s.Add(
        *(
            integral(integral(integral((v - u) * terms(a, b, c, u, v, w), u, e, v), v, e, w), w, e, o)
            for a, b, c in permutations(range(3))
        )
    )


def bilinear_ray_forcing(A, B, n):
    """Polarized QX without n(v_a.v_b), and complete polarized QV.

    Exact Fourier specialization of second_ray_forcing: partial_i is i*p_i.
    These are mixed second derivatives, not twice those derivatives.
    """
    Ga, Gb, ha, hb = A.G, B.G, A.hp, B.hp
    pa, pb, va, vb = A.p, B.p, A.v, B.v
    nn = lambda M: n.dot(M * n)
    QX = -(Ga * vb + Gb * va) / 2 - (Ga * Gb + Gb * Ga) * n / 8
    QX += s.I * (pb.dot(A.X) * Gb + pa.dot(B.X) * Ga) * n / 2
    raw = s.I * (pb.dot(A.X) * B.driver + pa.dot(B.X) * A.driver)
    raw -= (ha * vb + hb * va) / 2
    raw += (va * nn(hb) + vb * nn(ha)) / 2
    raw -= (comm(Ga, hb) + comm(Gb, ha)) * n / 8
    raw += s.I * (pa.dot(vb) * Ga + pb.dot(va) * Gb) * n / 2
    raw += s.I * (pa.dot(n) * Ga * vb + pb.dot(n) * Gb * va) / 2
    raw -= s.I * (pb * va.dot(Gb * n) + pa * vb.dot(Ga * n))
    raw += s.I * (pb.dot(Ga * n) * Gb + pa.dot(Gb * n) * Ga) * n / 4
    raw += s.I * (pb.dot(n) * comm(Ga, Gb) + pa.dot(n) * comm(Gb, Ga)) * n / 8
    raw -= s.I * (Ga * pb * nn(Gb) + Gb * pa * nn(Ga)) / 4
    return QX, (s.eye(3) - n * n.T) * raw


def polarized_bulk(states, n):
    """Complete cubic bulk except the explicitly restored R.n |v|^2 term."""
    nn = lambda M: n.dot(M * n)
    # All three local phases are common. Remove scalar identity pieces before
    # commutation, and sum the phase-free matrices before multiplying the phase.
    h = [A.G / A.phase - s.trace(A.G / A.phase) * s.eye(3) / 3 for A in states]
    hp = [A.hp / A.phase for A in states]
    local_matrix = sum((comm(h[a], comm(h[b], hp[c])) for a, b, c in permutations(range(3))), s.zeros(3))
    out = nn(local_matrix) * s.prod(A.phase for A in states) / 48
    for a, b, c in permutations(range(3)):
        A, B, C = states[a], states[b], states[c]
        out -= C.p.dot(A.X) * C.p.dot(B.X) * nn(C.hp) / 4
        out -= s.I * C.p.dot(B.X) * A.v.dot(C.hp * n)
        out += A.v.dot(B.hp * C.v) / 2 - A.v.dot(B.v) * nn(C.hp) / 2
    for a in range(3):
        b, c = tuple(i for i in range(3) if i != a)
        QX, QV = bilinear_ray_forcing(states[b], states[c], n)
        out -= states[a].R.dot(QX) + states[a].K.dot(QV)
    return out


def compose_mixed_exponents(U1, U2, U3, F1e, F2e, total_F1prime_e, d1, d2, D1, D2, o):
    """Polarized physical endpoints/clock, with the mixed terms counted once."""
    third = U3
    for a in range(3):
        b, c = tuple(i for i in range(3) if i != a)
        third -= F1e[a] * d2[b, c] + F2e[b, c] * d1[a] + total_F1prime_e[a] * d1[b] * d1[c]
        third += D1[a] * U2[b, c].diff(o) + D2[b, c] * U1[a].diff(o) + D1[b] * D1[c] * U1[a].diff(o, 2)
    for a, b, c in permutations(range(3)):
        third -= D1[a] * (F1e[b].diff(o) * d1[c] + F1e[b] * d1[c].diff(o))
    second = {}
    for b, c in combinations(range(3), 2):
        second[b, c] = U2[b, c] - F1e[b] * d1[c] - F1e[c] * d1[b] + D1[b] * U1[c].diff(o) + D1[c] * U1[b].diff(o)
    return second, third


def first_only_exponent_triple(seeds, n, e, o, initial):
    assert len(seeds) == 3
    n = s.Matrix(n)
    e, o, initial = map(s.sympify, (e, o, initial))
    t, u, v, w = (s.Dummy(x, real=True) for x in ("ray", "first", "middle", "last"))
    O = s.Dummy("observer", positive=True)
    p = [s.Matrix(seed.momentum) for seed in seeds]
    for seed in seeds:
        H = s.Matrix(seed.spatial)
        assert H.shape == (3, 3) and H == H.T and s.trace(H) == 0
        assert s.diff(s.sympify(seed.zeta), seed.time) == 0
    phase = lambda a, T: s.exp(s.I * p[a].dot(n) * (O - T))
    states = []
    atoms = []
    a_ray = []
    a_observer = []
    for a, seed in enumerate(seeds):
        h = seed.at(t)[1]
        hp = seed.prime(t)
        G = (2 * s.sympify(seed.zeta) * s.eye(3) + h) * phase(a, t)
        Hprime = hp * phase(a, t)
        atom = ray_polynomials(p[a], n, seed.zeta, seed.at(u)[1], seed.prime(u), t, u)
        atoms.append(atom)
        X = integral(atom.X_upper * phase(a, u), u, t, O)
        velocity = -integral(atom.v_upper * phase(a, u), u, t, O)
        R = integral(atom.R_lower * phase(a, u), u, e, t)
        K = integral(atom.K_lower * phase(a, u), u, e, t)
        driver = linear_direction_driver(n, Hprime, [s.I * p[a][i] * G for i in range(3)])
        states.append(RayFirst(p[a], phase(a, t), G, Hprime, driver, X, velocity, R, K))
        a_ray.append(s.sympify(seed.lapse).subs(seed.time, t) * phase(a, t))
        a_observer.append(s.sympify(seed.lapse).subs(seed.time, u))
    F1 = [n.dot(A.hp * n) / 2 for A in states]
    U1 = [integral(F, t, e, O) for F in F1]
    F2 = {}
    U2 = {}
    frame2 = {}
    d2 = {}
    D2 = {}
    d1 = [-e * a.subs(t, e) for a in a_ray]
    D1 = [-(initial**2 * a.subs(u, initial) + integral(u * a, u, initial, O)) / O for a in a_observer]
    for b, c in combinations(range(3), 2):
        B, C = states[b], states[c]
        F2[b, c] = (
            s.I * (C.p.dot(B.X) * n.dot(C.hp * n) + B.p.dot(C.X) * n.dot(B.hp * n)) / 2 - B.v.dot(C.hp * n) - C.v.dot(B.hp * n)
        )
        U2[b, c] = integral(F2[b, c], t, e, O)
        hb, hc = seeds[b].at(u)[1], seeds[c].at(u)[1]
        frame2[b, c] = integral((comm(hb, hc.diff(u)) + comm(hc, hb.diff(u))) / 8, u, initial, O)
        ae_b, ae_c = a_ray[b].subs(t, e), a_ray[c].subs(t, e)
        grad_b = s.I * p[b] * ae_b
        grad_c = s.I * p[c] * ae_c
        d2[b, c] = e * (
            ae_b * ae_c
            - B.X.subs(t, e).dot(grad_c)
            - C.X.subs(t, e).dot(grad_b)
            - d1[b] * a_ray[c].diff(t).subs(t, e)
            - d1[c] * a_ray[b].diff(t).subs(t, e)
        )
        ab, ac = a_observer[b], a_observer[c]
        D2[b, c] = (
            initial**3 * (ab * ac.diff(u) + ac * ab.diff(u)).subs(u, initial)
            - integral(u * ab * ac, u, initial, O)
            - O * (ab.subs(u, O) * D1[c] + ac.subs(u, O) * D1[b])
            - D1[b] * D1[c]
        ) / O
    third = integral(polarized_bulk(states, n), t, e, O)
    for a in range(3):
        b, c = tuple(i for i in range(3) if i != a)
        third -= states[a].K.subs(t, O).dot(frame2[b, c] * n)
    second, third = compose_mixed_exponents(
        U1,
        U2,
        third,
        [F.subs(t, e) for F in F1],
        {k: F.subs(t, e) for k, F in F2.items()},
        [F.diff(t).subs(t, e) for F in F1],
        d1,
        d2,
        D1,
        D2,
        O,
    )

    def radial(a, b, c, U, V, W):
        return (
            n.dot(atoms[a].R_lower.subs(u, U))
            * atoms[b].v_upper.subs(u, V).dot(atoms[c].v_upper.subs(u, W))
            * phase(a, U)
            * phase(b, V)
            * phase(c, W)
        )

    restored = ordered_radial_transport(radial, u, v, w, e, O)
    return DirectExponentTriple(
        tuple(x.subs(O, o) for x in U1),
        {k: x.subs(O, o) for k, x in second.items()},
        (third + restored).subs(O, o),
        restored.subs(O, o),
        n,
    )


def normalized_cubic_raw(data):
    """Complete third map before FINAL monopole centering only."""
    n = data.n
    average = lambda expr: contract_phase_integrals(expr, n).expression
    centered = tuple(expr - average(expr) for expr in data.first)
    raw = -data.third - s.prod(centered)
    for a in range(3):
        b, c = tuple(i for i in range(3) if i != a)
        pair = data.second[b, c]
        raw += centered[a] * (pair - average(pair)) + average(centered[b] * centered[c]) * centered[a]
    return raw


def measurement_phase_plan(data):
    """Unexpanded sectors for combining with adjoint terms BEFORE moments.

    Prefer this to expanding a dense arbitrary-triple Fr expression in isolation.
    All remaining sky operations are the explicit finite moment formula.
    """
    return phase_sectors(normalized_cubic_raw(data), data.n)


def contract_direct_measurement(data, C, maximum_rank=13):
    """Full normalized quadrupole J_111, optionally materialized in Fr.

    Only final monopole centering is omitted: the weight has zero average.
    No early dipole removal. Frame is already included; do NOT add it again.
    Dense noncoplanar Fr materialization can be expensive: measurement_phase_plan
    preserves the exact expression for merging with the other stationary blocks.
    """
    C = s.Matrix(C)
    assert C.shape == (3, 3) and C == C.T and s.trace(C) == 0
    weight = s.Rational(15, 2) * data.n.dot(C * data.n)
    return contract_phase_integrals(normalized_cubic_raw(data), data.n, weight, maximum_rank)
