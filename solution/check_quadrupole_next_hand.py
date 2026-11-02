"""Independent O(k0^6) check of C2^(11) and C2^(22)/4, using hand map formulas.

No observable generator or project module is imported. Fields are obtained from
an elementary curvature expansion through t^4; the ray, emission and clock terms
are explicit hand-expanded formulas. Only polynomial integrations and exact
isotropic moments are used. This does NOT compute the cubic C2^(13) correction.

Optional --compare reads an already saved generator OUTPUT, after the independent
calculation, to compare all six quadrupole-tensor components. It never runs or
imports the generator. The pickle must be the trusted local project output.

Run: .venv/bin/python -B solution/check_quadrupole_next_hand.py [--compare]
"""

import argparse
import pickle
from functools import lru_cache
from pathlib import Path

import sympy as s

b, c, t, o, e, ti = s.symbols("b c t o e ti", real=True)
nx, ny, nz = s.symbols("nx ny nz", real=True)
n = s.Matrix([nx, ny, nz])
identity = s.eye(3)
p, q = s.Matrix([b, 0, c]), s.Matrix([-b, 0, c])
K, k2, dot = p + q, b**2 + c**2, p.dot(q)
P, Q = p * p.T - identity * k2 / 3, q * q.T - identity * k2 / 3


def tf(matrix):
    return matrix - identity * s.trace(matrix) / 3


def nn(matrix):
    return (n.T * matrix * n)[0]


def connection(metric, momentum):
    return [
        [
            [s.I * (momentum[i] * metric[a, j] + momentum[j] * metric[a, i] - momentum[a] * metric[i, j]) / 2 for j in range(3)]
            for i in range(3)
        ]
        for a in range(3)
    ]


def linear_ricci(christoffel, momentum):
    return s.Matrix(
        3,
        3,
        lambda i, j: sum(s.I * momentum[a] * christoffel[a][i][j] - s.I * momentum[j] * christoffel[a][i][a] for a in range(3)),
    )


def gradient_fields():
    """Return A2=t^2 U+t^4 V from the local equations, without SVT inverses."""
    h, nu, hd, nd = s.symbols("h nu hd nd")
    gp, gq = 2 * identity + h * P, 2 * identity + h * Q
    pair = 4 * identity + 2 * h * (P + Q) + h**2 * (P * Q + Q * P) / 2
    cp, cq = connection(gp, p), connection(gq, q)
    cl = connection(pair, K)
    cpair = [
        [
            [s.expand(cl[a][i][j] - sum(gp[a, l] * cq[l][i][j] + gq[a, l] * cp[l][i][j] for l in range(3))) for j in range(3)]
            for i in range(3)
        ]
        for a in range(3)
    ]
    ricci = (
        linear_ricci(cpair, K)
        + s.Matrix(
            3,
            3,
            lambda i, j: sum(
                cp[a][a][l] * cq[l][i][j] + cq[a][a][l] * cp[l][i][j] - cp[a][j][l] * cq[l][i][a] - cq[a][j][l] * cp[l][i][a]
                for a in range(3)
                for l in range(3)
            ),
        )
        - gp * linear_ricci(cq, q)
        - gq * linear_ricci(cp, p)
    )
    hessian = s.Matrix(3, 3, lambda i, j: -s.I * nu * sum(cp[a][i][j] * q[a] + cq[a][i][j] * p[a] for a in range(3))) + nu * (
        gp * q * q.T + gq * p * p.T
    )
    source = tf(nd * hd * (P + Q) + ricci - 2 * hessian - 2 * nu**2 * (p * q.T + q * p.T))
    early = {h: -(t**2) / 3, nu: -k2 * t**2 / 9, hd: -2 * t / 3, nd: -2 * k2 * t / 9}
    source = source.applyfunc(lambda z: s.expand(z.subs(early)))
    U = -source.subs(t, 0) / 3
    assert (U - tf(2 * (p * p.T + q * q.T) + p * q.T + q * p.T) / 3).applyfunc(s.expand) == s.zeros(3)
    a2 = (4 * k2 - dot) / 9
    rlinear = tf((K.dot(K) * U - K * (K.T * U) - (U * K) * K.T) / 2)
    V = -(source.applyfunc(lambda z: z.coeff(t, 2)) + rlinear + 2 * tf(K * K.T) * a2) / 10
    V = V.applyfunc(s.factor)
    # The lapse t^4 coefficient is not needed in the quadratic map, but is checked.
    lapse_source = -s.trace(hessian) + hd**2 * s.trace(P * Q) / 2 + 4 * nu**2 * dot + 12 * nu**2 / t**2
    lapse_c2 = s.expand(lapse_source.subs(early)).coeff(t, 2)
    lapse_b2 = s.factor(-(K.dot(K) * a2 + lapse_c2) / 30)
    assert s.expand(lapse_b2 + k2 * (19 * k2 + 6 * dot) / 405) == 0
    # Independent checks of the constraints at the orders supplied by these data.
    profiles = {
        h: -(t**2) / 3 + k2 * t**4 / 60,
        nu: -k2 * t**2 / 9 + k2**2 * t**4 / 270,
        hd: -2 * t / 3 + k2 * t**3 / 15,
        nd: -2 * k2 * t / 9 + 2 * k2**2 * t**3 / 135,
    }
    metric2, lapse2 = t**2 * U + t**4 * V, t**2 * a2 + t**4 * lapse_b2
    ham_source = (
        s.trace(ricci) - 4 * s.trace(hessian) + 4 * nu**2 * dot - hd**2 * s.trace(P * Q) / 2 + 12 * nd**2 - 24 * nu**2 / t**2
    )
    ham = (
        -(K.T * metric2 * K)[0]
        + 4 * K.dot(K) * lapse2
        + 12 * s.diff(lapse2, t) / t
        + 12 * lapse2 / t**2
        + ham_source.subs(profiles)
    )
    for degree in (0, 2):
        assert s.expand(ham).coeff(t, degree) == 0
    for i in range(3):
        div = (
            hd
            / 2
            * sum(
                cp[j][j][a] * Q[a, i] + cq[j][j][a] * P[a, i] - cp[a][j][i] * Q[j, a] - cq[a][j][i] * P[j, a]
                for j in range(3)
                for a in range(3)
            )
        )
        mom_source = div + s.I * hd * nu * (P * q + Q * p)[i] + 2 * s.I * nd * nu * K[i]
        momentum = (
            s.I * (metric2.diff(t) * K)[i] / 2 - 2 * s.I * K[i] * (s.diff(lapse2, t) - lapse2 / t) + mom_source.subs(profiles)
        )
        for degree in (1, 3):
            assert s.expand(momentum).coeff(t, degree) == 0
    print("[ok] local t^2/t^4 gradient data and Hamiltonian/momentum constraints; no inverse total momenta")
    return U, V


def time_integral(expr):
    return s.expand(
        sum(
            coef * (o ** (power[0] + 1) - e ** (power[0] + 1)) / (power[0] + 1)
            for power, coef in s.Poly(s.expand(expr), t).terms()
        )
    )


@lru_cache(None)
def sphere_moment(a, bb, cc):
    if a % 2 or bb % 2 or cc % 2:
        return s.S.Zero
    return s.factorial2(a - 1) * s.factorial2(bb - 1) * s.factorial2(cc - 1) / s.factorial2(a + bb + cc + 1)


def sky(expr):
    return s.expand(sum(coef * sphere_moment(*powers) for powers, coef in s.Poly(s.expand(expr), nx, ny, nz).terms()))


def shell_integral(expr):
    mu, k = s.symbols("mu k", real=True)
    out = 0
    for (a, bb), coef in s.Poly(s.expand(expr), b, c).terms():
        assert a % 2 == 0 and bb % 2 == 0
        out += coef * (k**2 * (1 - mu) / 2) ** (a // 2) * (k**2 * (1 + mu) / 2) ** (bb // 2)
    return s.factor(s.integrate(s.expand(out), (mu, -1, 1)))


def hand_map(U, V):
    """Even-in-n part at degrees 2 and 4; odd degrees cannot enter ell=2."""
    lam, pn, qn = o - t, p.dot(n), q.dot(n)
    sigp, sigq = -t * P / 3, -t * Q / 3
    sigp4, sigq4 = k2 * t**3 * P / 30, k2 * t**3 * Q / 30
    sig2 = t * U - 2 * t * (P + Q) / 3
    sig24 = 2 * t**3 * V + t**3 * (P * Q + Q * P) / 9 + k2 * t**3 * (P + Q) / 15
    field = nn(sig24) - K.dot(n) ** 2 * lam**2 * nn(sig2) / 2
    position = lam**2 * (nn(sigq) * (qn**2 + pn * qn - dot / 2) + nn(sigp) * (pn**2 + pn * qn - dot / 2))

    def direction(pn, momentum, mat, other, other4, otherpn):
        vp = 2 * pn * n - momentum
        ep2 = (o**2 / 6 - t**2 / 3) * mat * n - (o**2 - t**2) * nn(mat) * n / 6 - lam**2 * pn * vp / 2
        # Check this direction against the explicitly integrated first-order ray.
        xp2 = (
            lam * o**2 * mat * n / 6
            - lam**2 * (t + 2 * o) * mat * n / 9
            + lam**2 * (t + 2 * o) * nn(mat) * n / 18
            + lam**3 * pn * vp / 6
        )
        assert (xp2.diff(t) - ep2).applyfunc(s.expand) == s.zeros(3, 1)
        return (
            -2 * nn(other4)
            - 2 * (n.T * other * ep2)[0]
            + 2 * lam**2 * otherpn * (n.T * other * vp)[0]
            + lam**2 * otherpn**2 * nn(other)
        )

    dirs = direction(pn, p, P, sigq, sigq4, qn) + direction(qn, q, Q, sigp, sigp4, pn)
    tau2 = (o**2 - e**2) / 6
    leading = -s.expand((o**2 - e**2) * nn(U) / 2)
    parts = {
        "field": -time_integral(field),
        "position": -time_integral(position),
        "direction": -time_integral(dirs),
        "endpoint": -k2 * e**4 * nn(P + Q) / 27,
        "clock": k2 * (o**4 + 3 * ti**4) * nn(P + Q) / 108,
        "normalization": tau2**2 * nn(P) * nn(Q),
    }
    # Each term may be centred separately. Constants do not enter the quadrupole.
    parts = {name: s.expand(value - sky(value)) for name, value in parts.items()}
    print("[ok] hand-expanded field, ray, endpoint, clock and random-monopole terms assembled")
    return leading, parts


def check_linear():
    tau2, lam, pn = (o**2 - e**2) / 6, o - t, p.dot(n)
    linear2 = tau2 * nn(P)
    linear4 = -time_integral(k2 * t**3 * nn(P) / 30 + t * lam**2 * pn**2 * nn(P) / 6)
    result = s.factor(8 * s.pi * sky(linear2 * linear4) / 5)
    # Separate derivation from the linear transfer integral: w2=2/15-11 y^2/315.
    k = s.Symbol("k", real=True)
    transfer2 = -(k**2) * (o**2 - e**2) / 45
    transfer4 = time_integral(k**4 * (t**3 / 225 + 11 * t * (o - t) ** 2 / 945))
    expected = s.expand(8 * s.pi * transfer2 * transfer4)
    assert s.expand(result - expected.subs(k**6, k2**3)) == 0
    print("[ok] C11 k0^6 agrees independently with the linear transfer integral")
    print("C11 k0^6 / sigma^2:", s.factor(expected))


def compare_saved_output(hand4):
    """Compare to a saved OUTPUT only; no import of any generator or its helpers."""
    path = Path(__file__).resolve().parent / "results" / "superhorizon_kernels_2.pkl"
    with path.open("rb") as stream:
        data = pickle.load(stream)
    gens = s.symbols(data["gens"])
    G = dict(zip(data["gens"], gens, strict=True))
    directions = [G["n1"], G["n2"], G["n3"]]
    mu, z = s.symbols("mu z", real=True)
    substitutions = {
        G["p"]: 1,
        G["b1"]: mu,
        G["b2"]: z,
        G["ip"]: 1,
        G["u"]: 1 / (2 * (1 + mu)),
        G["v"]: 1 / (2 * (1 - mu)),
        G["its"]: 1 / G["ts"],
        G["itI"]: 1 / G["tI"],
        G["itob"]: 1 / G["tob"],
    }

    def output_moment(extra):
        terms = []
        for mon, coefficient in data["kernels"][(1, 2)].items():
            if mon[gens.index(G["kap"])] != 4:
                continue
            moment = sphere_moment(*(mon[gens.index(nn)] + degree for nn, degree in zip(directions, extra, strict=True)))
            if moment:
                terms.append(
                    s.sympify(coefficient)
                    * moment
                    * s.Mul(*(g**power for g, power in zip(gens, mon, strict=True) if g not in directions + [G["kap"]]))
                )
        poly = s.Poly(s.expand(s.Add(*terms).subs(substitutions, simultaneous=True)), z)
        return s.factor(sum(coef * (1 - mu**2) ** (power // 2) * z ** (power % 2) for (power,), coef in poly.terms()))

    mean = output_moment((0, 0, 0))
    moments = s.Matrix(3, 3, lambda i, j: output_moment(tuple(int(a == i) + int(a == j) for a in range(3))))
    output_tensor = s.Rational(15, 2) * (moments - identity * mean / 3)
    own_tensor = s.Matrix(3, 3, lambda i, j: s.Rational(15, 2) * (sky(hand4 * n[i] * n[j]) - (sky(hand4) / 3 if i == j else 0)))
    rotation = s.Matrix([[b, 0, c], [-c, 0, b], [0, -1, 0]])
    rotated = rotation.T * output_tensor.subs({mu: c**2 - b**2, z: 2 * b * c}) * rotation
    own_tensor = own_tensor.subs({o: G["tob"], e: G["ts"], ti: G["tI"]})
    for component in rotated - own_tensor:
        numerator = s.together(component).as_numer_denom()[0]
        remainder = s.rem(s.Poly(s.expand(numerator), c), s.Poly(c**2 + b**2 - 1, c)).as_expr()
        assert s.expand(remainder) == 0, s.factor(remainder)
    print("[ok] independent comparison to saved OUTPUT: all quadrupole-tensor components agree")


def main(compare=False):
    U, V = gradient_fields()
    leading, parts = hand_map(U, V)
    k = s.Symbol("k", real=True)
    lead_spectrum = s.factor(s.pi * shell_integral(sky(leading**2)) / 5)
    assert s.expand(lead_spectrum - 68 * s.pi * k**4 * (o**2 - e**2) ** 2 / 6075) == 0
    print("C22/4 k0^4 / sigma^4:", lead_spectrum)
    contributions = {name: s.factor(2 * s.pi * shell_integral(sky(leading * part)) / 5) for name, part in parts.items()}
    for name, value in contributions.items():
        print(f"C22/4 k0^6 {name} / sigma^4: {value}")
    total = s.factor(sum(contributions.values()))
    target = (
        -8
        * s.pi
        * k**6
        * (o**2 - e**2)
        * (105 * ti**4 + 719 * o**4 - 1920 * o**2 * e**2 + 2560 * o * e**3 - 1464 * e**4)
        / 637875
    )
    assert s.expand(total - target) == 0
    assert contributions["normalization"] == 0
    for name in ("position", "direction", "endpoint", "clock"):
        assert contributions[name] != 0
    print("[ok] negative controls: each ray/endpoint/clock omission changes the spectrum")
    print("TOTAL C22/4 k0^6 / sigma^4:", total)
    check_linear()
    if compare:
        compare_saved_output(s.expand(sum(parts.values())))
    print("INDEPENDENT HAND-MAP QUADRUPOLE CHECK PASSED (C11 and C22 only; not C13)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--compare", action="store_true", help="compare to trusted saved kernel output, without running the generator"
    )
    main(parser.parse_args().compare)
