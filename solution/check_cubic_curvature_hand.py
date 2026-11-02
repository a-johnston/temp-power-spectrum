"""Independent local-curvature check of the k0^6 cubic quadrupole.

No project modules, field engine, ray generator, or saved output are imported.
All arithmetic is exact. The nilpotent labels extract symmetric Fourier kernels;
the local covariant gradient expansion is derived in report/critical_review.md.
Run: .venv/bin/python -B solution/check_cubic_curvature_hand.py
"""

from functools import lru_cache

import sympy as s

mu, v = s.symbols("mu v", real=True)
o, e, ti = s.symbols("o e ti", real=True)
nx, ny, nz = s.symbols("nx ny nz", real=True)
n = (nx, ny, nz)
I = s.I
momenta = ((1, 0, 0), (mu, v, 0), (-mu, -v, 0))
K = [tuple(sum(momenta[a][j] for a in range(3) if mask & (1 << a)) for j in range(3)) for mask in range(8)]


def simp(x):
    return s.Poly(s.expand(x), v).rem(s.Poly(v * v + mu * mu - 1, v)).as_expr().expand()


class J:
    """Commutative algebra with z_a^2=0, using no factorial rescaling."""

    def __init__(self, x=0):
        self.c = x if isinstance(x, dict) else ({0: s.sympify(x)} if x != 0 else {})

    def __add__(self, other):
        other = other if isinstance(other, J) else J(other)
        out = dict(self.c)
        for m, c in other.c.items():
            out[m] = simp(out.get(m, 0) + c)
        return J({m: c for m, c in out.items() if c != 0})

    __radd__ = __add__

    def __neg__(self):
        return J({m: -c for m, c in self.c.items()})

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return -self + other

    def __mul__(self, other):
        other = other if isinstance(other, J) else J(other)
        out = {}
        for a, x in self.c.items():
            for b, y in other.c.items():
                if not a & b:
                    out[a | b] = out.get(a | b, 0) + x * y
        return J({m: simp(c) for m, c in out.items() if c != 0})

    __rmul__ = __mul__

    def __truediv__(self, x):
        return self * (s.S.One / x)

    def d(self, j):
        return J({m: simp(I * K[m][j] * c) for m, c in self.c.items() if K[m][j] != 0})

    def coeff(self, mask):
        return self.c.get(mask, s.S.Zero)


def expg(a):
    return J({m: s.Integer(a) ** m.bit_count() for m in range(8)})


def matrix(f):
    return [[f(i, j) for j in range(3)] for i in range(3)]


def tf(A):
    tr = sum(A[j][j] for j in range(3))
    return matrix(lambda i, j: A[i][j] - (tr / 3 if i == j else 0))


def nn(A):
    return sum(n[i] * n[j] * A[i][j] for i in range(3) for j in range(3))


@lru_cache(None)
def moment(a, b, c):
    if a % 2 or b % 2 or c % 2:
        return s.S.Zero
    return s.factorial2(a - 1) * s.factorial2(b - 1) * s.factorial2(c - 1) / s.factorial2(a + b + c + 1)


def sky(x):
    return simp(sum(c * moment(*p) for p, c in s.Poly(s.expand(x), *n).terms()))


def qaverage(x):
    x = simp(x)
    assert not x.has(v), x
    return s.integrate(x, (mu, -1, 1)) / 2


def main():
    g = J({1: s.S.One, 2: s.S.One, 4: s.S.One})
    dg = [g.d(i) for i in range(3)]
    lapg = sum(g.d(i).d(i) for i in range(3))
    dg2 = sum(x * x for x in dg)
    # gamma0=e^(2g) delta. C is the covariant tensor gamma0*A.
    C = tf(matrix(lambda i, j: (g.d(i).d(j) - dg[i] * dg[j]) / 3))
    a0 = (2 * lapg + dg2) / 18
    Gamma = [
        [[(dg[j] if a == i else 0) + (dg[i] if a == j else 0) - (dg[a] if i == j else 0) for j in range(3)] for i in range(3)]
        for a in range(3)
    ]
    DC = [
        [
            [C[i][j].d(l) - sum(Gamma[a][l][i] * C[a][j] + Gamma[a][l][j] * C[i][a] for a in range(3)) for j in range(3)]
            for i in range(3)
        ]
        for l in range(3)
    ]
    DDC = [
        [
            [
                [
                    DC[l][i][j].d(k)
                    - sum(
                        Gamma[a][k][l] * DC[a][i][j] + Gamma[a][k][i] * DC[l][a][j] + Gamma[a][k][j] * DC[l][i][a]
                        for a in range(3)
                    )
                    for j in range(3)
                ]
                for i in range(3)
            ]
            for l in range(3)
        ]
        for k in range(3)
    ]
    em2, em4 = expg(-2), expg(-4)
    A = matrix(lambda i, j: em2 * C[i][j])
    LapA = matrix(lambda i, j: em4 * sum(DDC[k][k][i][j] for k in range(3)))
    A2 = matrix(lambda i, j: em4 * sum(C[i][a] * C[a][j] for a in range(3)))
    aA = matrix(lambda i, j: em4 * a0 * C[i][j])
    B = tf(matrix(lambda i, j: LapA[i][j] / 20 + aA[i][j] / 5 + 3 * A2[i][j] / 5))
    # Independent verification of B from the linear variation of the connection,
    # not the Ricci-flow identity used to simplify it in the derivation.
    dGamma = [[[em2 * (DC[i][a][j] + DC[j][a][i] - DC[a][i][j]) / 2 for j in range(3)] for i in range(3)] for a in range(3)]

    def cov_dGamma(k, a, i, j):
        return dGamma[a][i][j].d(k) + sum(
            Gamma[a][k][b] * dGamma[b][i][j] - Gamma[b][k][i] * dGamma[a][b][j] - Gamma[b][k][j] * dGamma[a][i][b]
            for b in range(3)
        )

    a = em2 * a0
    Ric = matrix(lambda i, j: -3 * A[i][j] - (12 * a if i == j else 0))
    dRic = matrix(
        lambda i, j: (
            em2 * sum(cov_dGamma(k, k, i, j) - cov_dGamma(j, k, i, k) for k in range(3))
            - sum(A[i][k] * Ric[k][j] for k in range(3))
        )
    )
    Hessa = matrix(lambda i, j: em2 * (a.d(i).d(j) - sum(Gamma[k][i][j] * a.d(k) for k in range(3))))
    Bdirect = tf(matrix(lambda i, j: -(dRic[i][j] - 2 * Hessa[i][j] + 4 * aA[i][j]) / 10))
    HessaTF, A2TF = tf(Hessa), tf(A2)
    for i in range(3):
        for j in range(3):
            assert not (Bdirect[i][j] - B[i][j]).c
            # Vanishing Cotton tensor of the conformally flat primordial metric.
            assert not (LapA[i][j] - 3 * HessaTF[i][j] + 9 * A2TF[i][j] + 18 * aA[i][j]).c
    print(
        "[ok] local fourth-gradient field and Cotton identities from direct covariant derivatives, through cubic amplitude",
        flush=True,
    )
    DnnA = em4 * sum(
        n[k] * n[l] * n[i] * n[j] * DDC[k][l][i][j] for k in range(3) for l in range(3) for i in range(3) for j in range(3)
    )
    Ann2 = em4 * nn(C) * nn(C)
    Delta = o * o - e * e
    L = o**4 / 12 - o * o * e * e / 2 + 2 * o * e**3 / 3 - e**4 / 4  # integral_e^o t(o-t)^2 dt
    pieces = {
        "spatial Hessian": (-L / 2, DnnA),
        "fourth-gradient field": (-(o**4 - e**4) / 2, nn(B)),
        "direction A^2": (-(Delta**2) / 4, nn(A2)),
        "direction+normalization Ann^2": (3 * Delta**2 / 8, Ann2),
        "emission": (-(e**4), nn(aA)),
        "clock": ((o**4 + 3 * ti**4) / 4, nn(aA)),
    }
    leading = nn(A) * (-Delta / 2)
    # The whole even-sky map reduces, in the quadrupole, to only three
    # fourth-gradient tensors. Verify all tensor components for every mask.
    compact4 = (
        (-11 * L / 14 - 3 * (o**4 - e**4) / 40) * nn(HessaTF)
        + (27 * L / 14 - 3 * (o**4 - e**4) / 40 - Delta**2 / 28) * nn(A2TF)
        + (27 * L / 7 + 7 * (o**4 - e**4) / 20 + (o**4 + 3 * ti**4) / 4 - e**4) * nn(aA)
    )
    for mask in range(1, 8):
        for i in range(3):
            for j in range(i, 3):
                test = n[i] * n[j] - (s.Rational(1, 3) if i == j else 0)
                # Apply the sky functional before introducing symbolic times.
                projected_full = sum(time * sky(test * value.coeff(mask)) for time, value in pieces.values())
                assert s.expand(projected_full - sky(test * compact4.coeff(mask))) == 0
    print("[ok] compact three-tensor quadrupole equals full map, every component through cubic amplitude", flush=True)
    # Quadrupole projection against the external k=(1,0,0).
    P2 = (3 * nx * nx - 1) / 2
    projected = {}
    for name, (time, value) in pieces.items():
        projected[name] = s.factor(time * qaverage(sky(P2 * value.coeff(7))))
        print("cubic <P2 K3>,", name, ":", projected[name], flush=True)
    cubic2 = s.factor(qaverage(sky(P2 * leading.coeff(7))))
    cubic4 = s.factor(sum(projected.values()))
    linear2 = leading.coeff(1)
    linear4 = sum(time * value.coeff(1) for time, value in pieces.values())
    # K1_ell=5 <P2 K1> P2. C13/3=4pi/5 <K1_ell <K3_ell>_q>_n.
    # Equivalently 4pi <P2 K1> <P2 K3>. Both k^2*k^4 orderings retained.
    p12, p14 = sky(P2 * linear2), sky(P2 * linear4)
    # Also check the linear and quadratic coefficients from this same geometric
    # map, independently of the earlier hand script and of the fast generator.
    c11next = s.factor(8 * s.pi * p12 * p14)
    expected11 = -4 * s.pi * (o - e) ** 2 * (o + e) * (59 * o**3 + 59 * o**2 * e - 106 * o * e**2 + 114 * e**3) / 637875
    assert s.expand(c11next - expected11) == 0
    k22 = leading.coeff(3)
    k24 = sum(time * value.coeff(3) for time, value in pieces.values())
    c22lead = s.factor(2 * s.pi * qaverage(sky(k22 * k22)) / 5)
    c22next = s.factor(4 * s.pi * qaverage(sky(k22 * k24)) / 5)
    expected22 = -8 * s.pi * Delta * (105 * ti**4 + 719 * o**4 - 1920 * o**2 * e**2 + 2560 * o * e**3 - 1464 * e**4) / 637875
    assert s.expand(c22lead - 68 * s.pi * Delta**2 / 6075) == 0
    assert s.expand(c22next - expected22) == 0
    print("[ok] independent C11 and C22/4 coefficients through k0^6", flush=True)
    c13lead = s.factor(4 * s.pi * p12 * cubic2)
    c13next = s.factor(4 * s.pi * (p12 * cubic4 + p14 * cubic2))
    expected = -2 * s.pi * Delta * (2835 * ti**4 + 9973 * o**4 - 25440 * o**2 * e**2 + 33920 * o * e**3 - 21288 * e**4) / 1913625
    print("independent C13/3 k^4 coefficient:", c13lead)
    print("independent C13/3 k^6 coefficient:", c13next)
    print("difference from report:", s.factor(c13next - expected))
    assert s.expand(c13lead - 16 * s.pi * Delta**2 / 2025) == 0
    assert s.expand(c13next - expected) == 0
    betaH = -11 * L / 14 - 3 * (o**4 - e**4) / 40
    betaa = 27 * L / 7 + 7 * (o**4 - e**4) / 20 + (o**4 + 3 * ti**4) / 4 - e**4
    assert s.expand(c11next - 16 * s.pi * Delta * betaH / 6075) == 0
    assert s.expand(c22next - 32 * s.pi * Delta * (30 * betaH - betaa) / 18225) == 0
    assert s.expand(c13next - 8 * s.pi * Delta * (110 * betaH - 9 * betaa) / 18225) == 0
    assert s.expand(c22next + c13next - 8 * s.pi * Delta * (230 * betaH - 13 * betaa) / 18225) == 0
    # The third tensor does not contribute to the quartic spectrum at this order.
    assert qaverage(sky(k22 * nn(A2TF).coeff(3))) == 0
    assert qaverage(sky(P2 * nn(A2TF).coeff(7))) == 0
    print("[ok] compact k0^6 spectrum depends on betaH and betaa only; the A^2 tensor drops out after Wick contraction")
    # Negative controls: remove each nonzero cubic contribution separately.
    for name in ("spatial Hessian", "fourth-gradient field", "emission", "clock"):
        assert s.expand(4 * s.pi * p12 * projected[name]) != 0
    assert s.expand(4 * s.pi * p14 * cubic2) != 0  # do not omit the K1^[4] K3^[2] ordering
    assert projected["direction A^2"] == 0
    assert projected["direction+normalization Ann^2"] == 0
    print("[ok] negative controls detect omitted Hessian/field/endpoint/clock and linear k^4 term")
    print("[ok] complete cubic k0^6 coefficient, with symbolic e,o,i; no generator/cache imports")


if __name__ == "__main__":
    main()
