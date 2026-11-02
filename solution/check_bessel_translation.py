"""Exact check of the separated (translation) form of the mixed-time sky integrals.

For V = a khat + b qhat with khat = z, khat.qhat = mu, the Funk-Hecke form
    int dOmega P_l(z.n) e^{i n.V} = 4 pi i^l j_l(|V|) P_l(z.Vhat),   |V| = sqrt(a^2 + b^2 + 2 a b mu),
contains the mixed-time square root. Expanding each plane wave separately gives
    = 4 pi sum_{L1,L2} (2L1+1)(2L2+1) i^{L1+L2} j_{L1}(a) j_{L2}(b) (l L1 L2; 0 0 0)^2 P_{L2}(mu),
with single-argument Bessel functions and mu only in the polynomial P_{L2}(mu).
Both sides are entire in (a, b); the identity is checked exactly order by order through total
degree N in (a, b), with symbolic mu, for l = 0..4. A negative control drops the 3j weight.
"""

import sympy as sp
from sympy.physics.wigner import wigner_3j

a, b, mu, x = sp.symbols("a b mu x", real=True)
N = 10


def jl_series(L, arg, order):
    """Series of j_L(arg) through arg^order: j_L(x) = sum_n (-1)^n x^{L+2n} / (2^n n! (2L+2n+1)!!)."""
    out = 0
    n = 0
    while L + 2 * n <= order:
        out += (-1) ** n * arg ** (L + 2 * n) / (2**n * sp.factorial(n) * sp.factorial2(2 * L + 2 * n + 1))
        n += 1
    return out


def truncate(e, order):
    e = sp.expand(e)
    P = sp.Poly(e, a, b)
    return sp.Add(*[c * a**i * b**j for (i, j), c in P.terms() if i + j <= order])


def lhs(l):
    # j_l(R) P_l(z/R) with R^2 = a^2+b^2+2ab mu, z = a + b mu, written without square roots:
    # j_l(R) = R^l * f(R^2) and R^l P_l(z/R) is a polynomial in z and R^2.
    R2 = a**2 + b**2 + 2 * a * b * mu
    z = a + b * mu
    f = sum((-1) ** n * R2**n / (2**n * sp.factorial(n) * sp.factorial2(2 * l + 2 * n + 1)) for n in range(N // 2 + 1))
    t = sp.Symbol("t")
    Pl = sp.legendre(l, t)
    solid = sp.expand(sp.Poly(Pl, t).as_expr())
    # R^l P_l(z/R) = sum_k c_k z^k R^{l-k}, with only l-k even
    sh = sum(c * z**k * R2 ** ((l - k) // 2) for (k,), c in sp.Poly(solid, t).terms())
    return truncate(4 * sp.pi * sp.I**l * f * sh, N)


def rhs(l, weight=True):
    out = 0
    for L1 in range(N + 1):
        for L2 in range(N + 1 - L1):
            w = wigner_3j(l, L1, L2, 0, 0, 0)
            if w == 0:
                continue
            w2 = w**2 if weight else 1
            out += (
                (2 * L1 + 1)
                * (2 * L2 + 1)
                * sp.I ** (L1 + L2)
                * jl_series(L1, a, N)
                * jl_series(L2, b, N)
                * w2
                * sp.legendre(L2, mu)
            )
    return truncate(4 * sp.pi * out, N)


if __name__ == "__main__":
    for l in range(5):
        assert sp.expand(lhs(l) - rhs(l)) == 0, l
    print(f"[ok] translation form = Funk-Hecke form through total degree {N} in (a, b), l = 0..4, symbolic mu")
    assert sp.expand(lhs(2) - rhs(2, weight=False)) != 0
    print("[ok] negative control: dropping the 3j weight fails")
    print("BESSEL TRANSLATION CHECK PASSED (exact)")
