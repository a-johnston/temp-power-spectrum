"""Exact check: the redshift driver in fully orthonormal photon directions.

With gamma = e^{2 zeta} e^H and the orthonormal direction v = e^{zeta} e^{H/2} e (v^T v = 1),
    sigma_ij e^i e^j = v^T B v,   B = (1/2) e^{-H/2} (e^H)' e^{-H/2},
and, for arbitrary non-commuting symmetric H(t) = eps H1 + eps^2 H2/2 + eps^3 H3/6,
    B = (1/2) H' + (1/48) [H, [H, H']] + O(eps^4):
no quadratic matrix product appears, and the e^{2 zeta} factor drops out exactly. On the
monochromatic shell H1(t, X) = f(t) B(X), so [H1, H1'] = 0 and B = H'/2 + O(eps^4).
"""

import sympy as sp

eps = sp.Symbol("epsilon")
A1, A2, A3, D1, D2, D3 = sp.symbols("A1 A2 A3 D1 D2 D3", commutative=False)  # H_r and H_r'
H = eps * A1 + eps**2 * A2 / 2 + eps**3 * A3 / 6
Hp = eps * D1 + eps**2 * D2 / 2 + eps**3 * D3 / 6


def trunc(e, order=3):
    e = sp.expand(e)
    return sp.Add(*[e.coeff(eps, k) * eps**k for k in range(order + 1)])


def expm(X, n=4):
    return trunc(sum(X**k / sp.factorial(k) for k in range(n)))


def dexpm(X, Xp, n=4):
    """d/dt e^X = sum_k 1/k! sum_j X^j X' X^{k-1-j}."""
    return trunc(sum(sum(X**j * Xp * X ** (k - 1 - j) for j in range(k)) / sp.factorial(k) for k in range(1, n)))


B = trunc(expm(-H / 2) * dexpm(H, Hp) * expm(-H / 2) / 2)
comm = lambda X, Y: X * Y - Y * X
target = trunc(Hp / 2 + comm(H, comm(H, Hp)) / 48)
assert sp.expand(B - target) == 0
print("[ok] B = (1/2) e^{-H/2} (e^H)' e^{-H/2} = H'/2 + [H,[H,H']]/48 + O(eps^4) for non-commuting H")
shell = {D1: sp.Symbol("fp") * sp.Symbol("Bm", commutative=False), A1: sp.Symbol("f") * sp.Symbol("Bm", commutative=False)}
assert sp.expand(comm(H, comm(H, Hp)).subs(shell)).coeff(eps, 3) == 0
print("[ok] shell: H1 = f(t) B(X) makes the cubic commutator vanish, so B = H'/2 + O(eps^4)")
print("ORTHONORMAL DRIVER CHECK PASSED (exact symbolic algebra)")
