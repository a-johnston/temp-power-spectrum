"""Third-order dilation Ward identity for the fields (two exactly constant long modes).

With zeta_L = z2 + z3 (zero momentum), exp(-zeta_L k.d_k) acting on the short-mode solution gives
    K_2(k,0) = -k.d_k K_1,   K_2(0,0) = 0,   K_3(k,0,0) = (k.d_k)^2 K_1.
The prediction is inserted into ALL equations at order z1 z2 z3 and must cancel exactly.
"""

import sympy as sp
from engine import Ctx, Field, equations, t

p = sp.Symbol("p", positive=True)
kk = sp.Symbol("kk", positive=True)
s3 = sp.sqrt(3)
Ctx.maxdeg = 3
Ctx.momenta = {1: (p, 0, 0), 2: (0, 0, 0), 3: (0, 0, 0)}
x = kk * t / s3
nu1 = -(sp.sin(x) / x - sp.cos(x))
Q = -6 * (2 * sp.sin(x) / x - sp.cos(x) - 1) / kk**2
H1 = sp.Matrix(3, 3, lambda i, j: Q * ((kk if i == 0 else 0) * (kk if j == 0 else 0) - (kk**2 if i == j else 0) / 3))
D = lambda e: -kk * sp.diff(e, kk)
at = lambda e: e.subs(kk, p)
nu = Field(
    {frozenset([1]): at(nu1), frozenset([1, 2]): at(D(nu1)), frozenset([1, 3]): at(D(nu1)), frozenset([1, 2, 3]): at(D(D(nu1)))}
)
zeta = Field({frozenset([1]): 1, frozenset([2]): 1, frozenset([3]): 1})
H = [
    [
        Field(
            {
                frozenset([1]): at(H1[i, j]),
                frozenset([1, 2]): at(D(H1[i, j])),
                frozenset([1, 3]): at(D(H1[i, j])),
                frozenset([1, 2, 3]): at(D(D(H1[i, j]))),
            }
        )
        for j in range(3)
    ]
    for i in range(3)
]
E = equations(nu, zeta, H)
lab = [1, 2, 3]


def z(e, name):
    r = sp.simplify(sp.expand(sp.sympify(e)))
    assert r == 0, (name, r)


z(E["N"].coeff(lab), "N")
for i in range(3):
    z(E["M"][i].coeff(lab), f"M{i}")
    for j in range(3):
        z(E["T"][i][j].coeff(lab), f"T{i}{j}")
z(E["Ham"].coeff(lab), "Ham")
lead = sp.series(at(D(D(nu1))), t, 0, 3).removeO()
c3 = -sp.Rational(2, 9) * (2 * p**2 + 5 * 0)  # c3 at (k, q, -q) with q = 0
assert sp.simplify(lead - c3 * t**2) == 0, lead  # t^2 coefficient = c3, fixed by (Ham)
print("[ok] third-order dilation Ward identity: K_3(k,0,0) = (k.d_k)^2 K_1 solves every equation; t^2 coefficient = c3")
