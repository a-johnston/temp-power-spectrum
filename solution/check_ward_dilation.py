"""Dilation Ward identity for the second-order fields (squeezed limit, exact).

A constant zeta_L is removed by the relabelling X -> e^{zeta_L} X, so the exact solution for
zeta_s + zeta_L is the solution for zeta_s(e^{-zeta_L} X) in rescaled labels.  For the kernels:
    K_2(k, q = 0) = -k . d_k K_1(k)          for nu and for every component of H.
Here the long mode is exactly constant (q = 0: its nu, H vanish; only zeta_L = 1).  We insert the
predicted K_2 into ALL equations (N), (M), (T), (Ham) at order z_s z_L and verify exact
cancellation in closed form; the prediction also lies in the regular class (nu_2 = 2k^2 t^2/9
+ O(t^4) = c2 t^2 with c2 = [2(k^2 + 0) - 0]/9, H_2(0) = 0), so by uniqueness it IS the solution.
"""
import sympy as sp
from engine import Field, Ctx, equations, t

p = sp.Symbol('p', positive=True)
s3 = sp.sqrt(3)
Ctx.maxdeg = 2
Ctx.momenta = {1: (p, 0, 0), 2: (0, 0, 0)}


def K1(kk):
    x = kk * t / s3
    f = sp.sin(x) / x - sp.cos(x)
    s = 2 * sp.sin(x) / x - sp.cos(x)
    Q = -6 * (s - 1) / kk**2
    nu = -f
    H = sp.Matrix(3, 3, lambda i, j: Q * ((kk if i == 0 else 0) * (kk if j == 0 else 0) - (kk**2 if i == j else 0) / 3))
    return nu, H


kk = sp.Symbol('kk', positive=True)
nu1, H1 = K1(kk)
dil = lambda e: (-kk * sp.diff(e, kk)).subs(kk, p)
nu12 = dil(nu1)
H12 = H1.applyfunc(dil)
nu1p, H1p = nu1.subs(kk, p), H1.subs(kk, p)

nu = Field({frozenset([1]): nu1p, frozenset([1, 2]): nu12})
zeta = Field({frozenset([1]): 1, frozenset([2]): 1})
H = [[Field({frozenset([1]): H1p[i, j], frozenset([1, 2]): H12[i, j]}) for j in range(3)] for i in range(3)]
E = equations(nu, zeta, H)
lab = [1, 2]
def z(e, name):
    r = sp.simplify(sp.expand(e))
    assert r == 0, (name, r)
z(E['N'].coeff(lab), 'N')
for i in range(3):
    z(E['M'][i].coeff(lab), f'M{i}')
    for j in range(3):
        z(E['T'][i][j].coeff(lab), f'T{i}{j}')
z(E['Ham'].coeff(lab), 'Ham')
lead = sp.series(nu12, t, 0, 3).removeO()
assert sp.simplify(lead - 2 * p**2 * t**2 / 9) == 0
assert all(sp.limit(e, t, 0) == 0 for e in H12)
print('[ok] dilation Ward identity: K_2(k, 0) = -k.d_k K_1(k) solves every equation and the regular class')
