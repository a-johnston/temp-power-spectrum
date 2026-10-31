"""Engine check: the known regular first-order solution satisfies (N), (M), (T), (Ham) exactly,
for a wave with GENERIC direction (k1,k2,k3)."""
import sympy as sp
from engine import Field, Ctx, equations, t

k1, k2, k3 = sp.symbols('k1 k2 k3', real=True)
k = sp.Symbol('k', positive=True)
Ctx.maxdeg = 1
Ctx.momenta = {1: (k1, k2, k3)}
kv = (k1, k2, k3)
x = k * t / sp.sqrt(3)


def first_order(label, kvec, kmag):
    xx = kmag * t / sp.sqrt(3)
    f = sp.sin(xx) / xx - sp.cos(xx)
    s = 2 * sp.sin(xx) / xx - sp.cos(xx)
    nu = Field({frozenset([label]): -f})
    zeta = Field({frozenset([label]): sp.Integer(1)})
    Q = -6 * (s - 1) / kmag**2
    H = [[Field({frozenset([label]): Q * (kvec[i] * kvec[j] - (kmag**2 if i == j else 0) / 3)})
          for j in range(3)] for i in range(3)]
    return nu, zeta, H


if __name__ == '__main__':
    nu, zeta, H = first_order(1, kv, k)
    E = equations(nu, zeta, H)
    on = {k1**2: k**2 - k2**2 - k3**2}
    def chk(F, lab):
        c = sp.expand(F.coeff([1]))
        c = sp.simplify(c.subs(k, sp.sqrt(k1**2 + k2**2 + k3**2)))
        assert c == 0, (lab, c)
    chk(E['N'], 'N')
    for i in range(3):
        chk(E['M'][i], f'M{i}')
        for j in range(3):
            chk(E['T'][i][j], f'T{i}{j}')
    chk(E['Ham'], 'Ham')
    for F in [E['N'], E['Ham']] + E['M']:
        assert sp.simplify(F.coeff([])) == 0
    print('[ok] first-order regular solution satisfies (N),(M),(T),(Ham) for generic direction')
