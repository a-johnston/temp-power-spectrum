"""The short derivation of the K=0 pair fields: sources of (Ham) and (T), and selection of the regular solution.

q = 1, x = t/sqrt3. At K=0, (Ham) is first order in n2 and (T) is second order in a:
    (t n2)' = -t^2 c_H / 12,     (t^2 a')' = -t^2 S_T.
Regularity (nu, H even and analytic, vanishing at t=0) removes the homogeneous solutions n2 ~ 1/t, a ~ 1/t and a ~ const.
"""

import k0_fields as kf
import sympy as sp
from plane_reduction import f, nu, reduced_equations, t, y, zeta

X = kf.X
s3 = kf.s3


def pair_source(key):
    eq = reduced_equations()
    z1, z2, eps = sp.symbols("z1 z2 eps")
    E = sp.exp(sp.I * y)
    g = z1 * E + z2 / E
    e = eq[key].subs({f: eps * kf.to_t(kf.h1) * g, nu: eps * kf.to_t(kf.n1) * g, zeta: eps * g}).doit()
    c = sp.expand(sp.diff(e, eps, 2).subs(eps, 0) / 2).coeff(z1, 1).coeff(z2, 1)
    return sp.simplify(sp.expand(sp.expand_trig(sp.expand(c.subs(t, s3 * X)))))


if __name__ == "__main__":
    print("reduced equation forms:")
    c_H = pair_source("Ham")
    print("c_H =", c_H)
    # (Ham) at K=0, order eps^2: 12 n2'/t + 12 n2/t^2 + c_H = 0, so (t n2)' = -t^2 c_H/12; in x: (x n2)_x = -x^2 c_H / 4
    res = sp.simplify(sp.expand(sp.expand_trig(sp.diff(X * kf.n2, X) + X**2 * c_H / 4)))
    print("(x n2)_x + x^2 c_H/4 =", res)
    assert res == 0
    print("x^2 c_H / 4 =", sp.simplify(sp.expand(sp.expand_trig(X**2 * c_H / 4))))
    # small-x data
    print("c_H series:", sp.series(c_H, X, 0, 4))
    print("PASS")
