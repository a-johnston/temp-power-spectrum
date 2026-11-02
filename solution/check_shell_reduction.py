"""Exact symmetry reductions of the second-order sources (all finite t, all mu).

Frame of second_order.py: k1 = (p, 0, 0), k2 = (b1, b2, 0), K = k1 + k2.

1. Reflection z -> -z (both momenta in the xy plane): c_M,z = 0 and c_T,xz = c_T,yz = 0 for general
   p, q, b1, b2. Hence the vector response lies in the plane and the TT response has a single
   polarisation, v v - w w with v = unit(k1 - k2) projected, w = z.
2. Equal-magnitude shell |k1| = |k2|: K x c_M = 0 identically in t and mu, so P_perp c_M = 0 and the
   second-order vector response vanishes for K != 0. Together with 1 and the scalar projector,
       A_2 = F_S (u u - delta/3) + F_T (v v - w w),   u = K/|K|.
   A negative control (|k1| != |k2|) confirms that the vector source is generically nonzero.
"""

import second_order as so
import sympy as sp

p, q, b1, b2 = so.p, so.q, so.b1, so.b2
t = [s for s in so.sources()["N"].free_symbols if s.name == "t"][0]
mu = sp.Symbol("mu", real=True)
s3 = sp.sqrt(3)


def check():
    src = so.sources()
    cM, cT = src["M"], src["T"]
    # 1. reflection symmetry: general momenta
    for e in (cM[2], cT[0, 2], cT[1, 2]):
        assert sp.simplify(e) == 0
    print("[ok] reflection: c_M,z = c_T,xz = c_T,yz = 0 for general momenta (single TT polarisation)")
    # 2. shell: K x c_M = 0  (p = q = 1, b1 = mu, b2^2 = 1 - mu^2), all t and mu
    shell = {p: 1, q: 1, b1: mu}
    KV = (1 + mu, b2, 0)
    cross = sp.expand((KV[0] * cM[1] - KV[1] * cM[0]).subs(shell))
    Pz = sp.Symbol("Pz")
    e = sp.expand(cross.rewrite(sp.exp))
    rep = {}
    for f in e.atoms(sp.exp):
        a = sp.nsimplify(sp.expand(f.args[0] / (sp.I * t / s3)))
        assert a.is_integer, f
        rep[f] = Pz ** int(a)
    e = sp.expand(e.xreplace(rep))
    poly = sp.Poly(e, b2)
    e = sum(c * (b2 ** (m % 2)) * (1 - mu**2) ** (m // 2) for (m,), c in poly.terms())
    assert sp.cancel(sp.together(sp.expand(e))) == 0
    print("[ok] shell: K x c_M = 0 identically in t and mu  =>  P_perp c_M = 0, no second-order vector")
    # negative control: |k2| = 2 |k1|
    m0 = sp.Rational(1, 3)
    off = {p: 1, q: 2, b1: 2 * m0, b2: 2 * sp.sqrt(1 - m0**2), t: sp.Rational(7, 5)}
    KVoff = (1 + 2 * m0, 2 * sp.sqrt(1 - m0**2))
    val = (KVoff[0] * cM[1] - KVoff[1] * cM[0]).subs(off)
    assert abs(complex(sp.N(val))) > 1e-6
    print("[ok] negative control: off the shell the vector source is nonzero")
    print("SHELL REDUCTION CHECKS PASSED (exact symbolic algebra)")


if __name__ == "__main__":
    check()
