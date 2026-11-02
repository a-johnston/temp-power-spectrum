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


def check_tensor_source_against_literature_form():
    """Ties the independent Kohri-Terada comparison (check_tensor_newtonian.py) to THIS repo's source:
    the TT part of second_order.py's c_T on the shell equals -(1/2) S_T (v v - w w), with S_T the
    independently derived shell tensor source, identically in t and mu."""
    src = so.sources()
    cT = src["T"]
    shell = {p: 1, q: 1, b1: mu}
    KV = sp.Matrix([1 + mu, b2, 0])
    K2 = 2 + 2 * mu
    Pm = sp.eye(3) - KV * KV.T / K2
    C = cT.applyfunc(lambda e: e.subs(shell))
    PMP = Pm * C * Pm
    TT = PMP - Pm * PMP.trace() / 2
    vv = sp.Matrix([1 - mu, -b2, 0])  # k1 - k2, |k1 - k2|^2 = 2 - 2 mu
    eplus = vv * vv.T / (2 - 2 * mu) - sp.Matrix(3, 3, lambda i, j: 1 if i == j == 2 else 0)
    amp = sum(TT[i, j] * eplus[i, j] for i in range(3) for j in range(3)) / 2  # TT = amp * e_plus
    # independent source (report/tensor_newtonian_crosscheck.md, Sec. 5), k0 = 1
    x = t / s3
    s_ = 2 * sp.sin(x) / x - sp.cos(x)
    B = 3 * (s_ - 1)
    h = -2 * B
    n = sp.cos(x) - sp.sin(x) / x
    ST = (
        (1 - mu)
        / 36
        * (
            h**2 * (9 * mu**2 + 12 * mu + 1)
            - 6 * h * (6 * mu * n + 3 * mu + 2)
            - 36 * sp.diff(h, t) * sp.diff(n, t)
            - 72 * n**2
            + 36
        )
    )
    diff = sp.expand(amp + ST / 2)
    Pz = sp.Symbol("Pz")
    e = sp.expand(diff.rewrite(sp.exp))
    rep = {}
    for f in e.atoms(sp.exp):
        a = sp.nsimplify(sp.expand(f.args[0] / (sp.I * t / s3)))
        assert a.is_integer, f
        rep[f] = Pz ** int(a)
    e = sp.expand(sp.together(e.xreplace(rep)))
    num = sp.expand(sp.fraction(sp.together(e))[0])
    poly = sp.Poly(num, b2)
    red = sum(c * (b2 ** (m % 2)) * (1 - mu**2) ** (m // 2) for (m,), c in poly.terms())
    assert sp.cancel(sp.expand(red)) == 0
    print("[ok] repo c_T (second_order.py) has TT part -(1/2) S_T (v v - w w): matches the independent KT-comparison source")


if __name__ == "__main__":
    check()
    check_tensor_source_against_literature_form()
    print("SHELL REDUCTION CHECKS PASSED (exact symbolic algebra)")
