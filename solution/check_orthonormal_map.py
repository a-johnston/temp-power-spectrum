"""Exact check: the second-order map in orthonormal photon directions v = e^zeta e^{H/2} e.

Plane-wave mode a: phase ph_a = e^{i k_a . n (t_o - t)} on the background ray, primordial amplitude z_a,
first-order metric H_a(t) (symmetric, trace-free) with derivative H_a'(t).

1. First-order direction. With v_a = e_a - ph_a (z_a n + H_a n / 2) and the report's ray source
   A_a = H_a' n - Gamma_a(n,n) - n sigma_a,nn, the total derivative along the ray is
       d v_a / dt = ph_a [ (1/2) P_perp H_a' n + P_perp (i z_a k_a) + (i/2) (k_a H_a,nn - (k_a.n) H_a n) ],
   which is transverse (n . dv_a/dt = 0), and v_a(t_o) = 0 (symmetric frame).
2. Second-order integrand. With the report's sigma_ab = (1/2) H_ab' + (1/4)(H_a H_b + H_b H_a)' + z_a H_b' + z_b H_a',
       sigma_ab,nn - 2 sigma_b(n, e_a) - 2 sigma_a(n, e_b) = (1/2) H_ab',nn - H_b'(n, v_a) - H_a'(n, v_b):
   the local products cancel identically against e - v (generic symmetric matrices).
"""

import sympy as sp

n = sp.Matrix(sp.symbols("n1:4", real=True))
unit = {n[2]: sp.sqrt(1 - n[0] ** 2 - n[1] ** 2)}


def sym(prefix):
    e = {(i, j): sp.Symbol(f"{prefix}{min(i, j)}{max(i, j)}") for i in range(3) for j in range(3)}
    return sp.Matrix(3, 3, lambda i, j: e[(i, j)])


def check_direction():
    k = sp.Matrix(sp.symbols("k1:4", real=True))
    z = sp.Symbol("z")
    Hm, Hp = sym("H"), sym("Hp")
    kn = (k.T * n)[0]
    nn = lambda M: (n.T * M * n)[0]
    Gam = sp.I * z * (2 * kn * n - k) + sp.I / 2 * (2 * kn * Hm * n - k * nn(Hm))
    A = Hp * n - Gam - n * nn(Hp) / 2
    # d/dt [ph (z n + H n / 2)] = ph [ -i (k.n)(z n + H n/2) + H' n / 2 ]
    dv = A - (-sp.I * kn * (z * n + Hm * n / 2) + Hp * n / 2)
    Pp = sp.eye(3) - n * n.T
    target = Pp * Hp * n / 2 + Pp * (sp.I * z * k) + sp.I / 2 * (k * nn(Hm) - kn * Hm * n)
    assert sp.simplify((dv - target).subs(unit)) == sp.zeros(3, 1)
    assert sp.simplify((n.T * target)[0].subs(unit)) == 0
    print("[ok] first-order orthonormal direction: dv/dt = P_perp(H'n/2 + i z k) + (i/2)(k H_nn - (k.n) H n), transverse")


def check_integrand():
    Ha, Hb, Hap, Hbp, Habp = sym("Ha"), sym("Hb"), sym("Hap"), sym("Hbp"), sym("Habp")
    za, zb = sp.symbols("za zb")
    va = sp.Matrix(sp.symbols("va1:4"))
    vb = sp.Matrix(sp.symbols("vb1:4"))
    nn = lambda M: (n.T * M * n)[0]
    bil = lambda M, x, y: (x.T * M * y)[0]
    sig_ab = Habp / 2 + (Hap * Hb + Ha * Hbp + Hbp * Ha + Hb * Hap) / 4 + za * Hbp + zb * Hap
    sa, sb = Hap / 2, Hbp / 2
    ea = va + za * n + Ha * n / 2
    eb = vb + zb * n + Hb * n / 2
    old = nn(sig_ab) - 2 * bil(sb, n, ea) - 2 * bil(sa, n, eb)
    new = nn(Habp) / 2 - bil(Hbp, n, va) - bil(Hap, n, vb)
    assert sp.expand(old - new) == 0
    print("[ok] second-order integrand: sigma_ab,nn - 2 sigma(n,e) terms = H_ab'_nn/2 - H'(n, v) terms (local products cancel)")


if __name__ == "__main__":
    check_direction()
    check_integrand()
    print("ORTHONORMAL MAP CHECKS PASSED (exact symbolic algebra)")
