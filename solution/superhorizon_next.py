"""Next order of the superhorizon expansion: map kernels through O(kappa^4) from the fast generator.

    python superhorizon_next.py kernels2      K1(k), K2(k,q) through O(kappa^4)      (saved to results/)
    python superhorizon_next.py kernels3      K3(k,q,-q) through O(kappa^4)          (saved to results/)
    python superhorizon_next.py spectrum      C_2 through O(k0^6) at fourth order from the saved kernels

All steps are exact (sparse polynomials over Q(i), then SymPy). The kernels are checked for
locality: beyond leading gradient order the FIELD kernels contain the non-local factors
u = 1/|k+q|^2, v = 1/|k-q|^2, and the observable must not.
"""

import os
import pickle
import sys
import time

KMAX = 4
MODE = sys.argv[1] if len(sys.argv) > 1 else "kernels2"
sys.argv = [sys.argv[0], "SH", str(KMAX)]
import observable_fast as of  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "results")
T0 = time.time()


def save(name, Th, labels_list):
    data = {}
    for labels in labels_list:
        poly = Th.coeff(labels)
        data[tuple(labels)] = {m: str(c) for m, c in poly.items()}  # exact Gaussian rationals as strings
    path = os.path.join(OUT, f"superhorizon_kernels_{name}.pkl")
    with open(path, "wb") as fh:
        pickle.dump({"gens": [str(g) for g in of.R.gens], "kernels": data, "KMAX": KMAX}, fh)
    return path


def locality(Th, labels):
    poly = Th.coeff(labels)
    iu, iv = of.R.gens.index(of.U), of.R.gens.index(of.V)
    bad = [m for m in poly if m[iu] or m[iv]]
    return len(poly), len(bad)


if __name__ == "__main__":
    if MODE == "kernels2":
        nu, zeta, H = of.build_fields([1, 2])
        of.LOG("fields ready")
        Th = of.observable(nu, zeta, H)
        for lab in ([1], [1, 2]):
            n, b = locality(Th, lab)
            print(f"  kernel {lab}: {n} terms, {b} with non-local 1/|k+-q|^2", flush=True)
        print("  saved", save("2", Th, ([1], [1, 2])))
    elif MODE == "kernels3":
        nu, zeta, H = of.build_fields([1, 2, 3])
        of.LOG("fields ready")
        Th = of.observable(nu, zeta, H)
        for lab in ([1], [1, 2, 3]):
            n, b = locality(Th, lab)
            print(f"  kernel {lab}: {n} terms, {b} with non-local 1/|k+-q|^2", flush=True)
        print("  saved", save("3", Th, ([1], [1, 2, 3])))
    elif MODE == "spectrum":
        pass
    else:
        raise SystemExit(f"unknown mode {MODE}")
    print(f"  done in {time.time() - T0:.0f}s")


# ------------------------------------------------------------------ spectrum from saved kernels
def load(name):
    import sympy as sp

    with open(os.path.join(OUT, f"superhorizon_kernels_{name}.pkl"), "rb") as fh:
        d = pickle.load(fh)
    gens = sp.symbols(d["gens"])
    G = dict(zip(d["gens"], gens, strict=True))
    p, b1, b2 = G["p"], G["b1"], G["b2"]
    sub = {
        G["ip"]: 1 / p, G["u"]: 1 / ((p + b1) ** 2 + b2**2), G["v"]: 1 / ((p - b1) ** 2 + b2**2),
        G["its"]: 1 / G["ts"], G["itI"]: 1 / G["tI"], G["itob"]: 1 / G["tob"],
    }  # fmt: skip
    out = {}
    for lab, terms in d["kernels"].items():
        e = sp.Add(*[sp.sympify(c) * sp.Mul(*[g**x for g, x in zip(gens, m, strict=True)]) for m, c in terms.items()])
        out[lab] = sp.expand(sp.cancel(sp.together(e.xreplace(sub))))
    return out, G


def sphere_moment(a, b, c):
    import sympy as sp

    if a % 2 or b % 2 or c % 2:
        return sp.S.Zero
    return (
        2
        * sp.gamma(sp.Rational(a + 1, 2))
        * sp.gamma(sp.Rational(b + 1, 2))
        * sp.gamma(sp.Rational(c + 1, 2))
        / sp.gamma(sp.Rational(a + b + c + 3, 2))
    )


def moments2(f, G):
    """M0 = int f dOmega and M_ij = int f n_i n_j dOmega for a polynomial f(n)."""
    import sympy as sp

    n = (G["n1"], G["n2"], G["n3"])
    P = sp.Poly(sp.expand(f), *n)
    M0 = 0
    M = [[0] * 3 for _ in range(3)]
    for mon, c in P.terms():
        M0 += c * sphere_moment(*mon)
        for i in range(3):
            for j in range(3):
                e = list(mon)
                e[i] += 1
                e[j] += 1
                M[i][j] += c * sphere_moment(*e)
    return M0, M


def S2(f, g, G):
    """sum_m f_2m g_2m^* = (5/4pi) int int f(n) g(n')^* P_2(n.n')."""
    import sympy as sp

    f0, F = moments2(f, G)
    g0, Gm = moments2(sp.expand(g).subs(sp.I, -sp.I), G)  # all symbols are real
    return sp.expand(
        5 / (4 * sp.pi) * (sp.Rational(3, 2) * sum(F[i][j] * Gm[i][j] for i in range(3) for j in range(3)) - f0 * g0 / 2)
    )


def truncate_kappa(e, kap, order):
    import sympy as sp

    e = sp.expand(e)
    return sp.Add(*[e.coeff(kap, k) * kap**k for k in range(order + 1)])


def spectrum():
    import sympy as sp

    K2s, G = load("2")
    kap, p, b1, b2, ts, tob = G["kap"], G["p"], G["b1"], G["b2"], G["ts"], G["tob"]
    mu, sq = sp.symbols("mu sq", real=True)
    K1 = K2s[(1,)]
    C11 = truncate_kappa(S2(K1, K1, G) / 5, kap, 6)  # sigma^2 = 1
    # independent check of the first-order kernel: small-k expansion of T_2 from its integral definition
    k, t, x = sp.symbols("k t x", positive=True)
    y = sp.Symbol("y", positive=True)
    j2 = sp.expand_func(sp.jn(2, y))
    w2 = sp.series(j2 / 3 + sp.diff(j2, y, 2), y, 0, 6).removeO()
    sprime = sp.series(sp.diff(2 * sp.sin(x) / x - sp.cos(x), x), x, 0, 6).removeO()
    T2 = sp.integrate(
        sp.expand(-3 * sprime.subs(x, k * t / sp.sqrt(3)) * w2.subs(y, k * (tob - t)) * k / sp.sqrt(3)), (t, ts, tob)
    )
    C11_T = truncate_kappa(sp.expand(4 * sp.pi * T2**2).subs(k, kap * p), kap, 6)
    assert sp.simplify(C11 - C11_T) == 0, sp.factor(C11 - C11_T)
    print("[ok] generator first-order quadrupole through O(k0^6) equals 4 pi T_2^2 from the transfer integral")
    K2 = K2s[(1, 2)]
    shell = lambda e: sp.expand(sp.expand(e.subs({b1: p * mu, b2: p * sq})).subs(sq**2, 1 - mu**2))
    dens = shell(S2(K2, K2, G))
    assert sp.expand(dens.coeff(sq, 1)) == 0, "odd powers of sqrt(1-mu^2) must cancel (reflection symmetry)"
    dens = dens.coeff(sq, 0)
    quarterC22 = truncate_kappa(sp.integrate(dens, (mu, -1, 1)) / 5 / 4, kap, 6)
    third = None
    if os.path.exists(os.path.join(OUT, "superhorizon_kernels_3.pkl")):
        K3s, _ = load("3")  # same symbol names, hence identical SymPy symbols
        assert sp.expand(K3s[(1,)] - K1) == 0, "first-order kernels of the two runs disagree"
        K3 = K3s[(1, 2, 3)]
        # (1/3) Re C13 = sigma^4/(2l+1) * (1/2) int dmu Re S_l[K1(k), K3(k,q,-q)]
        d13 = shell(S2(K1, K3, G))
        assert sp.expand(d13.coeff(sq, 1)) == 0, "odd powers of sqrt(1-mu^2) must cancel"
        d13 = sp.re(sp.expand(d13.coeff(sq, 0)))
        third = truncate_kappa(sp.integrate(d13, (mu, -1, 1)) / 2 / 5, kap, 6)
    return C11, quarterC22, third, G


if __name__ == "__main__" and MODE == "spectrum":
    import sympy as sp

    C11, qC22, tC13, G = spectrum()
    kap = G["kap"]
    print("  C11  kappa^4:", sp.factor(C11.coeff(kap, 4)))
    print("  C11  kappa^6:", sp.factor(C11.coeff(kap, 6)))
    print("  C22/4 kappa^4:", sp.factor(qC22.coeff(kap, 4)), "  ratio", sp.simplify(qC22.coeff(kap, 4) / C11.coeff(kap, 4)))
    print("  C22/4 kappa^6:", sp.factor(qC22.coeff(kap, 6)))
    if tC13 is not None:
        L4 = C11.coeff(kap, 4)
        print("  C13/3 kappa^4:", sp.factor(tC13.coeff(kap, 4)), "  ratio", sp.simplify(tC13.coeff(kap, 4) / L4))
        print("  C13/3 kappa^6:", sp.factor(tC13.coeff(kap, 6)))
        tot4 = sp.simplify((qC22.coeff(kap, 4) + tC13.coeff(kap, 4)) / L4)
        assert tot4 == sp.Rational(29, 3), tot4
        print("[ok] leading fourth-order factor 29/3 reproduced from the kappa^4 kernels")
        print("  complete eps^4 k0^6 coefficient (1/4 C22 + 1/3 Re C13):", sp.factor(qC22.coeff(kap, 6) + tC13.coeff(kap, 6)))
    else:
        print("  (cubic kernels not yet available: run kernels3)")
