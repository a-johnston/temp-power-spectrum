"""Second-order field kernels in closed form (entire functions Si, Cin, Ephi), general momenta,
and an exact closed-form check of EVERY component of (N), (M), (T), (Ham).

Frame: k1 = (p,0,0), k2 = (b1,b2,0), with b1 eliminated in favour of K = |k1+k2|:
    b1 = (K^2 - p^2 - q^2)/(2p),   b2^2 = q^2 - b1^2.
Kernels (z1 z2 coefficients, no factorials):
  nu_2  = (9 c2/K^2) u1(x) - (1/(3c)) [u2(x) P1(t) - u1(x) P2(t)],  x = K t/sqrt3, c = K/sqrt3,
          P_a = int_0^t u_a(c s) c_N(s) ds,   u1 = sin x/x - cos x,  u2 = -cos x/x - sin x
  h_S   = (3/(2K^2)) [4K^2 nu + 12 nu'/t + 12 nu/t^2 + c_Ham]           (Hamiltonian constraint;
          algebraic -- no integral)
  w     = (2i/K^2) int_0^t P_perp c_M ds,   H^V = K w + w K
  h     = -(2/(K t)) [sin(Kt) int_0^t s cos(Ks) TT c_T ds - cos(Kt) int_0^t s sin(Ks) TT c_T ds]
  H_2   = (K K/K^2 - delta/3) h_S + H^V + h
"""

import sys
import time

import second_order as so
import sympy as sp
from closed_forms import Cin, Freq, check_antiderivative, integrate0, split_frequencies, to_elementary
from engine import I, t

p, q, b1, b2 = so.p, so.q, so.b1, so.b2
K = sp.Symbol("K", positive=True)
s = sp.Symbol("s", positive=True)
s3 = sp.sqrt(3)
B1 = (K**2 - p**2 - q**2) / (2 * p)
SYMS = (p, q, K)


def red(e):
    """Exact reduction: b1 -> (K^2-p^2-q^2)/(2p), b2^2 -> q^2 - b1^2, then cancel."""
    e = sp.together(sp.sympify(e).subs(b1, B1))
    n, d = sp.fraction(e)

    def rb(x):
        x = sp.expand(x)
        if not x.has(b2):
            return x
        P = sp.Poly(x, b2)
        return sp.expand(sum(c * (q**2 - B1**2) ** (m // 2) * b2 ** (m % 2) for (m,), c in P.terms()))

    return sp.cancel(rb(n) / rb(d))


_T0 = time.time()


def LOG(msg):
    print(f"  [{time.time() - _T0:7.1f}s] {msg}", flush=True)


def build():
    src = so.sources()
    unk = src["unk"]
    fnl = [unk["n"]] + [e for e in set(unk["A"]) if isinstance(e, sp.Function)]
    zero = {f: 0 for f in fnl}
    full = src["Ham_full"]
    cH = full.subs({d: 0 for d in full.atoms(sp.Derivative)}).subs(zero)
    linH = sp.expand(full - cH)
    sub = {b1: B1}
    cN = sp.expand(src["N"].subs(sub))
    cM = [sp.expand(c.subs(sub)) for c in src["M"]]
    cT = src["T"].applyfunc(lambda e: sp.expand(e.subs(sub)))
    cH = sp.expand(cH.subs(sub))
    return src, unk, linH, cN, cM, cT, cH


def solve():
    src, unk, linH, cN, cM, cT, cH = build()
    LOG("sources ready")
    c = K / s3
    x = c * t
    u1 = lambda X: sp.sin(X) / X - sp.cos(X)
    u2 = lambda X: -sp.cos(X) / X - sp.sin(X)
    c2 = red((2 * (p**2 + q**2) - p * b1) / 9)
    cNs = cN.subs(t, s)
    B_1 = sp.expand(u1(c * s) * cNs)
    B_2 = sp.expand(u2(c * s) * cNs)
    P1 = integrate0(B_1, s, t, SYMS, red)
    LOG("P1 integrated")
    P2 = integrate0(B_2, s, t, SYMS, red)
    LOG("P2 integrated")
    check_antiderivative(P1, B_1, s, t, SYMS, red)
    check_antiderivative(P2, B_2, s, t, SYMS, red)
    LOG("[ok] P1, P2 verified by exact differentiation")
    nu = 9 * c2 / K**2 * u1(x) - (u2(x) * P1 - u1(x) * P2) / (3 * c)
    # vector
    KV = [p + B1, b2, 0]
    P = sp.Matrix(3, 3, lambda i, j: (1 if i == j else 0) - KV[i] * KV[j] / K**2)
    cMv = sp.Matrix(cM)
    Pc = (P * cMv).applyfunc(lambda e: sp.expand(e))
    w = []
    for i in range(3):
        if red(Pc[i]) == 0:
            w.append(0)
            continue
        Bi = sp.expand(Pc[i].subs(t, s))
        Wi = integrate0(Bi, s, t, SYMS, red)
        check_antiderivative(Wi, Bi, s, t, SYMS, red)
        w.append(2 * I / K**2 * Wi)
    LOG("[ok] vector kernel integrated and verified")
    # tensor
    PMP = P * cT * P
    TTc = PMP - P * PMP.trace() / 2
    h = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            e = sp.expand(TTc[i, j])
            if red(e) == 0:
                continue
            es = e.subs(t, s)
            Bc = sp.expand(s * sp.cos(K * s) * es)
            Bs = sp.expand(s * sp.sin(K * s) * es)
            Pc_ = integrate0(Bc, s, t, SYMS, red)
            Ps_ = integrate0(Bs, s, t, SYMS, red)
            check_antiderivative(Pc_, Bc, s, t, SYMS, red)
            check_antiderivative(Ps_, Bs, s, t, SYMS, red)
            h[i, j] = h[j, i] = -2 / (K * t) * (sp.sin(K * t) * Pc_ - sp.cos(K * t) * Ps_)
    LOG("[ok] tensor kernel integrated and verified")
    # scalar from the Hamiltonian constraint: linH = Lnu-part + Lh-part; solve for h_S
    hS_sym = sp.Symbol("hS")
    Kh = [e / K for e in KV]
    Aprobe = sp.Matrix(3, 3, lambda i, j: (Kh[i] * Kh[j] - sp.Rational(1 if i == j else 0, 3)) * hS_sym)
    fns = {unk["n"]: sp.Function("NU")(t)}
    for i in range(3):
        for j in range(i, 3):
            if isinstance(unk["A"][i, j], sp.Function):
                fns[unk["A"][i, j]] = Aprobe[i, j]
    lin = sp.expand(linH.subs(fns).doit())
    lin = red(lin)
    coef_h = red(sp.diff(lin, hS_sym))
    rest = sp.expand(lin - coef_h * hS_sym)
    assert red(coef_h + sp.Rational(2, 3) * K**2) == 0, ("Ham coefficient of h_S", coef_h)
    NU = sp.Function("NU")(t)
    exp_rest = 4 * K**2 * NU + 12 * NU.diff(t) / t + 12 * NU / t**2
    assert red(rest - exp_rest) == 0, ("Ham lapse part", rest)
    LOG("[ok] Hamiltonian linear part = -(2/3)K^2 h_S + 4K^2 nu + 12 nu'/t + 12 nu/t^2")
    hS = 3 / (2 * K**2) * (4 * K**2 * nu + 12 * sp.diff(nu, t) / t + 12 * nu / t**2 + cH)
    A = sp.Matrix(
        3, 3, lambda i, j: (Kh[i] * Kh[j] - sp.Rational(1 if i == j else 0, 3)) * hS + KV[i] * w[j] + w[i] * KV[j] + h[i, j]
    )
    return dict(nu=nu, A=A, hS=hS, w=w, h=h, P1=P1, P2=P2, unk=unk, src=src, linH=linH, cH=cH, cN=cN, cM=cM, cT=cT)


def residual_zero(e, label):
    """Exact zero test for expressions in elementary functions, Si, Cin, Ephi of frequency args."""
    e = to_elementary(sp.expand(e))

    # canonical arguments: Si/Cin of algebraically equal arguments must be the same atom
    # (Si odd and Cin even are normalised by their eval rules on the expanded argument)
    def canon(f):
        a = sp.expand(f.args[0])
        for v in (K, p, q):
            c = a.coeff(v)
            if c != 0:
                if (c / sp.Abs(c)).is_negative or sp.sign(c) == -1:
                    # Si odd, Cin even: fix the sign of the argument deterministically
                    return -sp.Si(-a) if isinstance(f, sp.Si) else Cin(-a)
                break
        return f.func(a)

    e = e.replace(lambda f: isinstance(f, (sp.Si, Cin)), canon)
    # Si/Cin of distinct canonical arguments are treated as independent transcendentals;
    # they are replaced by symbols BEFORE rewriting trig functions as exponentials
    specials = sorted(e.atoms(sp.Si, Cin), key=str)
    syms = [sp.Symbol(f"_F{i}") for i in range(len(specials))]
    e = e.xreplace(dict(zip(specials, syms, strict=True)))
    e = sp.expand(sp.expand(e).rewrite(sp.exp))
    groups = sp.collect(e, syms, evaluate=False) if syms else {1: e}
    fr = Freq(SYMS, t)
    for key, cf in groups.items():
        parts = split_frequencies(cf, fr, t)
        for fk, lau in parts.items():
            for m, cc in lau.items():
                assert red(cc) == 0, (label, key, fk, m)
    print(f"[ok] {label}")


if __name__ == "__main__":
    import os
    import pickle

    cache = os.path.join(os.path.dirname(__file__) or ".", "results", ".closed_second_order.pkl")
    if os.path.exists(cache) and "--fresh" not in sys.argv:
        with open(cache, "rb") as fh:
            sol = pickle.load(fh)
        print("  (loaded cached closed forms; run with --fresh to rebuild and re-verify antiderivatives)")
    else:
        sol = solve()
        with open(cache, "wb") as fh:
            pickle.dump(sol, fh)
    with open("results/second_order_closed_forms.txt", "w") as fh:
        fh.write("# Exact second-order kernels in closed form (Si, Cin, Ephi entire). Frame and symbols: see\n")
        fh.write("# closed_second_order.py.  Cin(x)=int_0^x (1-cos u)/u du; Ephi(m,z)=int_0^1 u^m e^{zu} du.\n")
        fh.write(f"nu_2 = {sol['nu']}\n")
        fh.write(f"h_S = {sol['hS']}\n")
        for i in range(3):
            fh.write(f"w[{i}] = {sol['w'][i]}\n")
        for i in range(3):
            for j in range(i, 3):
                fh.write(f"h_TT[{i},{j}] = {sol['h'][i, j]}\n")
    print(
        "[ok] closed forms written to results/second_order_closed_forms.txt "
        "(every antiderivative verified by exact differentiation when built)"
    )
    # Optional, slow (hours): substitute the closed forms back into the equations.  This is
    # redundant with the antiderivative checks + verified Green functions + the exact series
    # tests through O(t^8); use --check N (lapse equation only) or --check all.
    mode = sys.argv[sys.argv.index("--check") + 1] if "--check" in sys.argv else None
    if mode:
        unk = sol["unk"]
        LN, LM, LT = so.expected_linear(unk)
        fns = {unk["n"]: sol["nu"]}
        for i in range(3):
            for j in range(i, 3):
                if isinstance(unk["A"][i, j], sp.Function):
                    fns[unk["A"][i, j]] = sol["A"][i, j]

        def app(e):
            return sp.expand(e.subs({b1: B1}).subs(fns).doit())

        # the Hamiltonian constraint holds by construction of h_S; the rest are independent checks
        residual_zero(app(LN) + sol["cN"], "closed form: (N) lapse equation")
        if mode == "N":
            sys.exit(0)
        for i in range(3):
            residual_zero(app(LM[i]) + sol["cM"][i], f"closed form: (M) momentum constraint, component {i}")
        for i in range(3):
            for j in range(i, 3):
                residual_zero(app(LT[i, j]) + sol["cT"][i, j], f"closed form: (T) trace-free evolution [{i}{j}]")
        print("SECOND-ORDER CLOSED FORMS VERIFIED IN EVERY EQUATION")
