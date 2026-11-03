"""Subhorizon asymptotics of the first-order transfer function T_l (exact derivation).

Limit: k0 -> infinity at fixed physical times, i.e. x_o -> infinity with x* = x_o / R, R > 1 fixed.

Start from the entire-function representation (Prop. Tl, transfer_entire.py). The exact recursion
    J_m(w, t) = -t^(1-m) g_m(t) / (m-1) + i w J_{m-1}(w, t) / (m-1),   g_m = e^{iwt} - sum_{j<m} (iwt)^j / j!,
reduces every J_m to J_1 plus elementary terms, and
    J_1(w, t) = -Ein(-i w t) = -gamma - ln(|w| t) + i (pi/2) sgn(w) - E_1(-i w t)
is exact (E_1 the exponential integral).  The only asymptotic input is the standard expansion
    E_1(z) ~ e^{-z}/z sum_k (-1)^k k! / z^k,   |z| -> infinity,  |arg z| < 3 pi / 2.
Substituting, the result is a finite sum of phases e^{i nu x_o} times Laurent series in x_o (and
logarithms, which must cancel).  The script collects the phases, checks that the non-oscillating part
and all logarithms cancel to the order kept, and prints the leading terms.  Results asserted exactly: every growing, constant and log(x_o) term cancels; the leading term is
    T_l = -[cos x* j_l(y) + sqrt3 sin x* j_l'(y)] + O(x_o^-2),   y = sqrt3 (x_o - x*) = k0 (t_o - t*);
and (order >= 2) the observer first enters at O(x_o^-2), at the frequency of its own acoustic
oscillation, e^{+-i x_o}: exactly 2 cos(x_o)/x_o^2 for l = 2 (w_2(0) = 2/15: the local fluid shear),
and with transcendental constants (pi, log(2+sqrt3)) for l > 2, which are printed.  A high-precision evaluation of Prop. Tl at large x_o is printed as a spot
check of the expansion (not as evidence).

    python transfer_subhorizon.py [l] [order]
"""

import sys

import mpmath as mp
import sympy as sp
import transfer_entire as te

L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
ORDER = int(sys.argv[2]) if len(sys.argv) > 2 else 3  # keep terms through x_o^(-ORDER)
xo, R = sp.symbols("x_o R", positive=True)
Ssym, Dsym = sp.symbols("S D", positive=True)  # x* = S x_o and x_o - x* = D x_o, with S = 1/R, D = 1 - 1/R
I, s3 = sp.I, sp.sqrt(3)
EG = sp.EulerGamma
KMAX = ORDER + 2  # terms of the E_1 series: the phase coefficients are O(1) at most
# Phase generators e^{i x_o}, e^{i sqrt3 x_o}, e^{i x*}, e^{i sqrt3 x*} and their inverses as independent
# symbols (so no generator is ever folded into a denominator); log(x_o) as lg.
GEN = sp.symbols("P1 P2 P3 P4")
INV = sp.symbols("Q1 Q2 Q3 Q4")
lg = sp.Symbol("lg")
TS = {"o": (xo, sp.S.Zero), "s": (Ssym * xo, sp.log(Ssym)), "d": (Dsym * xo, sp.log(Dsym))}


def gpow(i, n):
    return GEN[i] ** n if n >= 0 else INV[i] ** (-n)


def phase(w, an):
    """e^{i w t} for w = alpha + beta sqrt3 (integers) and t of kind an."""
    w = sp.expand(w)
    beta = int(w.coeff(s3))
    alpha = int(sp.expand(w - beta * s3))
    if an == "o":
        return gpow(0, alpha) * gpow(1, beta)
    if an == "s":
        return gpow(2, alpha) * gpow(3, beta)
    return gpow(0, alpha) * gpow(1, beta) * gpow(2, -alpha) * gpow(3, -beta)


def J1(w, an):
    t, logfac = TS[an]
    zc = -I * w  # z = zc * t; rationalise the numerical factor so no surd sits in a denominator
    ser = sum((-1) ** k * sp.factorial(k) * sp.radsimp(sp.expand(zc ** (-(k + 1)))) * t ** (-(k + 1)) for k in range(KMAX))
    return -EG - sp.log(abs(w)) - (lg + logfac) + I * sp.pi / 2 * sp.sign(w) - phase(w, an) * ser


def J(m, w, an):
    t = TS[an][0]
    if m == 1:
        return J1(w, an)
    g = phase(w, an) - sum((I * w * t) ** j / sp.factorial(j) for j in range(m))
    return -(t ** (1 - m)) * g / (m - 1) + I * w * J(m - 1, w, an) / (m - 1)


def representation():
    _, coeffs = te.coefficients(L)
    expr = 0
    for sg in (+1, -1):
        a1, a2, b = coeffs[sg]
        a1, a2 = a1.subs(te.xo, xo), a2.subs(te.xo, xo)
        w = 1 - sg * s3
        expr += gpow(1, sg) * (a1 * (J(1, w, "o") - J(1, w, "s")) + a2 * (J(2, w, "o") - J(2, w, "s")))
        expr += GEN[0] * sum(b[m].subs(te.xo, xo) * J(m, -w, "d") for m in b)
    return expr


def split(T):
    """{(phase exponents (4), log power, x_o power): coefficient}."""
    out = {}
    for term in sp.Add.make_args(sp.expand(T)):
        pd = term.as_powers_dict()
        ph = tuple(int(pd.get(GEN[i], 0)) - int(pd.get(INV[i], 0)) for i in range(4))
        key = (ph, int(pd.get(lg, 0)), int(pd.get(xo, 0)))
        rest = term / sp.Mul(*[GEN[i] ** pd.get(GEN[i], 0) * INV[i] ** pd.get(INV[i], 0) for i in range(4)])
        rest = rest / (lg ** pd.get(lg, 0) * xo ** pd.get(xo, 0))
        assert not rest.has(xo, lg, *GEN, *INV), rest
        out[key] = out.get(key, 0) + rest
    return out


def main():
    expr = sp.expand(representation())
    sw = {I: -I, **{GEN[i]: INV[i] for i in range(4)}, **{INV[i]: GEN[i] for i in range(4)}}
    conj = expr.xreplace(sw)
    T = sp.expand(-sp.Rational(3, 2) * (expr + conj) / 2)
    parts = split(T)
    phys = {Ssym: 1 / R, Dsym: 1 - 1 / R}
    res, growing, logs = {}, [], []
    for (ph, lp, xp), c in parts.items():
        c = sp.expand(sp.expand_log(sp.expand(c).subs(phys), force=True))
        if c == 0:
            continue
        if xp > 0 or (xp == 0 and ph == (0, 0, 0, 0)):
            growing.append(((ph, lp, xp), c))
            continue
        if lp:
            logs.append(((ph, lp, xp), c))
        if -xp <= ORDER:
            res[(ph, lp, xp)] = c
    assert not growing, f"non-decaying terms survive: {growing[:4]}"
    print(f"[ok] all growing and constant terms cancel exactly; log(x_o) terms present: {len(logs)}")
    out = {}
    for (ph, lp, xp), c in res.items():
        nu = sp.nsimplify(ph[0] + ph[1] * s3 + (ph[2] + ph[3] * s3) / R)
        out.setdefault(nu, 0)
        out[nu] += c * xo**xp * sp.log(xo) ** lp
    out = {k: sp.collect(sp.expand(v), xo) for k, v in out.items() if sp.expand(v) != 0}
    for nu in sorted(out, key=lambda q: float(sp.N(q.subs(R, 3)))):
        print(f"  e^(i x_o ({nu})):  {out[nu]}")
    return out


def assemble(res, xo_v, R_v):
    tot = 0
    for nu, c in res.items():
        tot += complex(sp.N((sp.exp(I * nu * xo) * c).subs({xo: xo_v, R: R_v}), 30))
    return tot


def check_acoustic_form(res):
    """Leading order equals -[cos x* j_l(y) + sqrt3 sin x* j_l'(y)], y = sqrt3 (x_o - x*), with
    j_l(y) ~ sin(y - l pi/2)/y and j_l'(y) ~ cos(y - l pi/2)/y."""
    xs = xo / R
    y = s3 * (xo - xs)
    cand = -(sp.cos(xs) * sp.sin(y - L * sp.pi / 2) + s3 * sp.sin(xs) * sp.cos(y - L * sp.pi / 2)) / y
    lead = sum(sp.exp(I * nu * xo) * c.coeff(xo, -1) / xo for nu, c in res.items())
    diff = sp.simplify(sp.expand((lead - cand).rewrite(sp.exp)))
    assert diff == 0, diff
    print(f"[ok] l={L}: leading subhorizon term = -[cos x* j_l(y) + sqrt3 sin x* j_l'(y)] with y = k0 (t_o - t*)")
    if ORDER >= 2:
        obs = {nu: sp.expand(c).coeff(xo, -2) * xo**-2 for nu, c in res.items() if nu in (1, -1)}
        assert set(obs) == {1, -1}, obs
        if L == 2:
            assert {k: sp.expand(v) for k, v in obs.items()} == {1: xo**-2, -1: xo**-2}, obs
            print("[ok] l=2: observer term at O(x_o^-2) is exactly 2 cos(x_o)/x_o^2")
        else:
            c = sp.expand(sp.exp(I * xo) * obs[1] + sp.exp(-I * xo) * obs[-1])
            c = sp.simplify(sp.expand(c.rewrite(sp.cos)))
            print(f"     l={L}: observer term at O(x_o^-2): {c}")


def exact_T(xo_v, xs_v):
    """Prop. Tl evaluated at high precision (the elementary parts of J_m cancel strongly at large x)."""
    _, coeffs = te.coefficients(L)
    with mp.workdps(40 + int(xo_v / 4)):
        X, Xs = mp.mpf(xo_v), mp.mpf(str(sp.N(xs_v, mp.mp.dps)))
        num = lambda e: mp.mpc(*[mp.mpf(str(sp.N(part, mp.mp.dps))) for part in (sp.re(e), sp.im(e))])

        def Jn(m, w, tt):
            if m == 1:
                z = abs(w) * tt
                return -(mp.euler + mp.log(z) - mp.ci(z)) + 1j * mp.sign(w) * mp.si(z)
            EM = mp.expj(w * tt) - sum((1j * w * tt) ** j / mp.factorial(j) for j in range(m))
            return tt ** (1 - m) * EM / (1 - m) + 1j * w * Jn(m - 1, w, tt) / (m - 1)

        val = mp.mpc(0)
        for sg in (+1, -1):
            a1, a2, b = coeffs[sg]
            w = mp.mpf(1) - sg * mp.sqrt(3)
            f = lambda e: num(e.subs(te.xo, sp.Rational(str(xo_v))))
            val += mp.expj(sg * mp.sqrt(3) * X) * (f(a1) * (Jn(1, w, X) - Jn(1, w, Xs)) + f(a2) * (Jn(2, w, X) - Jn(2, w, Xs)))
            val += mp.expj(X) * sum(f(b[m]) * Jn(m, -w, X - Xs) for m in b)
        return -1.5 * val.real


if __name__ == "__main__":
    res = main()
    check_acoustic_form(res)
    _, coeffs = te.coefficients(L)
    for xo_v, R_v in ((40, 3), (120, 3), (360, 3), (200, 1.5)):
        exact = exact_T(xo_v, sp.Rational(str(xo_v)) / sp.Rational(str(R_v)))
        approx = assemble(res, xo_v, sp.Rational(str(R_v)))
        print(
            f"  x_o={xo_v:4d} R={R_v}: exact {mp.nstr(exact, 12)}   asymptotic {approx.real:.12g}  (imag {abs(approx.imag):.1e})"
            f"   rel. diff {abs(approx.real - float(exact)) / abs(float(exact)):.1e}"
        )
