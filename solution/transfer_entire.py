"""Entire-function representation of the first-order transfer function T_l.

    T_l(x*, xo) = -3 int_{x*}^{xo} s'(x) w_l(sqrt3 (xo - x)) dx,
    s(x) = 2 sin x/x - cos x,  w_l(y) = j_l(y)/3 + j_l''(y).

Write s' = Re S(x) and w_l = Re W_l(y) with the complex (Hankel) parts
    S(x)   = e^{ix} Sr(x),      Sr(x) = -i + 2/x + 2i/x^2,
    W_l(y) = e^{iy} Vr(y),      W_l = (1/3 + d_y^2) h_l^(1)(y),   Vr(y) = sum_{n=1}^{l+3} v_n y^{-n}.
Then
    T_l = -(3/2) Re sum_{sigma=+-1} int_{x*}^{xo} S(x) W_l^sigma(sqrt3(xo-x)) dx,
with W^+ = W_l, W^- = conj W_l.  Each rational factor has poles only at x = 0 (order 2) and
x = xo (order l+3) and vanishes at infinity, so it equals the sum of its two principal parts.
Integrating e^{i w x} x^{-m} with the subtracted kernels of the report,
    J_m(w, t) = int_0^t s^{-m} [ e^{i w s} - sum_{j<m} (i w s)^j / j! ] ds    (entire in w t),
leaves elementary subtraction terms whose total is the principal part of the REAL, regular
integrand at x = 0 and at x = xo; it therefore vanishes under Re.  Hence (no logs, no Ci, no
elementary remainder)

    T_l = -(3/2) Re sum_sigma { e^{i sigma sqrt3 xo} sum_{m=1,2} a_{sigma m} [J_m(w_sigma, xo) - J_m(w_sigma, x*)]
                                + e^{i xo} sum_{m=1}^{l+3} b_{sigma m} J_m(-w_sigma, xo - x*) },
    w_sigma = 1 - sigma sqrt3,
    a_{sigma 2} = 2i V_sigma(sqrt3 xo),   a_{sigma 1} = 2 V_sigma(sqrt3 xo) - 2 sqrt3 i V_sigma'(sqrt3 xo),
    b_{sigma m} = [u^{-m}] Sr(xo - u) V_sigma(sqrt3 u)        (Laurent coefficient at u = 0).

EXACT CHECKS (no quadrature as evidence):
  1. d/dx* of the representation equals +3 s'(x*) w_l(sqrt3(xo-x*)) identically (exponential
     polynomial identity, proved by mapping e^{i x*}, e^{i sqrt3 x*}, ... to independent Laurent
     generators), and the representation vanishes at x* = xo.  Together these prove it.
     The construction itself is valid for every l; the identity is checked for each l given.
A high-precision comparison with the earlier Ci/Si closed form is printed as a supplementary
spot check only.
"""

import sys

import sympy as sp

x, u, xo, xs, t = sp.symbols("x u x_o x_star t", positive=True)
y = sp.Symbol("y", positive=True)
s3 = sp.sqrt(3)
I = sp.I


def hankel_rational(n):
    """h_n^(1)(y) = e^{iy} * this rational function of y."""
    return (
        (-I) ** (n + 1)
        / y
        * sum((I / (2 * y)) ** m * sp.factorial(n + m) / (sp.factorial(m) * sp.factorial(n - m)) for m in range(n + 1))
    )


def V_rational(l):
    """W_l = (1/3 + d^2) h_l^(1) = e^{iy} V(y)."""
    Q = hankel_rational(l)
    return sp.expand(sp.Rational(1, 3) * Q - Q + 2 * I * sp.diff(Q, y) + sp.diff(Q, y, 2))


Sr = -I + 2 / x + 2 * I / x**2


def conj_rational(expr):
    return expr.subs(I, -I)  # coefficients are in Q(i); y, x real


def coefficients(l):
    V = {+1: V_rational(l)}
    V[-1] = conj_rational(V[+1])
    out = {}
    for sg in (+1, -1):
        Vs = V[sg]
        a2 = 2 * I * Vs.subs(y, s3 * xo)
        a1 = 2 * Vs.subs(y, s3 * xo) - 2 * s3 * I * sp.diff(Vs, y).subs(y, s3 * xo)
        prod = sp.expand(sp.series(Sr.subs(x, xo - u) * Vs.subs(y, s3 * u), u, 0, 1).removeO())
        b = {m: sp.factor(sp.expand(prod.coeff(u, -m))) for m in range(1, l + 4)}
        out[sg] = (sp.factor(a1), sp.factor(a2), b)
    return V, out


J = sp.Function("J")  # J(m, w, t)


def representation(l, coeffs):
    expr = 0
    for sg in (+1, -1):
        a1, a2, b = coeffs[sg]
        w = 1 - sg * s3
        expr += sp.exp(I * sg * s3 * xo) * (a1 * (J(1, w, xo) - J(1, w, xs)) + a2 * (J(2, w, xo) - J(2, w, xs)))
        expr += sp.exp(I * xo) * sum(b[m] * J(m, -w, xo - xs) for m in b)
    return expr  # T_l = -(3/2) Re expr


def dJdt(m, w, tt):
    return tt ** (-m) * (sp.exp(I * w * tt) - sum((I * w * tt) ** j / sp.factorial(j) for j in range(m)))


def laurent_generators(expr, syms):
    """Map sin/cos/exp with arguments sum (a + b sqrt3) s to independent Laurent generators."""
    gens = {s_: (sp.Symbol(f"P_{s_}"), sp.Symbol(f"Q_{s_}")) for s_ in syms}
    expr = sp.expand(expr.rewrite(sp.exp))
    rep = {}
    for f in expr.atoms(sp.exp):
        a = sp.expand(f.args[0] / I)
        Z = 1
        rest = a
        for s_ in syms:
            c = sp.expand(a.coeff(s_))
            rest = sp.expand(rest - c * s_)
            bq = c.coeff(s3)
            ap = sp.expand(c - bq * s3)
            assert ap.is_integer and bq.is_integer, f
            Pz, Qz = gens[s_]
            Z *= Pz ** int(ap) * Qz ** int(bq)
        assert rest == 0, f
        rep[f] = Z
    return expr.xreplace(rep)


def check(l, coeffs):
    rep = representation(l, coeffs)
    # derivative with respect to x*: J(m,w,xo)-terms are constant; J(m,w,x*) and J(m,-w,xo-x*)
    d = 0
    for sg in (+1, -1):
        a1, a2, b = coeffs[sg]
        w = 1 - sg * s3
        d += sp.exp(I * sg * s3 * xo) * (-a1 * dJdt(1, w, xs) - a2 * dJdt(2, w, xs))
        d += sp.exp(I * xo) * sum(-b[m] * dJdt(m, -w, xo - xs) for m in b)
    dT = -sp.Rational(3, 4) * (d + d.subs(I, -I))  # -(3/2) Re d, all symbols real
    jl = sp.expand_func(sp.jn(l, y))
    wl = jl / 3 + sp.diff(jl, y, 2)
    sprime = sp.diff(2 * sp.sin(x) / x - sp.cos(x), x)
    target = 3 * sprime.subs(x, xs) * wl.subs(y, s3 * (xo - xs))  # d/dx* of -3 int_{x*}
    resid = laurent_generators(dT - target, [xs, xo])
    assert sp.cancel(sp.together(resid)) == 0, f"derivative check failed for l={l}"
    return rep


def numeric(l, coeffs, xs_v, xo_v):
    import mpmath as mp

    def Jn(m, w, tt):
        # exact recursion of the report: J_1 = -Cin(|w| t) + i sgn(w) Si(|w| t),
        # J_M = t^{1-M} E_M / (1-M) + i w J_{M-1} / (M-1)
        if m == 1:
            z = abs(w) * tt
            cin = mp.euler + mp.log(z) - mp.ci(z)
            return -cin + 1j * mp.sign(w) * mp.si(z)
        EM = mp.e ** (1j * w * tt) - sum((1j * w * tt) ** j / mp.factorial(j) for j in range(m))
        return tt ** (1 - m) * EM / (1 - m) + 1j * w * Jn(m - 1, w, tt) / (m - 1)

    val = 0
    for sg in (+1, -1):
        a1, a2, b = coeffs[sg]
        w = float(1 - sg * s3)
        f = lambda e: complex(sp.N(e.subs(xo, xo_v), 30))
        val += mp.e ** (1j * sg * float(s3) * xo_v) * (
            f(a1) * (Jn(1, w, xo_v) - Jn(1, w, xs_v)) + f(a2) * (Jn(2, w, xo_v) - Jn(2, w, xs_v))
        )
        val += mp.e ** (1j * xo_v) * sum(f(b[m]) * Jn(m, -w, xo_v - xs_v) for m in b)
    return -1.5 * mp.re(val)


if __name__ == "__main__":
    ls = [int(a) for a in sys.argv[1:]] or [2, 3]
    for l in ls:
        V, coeffs = coefficients(l)
        check(l, coeffs)
        print(f"[ok] l={l}: entire-J representation proved (exact x*-derivative identity; vanishes at x*=xo)")
        print(f"     V_l(y) = {sp.factor(V[1])}")
        for sg in (+1, -1):
            a1, a2, b = coeffs[sg]
            print(f"     sigma={sg:+d}: a1={a1}\n                a2={a2}")
            for m in b:
                print(f"                b{m}={b[m]}")
        if l == 2:
            import mpmath as mp

            mp.mp.dps = 20
            for a_, b_ in [(0.3, 4.1), (2.0, 17.5)]:
                jl = lambda q, l=l: mp.sqrt(mp.pi / (2 * q)) * mp.besselj(l + 0.5, q)

                def integrand(xx, b_=b_, jl=jl):
                    yy = mp.sqrt(3) * (b_ - xx)
                    return -3 * (2 * mp.cos(xx) / xx - 2 * mp.sin(xx) / xx**2 + mp.sin(xx)) * (jl(yy) / 3 + mp.diff(jl, yy, 2))

                q = mp.quad(integrand, [a_, b_ - mp.mpf("1e-15")])
                print(
                    f"     supplementary spot check x*={a_}, xo={b_}: representation {mp.nstr(numeric(l, coeffs, a_, b_), 12)}  quadrature {mp.nstr(q, 12)}"
                )
