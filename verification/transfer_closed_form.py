"""Closed form of the first-order transfer function T_l(k) and of C_l^(11).

    T_l(k; t*, to) = -3 zeta-normalised  int_{x*}^{xo} dx s'(x) w_l(sqrt3 (xo - x)),
    s(x) = 2 sin x / x - cos x,   w_l(y) = j_l(y)/3 + j_l''(y),   x = k t / sqrt3.

Method (exact): write the integrand in u = xo - x as a finite sum
  sum_j c_j(u) exp(i beta_j u),  c_j rational in u,
partial-fraction each c_j, integrate term by term with the recursion
  int e^{i b v} v^{-m} dv = -e^{i b v} v^{1-m}/(m-1) + (i b/(m-1)) int e^{i b v} v^{1-m} dv,
  int e^{i b v} v^{-1} dv = Ci(|b| v) + i sgn(b) Si(|b| v),
CHECK the antiderivative by exact differentiation, and take the u -> 0 (observer) endpoint
as an exact limit (log and pole terms must cancel; asserted).

Monochromatic spectrum  Delta_zeta^2(k) = A_zeta delta_D(ln(k/k0))  (A_zeta = <zeta^2> held fixed):
  C_l^(11) = 4 pi A_zeta T_l(k0)^2.
"""
import sys
import sympy as sp

u, v = sp.symbols('u v', positive=True)
xo, xs, U = sp.symbols('x_o x_star U', positive=True)
y = sp.Symbol('y', positive=True)
s3 = sp.sqrt(3)


def antideriv(beta, m, var):
    """Antiderivative of exp(I beta var) * var**(-m) (beta real nonzero, var>0)."""
    e = sp.exp(sp.I * beta * var)
    if m <= 0:
        return sp.integrate(e * var**(-m), var)
    if m == 1:
        b = sp.Abs(beta)
        sg = sp.sign(beta)
        return sp.Ci(b * var) + sp.I * sg * sp.Si(b * var)
    return -e * var**(1 - m) / (m - 1) + sp.I * beta / (m - 1) * antideriv(beta, m - 1, var)


def split_exp(expr, var):
    """Trig expression in var (frequencies 1 and sqrt3) -> {beta: c_beta(var)} with
    expr = sum_beta c_beta(var) exp(I beta var), c_beta rational in var."""
    P, Q = sp.symbols('P Q')  # P = e^{i var}, Q = e^{i sqrt3 var}
    e = sp.expand(sp.expand_trig(expr))
    rep = {}
    for f in e.atoms(sp.sin, sp.cos):
        a = sp.expand(f.args[0])
        if not a.has(var):
            continue
        c = sp.expand(a.coeff(var))
        assert sp.expand(a - c * var) == 0, f
        bq = c.coeff(s3)
        ap = sp.expand(c - bq * s3)
        assert ap.is_integer and bq.is_integer, f
        Z = P**int(ap) * Q**int(bq)
        rep[f] = (Z + 1 / Z) / 2 if isinstance(f, sp.cos) else (Z - 1 / Z) / (2 * sp.I)
    e = sp.together(e.xreplace(rep))
    num, den = sp.fraction(e)
    dpq = sp.Poly(den, P, Q)
    assert len(dpq.terms()) == 1, "denominator must be monomial in P, Q"
    (mon, dc), = dpq.terms()
    out = {}
    for (i, j), cnum in sp.Poly(sp.expand(num), P, Q).terms():
        beta = sp.nsimplify((i - mon[0]) + (j - mon[1]) * s3)
        out[beta] = out.get(beta, 0) + cnum / dc
    return out


def transfer(l):
    j = sp.expand_func(sp.jn(l, y))
    w = sp.simplify(j / 3 + sp.diff(j, y, 2))
    sp_ = lambda X: 2 * sp.cos(X) / X - 2 * sp.sin(X) / X**2 + sp.sin(X)   # s'(x)
    integrand_u = -3 * sp_(xo - u) * w.subs(y, s3 * u)                      # dx = -du, limits flip
    parts = split_exp(integrand_u, u)
    G_sing = 0   # pieces singular at u=0 (powers of u), integrated in u
    G_far = 0    # pieces with (xo-u) poles, integrated in w = xo - u
    for beta, cf in parts.items():
        pf = sp.apart(sp.together(cf), u)
        for term in sp.Add.make_args(sp.expand(pf)):
            # classify: c * u^-m  or  c * (u - xo)^-p  or polynomial in u
            num, den = sp.fraction(sp.together(term))
            if den.has(xo) and sp.Poly(den, u).degree() > 0 and sp.simplify(den.subs(u, xo)) == 0:
                p = sp.Poly(den, u).degree()
                cst = sp.simplify(num / sp.LC(sp.Poly(den, u)))
                # term = cst/(u - xo)^p = cst (-1)^p / w^p ; e^{i beta u} = e^{i beta xo} e^{-i beta w}; du = -dw
                G_far += -cst * (-1)**p * sp.exp(sp.I * beta * xo) * antideriv(-beta, p, v).subs(v, xo - u) \
                    if beta != 0 else -cst * (-1)**p * sp.integrate((xo - u)**(-p) * 0 + v**(-p), v).subs(v, xo - u)
            else:
                m = -sp.degree(term, u) if not term.has(u) or sp.Poly(sp.fraction(term)[1], u).degree() > 0 else 0
                # term = c * u^n, n integer (possibly negative)
                powd = sp.Poly(sp.fraction(sp.together(term))[1], u).degree()
                pown = sp.Poly(sp.fraction(sp.together(term))[0], u).degree()
                n = pown - powd
                cst = sp.simplify(term / u**n)
                assert not cst.has(u), term
                if beta == 0:
                    G_sing += cst * (u**(n + 1) / (n + 1) if n != -1 else sp.log(u))
                else:
                    G_sing += cst * antideriv(beta, -n, u)
    G = G_sing + G_far
    print("  [stage] antiderivative built", flush=True)
    # exact derivative check of the antiderivative
    resid = split_exp(sp.diff(G, u).rewrite(sp.cos) - integrand_u, u)
    for beta, cf in resid.items():
        assert sp.cancel(sp.together(cf)) == 0, f"antiderivative check failed l={l}, beta={beta}: {cf}"
    print("  [stage] derivative check passed", flush=True)
    # u -> 0 limit: Ci(c u) = EulerGamma + log c + log u + O(u^2), Si(cu) = O(u)
    lu = sp.Symbol('lu')
    Gs = G_sing.replace(sp.Ci, lambda a: sp.EulerGamma + sp.log(a / u) + lu + (sp.Ci(a) - sp.EulerGamma - sp.log(a)))
    Gs = Gs.replace(sp.Si, lambda a: sp.Si(a))
    ser = sp.series(Gs, u, 0, 1).removeO()
    ser = sp.expand(ser)
    assert sp.simplify(ser.coeff(lu)) == 0, "log(u) does not cancel"
    lim0 = sp.limit(sp.expand(ser.subs(lu, 0)), u, 0)
    assert lim0.is_finite is not False
    lim0 = lim0 + G_far.subs(u, 0)
    T = sp.expand(sp.expand_complex((G.subs(u, U) - lim0).rewrite(sp.cos)))
    T = T.subs(U, xo - xs)
    im = sp.im(T)
    assert sp.cancel(sp.together(sp.expand(im))) == 0, f"T_l not real: {im}"
    return sp.collect(sp.expand(sp.re(T)), list(T.atoms(sp.Si, sp.Ci)))


if __name__ == "__main__":
    ls = [int(a) for a in sys.argv[1:]] or [2]
    import mpmath as mp
    for l in ls:
        T = transfer(l)
        print(f"[ok] l={l}: antiderivative verified by exact differentiation; observer-endpoint log/poles cancel")
        print(f"T_{l} =", T)
        # supplementary (NOT part of the derivation): high-precision spot check vs quadrature
        fT = sp.lambdify((xs, xo), T, 'mpmath')
        def integrand(xx, xo_):
            yy = mp.sqrt(3) * (xo_ - xx)
            j = lambda q: mp.sqrt(mp.pi / (2 * q)) * mp.besselj(l + 0.5, q)
            w = j(yy) / 3 + mp.diff(j, yy, 2)
            sp_ = 2 * mp.cos(xx) / xx - 2 * mp.sin(xx) / xx**2 + mp.sin(xx)
            return -3 * sp_ * w
        mp.mp.dps = 30
        for (a, b) in [(mp.mpf('0.3'), mp.mpf('4.1')), (mp.mpf('2'), mp.mpf('17.5'))]:
            q = mp.quad(lambda xx: integrand(xx, b), [a, b - mp.mpf('1e-25')])
            val = fT(a, b)
            print(f"    supplementary spot check x*={a}, xo={b}: closed form {mp.nstr(mp.re(val), 18)}  quad {mp.nstr(q, 18)}")
