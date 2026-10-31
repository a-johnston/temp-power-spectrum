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
        a0 = sp.expand(a - c * var)            # var-independent phase
        assert not a0.has(var), f
        bq = c.coeff(s3)
        ap = sp.expand(c - bq * s3)
        assert ap.is_integer and bq.is_integer, f
        Z = P**int(ap) * Q**int(bq)
        Cz, Sz = (Z + 1 / Z) / 2, (Z - 1 / Z) / (2 * sp.I)
        if isinstance(f, sp.cos):
            rep[f] = sp.cos(a0) * Cz - sp.sin(a0) * Sz
        else:
            rep[f] = sp.sin(a0) * Cz + sp.cos(a0) * Sz
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


def laurentize(expr, syms):
    """Map every sin/cos whose argument is sum_s (a_s + b_s sqrt3) s, a_s, b_s integers, to
    Laurent monomials in P_s = e^{i s}, Q_s = e^{i sqrt3 s}.  An expression built from such
    trig functions and rational functions of syms vanishes identically iff its image cancels
    to zero as a rational function (the map is a ring homomorphism into Q(sqrt3, i)(syms, P, Q))."""
    gens = {s_: (sp.Symbol(f'P_{s_}'), sp.Symbol(f'Q_{s_}')) for s_ in syms}
    rep = {}
    for f in expr.atoms(sp.sin, sp.cos):
        a = sp.expand(f.args[0])
        Z = 1
        rest = a
        for s_ in syms:
            c = sp.expand(a.coeff(s_))
            rest = sp.expand(rest - c * s_)
            bq = c.coeff(s3)
            ap = sp.expand(c - bq * s3)
            assert ap.is_integer and bq.is_integer, f
            Pz, Qz = gens[s_]
            Z *= Pz**int(ap) * Qz**int(bq)
        assert (rest / sp.pi).is_rational, f      # only exact phases such as -pi/2 allowed
        Z *= sp.exp(sp.I * rest)                  # evaluates to an exact algebraic number
        rep[f] = (Z + 1 / Z) / 2 if isinstance(f, sp.cos) else (Z - 1 / Z) / (2 * sp.I)
    return expr.xreplace(rep)


def partial_fractions(cf, var, pole):
    """Exact partial fractions of a rational function of var whose only finite poles are
    var = 0 and var = pole.  Returns list of (kind, power, coeff) with kind in
    {'u': coeff*var**power (power may be negative), 'far': coeff*(var-pole)**(-power)}.
    The decomposition is re-assembled and checked against cf exactly."""
    cf = sp.cancel(sp.together(cf))
    num, den = sp.fraction(cf)
    dpoly = sp.Poly(den, var)
    m = 0
    while dpoly.eval(0) == 0:
        dpoly = sp.Poly(sp.cancel(dpoly.as_expr() / var), var)
        m += 1
    p = 0
    while sp.simplify(dpoly.as_expr().subs(var, pole)) == 0:
        dpoly = sp.Poly(sp.cancel(dpoly.as_expr() / (var - pole)), var)
        p += 1
    assert dpoly.degree() == 0, f"unexpected poles: {dpoly}"
    out = []
    recon = 0
    for (order, centre, kind) in ((m, 0, 'u'), (p, pole, 'far')):
        if order == 0:
            continue
        g = sp.cancel(cf * (var - centre)**order)
        d = g
        for n in range(order):
            c = sp.cancel(d.subs(var, centre) / sp.factorial(n))
            if c != 0:
                if kind == 'u':
                    out.append(('u', n - order, c))
                else:
                    out.append(('far', order - n, c))
                recon += c * (var - centre)**(n - order)
            d = sp.diff(d, var)
    rest = sp.cancel(cf - recon)
    assert sp.fraction(rest)[1].has(var) is False, f"remainder not polynomial: {rest}"
    for (deg,), c in sp.Poly(rest, var).terms():
        out.append(('u', deg, c))
        recon += c * var**deg
    assert sp.cancel(recon - cf) == 0
    return out


def transfer(l):
    j = sp.expand_func(sp.jn(l, y))
    w = sp.simplify(j / 3 + sp.diff(j, y, 2))
    sp_ = lambda X: 2 * sp.cos(X) / X - 2 * sp.sin(X) / X**2 + sp.sin(X)   # s'(x)
    integrand_u = -3 * sp_(xo - u) * w.subs(y, s3 * u)                      # dx = -du, limits flip
    parts = split_exp(integrand_u, u)
    Sx, Cx = sp.symbols('S_xo C_xo')     # sin(xo), cos(xo) treated as constants during algebra
    back = {Sx: sp.sin(xo), Cx: sp.cos(xo)}
    G_sing = 0   # pieces with powers of u (singular at the observer u=0), integrated in u
    G_far = 0    # pieces with (u - xo) poles (emission side), integrated in v = xo - u
    for beta, cf in parts.items():
        cf = cf.xreplace({sp.sin(xo): Sx, sp.cos(xo): Cx})
        for kind, pw, cst in partial_fractions(cf, u, xo):
            cst = cst.xreplace(back)
            if kind == 'u':
                n = pw                                   # cst * u**n
                if beta == 0:
                    G_sing += cst * (u**(n + 1) / (n + 1) if n != -1 else sp.log(u))
                else:
                    G_sing += cst * antideriv(beta, -n, u)
            else:
                # cst/(u - xo)^pw = cst (-1)^pw / v^pw, v = xo - u, du = -dv, e^{i b u} = e^{i b xo} e^{-i b v}
                if beta == 0:
                    prim = (v**(1 - pw) / (1 - pw)) if pw != 1 else sp.log(v)
                    G_far += -cst * (-1)**pw * prim.subs(v, xo - u)
                else:
                    G_far += -cst * (-1)**pw * sp.exp(sp.I * beta * xo) * antideriv(-beta, pw, v).subs(v, xo - u)
    G = G_sing + G_far
    print("  [stage] antiderivative built", flush=True)
    # exact derivative check of the antiderivative
    resid = laurentize(sp.expand(sp.diff(G, u).rewrite(sp.cos) - integrand_u), [u, xo])
    assert sp.cancel(sp.together(resid)) == 0, f"antiderivative check failed for l={l}"
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
    # Evaluate at u = U and parametrise xo = x* + U with x*, U > 0 so every Si/Ci argument is
    # manifestly positive (Si, Ci are then real); prove realness exactly via conjugation.
    T = (G.subs(u, U) - lim0).rewrite(sp.cos).subs(xo, xs + U)
    T = T.replace(sp.Si, lambda a: sp.Si(sp.expand(a))).replace(sp.Ci, lambda a: sp.Ci(sp.expand(a)))
    specials = sorted(T.atoms(sp.Si, sp.Ci), key=str)
    for f in specials:
        assert sp.factor(f.args[0]).is_positive, f"argument not manifestly positive: {f}"
    reals = {f: sp.Symbol(f'F{i}', real=True) for i, f in enumerate(specials)}
    TL = laurentize(sp.expand(T.xreplace(reals)), [xs, U])
    gens = [g_ for g_ in TL.free_symbols if g_.name.startswith(('P_', 'Q_'))]
    conj = TL.xreplace({g_: 1 / g_ for g_ in gens}).subs(sp.I, -sp.I)
    assert sp.cancel(sp.together(TL - conj)) == 0, "T_l is not real"
    Tre = sp.expand(sp.re(sp.expand_complex(sp.expand(T.xreplace(reals)))))
    Tre = Tre.xreplace({v_: k_ for k_, v_ in reals.items()})
    return Tre, specials


def present(Tre, specials):
    """Group by Si/Ci, return to (x*, xo), simplify each coefficient separately."""
    d = sp.collect(sp.expand(Tre), specials, evaluate=False)
    out = 0
    pieces = []
    for key in sorted(d, key=str):
        cf = sp.simplify(sp.trigsimp(sp.expand_trig(sp.expand(d[key].subs(U, xo - xs)))))
        arg_key = key.xreplace({U: xo - xs}) if key != 1 else 1
        if arg_key != 1:
            arg_key = arg_key.func(sp.factor(arg_key.args[0]))
        pieces.append((arg_key, cf))
        out += cf * arg_key
    return out, pieces


if __name__ == "__main__":
    ls = [int(a) for a in sys.argv[1:]] or [2]
    import mpmath as mp
    for l in ls:
        Tre, specials = transfer(l)
        T, pieces = present(Tre, specials)
        import os
        os.makedirs(os.path.join(os.path.dirname(__file__) or '.', 'results'), exist_ok=True)
        with open(os.path.join(os.path.dirname(__file__) or '.', 'results', f'T{l}_closed_form.txt'), 'w') as fh:
            fh.write(f"# T_{l}(x*, xo), x = k t / sqrt(3); exact, generated by transfer_closed_form.py\n")
            for key, cf in pieces:
                fh.write(f"[{key}] * (\n    {cf}\n)\n")
            fh.write("\n# LaTeX\n" + sp.latex(T) + "\n")
        print(f"    written results/T{l}_closed_form.txt ({len(pieces)} Si/Ci groups)")
        print(f"[ok] l={l}: antiderivative verified by exact differentiation; observer-endpoint log/poles cancel; result exactly real")
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
