"""Generator of the complete observed-temperature kernels Theta^(1), Theta^(2), Theta^(3).

This is an INDEPENDENT implementation of the measurement: it does not use the hand-expanded
formulas of the derivation.  It solves the exact definitions by fixed-point iteration in the
multilinear algebra (labels = plane waves; products of distinct labels only):
  * ray:      X'' = A(t, Xbar + dX, ebar + dX'),  dX(t_o) = 0,  dX'(t_o) = -(gamma_o^{-1/2} - 1) n,
              A = -2 S e - Gamma(e, e) + e sigma(e, e), fields Taylor-evaluated at the ray;
  * emission: ln(t_e/tbar_*) + nu(t_e, X(t_e)) = 0;
  * exponent: I = int_{t_e}^{t_o} sigma(e, e) dt;   T_obs/T_fluid(o) = e^{-I};
  * normalise: Theta_c = e^{-I}/<e^{-I}>_sky - 1;
  * clock:    int_{t_I}^{t_o} t e^{nu(t,0)} dt = (tbar_o^2 - tbar_I^2)/2,  ln(t_I/tbar_I) + nu(t_I, 0) = 0;
  * compose:  Theta = Theta_c(t_o = tbar_o + D).
Fields come from the triangular solver as exact Taylor series (third_order_consistency.solve_top).
Mode: all momenta carry a bookkeeping factor kappa; every quantity is an exact polynomial in kappa
truncated at kappa^KMAX (the kernels are analytic in the momenta at fixed finite times), and in
the sky direction n (sky means use exact spherical moments).

Checks (exact, order by order in kappa):
  W2  dilation Ward identity of the OBSERVABLE at second order: Theta_{sL} = -kappa d_kappa Theta_s
  W3  at third order: Theta_{sLL} = (kappa d_kappa)^2 Theta_s
  SH  superhorizon closed form at O(kappa^2) for a general pair (k, q):
      Theta = -tau^2 e^{-2 zeta} TF[zeta_ij - zeta_i zeta_j] n n,  tau^2 = (t_o^2 - tbar_*^2)/6.
"""

import sys
import time

import second_order as so
import sympy as sp
from engine import Ctx, Field, I, equations, fexp, t

KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 4
_argv, sys.argv = sys.argv, sys.argv[:1]
import third_order_consistency as tc  # noqa: E402  (exact series depth 6 >= KMAX)

sys.argv = _argv
kap = sp.Symbol("kappa", positive=True)
to, ts, tIb, tob = sp.symbols("t_o tbar_star tbar_I tbar_o", positive=True)
tp = sp.Symbol("tp", positive=True)
nv = sp.symbols("n1 n2 n3", real=True)
T0 = time.time()


def LOG(m):
    print(f"  [{time.time() - T0:7.1f}s] {m}", flush=True)


def trunc(e):
    e = sp.expand(e)
    if not e.has(kap):
        return e
    P = sp.Poly(e, kap)
    return sp.expand(sum(c * kap**m for (m,), c in P.terms() if m <= KMAX))


# ---------------------------------------------------------------- zero-momentum (long) modes
_lin = so.lin_fields


def lin_fields(label, kvec, kmag):
    if kmag == 0:
        return Field(), Field({frozenset([label]): sp.Integer(1)}), [[Field() for _ in range(3)] for _ in range(3)]
    return _lin(label, kvec, kmag)


so.lin_fields = lin_fields  # ty: ignore[invalid-assignment]  (deliberate patch: zero-momentum modes)


def sky_mean(e):
    """Exact mean over the unit sphere of a polynomial in n1, n2, n3."""
    e = sp.expand(e)
    if not e.has(*nv):
        return e
    P = sp.Poly(e, *nv)
    out = 0
    for (a, b, c), cf in P.terms():
        if a % 2 or b % 2 or c % 2:
            continue
        df = lambda m: sp.factorial2(m - 1) if m > 0 else 1
        out += cf * df(a) * df(b) * df(c) / sp.factorial2(a + b + c + 1)
    return sp.expand(out)


def fmap(F, f):
    return Field({k: f(v) for k, v in F.d.items()})


def integ_poly(e, lo, hi, var=t):
    """Exact definite integral of a polynomial in var."""
    e = sp.expand(e)
    if e == 0:
        return 0
    P = sp.Poly(e, var)
    out = 0
    for (m,), c in P.terms():
        out += c * (hi ** (m + 1) - lo ** (m + 1)) / (m + 1)
    return sp.expand(out)


def phase(S):
    """exp(i K_S . n (t_o - t)) truncated in kappa."""
    K = Ctx.mom(S)
    arg = I * sum(K[i] * nv[i] for i in range(3)) * (to - t)
    return trunc(sum(arg**m / sp.factorial(m) for m in range(KMAX + 1)))


def on_ray(F):
    return Field({S: trunc(v * phase(S)) for S, v in F.d.items()})


def at_origin(F, tt):
    return Field({S: trunc(v.subs(t, tt)) for S, v in F.d.items()})


def taylor_on_ray(F, dX):
    """F(t, Xbar + dX) = sum_n (dX.d)^n F / n!, multilinear (dX has degree >= 1)."""
    out = on_ray(F)
    terms = [(F, Field.const(1))]
    for n in range(1, Ctx.maxdeg + 1):
        new = []
        for G, w in terms:
            for i in range(3):
                new.append((G.d_(i), w * dX[i]))
        terms = new
        acc = Field()
        for G, w in terms:
            acc = acc + w * on_ray(G)
        out = out + acc / sp.factorial(n)
    return out


def shift_eval(Fr, delta, tbar):
    """f(tbar + delta) for a ray function f(t) (Field of time functions), delta of degree >= 1."""
    out = fmap(Fr, lambda v: v.subs(t, tbar))
    dn = Field.const(1)
    G = Fr
    for n in range(1, Ctx.maxdeg + 1):
        dn = dn * delta
        G = G.dt()
        out = out + dn * fmap(G, lambda v: v.subs(t, tbar)) / sp.factorial(n)
    return out


def solve_root(phi, tbar):
    """ln(1 + delta/tbar) + phi(tbar + delta) = 0 for delta (degree >= 1), multilinear."""
    delta = Field()
    for _ in range(Ctx.maxdeg + 1):
        u = delta / tbar
        log_rest = Field()
        un = Field.const(1)
        for m in range(1, Ctx.maxdeg + 1):
            un = un * u
            if m >= 2:
                log_rest = log_rest + un * (sp.Integer(-1) ** (m + 1) / m)
        delta = (shift_eval(phi, delta, tbar) + log_rest) * (-tbar)
        delta = fmap(delta, trunc)
    return delta


def observable(nu, zeta, H):
    """Theta(n) kernels (dict label-set -> polynomial in n, kappa, times) for the given fields."""
    E = equations(nu, zeta, H)
    gam, S, Gam = E["gam"], E["S"], E["Gam"]
    sig = [[gam[i][j].dt() / 2 for j in range(3)] for i in range(3)]
    # initial direction: e(t_o) = -e^{-zeta} exp(-H/2) n at the origin
    from engine import mexp, mscale

    Ginvhalf = mscale(mexp(mscale(H, sp.Rational(-1, 2))), fexp(zeta * (-1)))
    G0 = [[at_origin(Ginvhalf[i][j], to) for j in range(3)] for i in range(3)]
    n = [Field.const(c) for c in nv]
    e0 = [sum((G0[i][j] * (-nv[j]) for j in range(3)), Field()) + nv[i] for i in range(3)]
    e0 = [e0[i].without_const() for i in range(3)]
    dX = [e0[i] * (t - to) for i in range(3)]
    for it in range(Ctx.maxdeg):
        e = [n[i] * (-1) + dX[i].dt() for i in range(3)]
        Sr = [[taylor_on_ray(S[i][j], dX) for j in range(3)] for i in range(3)]
        Gr = [[[taylor_on_ray(Gam[a][i][j], dX) for j in range(3)] for i in range(3)] for a in range(3)]
        sr = [[taylor_on_ray(sig[i][j], dX) for j in range(3)] for i in range(3)]
        see = sum((sr[i][j] * e[i] * e[j] for i in range(3) for j in range(3)), Field())
        A = []
        for i in range(3):
            Ai = sum((Sr[i][j] * e[j] for j in range(3)), Field()) * (-2)
            Ai = Ai - sum((Gr[i][j][k] * e[j] * e[k] for j in range(3) for k in range(3)), Field())
            Ai = Ai + e[i] * see
            A.append(fmap(Ai, trunc))
        dX = [
            Field(
                {
                    Sx: trunc(
                        e0[i].d.get(Sx, 0) * (t - to)
                        + integ_poly((tp - t) * sp.sympify(A[i].d.get(Sx, 0)).subs(t, tp), t, to, tp)
                    )
                    for Sx in set(e0[i].d) | set(A[i].d)
                }
            )
            for i in range(3)
        ]
        LOG(f"ray iteration {it + 1}")
    e = [n[i] * (-1) + dX[i].dt() for i in range(3)]
    sr = [[taylor_on_ray(sig[i][j], dX) for j in range(3)] for i in range(3)]
    F = fmap(sum((sr[i][j] * e[i] * e[j] for i in range(3) for j in range(3)), Field()), trunc)
    phi = fmap(taylor_on_ray(nu, dX), trunc)
    delta = solve_root(phi, ts)
    Iexp = Field({Sx: integ_poly(v, ts, to) for Sx, v in F.d.items()})
    # int_{ts}^{ts+delta} F = sum_{n>=1} delta^n/n! F^{(n-1)}(ts)
    corr = Field()
    dn = Field.const(1)
    G = F
    for m in range(1, Ctx.maxdeg + 1):
        dn = dn * delta
        corr = corr + dn * fmap(G, lambda v: v.subs(t, ts)) / sp.factorial(m)
        G = G.dt()
    Iexp = fmap(Iexp - corr, trunc)
    LOG("exponent I built (ray, emission root, endpoint)")
    Tn = fexp(Iexp * (-1))
    mean = fmap(Tn, sky_mean)
    m1 = (mean - 1).without_const()
    inv = Field.const(1)
    term = Field.const(1)
    for _ in range(Ctx.maxdeg):
        term = term * m1 * (-1)
        inv = inv + term
    Thc = fmap((Tn * inv - 1).without_const(), trunc)
    # clock
    omega = Field({Sx: v for Sx, v in nu.d.items()})  # nu(t, 0): no phase at the origin
    tI_shift = solve_root(omega, tIb)
    eom = fexp(omega)
    Wint = Field({Sx: integ_poly(t * v, 0, tp) for Sx, v in (eom - 1).without_const().d.items()})
    Wf = fmap(Wint, lambda v: v.subs(tp, t))
    D = Field()
    for _ in range(Ctx.maxdeg + 1):
        R = tI_shift * (2 * tIb) + tI_shift * tI_shift - shift_eval(Wf, D, tob) * 2 + shift_eval(Wf, tI_shift, tIb) * 2
        D = fmap((R - D * D) / (2 * tob), trunc)
    Th = fmap(Thc, lambda v: v.subs(to, tob))
    dn = Field.const(1)
    G = Thc
    for m in range(1, Ctx.maxdeg + 1):
        dn = dn * D
        G = fmap(G, lambda v: sp.diff(v, to))
        Th = Th + dn * fmap(G, lambda v: v.subs(to, tob)) / sp.factorial(m)
    LOG("normalised, clock composed")
    return fmap(Th, trunc)


def build_fields(W, MAG, labels, depth):
    tc.W, tc.MAG = W, MAG
    # size the exact field series to the kappa order actually needed (kappa^m <-> t^m), with the
    # same +6 safety margin as the residual tests; every retained coefficient stays exact
    tc.DEPTH = max(depth, 2)
    tc.TCUT = tc.DEPTH + 6
    so.TAYLOR = tc.DEPTH + 6
    second = {}
    pairs = [(a, b) for i, a in enumerate(labels) for b in labels[i + 1 :]]
    for pr in pairs:
        ka, kb = W[pr[0]], W[pr[1]]
        k1k2 = sum(ka[i] * kb[i] for i in range(3))
        c2 = sp.expand((2 * (MAG[pr[0]] ** 2 + MAG[pr[1]] ** 2) - k1k2) / 9)
        nus, As, _, _ = tc.solve_top(list(pr), {}, c2)
        second[pr] = (nus, As)
    third = None
    if len(labels) == 3:
        ks = [W[a] for a in labels]
        sq = sum(MAG[a] ** 2 for a in labels)
        dots = sum(sum(ks[i][c] * ks[j][c] for c in range(3)) for i in range(3) for j in range(i + 1, 3))
        c3 = sp.expand(sp.Rational(2, 9) * (-2 * sq + dots))
        third = tc.solve_top(list(labels), second, c3)[:2]
    Ctx.maxdeg = len(labels)
    Ctx.momenta = dict(W)
    nu, zeta, H = tc.fields(labels, second, (tuple(labels),) + third if third else None)
    return nu, zeta, H


def build_fields_fast(W, MAG, labels, kmax):
    """Same fields as build_fields, solved with the fast exact engine (fast_series.py) at the
    configuration (k, q, -q) and converted to SymPy with momenta scaled by kappa.

    Exactness: the fast engine is exact; a t^m coefficient is homogeneous of degree m in the
    momenta, so keeping t^0..t^(kmax+2) keeps every kappa^(<= kmax) term of every kernel.
    """
    import fast_series as fs

    K12 = (so.p + so.b1) ** 2 + so.b2**2
    K13 = (so.p - so.b1) ** 2 + so.b2**2

    plain = {"p": so.p, "b1": so.b1, "b2": so.b2, "ip": 1 / so.p, "u": 1 / K12, "v": 1 / K13}

    def conv(series):
        # a t^k coefficient is homogeneous of degree k in the momenta, so scaling every momentum
        # by kappa multiplies it by exactly kappa^k; only k <= kmax survives the kappa cut
        out = 0
        for k, c in series.d.items():
            if k > kmax:
                continue
            e = c.as_expr()
            e = e.xreplace({s_: plain[s_.name] for s_ in e.free_symbols})
            out += e * kap**k * t**k
        return out

    second_fast = {}
    for pr in ((1, 2), (1, 3), (2, 3)):
        second_fast[pr] = fs.solve(list(pr), {})
        LOG(f"fast engine: pair {pr} solved")
    top = fs.solve([1, 2, 3], second_fast)
    LOG("fast engine: third order solved")
    second = {pr: (conv(n_), sp.Matrix(3, 3, lambda i, j, A_=A_: conv(A_[i][j]))) for pr, (n_, A_) in second_fast.items()}
    third = (conv(top[0]), sp.Matrix(3, 3, lambda i, j: conv(top[1][i][j])))
    tc.W, tc.MAG = W, MAG
    so.TAYLOR = kmax + 8
    Ctx.maxdeg = 3
    Ctx.momenta = dict(W)
    return tc.fields(labels, second, (tuple(labels),) + third)


if __name__ == "__main__":
    test = sys.argv[1] if len(sys.argv) > 1 else "W2"
    p, q, b1, b2 = so.p, so.q, so.b1, so.b2
    Ctx.trunc = None
    if test in ("W2", "W3"):
        labels = [1, 2] if test == "W2" else [1, 2, 3]
        W = {1: (kap * p, 0, 0), 2: (0, 0, 0), 3: (0, 0, 0)}
        MAG = {1: kap * p, 2: 0, 3: 0}
        nu, zeta, H = build_fields(W, MAG, labels, KMAX)
        LOG("fields ready")
        Ctx.trunc = trunc
        Th = observable(nu, zeta, H)
        base = Th.coeff([1])
        Dk = lambda e: sp.expand(kap * sp.diff(e, kap))
        if test == "W2":
            pred = -Dk(base)
            got = Th.coeff([1, 2])
            lab = "Theta_{sL} = -kappa d_kappa Theta_s"
        else:
            pred = Dk(Dk(base))
            got = Th.coeff([1, 2, 3])
            lab = "Theta_{sLL} = (kappa d_kappa)^2 Theta_s"
        diff = trunc(sp.expand(got - pred))
        assert diff == 0, ("Ward identity of the observable fails", sp.factor(diff))
        assert trunc(base) != 0
        print(f"[ok] {test}: {lab}  through O(kappa^{KMAX}) (exact; all ray, endpoint, clock, frame, normalisation terms)")
    elif test == "SH":
        labels = [1, 2]
        W = {1: (kap * p, 0, 0), 2: (kap * b1, kap * b2, 0)}
        MAG = {1: kap * p, 2: kap * q}
        nu, zeta, H = build_fields(W, MAG, labels, KMAX)
        LOG("fields ready")
        Ctx.trunc = trunc
        Th = observable(nu, zeta, H)
        tau2 = (tob**2 - ts**2) / 6
        k1 = sp.Matrix([kap * p, 0, 0])
        k2 = sp.Matrix([kap * b1, kap * b2, 0])
        nn = sp.Matrix(nv)

        # multilinear coefficient of -tau^2 e^{-2 zeta} TF[zeta_ij - zeta_i zeta_j] n n at the origin
        def tf_nn(M):
            return (nn.T * M * nn)[0, 0] - M.trace() * (nn.T * nn)[0, 0] / 3

        zij = lambda k: -k * k.T
        pred1 = -tau2 * tf_nn(zij(k1))
        pred12 = -tau2 * (tf_nn(-2 * zij(k2)) + tf_nn(-2 * zij(k1)) - tf_nn(-(k1 * k2.T + k2 * k1.T)))
        r1 = trunc(sp.expand(Th.coeff([1]) - pred1))
        assert sp.expand(r1).coeff(kap, 2) == 0 and sp.expand(r1).coeff(kap, 1) == 0 and sp.expand(r1).coeff(kap, 0) == 0
        r12 = so.reduce_b2(sp.expand(sp.expand(Th.coeff([1, 2]) - pred12).coeff(kap, 2)))
        low = [so.reduce_b2(sp.expand(sp.expand(Th.coeff([1, 2])).coeff(kap, m))) for m in (0, 1)]
        assert r12 == 0 and low == [0, 0], (r12, low)
        print("[ok] SH: generator reproduces the superhorizon closed form at O(kappa^2), orders 1 and 2, general (k, q)")
    elif test == "SH3":
        labels = [1, 2, 3]
        W = {1: (kap * p, 0, 0), 2: (kap * b1, kap * b2, 0), 3: (-kap * b1, -kap * b2, 0)}
        MAG = {1: kap * p, 2: kap * q, 3: kap * q}
        nu, zeta, H = build_fields_fast(W, MAG, labels, KMAX)
        LOG("fields ready")
        Ctx.trunc = trunc
        Th = observable(nu, zeta, H)
        tau2 = (tob**2 - ts**2) / 6
        zs = sp.symbols("z1 z2 z3")
        ks = [sp.Matrix(W[a]) for a in labels]
        Z = sum(zs[a] for a in range(3))  # zeta at the origin
        Zi = sum((zs[a] * I * ks[a] for a in range(3)), sp.zeros(3, 1))
        Zij = sum((zs[a] * (-ks[a] * ks[a].T) for a in range(3)), sp.zeros(3, 3))
        M = sp.exp(-2 * Z) * (Zij - Zi * Zi.T)
        nn = sp.Matrix(nv)
        val = (nn.T * M * nn)[0, 0] - M.trace() * (nn.T * nn)[0, 0] / 3
        val = -tau2 * val
        pred = sp.diff(val, zs[0], zs[1], zs[2]).subs({z_: 0 for z_ in zs})
        got = sp.expand(Th.coeff([1, 2, 3]))
        r = so.reduce_b2(sp.expand(sp.expand(got - pred).coeff(kap, 2)))
        low = [so.reduce_b2(sp.expand(got.coeff(kap, m))) for m in (0, 1)]
        assert r == 0 and low == [0, 0], (r, low)
        print("[ok] SH3: generator reproduces the superhorizon closed form at O(kappa^2), third order at (k, q, -q)")
