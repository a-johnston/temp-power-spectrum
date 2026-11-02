"""Implementation test at THIRD order: the triangular solution satisfies all ten equations.

Configuration (k, q, -q) with k = (p,0,0), q = (b1,b2,0), b1^2 + b2^2 = q^2 (general momenta).
Every kernel is carried as an exact Taylor polynomial in t:
  1. first-order fields (Taylor polynomials of the closed forms);
  2. second-order pair kernels (k,q), (k,-q) [general K] and (q,-q) [exact K = 0 sector, no
     projectors], each solved by the triangular prescription of the derivation (Sec. 4.2);
  3. the third-order kernel at total momentum k by the same prescription, with
     c3 = -(2/9)(2k^2 + 5q^2);
  4. residuals of (N), (M), (T) [all components] and (Ham) at order z1 z2 z3.
Usage: python third_order_consistency.py DEPTH   (residuals checked through t^(DEPTH-2)).
"""

import sys
import time

import second_order as so
import sympy as sp
from engine import Ctx, Field, I, equations, t

DEPTH = int(sys.argv[1]) if len(sys.argv) > 1 else 6
so.TAYLOR = DEPTH + 6
p, q, b1, b2 = so.p, so.q, so.b1, so.b2
W = {1: (p, 0, 0), 2: (b1, b2, 0), 3: (-b1, -b2, 0)}
MAG = {1: p, 2: q, 3: q}
T0 = time.time()


def LOG(m):
    print(f"  [{time.time() - T0:7.1f}s] {m}", flush=True)


def rred(e):
    n, d = sp.fraction(sp.cancel(sp.together(e)))
    return sp.cancel(so.reduce_b2(n) / so.reduce_b2(d))


TCUT = DEPTH + 6  # highest power of t retained inside products


def trunc_t(e):
    """Drop powers of t above TCUT; keep all lower and negative powers.

    Exactness: every residual coefficient checked is of order <= DEPTH - 2.  A coefficient of
    t^m in a product depends only on factor terms of degree <= m + (degree lost later); the
    equations lower degrees by at most 3 (two t-derivatives and a 1/t^2, or one of each), so
    terms above DEPTH + 6 can never reach a checked coefficient.  Negative powers are kept, so
    the 'no negative power survives' assertion is unaffected.
    """
    e = sp.expand(e)
    if not e.has(t):
        return e
    out = 0
    for term in sp.Add.make_args(e):
        c, k = term.as_coeff_exponent(t)
        if c.has(t) or not k.is_integer or k <= TCUT:  # non-monomial terms are kept untouched
            out += term
    return out


def coeffs(e, n=DEPTH):
    s = sp.expand(e)
    out = [0] * (n + 1)
    if s == 0:
        return out
    P = sp.Poly(sp.expand(s * t**10), t)
    for (m,), c in P.terms():
        m -= 10
        if m < 0:
            assert rred(c) == 0, ("negative power", m)
        elif m <= n:
            out[m] = rred(c)
    return out


def sym_matrix(name):
    A = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            if (i, j) == (2, 2):
                continue
            A[i, j] = A[j, i] = sp.Function(f"{name}{i}{j}")(t)
    A[2, 2] = -A[0, 0] - A[1, 1]
    return A


def linear_ops(KV, nn, A):
    Ksq = sum(c**2 for c in KV)
    Kv = sp.Matrix(KV)
    LN = 3 * nn.diff(t, 2) + Ksq * nn - 6 * nn / t**2
    LM = [I / 2 * sum(KV[j] * A[j, i].diff(t) for j in range(3)) - 2 * I * KV[i] * (nn.diff(t) - nn / t) for i in range(3)]
    R1 = sp.Matrix(
        3, 3, lambda i, j: (Ksq * A[i, j] - sum(KV[a] * KV[i] * A[a, j] + KV[a] * KV[j] * A[a, i] for a in range(3))) / 2
    )
    R1 = R1 - sp.eye(3) * R1.trace() / 3
    LT = A.diff(t, 2) / 2 + A.diff(t) / t + R1 + 2 * (Kv * Kv.T - sp.eye(3) * Ksq / 3) * nn
    return LN, LM, LT


def fields(labels, second=None, top=None):
    nu, zeta = Field(), Field()
    H = [[Field() for _ in range(3)] for _ in range(3)]
    for a in labels:
        n_, z_, H_ = so.lin_fields(a, W[a], MAG[a])
        nu, zeta = nu + n_, zeta + z_
        H = [[H[i][j] + H_[i][j] for j in range(3)] for i in range(3)]
    for pair, (nk, Ak) in (second or {}).items():
        nu = nu + Field({frozenset(pair): nk})
        H = [[H[i][j] + Field({frozenset(pair): Ak[i, j]}) for j in range(3)] for i in range(3)]
    if top is not None:
        lab, nk, Ak = top
        nu = nu + Field({frozenset(lab): nk})
        H = [[H[i][j] + Field({frozenset(lab): Ak[i, j]}) for j in range(3)] for i in range(3)]
    return nu, zeta, H


def solve_top(labels, second, c_r):
    """Solve the top-order kernel for the given wave labels; return (nu_series, A_series, E, unk)."""
    Ctx.maxdeg = len(labels)
    Ctx.momenta = {a: W[a] for a in labels}
    nn = sp.Function("NTOP")(t)
    A = sym_matrix("ATOP")
    nu, zeta, H = fields(labels, second, (labels, nn, A))
    old, Ctx.trunc = Ctx.trunc, trunc_t
    try:
        E = equations(nu, zeta, H)
    finally:
        Ctx.trunc = old
    lab = frozenset(labels)
    fnl = [nn] + [e for e in set(A) if isinstance(e, sp.Function)]
    z = {f: 0 for f in fnl}

    def src(expr):
        expr = sp.sympify(expr)
        return expr.subs({d: 0 for d in expr.atoms(sp.Derivative) if d.expr in fnl}).subs(z)

    KV = [sp.expand(sum(W[a][i] for a in labels)) for i in range(3)]
    Ksq = so.reduce_b2(sum(c**2 for c in KV))
    cN = coeffs(src(E["N"].coeff(lab)))
    cM = [coeffs(src(E["M"][i].coeff(lab))) for i in range(3)]
    cT = [[coeffs(src(E["T"][i][j].coeff(lab))) for j in range(3)] for i in range(3)]
    LOG(f"sources for {tuple(labels)} as series")
    a = [0] * (DEPTH + 3)
    a[2] = c_r
    for m in range(4, DEPTH + 3):
        a[m] = rred((-cN[m - 2] - Ksq * a[m - 2]) / (3 * (m - 2) * (m + 1)))
    assert all(c == 0 for c in cN[:2])
    nus = sum(a[m] * t**m for m in range(DEPTH + 3))
    As = sp.zeros(3, 3)
    if Ksq != 0:
        P = sp.Matrix(3, 3, lambda i, j: (1 if i == j else 0) - KV[i] * KV[j] / Ksq)
        hSd = sp.expand(6 * (nus.diff(t) - nus / t))
        hS, w = 0, [0, 0, 0]
        for m in range(DEPTH + 1):
            KcM = sum(KV[i] * cM[i][m] for i in range(3))
            hS += rred(hSd.coeff(t, m) + 3 * I * KcM / Ksq) * t ** (m + 1) / (m + 1)
            for i in range(3):
                w[i] += rred(2 * I * sum(P[i, j] * cM[j][m] for j in range(3)) / Ksq) * t ** (m + 1) / (m + 1)
        hm = [sp.zeros(3, 3) for _ in range(DEPTH + 3)]
        h = sp.zeros(3, 3)
        for m in range(2, DEPTH + 3):
            Cm = sp.Matrix(3, 3, lambda i, j, m=m: cT[i][j][m - 2])
            PMP = P * Cm * P
            hm[m] = ((-2 * (PMP - P * PMP.trace() / 2) - Ksq * hm[m - 2]) / (m * (m + 1))).applyfunc(rred)
            h += hm[m] * t**m
        As = sp.Matrix(
            3,
            3,
            lambda i, j: (KV[i] * KV[j] / Ksq - sp.Rational(1 if i == j else 0, 3)) * hS + KV[i] * w[j] + w[i] * KV[j] + h[i, j],
        )
    else:
        assert all(c == 0 for comp in cM for c in comp), "momentum source must vanish at K = 0"
        for m in range(2, DEPTH + 3):
            Cm = sp.Matrix(3, 3, lambda i, j, m=m: cT[i][j][m - 2])
            As += (-2 * Cm / (m * (m + 1))) * t**m
    As = As.applyfunc(lambda e: sum(rred(sp.expand(e).coeff(t, m)) * t**m for m in range(DEPTH + 3)))
    return nus, As, E, (nn, A, KV, fnl, src, lab)


def residuals(nus, As, E, info, label):
    nn, A, KV, fnl, src, lab = info
    fns = {nn: nus}
    for i in range(3):
        for j in range(i, 3):
            if isinstance(A[i, j], sp.Function):
                fns[A[i, j]] = As[i, j]
    LN, LM, LT = linear_ops(KV, nn, A)
    full_H = sp.sympify(E["Ham"].coeff(lab))
    srcH = src(full_H)
    linH = sp.expand(full_H - srcH)
    items = (
        [("N", LN, E["N"].coeff(lab))]
        + [(f"M{i}", LM[i], E["M"][i].coeff(lab)) for i in range(3)]
        + [(f"T{i}{j}", LT[i, j], E["T"][i][j].coeff(lab)) for i in range(3) for j in range(3)]
        + [("Ham", linH, full_H)]
    )
    for key, L, full in items:
        lin = coeffs(sp.expand(L.subs(fns).doit()), DEPTH - 2)
        s_ = coeffs(src(full), DEPTH - 2)
        bad = [m for m in range(DEPTH - 1) if rred(lin[m] + s_[m]) != 0]
        assert not bad, (label, key, bad[:3])
    LOG(f"[ok] {label}: every component of (N), (M), (T), (Ham) vanishes through O(t^{DEPTH - 2})")


if __name__ == "__main__":
    second = {}
    c2 = {(1, 2): (2 * (p**2 + q**2) - p * b1) / 9, (1, 3): (2 * (p**2 + q**2) + p * b1) / 9, (2, 3): sp.Rational(5, 9) * q**2}
    for pair in ((1, 2), (1, 3), (2, 3)):
        nus, As, E, info = solve_top(list(pair), {}, sp.expand(c2[pair]))
        residuals(nus, As, E, info, f"second order {pair}")
        second[pair] = (nus, As)
    c3 = -sp.Rational(2, 9) * (2 * p**2 + 5 * q**2)
    nus, As, E, info = solve_top([1, 2, 3], second, c3)
    residuals(nus, As, E, info, "third order (k, q, -q)")
    print("THIRD-ORDER IMPLEMENTATION TEST PASSED")
