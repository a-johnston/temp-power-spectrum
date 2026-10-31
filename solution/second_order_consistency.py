"""Consistency of the triangular second-order solution with ALL field equations.

Solution prescription (derivations/UNIFIED_DERIVATION.md, Sec. 4):
  nu_2 : 3 nu'' + K^2 nu - 6 nu/t^2 = -c_N,   nu_2 = c2 t^2 + O(t^4),
         c2 = [2(k1^2 + k2^2) - k1.k2]/9   (t^2 coefficient of -R[e^{2 zeta} delta]/36, exact in zeta)
  K != 0:
    scalar  h_S' = 6 (nu' - nu/t) + 3 i (K.c_M)/K^2                (momentum constraint, longitudinal)
    vector  w'   = (2 i/K^2) P_perp c_M,  H^V = K w + w K          (momentum constraint, transverse)
    tensor  h'' + 2h'/t + K^2 h = -2 TT[c_T]                       (trace-free evolution, TT part)
    H_2 = (K K/K^2 - delta/3) h_S + H^V + h,  all vanishing at t = 0
  K = 0 (k2 = -k1): no projectors;  (1/2)(A'' + 2A'/t) = -c_T,  A(0) = 0, and c_M must vanish.

Checked here, as exact Taylor coefficients in t with general momenta (symbolic p, q, b1, b2):
the scalar and vector parts of (T) and the Hamiltonian constraint -- none of which were used to
build the solution -- vanish, and c_M = 0 identically at K = 0.  The residuals obey the
homogeneous constraint-propagation system, whose regular solution is fixed by finitely many
leading Taylor coefficients; see the derivation, Sec. 4.4.
"""
import sys
import sympy as sp
from engine import t, I
import second_order as so

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10
p, q, b1, b2 = so.p, so.q, so.b1, so.b2
# Linear time functions -> exact Taylor polynomials to t^(NMAX+6): every coefficient through
# t^NMAX of every derived quantity is then exact.
so.TAYLOR = NMAX + 6


def ser(e, n=NMAX):
    """Exact Taylor coefficients [t^0..t^n] of a Laurent polynomial in t; negative powers must cancel."""
    s = sp.expand(e)
    if s == 0:
        return [0] * (n + 1)
    lo = sp.Poly(sp.expand(s * t**8), t)
    out = [0] * (n + 1)
    for (m,), c in lo.terms():
        m -= 8
        if m < 0:
            assert rred(c) == 0, ('negative power survives', m)
        elif m <= n:
            out[m] = rred(c)
    return out


def rred(e):
    """Exact reduction of a rational function using b2^2 = q^2 - b1^2 in numerator and denominator."""
    n, d = sp.fraction(sp.cancel(sp.together(e)))
    return sp.cancel(so.reduce_b2(n) / so.reduce_b2(d))


def split(expr, fnlist):
    z = {f: 0 for f in fnlist}
    src = expr.subs({d: 0 for d in expr.atoms(sp.Derivative)}).subs(z)
    return sp.expand(expr - src), src


def run(general=True):
    src = so.sources()
    unk = src['unk']
    fnl = [unk['n']] + [e for e in set(unk['A']) if isinstance(e, sp.Function)]
    sub = {} if general else {b1: -p, b2: 0, q: p}
    KV = [sp.expand(sp.sympify(c).subs(sub)) for c in so.KV]
    Ksq = so.reduce_b2(sp.expand(so.Ksq.subs(sub)))
    c2 = sp.expand(((2 * (p**2 + q**2) - p * b1) / 9).subs(sub))
    cN = ser(src['N'].subs(sub))
    cM = [ser(c.subs(sub)) for c in src['M']]
    cT = [[ser(src['T'][i, j].subs(sub)) for j in range(3)] for i in range(3)]
    assert all(c == 0 for c in cN[:2]), 'c_N must start at O(t^2)'
    a = [0] * (NMAX + 3)
    a[2] = c2
    for m in range(4, NMAX + 3):
        a[m] = sp.cancel((-cN[m - 2] - Ksq * a[m - 2]) / (3 * (m - 2) * (m + 1)))
    assert all(sp.cancel(a[m]) == 0 for m in range(3, NMAX + 3, 2)), 'odd powers in nu_2'
    nu = sum(a[m] * t**m for m in range(NMAX + 3))
    A = sp.zeros(3, 3)
    if general:
        P = sp.Matrix(3, 3, lambda i, j: (1 if i == j else 0) - KV[i] * KV[j] / Ksq)
        hSd = sp.expand(6 * (nu.diff(t) - nu / t))
        hS = 0
        w = [0, 0, 0]
        for m in range(NMAX + 1):
            KcM = sum(KV[i] * cM[i][m] for i in range(3))
            hS += sp.cancel(hSd.coeff(t, m) + 3 * I * KcM / Ksq) * t**(m + 1) / (m + 1)
            for i in range(3):
                comp = sum(P[i, j] * cM[j][m] for j in range(3))
                w[i] += sp.cancel(2 * I * comp / Ksq) * t**(m + 1) / (m + 1)
        def TTm(Mx):
            PMP = P * Mx * P
            return PMP - P * PMP.trace() / 2
        h = sp.zeros(3, 3)
        hm = [sp.zeros(3, 3) for _ in range(NMAX + 3)]
        for m in range(2, NMAX + 3):
            Cm = sp.Matrix(3, 3, lambda i, j: cT[i][j][m - 2])
            hm[m] = ((-2 * TTm(Cm) - Ksq * hm[m - 2]) / (m * (m + 1))).applyfunc(sp.cancel)
            h += hm[m] * t**m
        A = sp.Matrix(3, 3, lambda i, j: (KV[i] * KV[j] / Ksq - sp.Rational(1 if i == j else 0, 3)) * hS
                      + KV[i] * w[j] + w[i] * KV[j] + h[i, j])
    else:
        assert all(sp.cancel(c) == 0 for comp in cM for c in comp), 'c_M must vanish at K = 0'
        print('[ok] K = 0: momentum-constraint source vanishes identically', flush=True)
        for m in range(2, NMAX + 3):
            Cm = sp.Matrix(3, 3, lambda i, j: cT[i][j][m - 2])
            A += (-2 * Cm / (m * (m + 1))) * t**m
    A = A.applyfunc(lambda e: sum(rred(sp.expand(e).coeff(t, m)) * t**m for m in range(NMAX + 4)))
    assert all(rred(sp.expand(A.trace()).coeff(t, m)) == 0 for m in range(NMAX + 4))
    fns = {unk['n']: nu}
    for i in range(3):
        for j in range(i, 3):
            if isinstance(unk['A'][i, j], sp.Function):
                fns[unk['A'][i, j]] = A[i, j]
    LN, LM, LT = so.expected_linear(unk)
    linH, srcH = split(src['Ham_full'], fnl)
    def app(e):
        return sp.expand(e.subs(sub).subs(fns).doit())
    res = {'N': (app(LN), src['N'])}
    for i in range(3):
        res[f'M{i}'] = (app(LM[i]), src['M'][i])
        for j in range(3):
            res[f'T{i}{j}'] = (app(LT[i, j]), src['T'][i, j])
    res['Ham'] = (app(linH), srcH)
    return sub, res


if __name__ == '__main__':
    for general in (True, False):
        lab = 'general K' if general else 'K = 0'
        sub, res = run(general)
        for key, (lin, srce) in res.items():
            L = ser(lin, NMAX - 2)
            Sx = ser(sp.expand(srce.subs(sub)), NMAX - 2)
            bad = [m for m in range(NMAX - 1) if rred(L[m] + Sx[m]) != 0]
            assert not bad, (lab, key, bad[:3])
        print(f'[ok] {lab}: every component of (N), (M), (T) and (Ham) vanishes through O(t^{NMAX - 2})',
              flush=True)
