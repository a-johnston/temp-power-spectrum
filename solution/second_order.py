"""Second-order sources of the unified temperature-frame system, exact, general momenta.

Frame: k1 = (p, 0, 0), k2 = (b1, b2, 0) with b1^2 + b2^2 = q^2.  This is a rotation of a generic
pair, so scalar kernels are general; tensor kernels are rotated covariantly afterwards.
At O(eps^2) every equation reads  L[K_2](t) + c(t) = 0  for the coefficient of z1 z2 e^{i K.X},
K = k1 + k2.  This module returns c for (N), (M), (T), (Ham) and checks that the linear part
acting on the unknown second-order kernels is the expected operator:

  L_N  = 3 nu'' + K^2 nu - 6 nu / t^2
  L_M  = (i/2) K_j A'_{ji} - 2 i K_i (nu' - nu/t)
  L_T  = (1/2)(A'' + 2 A'/t) + R1TF[A] + 2 (K_i K_j - delta K^2/3) nu,
         R1[A]_ij = (1/2)(K^2 A_ij - K_a K_i A_aj - K_a K_j A_ai)    (A traceless)
"""
import sympy as sp
from engine import Field, Ctx, equations, t, I

p, q = sp.symbols('p q', positive=True)
b1, b2 = sp.symbols('b1 b2', real=True)
s3 = sp.sqrt(3)
K1 = (p, 0, 0)
K2 = (b1, b2, 0)
KV = tuple(K1[i] + K2[i] for i in range(3))
Ksq = sp.expand(sum(c**2 for c in KV))
on_shell = {b2**2: q**2 - b1**2}


TAYLOR = None  # if an int M: use exact Taylor polynomials of the linear time functions to t^M


def _taylor(expr_of_x, x, M):
    return sp.expand(sp.series(expr_of_x, x, 0, M + 1).removeO())


def lin_fields(label, kvec, kmag):
    x = kmag * t / s3
    if TAYLOR is None:
        f = sp.sin(x) / x - sp.cos(x)
        s = 2 * sp.sin(x) / x - sp.cos(x)
    else:
        X_ = sp.Symbol('X_')
        f = _taylor(sp.sin(X_) / X_ - sp.cos(X_), X_, TAYLOR).subs(X_, x)
        s = _taylor(2 * sp.sin(X_) / X_ - sp.cos(X_), X_, TAYLOR).subs(X_, x)
    Q = -6 * (s - 1) / kmag**2
    nu = Field({frozenset([label]): -f})
    zeta = Field({frozenset([label]): sp.Integer(1)})
    H = [[Field({frozenset([label]): Q * (kvec[i] * kvec[j] - (kmag**2 if i == j else 0) / 3)})
          for j in range(3)] for i in range(3)]
    return nu, zeta, H


def reduce_b2(e):
    """Eliminate even powers of b2 using b2^2 = q^2 - b1^2 (exact)."""
    e = sp.expand(e)
    poly = sp.Poly(e, b2)
    out = 0
    for (m,), c in poly.terms():
        out += c * (q**2 - b1**2)**(m // 2) * b2**(m % 2)
    return sp.expand(out)


def build(unknowns=True):
    Ctx.maxdeg = 2
    Ctx.momenta = {1: K1, 2: K2}
    n1, z1, H1 = lin_fields(1, K1, p)
    n2, z2, H2 = lin_fields(2, K2, q)
    nu = n1 + n2
    zeta = z1 + z2
    H = [[H1[i][j] + H2[i][j] for j in range(3)] for i in range(3)]
    unk = {}
    if unknowns:
        nn = sp.Function('n2')(t)
        names = [[None] * 3 for _ in range(3)]
        A = sp.zeros(3, 3)
        for i in range(3):
            for j in range(i, 3):
                if (i, j) == (2, 2):
                    continue
                A[i, j] = A[j, i] = sp.Function(f'A{i}{j}')(t)
        A[2, 2] = -A[0, 0] - A[1, 1]
        nu = nu + Field({frozenset([1, 2]): nn})
        H = [[H[i][j] + Field({frozenset([1, 2]): A[i, j]}) for j in range(3)] for i in range(3)]
        unk = {'n': nn, 'A': A}
    return equations(nu, zeta, H), unk


def split(expr, unk):
    """expr = linear(unknowns) + source; return (linear, source)."""
    src = expr
    reps = {unk['n']: 0}
    for a in unk['A'].free_symbols | set(unk['A'].atoms(sp.Function)):
        if isinstance(a, sp.Function):
            reps[a] = 0
    src = expr.subs({d: 0 for d in expr.atoms(sp.Derivative)}).subs(reps)
    return sp.expand(expr - src), src


def expected_linear(unk):
    nn, A = unk['n'], unk['A']
    Kv = sp.Matrix(KV)
    LN = 3 * nn.diff(t, 2) + Ksq * nn - 6 * nn / t**2
    LM = [I / 2 * sum(KV[j] * A[j, i].diff(t) for j in range(3)) - 2 * I * KV[i] * (nn.diff(t) - nn / t)
          for i in range(3)]
    R1 = sp.Matrix(3, 3, lambda i, j: (Ksq * A[i, j] - sum(KV[a] * KV[i] * A[a, j] + KV[a] * KV[j] * A[a, i]
                                                         for a in range(3))) / 2)
    R1 = R1 - sp.eye(3) * R1.trace() / 3
    LT = A.diff(t, 2) / 2 + A.diff(t) / t + R1 + 2 * (Kv * Kv.T - sp.eye(3) * Ksq / 3) * nn
    return LN, LM, LT


def sources():
    """Exact second-order sources c_N, c_M (3-vector), c_T (3x3), c_H, as functions of t."""
    E, unk = build(unknowns=True)
    LN, LM, LT = expected_linear(unk)
    out = {}
    lin, src = split(E['N'].coeff([1, 2]), unk)
    assert sp.expand(lin - sp.expand(LN)) == 0, 'N operator'
    out['N'] = src
    out['M'] = []
    for i in range(3):
        lin, src = split(E['M'][i].coeff([1, 2]), unk)
        assert sp.expand(lin - sp.expand(LM[i])) == 0, ('M operator', i)
        out['M'].append(src)
    out['T'] = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            lin, src = split(E['T'][i][j].coeff([1, 2]), unk)
            assert reduce_b2(lin - LT[i, j]) == 0, ('T operator', i, j)
            out['T'][i, j] = src
    out['Ham_full'] = E['Ham'].coeff([1, 2])
    out['unk'] = unk
    return out


if __name__ == '__main__':
    import time
    t0 = time.time()
    out = sources()
    print(f'[ok] second-order linear operators (N),(M),(T) confirmed  [{time.time() - t0:.1f}s]')
    import os
    os.makedirs(os.path.join(os.path.dirname(__file__), 'results'), exist_ok=True)
    with open(os.path.join(os.path.dirname(__file__), 'results', 'second_order_sources.txt'), 'w') as fh:
        fh.write('# Exact second-order sources: L_K[K_2] + c = 0 for the z1 z2 coefficient.\n')
        fh.write('# Frame k1 = (p,0,0), k2 = (b1,b2,0), b1^2 + b2^2 = q^2; K = k1 + k2.\n')
        fh.write('# Linear operators: see module docstring of second_order.py.\n')
        fh.write(f'c_N = {out["N"]}\n')
        for i in range(3):
            fh.write(f'c_M[{i}] = {out["M"][i]}\n')
        for i in range(3):
            for j in range(i, 3):
                fh.write(f'c_T[{i},{j}] = {out["T"][i, j]}\n')
    print('  c_N has', sp.count_ops(out['N']), 'ops; sources written to results/second_order_sources.txt')
