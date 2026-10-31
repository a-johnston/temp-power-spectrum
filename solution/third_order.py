"""Third-order sources of the unified system at the configuration entering C^(13): (k, q, -q).

Waves: 1 -> k = (p, 0, 0), 2 -> q = (b1, b2, 0), 3 -> -q, total momentum K = k.
Second-order kernels enter as ARBITRARY functions of t:
    nu_2(1,2) = n12(t), H_2(1,2) = A12(t) (symmetric traceless), likewise (1,3), (2,3);
the pair (2,3) has zero total momentum (the forced homogeneous mode of the quadratic order),
so no projector is ever applied to it.  At O(eps^3) each equation reads
    L_K[K_3] + c3 = 0,
with the SAME linear operators as at second order (checked), and c3 is an explicit finite
expression in the first-order closed forms and the second-order kernels and their t-derivatives.
"""
import sympy as sp
from engine import Field, Ctx, equations, t, I
import second_order as so

p, q, b1, b2 = so.p, so.q, so.b1, so.b2
W = {1: (p, 0, 0), 2: (b1, b2, 0), 3: (-b1, -b2, 0)}
MAG = {1: p, 2: q, 3: q}


def sym_matrix(name):
    A = sp.zeros(3, 3)
    for i in range(3):
        for j in range(i, 3):
            if (i, j) == (2, 2):
                continue
            A[i, j] = A[j, i] = sp.Function(f'{name}{i}{j}')(t)
    A[2, 2] = -A[0, 0] - A[1, 1]
    return A


def build():
    Ctx.maxdeg = 3
    Ctx.momenta = dict(W)
    nu, zeta = Field(), Field()
    H = [[Field() for _ in range(3)] for _ in range(3)]
    for a in (1, 2, 3):
        n_, z_, H_ = so.lin_fields(a, W[a], MAG[a])
        nu, zeta = nu + n_, zeta + z_
        H = [[H[i][j] + H_[i][j] for j in range(3)] for i in range(3)]
    second = {}
    for pair in ((1, 2), (1, 3), (2, 3)):
        nm = ''.join(map(str, pair))
        nfun = sp.Function(f'n{nm}')(t)
        Am = sym_matrix(f'A{nm}_')
        second[pair] = (nfun, Am)
        nu = nu + Field({frozenset(pair): nfun})
        H = [[H[i][j] + Field({frozenset(pair): Am[i, j]}) for j in range(3)] for i in range(3)]
    n3 = sp.Function('n123')(t)
    A3 = sym_matrix('A123_')
    nu = nu + Field({frozenset((1, 2, 3)): n3})
    H = [[H[i][j] + Field({frozenset((1, 2, 3)): A3[i, j]}) for j in range(3)] for i in range(3)]
    return equations(nu, zeta, H), {'n': n3, 'A': A3}, second


def sources():
    import second_order as so2
    E, unk, second = build()
    # expected linear operator at total momentum K = k = (p,0,0): reuse the second-order builder
    saveKV, saveKsq = so2.KV, so2.Ksq
    so2.KV, so2.Ksq = (p, 0, 0), p**2
    LN, LM, LT = so2.expected_linear(unk)
    so2.KV, so2.Ksq = saveKV, saveKsq
    fnl = [unk['n']] + [e for e in set(unk['A']) if isinstance(e, sp.Function)]
    def split(expr):
        z = {f: 0 for f in fnl}
        src = expr.subs({d: 0 for d in expr.atoms(sp.Derivative) if d.expr in fnl}).subs(z)
        return sp.expand(expr - src), src
    lab = frozenset((1, 2, 3))
    out = {'second': second}
    lin, out['N'] = split(E['N'].coeff(lab))
    assert so.reduce_b2(lin - sp.expand(LN)) == 0, 'N operator (3rd order)'
    out['M'] = []
    for i in range(3):
        lin, s_ = split(E['M'][i].coeff(lab))
        assert so.reduce_b2(lin - sp.expand(LM[i])) == 0, ('M operator (3rd order)', i)
        out['M'].append(s_)
    out['T'] = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            lin, s_ = split(E['T'][i][j].coeff(lab))
            assert so.reduce_b2(lin - sp.expand(LT[i, j])) == 0, ('T operator (3rd order)', i, j)
            out['T'][i, j] = s_
    out['unk'] = unk
    return out


if __name__ == '__main__':
    import time, os
    t0 = time.time()
    out = sources()
    print(f'[ok] third-order linear operators at (k, q, -q) identical to second order  [{time.time() - t0:.0f}s]')
    os.makedirs(os.path.join(os.path.dirname(__file__), 'results'), exist_ok=True)
    with open(os.path.join(os.path.dirname(__file__), 'results', 'third_order_sources.txt'), 'w') as fh:
        fh.write('# Exact third-order sources c3 at (k, q, -q); L_K[K_3] + c3 = 0.\n')
        fh.write('# k = (p,0,0), q = (b1,b2,0), b1^2+b2^2 = q^2; second-order kernels n12,n13,n23, A12_ij, ...\n')
        fh.write(f'c3_N = {sp.srepr(out["N"]) if False else out["N"]}\n')
        for i in range(3):
            fh.write(f'c3_M[{i}] = {out["M"][i]}\n')
        for i in range(3):
            for j in range(i, 3):
                fh.write(f'c3_T[{i},{j}] = {out["T"][i, j]}\n')
    sizes = {'N': sp.count_ops(out['N']), 'M': [sp.count_ops(m) for m in out['M']]}
    print('  operation counts:', sizes)
