"""Fast exact series engine for the second- and third-order implementation tests.

Same mathematics as engine.py + third_order_consistency.py, different arithmetic:

* every coefficient is a truncated Laurent series in t, sum_k c_k t^k (k <= TCUT, negative powers
  kept), stored as a dict;
* each c_k is a sparse polynomial over the Gaussian rationals QQ_I in the symbols
      p, b1, b2   (k = (p,0,0), q = (b1,b2,0); q^2 is written b1^2 + b2^2, so no relation is
                    needed) and
      ip, u, v    (formal symbols standing for 1/p, 1/|k+q|^2, 1/|k-q|^2, the only denominators
                    that ever occur).
  All operations are ring operations (+, *, scalar division by integers), which commute with the
  homomorphism ip -> 1/p, u -> 1/|k+q|^2, v -> 1/|k-q|^2.  Zero tests apply that homomorphism
  exactly by clearing denominators (see `vanishes`).  No floating point, no gcd computations.

Exactness of truncation: residual coefficients are checked through t^(DEPTH-2); the equations
lower t-degrees by at most 3 (derivatives and 1/t^2), so terms above TCUT = DEPTH + 6 never reach a
checked coefficient.

The test is stronger than the SymPy version: after the triangular solution is built from the
sources, the FULL nonlinear equations (N), (M), (T), (Ham) are recomputed with the solved kernel
inserted, and every component of the top-order coefficient must vanish.
"""
import sys
import time
import sympy as sp
from sympy import QQ_I
from sympy.polys.rings import ring

R, P, B1, B2, IP, U, V = ring('p b1 b2 ip u v', QQ_I)
ONE, ZERO = R.one, R.zero
IU = R(QQ_I.from_sympy(sp.I))
DEPTH = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 4
TCUT = DEPTH + 6
T0 = time.time()


def LOG(m):
    print(f'  [{time.time() - T0:7.1f}s] {m}', flush=True)


def rat(a, b=1):
    return R(QQ_I.convert(sp.Rational(a, b)))


# ------------------------------------------------------------------ truncated Laurent series in t
class L:
    __slots__ = ('d',)

    def __init__(self, d=None):
        self.d = {k: c for k, c in (d or {}).items() if c != 0 and k <= TCUT}

    @staticmethod
    def c(x, k=0):
        x = x if hasattr(x, 'ring') else rat(x) if not isinstance(x, sp.Basic) else R(QQ_I.from_sympy(x))
        return L({k: x})

    def __add__(self, o):
        if not isinstance(o, L):
            o = L.c(o)
        d = dict(self.d)
        for k, c in o.d.items():
            d[k] = d.get(k, ZERO) + c
        return L(d)

    __radd__ = __add__

    def __neg__(self):
        return L({k: -c for k, c in self.d.items()})

    def __sub__(self, o):
        return self + (-(o if isinstance(o, L) else L.c(o)))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, L):
            o = L.c(o)
        d = {}
        for k1, c1 in self.d.items():
            for k2, c2 in o.d.items():
                k = k1 + k2
                if k > TCUT:
                    continue
                d[k] = d.get(k, ZERO) + c1 * c2
        return L(d)

    __rmul__ = __mul__

    def dt(self):
        return L({k - 1: c * k for k, c in self.d.items() if k != 0})

    def scale(self, s):
        return L({k: c * s for k, c in self.d.items()})

    def coeff(self, k):
        return self.d.get(k, ZERO)

    def is_zero(self):
        return not self.d


TINV = L({-1: ONE})


# ------------------------------------------------------------------ multilinear fields
class Ctx:
    maxdeg = 2
    momenta = {}          # label -> (kx, ky, kz) ring elements

    @classmethod
    def mom(cls, labels):
        v = [ZERO, ZERO, ZERO]
        for a in labels:
            for i in range(3):
                v[i] = v[i] + cls.momenta[a][i]
        return v


class F:
    __slots__ = ('d',)

    def __init__(self, d=None):
        self.d = {k: v for k, v in (d or {}).items() if not v.is_zero()}

    @staticmethod
    def const(x):
        return F({frozenset(): x if isinstance(x, L) else L.c(x)})

    def __add__(self, o):
        o = o if isinstance(o, F) else F.const(o)
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d[k] + v if k in d else v
        return F(d)

    __radd__ = __add__

    def __neg__(self):
        return F({k: -v for k, v in self.d.items()})

    def __sub__(self, o):
        return self + (-(o if isinstance(o, F) else F.const(o)))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, F):
            o = o if isinstance(o, L) else L.c(o)
            return F({k: v * o for k, v in self.d.items()})
        d = {}
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items():
                if k1 & k2:
                    continue
                k = k1 | k2
                if len(k) > Ctx.maxdeg:
                    continue
                pr = v1 * v2
                d[k] = d[k] + pr if k in d else pr
        return F(d)

    __rmul__ = __mul__

    def d_(self, i):
        out = {}
        for k, v in self.d.items():
            m = Ctx.mom(k)[i]
            if m != 0:
                out[k] = v.scale(IU * m)
        return F(out)

    def dt(self):
        return F({k: v.dt() for k, v in self.d.items()})

    def coeff(self, labels):
        return self.d.get(frozenset(labels), L())

    def without_const(self):
        return F({k: v for k, v in self.d.items() if k})


def fexp(X):
    X = X.without_const()
    out, term = F.const(1), F.const(1)
    for n in range(1, Ctx.maxdeg + 1):
        term = term * X * rat(1, n)
        out = out + term
    return out


def mz():
    return [[F() for _ in range(3)] for _ in range(3)]


def madd(A, B, s=1):
    return [[A[i][j] + B[i][j] * s for j in range(3)] for i in range(3)]


def mmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(3)), F()) for j in range(3)] for i in range(3)]


def mscale(A, c):
    return [[A[i][j] * c for j in range(3)] for i in range(3)]


def mexp(H):
    out = [[F.const(1 if i == j else 0) for j in range(3)] for i in range(3)]
    term = [row[:] for row in out]
    for n in range(1, Ctx.maxdeg + 1):
        term = mscale(mmul(term, H), rat(1, n))
        out = madd(out, term)
    return out


def tr(A):
    return A[0][0] + A[1][1] + A[2][2]


def tf(A):
    t3 = tr(A) * rat(1, 3)
    return [[A[i][j] - (t3 if i == j else F()) for j in range(3)] for i in range(3)]


def equations(nu, zeta, H):
    """Identical to engine.equations, in the fast arithmetic."""
    e2z, em2z = fexp(zeta * 2), fexp(zeta * (-2))
    M, Minv = mexp(H), mexp(mscale(H, -1))
    gam, gi = mscale(M, e2z), mscale(Minv, em2z)
    dg = [[[gam[i][j].d_(k) for k in range(3)] for j in range(3)] for i in range(3)]
    low = [[[(dg[b][j][i] + dg[b][i][j] - dg[i][j][b]) * rat(1, 2) for j in range(3)] for i in range(3)]
           for b in range(3)]
    Gam = [[[sum((gi[a][b] * low[b][i][j] for b in range(3)), F()) for j in range(3)] for i in range(3)]
           for a in range(3)]
    Ric = mz()
    for i in range(3):
        for j in range(i, 3):
            s = F()
            for a in range(3):
                s = s + Gam[a][i][j].d_(a) - Gam[a][i][a].d_(j)
                for b in range(3):
                    s = s + Gam[a][a][b] * Gam[b][i][j] - Gam[a][j][b] * Gam[b][i][a]
            Ric[i][j] = Ric[j][i] = s
    Rmix = mmul(gi, Ric)
    Rs = tr(Rmix)
    S = mscale(mmul(Minv, [[M[i][j].dt() for j in range(3)] for i in range(3)]), rat(1, 2))
    trSS = tr(mmul(S, S))
    dnu = [nu.d_(i) for i in range(3)]
    nud = nu.dt()
    DD = [[nu.d_(i).d_(j) - sum((Gam[k][i][j] * dnu[k] for k in range(3)), F()) for j in range(3)]
          for i in range(3)]
    lap = sum((gi[i][j] * DD[i][j] for i in range(3) for j in range(3)), F())
    g2 = sum((gi[i][j] * dnu[i] * dnu[j] for i in range(3) for j in range(3)), F())
    em2nu = fexp(nu * (-2))
    t2inv = TINV * TINV
    eqN = nud.dt() * 3 - lap - g2 * 2 + trSS + (em2nu - 1) * (t2inv * 3)
    eqM = []
    for i in range(3):
        div = F()
        for j in range(3):
            div = div + S[j][i].d_(j)
            for k in range(3):
                div = div + Gam[j][j][k] * S[k][i] - Gam[k][j][i] * S[j][k]
        sj = sum((S[j][i] * dnu[j] for j in range(3)), F())
        eqM.append(div - nud.d_(i) * 2 + dnu[i] * (TINV * 2) + sj * 2 + nud * dnu[i] * 2)
    dndn = [[dnu[i] * dnu[j] for j in range(3)] for i in range(3)]
    inner = madd(madd(madd([[S[i][j].dt() for j in range(3)] for i in range(3)],
                           mscale(S, nud * 2 + F.const(TINV * 2))), Rmix),
                 mmul(gi, madd(mscale(DD, -2), mscale(dndn, 2))))
    eqT = tf(inner)
    ot = F.const(TINV)
    eqH = Rs - lap * 4 - g2 * 2 - trSS + (ot + nud) * (ot + nud) * 6 - em2nu * (t2inv * 6)
    return {'N': eqN, 'M': eqM, 'T': eqT, 'Ham': eqH}


# ------------------------------------------------------------------ exact zero test
K12 = (P + B1)**2 + B2**2
K13 = (P - B1)**2 + B2**2


def vanishes(c):
    """Exact test that a ring element is zero under ip -> 1/p, u -> 1/K12, v -> 1/K13."""
    if c == 0:
        return True
    terms = c.terms()
    A = max(m[3] for m, _ in terms)
    B = max(m[4] for m, _ in terms)
    C = max(m[5] for m, _ in terms)
    g = ZERO
    for m, cf in terms:
        g += R({(m[0], m[1], m[2], 0, 0, 0): cf}) * P**(A - m[3]) * K12**(B - m[4]) * K13**(C - m[5])
    return g == 0


# ------------------------------------------------------------------ first-order Taylor fields
_tx = sp.Symbol('x')
F_SER = sp.Poly(sp.series(sp.sin(_tx) / _tx - sp.cos(_tx), _tx, 0, TCUT + 2).removeO(), _tx)
S1_SER = sp.Poly(sp.series(2 * sp.sin(_tx) / _tx - sp.cos(_tx) - 1, _tx, 0, TCUT + 2).removeO(), _tx)


def lin_fields(label, kvec):
    k2 = sum((c * c for c in kvec), ZERO)
    if k2 == 0:
        return F(), F({frozenset([label]): L.c(1)}), [[F() for _ in range(3)] for _ in range(3)]
    nu = {}
    Q = {}
    for (n,), c in F_SER.terms():          # x^n, n even; x^2 = k^2 t^2 / 3
        nu[n] = -R(QQ_I.convert(c)) * k2**(n // 2) * rat(1, 3**(n // 2))
    for (n,), c in S1_SER.terms():         # Q = -6 (s-1)/k^2
        Q[n] = R(QQ_I.convert(-6 * c)) * k2**(n // 2 - 1) * rat(1, 3**(n // 2))
    Ql = L(Q)
    nuF = F({frozenset([label]): L(nu)})
    zF = F({frozenset([label]): L.c(1)})
    H = [[F({frozenset([label]): Ql.scale(kvec[i] * kvec[j] - (k2 * rat(1, 3) if i == j else ZERO))})
          for j in range(3)] for i in range(3)]
    return nuF, zF, H


# ------------------------------------------------------------------ triangular solution
W = {1: (P, ZERO, ZERO), 2: (B1, B2, ZERO), 3: (-B1, -B2, ZERO)}
INVK2 = {frozenset([1, 2]): U, frozenset([1, 3]): V, frozenset([1, 2, 3]): IP * IP}


def build(labels, second, top=None):
    Ctx.maxdeg = len(labels)
    Ctx.momenta = {a: W[a] for a in labels}
    nu, zeta, H = F(), F(), [[F() for _ in range(3)] for _ in range(3)]
    for a in labels:
        n_, z_, H_ = lin_fields(a, W[a])
        nu, zeta = nu + n_, zeta + z_
        H = [[H[i][j] + H_[i][j] for j in range(3)] for i in range(3)]
    for pair, (nk, Ak) in second.items():
        if set(pair) <= set(labels):
            nu = nu + F({frozenset(pair): nk})
            H = [[H[i][j] + F({frozenset(pair): Ak[i][j]}) for j in range(3)] for i in range(3)]
    if top is not None:
        nk, Ak = top
        lab = frozenset(labels)
        nu = nu + F({lab: nk})
        H = [[H[i][j] + F({lab: Ak[i][j]}) for j in range(3)] for i in range(3)]
    return nu, zeta, H


def c_r(labels):
    ks = [W[a] for a in labels]
    dot = lambda a, b: sum((a[i] * b[i] for i in range(3)), ZERO)
    sq = sum((dot(k, k) for k in ks), ZERO)
    pairs = sum((dot(ks[i], ks[j]) for i in range(len(ks)) for j in range(i + 1, len(ks))), ZERO)
    if len(labels) == 2:
        return (sq * 2 - pairs) * rat(1, 9)
    c3 = (sq * (-2) + pairs) * rat(2, 9)
    if 'mutate-c3' in sys.argv:            # negative control: wrong free lapse coefficient
        c3 = c3 * rat(11, 10)
    return c3


def solve(labels, second):
    lab = frozenset(labels)
    E = equations(*build(labels, second))
    LOG(f'sources for {tuple(labels)}')
    cN = E['N'].coeff(lab)
    cM = [E['M'][i].coeff(lab) for i in range(3)]
    cT = [[E['T'][i][j].coeff(lab) for j in range(3)] for i in range(3)]
    KV = Ctx.mom(lab)
    Ksq = sum((c * c for c in KV), ZERO)
    assert all(vanishes(cN.coeff(m)) for m in range(-6, 2)), 'c_N must start at O(t^2)'
    a = {2: c_r(labels)}
    for m in range(4, TCUT + 1, 2):
        a[m] = (-cN.coeff(m - 2) - Ksq * a.get(m - 2, ZERO)) * rat(1, 3 * (m - 2) * (m + 1))
    nu = L(a)
    A = [[L() for _ in range(3)] for _ in range(3)]
    if Ksq != 0:
        iK2 = INVK2[lab]
        Pm = [[(ONE if i == j else ZERO) - KV[i] * KV[j] * iK2 for j in range(3)] for i in range(3)]
        hSd = (nu.dt() - nu * TINV) * 6
        hS, w = {}, [{}, {}, {}]
        for m in range(-1, TCUT):
            KcM = sum((KV[i] * cM[i].coeff(m) for i in range(3)), ZERO)
            val = hSd.coeff(m) + IU * KcM * iK2 * 3
            if val != 0:
                assert m >= 0, 'h_S integrand singular'
                hS[m + 1] = val * rat(1, m + 1)
            for i in range(3):
                comp = sum((Pm[i][j] * cM[j].coeff(m) for j in range(3)), ZERO)
                if comp != 0:
                    assert m >= 0
                    w[i][m + 1] = comp * IU * iK2 * 2 * rat(1, m + 1)
        hS = L(hS)
        w = [L(x) for x in w]
        hm = {}
        for m in range(2, TCUT + 1):
            Cm = [[cT[i][j].coeff(m - 2) for j in range(3)] for i in range(3)]
            PCP = [[sum((Pm[i][k] * Cm[k][l] * Pm[l][j] for k in range(3) for l in range(3)), ZERO)
                    for j in range(3)] for i in range(3)]
            trPCP = PCP[0][0] + PCP[1][1] + PCP[2][2]
            prev = hm.get(m - 2)
            hm[m] = [[(-(PCP[i][j] - Pm[i][j] * trPCP * rat(1, 2)) * 2
                       - (Ksq * prev[i][j] if prev else ZERO)) * rat(1, m * (m + 1))
                      for j in range(3)] for i in range(3)]
        for i in range(3):
            for j in range(3):
                h_ij = L({m: hm[m][i][j] for m in hm})
                A[i][j] = (hS.scale(KV[i] * KV[j] * iK2) - (hS * rat(1, 3) if i == j else L())
                           + w[j].scale(KV[i]) + w[i].scale(KV[j]) + h_ij)
    else:
        assert all(vanishes(cM[i].coeff(m)) for i in range(3) for m in range(-6, TCUT + 1)), \
            'momentum source must vanish at K = 0'
        LOG('[ok] K = 0: momentum-constraint source vanishes identically')
        for i in range(3):
            for j in range(3):
                A[i][j] = L({m: cT[i][j].coeff(m - 2) * rat(-2, m * (m + 1)) for m in range(2, TCUT + 1)})
    trA = A[0][0] + A[1][1] + A[2][2]
    assert all(vanishes(c) for c in trA.d.values()), 'H must be trace-free'
    return nu, A


def residual_check(labels, second, sol, label):
    lab = frozenset(labels)
    E = equations(*build(labels, second, sol))
    comps = [('N', E['N'])] + [(f'M{i}', E['M'][i]) for i in range(3)] \
        + [(f'T{i}{j}', E['T'][i][j]) for i in range(3) for j in range(3)] + [('Ham', E['Ham'])]
    for name, Fq in comps:
        c = Fq.coeff(lab)
        bad = [m for m in range(-8, DEPTH - 1) if not vanishes(c.coeff(m))]
        assert not bad, (label, name, bad[:3])
    LOG(f'[ok] {label}: every component of the FULL equations (N), (M), (T), (Ham) vanishes '
        f'through O(t^{DEPTH - 2})')


if __name__ == '__main__':
    second = {}
    for pair in ((1, 2), (1, 3), (2, 3)):
        sol = solve(list(pair), {})
        residual_check(list(pair), {}, sol, f'second order {pair}')
        if 'mutate-K2' in sys.argv and pair == (2, 3):   # negative control: corrupt the K = 0 pair
            nu2, A2 = sol
            sol = (nu2 + L({4: P * P * B1 * B1}), A2)
        second[pair] = sol
    sol = solve([1, 2, 3], second)
    residual_check([1, 2, 3], second, sol, 'third order (k, q, -q)')
    print('THIRD-ORDER IMPLEMENTATION TEST PASSED (fast exact engine)')
