"""Exact multilinear plane-wave engine for the temperature-frame field equations.

A field is a dict {frozenset(wave labels): coefficient}.  The seed is
    g(X) = sum_a z_a exp(i k_a . X),
and only multilinear monomials z_{a1}...z_{ar} (distinct labels) are kept.  With the
convention  Y = sum_r eps^r Y^(r)/r!,  Y^(r) = int K_r(k_1..k_r) g(k_1)..g(k_r) e^{i K.X},
the coefficient of z_{a1}..z_{ar} exp(i(k_{a1}+..).X) is exactly the symmetric kernel
K_r(k_{a1},..,k_{ar}) -- no factorials (the r! of the Taylor convention cancels the r!
orderings).  Spatial derivatives act as i * (total momentum of the monomial).

Everything is exact SymPy algebra.
"""
import itertools
import sympy as sp

t = sp.Symbol('t', positive=True)
I = sp.I


class Field:
    __slots__ = ('d',)

    def __init__(self, d=None):
        self.d = {k: v for k, v in (d or {}).items() if v != 0}

    @staticmethod
    def const(c):
        return Field({frozenset(): sp.sympify(c)})

    def __add__(self, o):
        o = o if isinstance(o, Field) else Field.const(o)
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, 0) + v
        return Field(d)

    __radd__ = __add__

    def __neg__(self):
        return Field({k: -v for k, v in self.d.items()})

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Field) else Field.const(o)))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Field):
            o = sp.sympify(o)
            return Field({k: v * o for k, v in self.d.items()})
        d = {}
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items():
                if k1 & k2:
                    continue
                k = k1 | k2
                if len(k) > Ctx.maxdeg:
                    continue
                d[k] = d.get(k, 0) + v1 * v2
        if Ctx.trunc is not None:
            d = {k: Ctx.trunc(v) for k, v in d.items()}
        return Field(d)

    __rmul__ = __mul__

    def __truediv__(self, c):
        return self * (1 / sp.sympify(c))

    def d_(self, i):
        return Field({k: I * Ctx.mom(k)[i] * v for k, v in self.d.items()})

    def dt(self):
        return Field({k: sp.diff(v, t) for k, v in self.d.items()})

    def apply(self, f):
        return Field({k: f(v) for k, v in self.d.items()})

    def coeff(self, labels):
        return self.d.get(frozenset(labels), 0)

    def without_const(self):
        return Field({k: v for k, v in self.d.items() if k})


class Ctx:
    maxdeg = 2
    momenta = {}
    trunc = None   # optional coefficient truncation (e.g. a power series cut), applied after products

    @classmethod
    def mom(cls, labels):
        v = [0, 0, 0]
        for a in labels:
            for i in range(3):
                v[i] += cls.momenta[a][i]
        return v


def fexp(F):
    """exp(F) for a field F without constant term, truncated multilinearly."""
    F = F.without_const()
    out, term = Field.const(1), Field.const(1)
    for n in range(1, Ctx.maxdeg + 1):
        term = term * F / n
        out = out + term
    return out


def mzero():
    return [[Field() for _ in range(3)] for _ in range(3)]


def meye():
    return [[Field.const(1 if i == j else 0) for j in range(3)] for i in range(3)]


def madd(A, B, s=1):
    return [[A[i][j] + B[i][j] * s for j in range(3)] for i in range(3)]


def mmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(3)), Field()) for j in range(3)] for i in range(3)]


def mscale(A, c):
    return [[A[i][j] * c for j in range(3)] for i in range(3)]


def mexp(H):
    out, term = meye(), meye()
    for n in range(1, Ctx.maxdeg + 1):
        term = mscale(mmul(term, H), sp.Rational(1, n))
        out = madd(out, term)
    return out


def mtrace(A):
    return A[0][0] + A[1][1] + A[2][2]


def tf(A):
    tr = mtrace(A) / 3
    return [[A[i][j] - (tr if i == j else 0) for j in range(3)] for i in range(3)]


def equations(nu, zeta, H):
    """Exact temperature-frame equations (unified form), background a = a1 t.

    gamma = e^{2 zeta} exp(H), tr H = 0;  N = a1 t e^{nu}.
    Returns dict with keys 'N', 'M' (list of 3), 'T' (3x3 trace-free), 'Ham'.
      (N)   3 nu'' - Delta nu - 2|D nu|^2 + tr(S S) + (3/t^2)(e^{-2 nu} - 1)
      (M)_i D_j S^j_i - 2 d_i nu' + (2/t) d_i nu + 2 S^j_i d_j nu + 2 nu' d_i nu
      (T)   TF[ S' + 2(1/t + nu') S + R^k_j - 2 D^k D_j nu + 2 D^k nu D_j nu ]
      (Ham) R - 4 Delta nu - 2|D nu|^2 - tr(S S) + 6 (1/t + nu')^2 - 6 e^{-2 nu}/t^2
    """
    e2z = fexp(zeta * 2)
    em2z = fexp(zeta * (-2))
    M = mexp(H)
    Minv = mexp(mscale(H, -1))
    gam = mscale(M, e2z)
    gi = mscale(Minv, em2z)
    dgam = [[[gam[i][j].d_(k) for k in range(3)] for j in range(3)] for i in range(3)]
    # Christoffel Gamma^a_{ij}
    low = [[[(dgam[b][j][i] + dgam[b][i][j] - dgam[i][j][b]) * sp.Rational(1, 2)
             for j in range(3)] for i in range(3)] for b in range(3)]
    Gam = [[[sum((gi[a][b] * low[b][i][j] for b in range(3)), Field()) for j in range(3)]
            for i in range(3)] for a in range(3)]
    Ric = mzero()
    for i in range(3):
        for j in range(i, 3):
            s = Field()
            for a in range(3):
                s = s + Gam[a][i][j].d_(a) - Gam[a][i][a].d_(j)
                for b in range(3):
                    s = s + Gam[a][a][b] * Gam[b][i][j] - Gam[a][j][b] * Gam[b][i][a]
            Ric[i][j] = s
            Ric[j][i] = s
    Rmix = mmul(gi, Ric)
    Rs = mtrace(Rmix)
    S = mscale(mmul(Minv, [[M[i][j].dt() for j in range(3)] for i in range(3)]), sp.Rational(1, 2))
    trSS = mtrace(mmul(S, S))
    dnu = [nu.d_(i) for i in range(3)]
    nudot = nu.dt()
    DDnu = [[nu.d_(i).d_(j) - sum((Gam[k][i][j] * dnu[k] for k in range(3)), Field())
             for j in range(3)] for i in range(3)]
    lap = sum((gi[i][j] * DDnu[i][j] for i in range(3) for j in range(3)), Field())
    grad2 = sum((gi[i][j] * dnu[i] * dnu[j] for i in range(3) for j in range(3)), Field())
    em2nu = fexp(nu * (-2))
    eqN = nudot.dt() * 3 - lap - grad2 * 2 + trSS + (em2nu - 1) * (3 / t**2)
    eqM = []
    for i in range(3):
        div = Field()
        for j in range(3):
            div = div + S[j][i].d_(j)
            for k in range(3):
                div = div + Gam[j][j][k] * S[k][i] - Gam[k][j][i] * S[j][k]
        sj = sum((S[j][i] * dnu[j] for j in range(3)), Field())
        eqM.append(div - nudot.d_(i) * 2 + dnu[i] * (2 / t) + sj * 2 + nudot * dnu[i] * 2)
    dndn = [[dnu[i] * dnu[j] for j in range(3)] for i in range(3)]
    inner = madd(madd(madd([[S[i][j].dt() for j in range(3)] for i in range(3)],
                           mscale(S, nudot * 2 + 2 / t)), Rmix),
                 mmul(gi, madd(mscale(DDnu, -2), mscale(dndn, 2))))
    eqT = tf(inner)
    one_t = Field.const(1 / t)
    eqH = Rs - lap * 4 - grad2 * 2 - trSS + (one_t + nudot) * (one_t + nudot) * 6 - em2nu * (6 / t**2)
    return {'N': eqN, 'M': eqM, 'T': eqT, 'Ham': eqH, 'S': S, 'gam': gam, 'gi': gi, 'Gam': Gam}


def tt_proj(Kv):
    """Projectors for total momentum Kv (list of 3 sympy exprs, |K| > 0)."""
    K2 = sum(c**2 for c in Kv)
    P = sp.Matrix(3, 3, lambda i, j: (1 if i == j else 0) - Kv[i] * Kv[j] / K2)
    def TT(A):
        A = sp.Matrix(A)
        PAP = P * A * P
        return PAP - P * (PAP.trace()) / 2
    return P, TT
