"""Fast exact port of the independent observable generator (observable_generator.py).

Same mathematics, operation for operation: the measured map is solved from its definitions
(ray equation, emission root, exponent, sky normalisation, observer clock, composition with the
random observation time). Only the arithmetic changes:

* a coefficient is a sparse polynomial over the Gaussian rationals QQ_I in
      kap (momentum scaling), t, tp (integration variable), to (coordinate observation time),
      ts, tI, tob (background emission, initial and observation times) and its, itI, itob
      (formal symbols for their inverses), n1, n2, n3 (sky direction), p, b1, b2 (momenta,
      k = (p,0,0), q = (b1,b2,0)), and ip, u, v (formal 1/p, 1/|k+q|^2, 1/|k-q|^2, as in
      fast_series.py);
* every product is truncated at kap^KMAX and the pairs (ts, its), (tI, itI), (tob, itob) are
  cancelled;
* the fields come from the fast exact series engine (fast_series.py); a t^m coefficient of a
  kernel is homogeneous of degree m in the momenta and is multiplied by kap^m.

Modes (python observable_fast.py MODE [KMAX]):
  SH    superhorizon closed form at O(kap^2), orders 1-2, general (k, q)     [cross-check vs slow]
  W2    observable dilation Ward identity, second order, through O(kap^KMAX)  [cross-check vs slow]
  W3    observable dilation Ward identity, third order, through O(kap^KMAX)   [cross-check vs slow]
  SH3   superhorizon closed form at O(kap^2), third order at (k, q, -q)       [new: slow version did not finish]
A third argument mutates the map. noezeta (drop e^{-zeta} in the initial direction) breaks dilation
covariance and must fail W2/W3; noray, nonorm and noend do NOT affect W2/W3, because a zero-momentum
long mode neither deflects the ray nor shifts the endpoint, and every piece of the map is separately
dilation covariant. The dilation Ward tests therefore check the fields and the covariance, not the
ray, endpoint or normalisation terms.
"""

import sys
import time

import sympy as sp
from sympy import QQ_I
from sympy.polys.rings import ring

MODE = sys.argv[1] if len(sys.argv) > 1 else "SH"
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2
MUTATE = sys.argv[3] if len(sys.argv) > 3 else ""  # controls: noray, nonorm, noend (invisible to W2/W3), noezeta (must fail)
_argv, sys.argv = sys.argv, sys.argv[:1]
import fast_series as fs  # noqa: E402  (fast_series parses argv at import)

sys.argv = _argv

R, KAP, T, TP, TO, TS, TI, TOB, ITS, ITI, ITOB, N1, N2, N3, P, B1, B2, IP, U, V = ring(
    "kap t tp to ts tI tob its itI itob n1 n2 n3 p b1 b2 ip u v", QQ_I
)
ONE, ZERO = R.one, R.zero
IU = R(QQ_I.from_sympy(sp.I))
NV = (N1, N2, N3)
T0 = time.time()
_IKAP = 0
_INV = ((R.gens.index(TS), R.gens.index(ITS)), (R.gens.index(TI), R.gens.index(ITI)), (R.gens.index(TOB), R.gens.index(ITOB)))
_IT, _ITP, _ITO = R.gens.index(T), R.gens.index(TP), R.gens.index(TO)
_IN = tuple(R.gens.index(x) for x in NV)


def LOG(m):
    print(f"  [{time.time() - T0:7.1f}s] {m}", flush=True)


def rat(a, b=1):
    return R(QQ_I.convert(sp.Rational(a, b)))


def norm(poly):
    """Truncate at kap^KMAX and cancel ts*its, tI*itI, tob*itob."""
    out = {}
    for mon, c in poly.items():
        if mon[_IKAP] > KMAX:
            continue
        m = list(mon)
        for a, b in _INV:
            k = min(m[a], m[b])
            if k:
                m[a] -= k
                m[b] -= k
        m = tuple(m)
        out[m] = out.get(m, QQ_I.zero) + c
    return R({m: c for m, c in out.items() if c})


def as_ring(x):
    if hasattr(x, "ring"):
        return x
    return R(QQ_I.convert(sp.sympify(x)))


# ------------------------------------------------------------------ multilinear fields
class Ctx:
    maxdeg = 2
    momenta = {}

    @classmethod
    def mom(cls, labels):
        v = [ZERO, ZERO, ZERO]
        for a in labels:
            for i in range(3):
                v[i] = v[i] + cls.momenta[a][i]
        return v


class Fld:
    __slots__ = ("d",)

    def __init__(self, d=None):
        self.d = {k: v for k, v in (d or {}).items() if v}

    @staticmethod
    def const(c):
        return Fld({frozenset(): as_ring(c)})

    def __add__(self, o):
        o = o if isinstance(o, Fld) else Fld.const(o)
        d = dict(self.d)
        for k, v in o.d.items():
            d[k] = d.get(k, ZERO) + v
        return Fld(d)

    __radd__ = __add__

    def __neg__(self):
        return Fld({k: -v for k, v in self.d.items()})

    def __sub__(self, o):
        return self + (-(o if isinstance(o, Fld) else Fld.const(o)))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Fld):
            o = as_ring(o)
            return Fld({k: norm(v * o) for k, v in self.d.items()})
        d = {}
        for k1, v1 in self.d.items():
            for k2, v2 in o.d.items():
                if k1 & k2:
                    continue
                k = k1 | k2
                if len(k) > Ctx.maxdeg:
                    continue
                d[k] = d.get(k, ZERO) + v1 * v2
        return Fld({k: norm(v) for k, v in d.items()})

    __rmul__ = __mul__

    def d_(self, i):
        return Fld({k: norm(IU * Ctx.mom(k)[i] * v) for k, v in self.d.items()})

    def dt(self):
        return Fld({k: v.diff(T) for k, v in self.d.items()})

    def map(self, f):
        return Fld({k: f(v) for k, v in self.d.items()})

    def coeff(self, labels):
        return self.d.get(frozenset(labels), ZERO)

    def without_const(self):
        return Fld({k: v for k, v in self.d.items() if k})


def fexp(F):
    F = F.without_const()
    out, term = Fld.const(1), Fld.const(1)
    for n in range(1, Ctx.maxdeg + 1):
        term = term * F * rat(1, n)
        out = out + term
    return out


def mzero():
    return [[Fld() for _ in range(3)] for _ in range(3)]


def madd(A, B, s=1):
    return [[A[i][j] + B[i][j] * s for j in range(3)] for i in range(3)]


def mmul(A, B):
    return [[sum((A[i][k] * B[k][j] for k in range(3)), Fld()) for j in range(3)] for i in range(3)]


def mscale(A, c):
    return [[A[i][j] * c for j in range(3)] for i in range(3)]


def mexp(H):
    out = [[Fld.const(1 if i == j else 0) for j in range(3)] for i in range(3)]
    term = [[Fld.const(1 if i == j else 0) for j in range(3)] for i in range(3)]
    for n in range(1, Ctx.maxdeg + 1):
        term = mscale(mmul(term, H), rat(1, n))
        out = madd(out, term)
    return out


def geometry(zeta, H):
    """gamma, S = gamma^{-1} sigma / ... and Christoffels: the parts of engine.equations the map uses."""
    e2z, em2z = fexp(zeta * 2), fexp(zeta * (-2))
    M, Minv = mexp(H), mexp(mscale(H, -1))
    gam, gi = mscale(M, e2z), mscale(Minv, em2z)
    dgam = [[[gam[i][j].d_(k) for k in range(3)] for j in range(3)] for i in range(3)]
    low = [[[(dgam[b][j][i] + dgam[b][i][j] - dgam[i][j][b]) * rat(1, 2) for j in range(3)] for i in range(3)] for b in range(3)]
    Gam = [[[sum((gi[a][b] * low[b][i][j] for b in range(3)), Fld()) for j in range(3)] for i in range(3)] for a in range(3)]
    S = mscale(mmul(Minv, [[M[i][j].dt() for j in range(3)] for i in range(3)]), rat(1, 2))
    return gam, S, Gam


# ------------------------------------------------------------------ polynomial operations
def subs(poly, idx, val):
    return norm(poly.compose(R.gens[idx], val))


def integ(poly, idx, lo, hi):
    """Exact definite integral in generator idx from lo to hi (ring elements)."""
    by = {}
    for mon, c in poly.items():
        m = mon[idx]
        mon0 = mon[:idx] + (0,) + mon[idx + 1 :]
        by.setdefault(m, {})[mon0] = c
    out = ZERO
    for m, dct in by.items():
        out = out + R(dct) * (hi ** (m + 1) - lo ** (m + 1)) * rat(1, m + 1)
    return norm(out)


def sky_mean(poly):
    df = lambda m: sp.factorial2(m - 1) if m > 0 else 1
    out = {}
    for mon, c in poly.items():
        a, b, cc = (mon[i] for i in _IN)
        if a % 2 or b % 2 or cc % 2:
            continue
        w = QQ_I.convert(sp.Rational(df(a) * df(b) * df(cc), sp.factorial2(a + b + cc + 1)))
        m = list(mon)
        for i in _IN:
            m[i] = 0
        m = tuple(m)
        out[m] = out.get(m, QQ_I.zero) + c * w
    return R({m: c for m, c in out.items() if c})


def phase(S):
    K = Ctx.mom(S)
    arg = IU * sum((K[i] * NV[i] for i in range(3)), ZERO) * (TO - T)
    out, term = ONE, ONE
    for m in range(1, KMAX + 1):
        term = norm(term * arg) * rat(1, m)
        out = out + term
    return norm(out)


def on_ray(F):
    return Fld({S: norm(v * phase(S)) for S, v in F.d.items()})


def taylor_on_ray(F, dX):
    out = on_ray(F)
    terms = [(F, Fld.const(1))]
    for n in range(1, Ctx.maxdeg + 1):
        new = []
        for G, w in terms:
            for i in range(3):
                new.append((G.d_(i), w * dX[i]))
        terms = new
        acc = Fld()
        for G, w in terms:
            acc = acc + w * on_ray(G)
        out = out + acc * rat(1, sp.factorial(n))
    return out


def shift_eval(Fr, delta, tbar):
    """f(tbar + delta) for a ray/time function f(t) (Fld of time functions); delta has degree >= 1."""
    out = Fr.map(lambda v: subs(v, _IT, tbar))
    dn = Fld.const(1)
    G = Fr
    for n in range(1, Ctx.maxdeg + 1):
        dn = dn * delta
        G = G.dt()
        out = out + dn * G.map(lambda v: subs(v, _IT, tbar)) * rat(1, sp.factorial(n))
    return out


def solve_root(phi, tbar, itbar):
    """ln(1 + delta/tbar) + phi(tbar + delta) = 0 for delta (degree >= 1), multilinear."""
    delta = Fld()
    for _ in range(Ctx.maxdeg + 1):
        u = delta * itbar
        log_rest = Fld()
        un = Fld.const(1)
        for m in range(1, Ctx.maxdeg + 1):
            un = un * u
            if m >= 2:
                log_rest = log_rest + un * rat((-1) ** (m + 1), m)
        delta = (shift_eval(phi, delta, tbar) + log_rest) * (-tbar)
    return delta


def observable(nu, zeta, H):
    gam, S, Gam = geometry(zeta, H)
    LOG("geometry (gamma, S, Christoffels)")
    sig = [[gam[i][j].dt() * rat(1, 2) for j in range(3)] for i in range(3)]
    Ginvhalf = mscale(mexp(mscale(H, rat(-1, 2))), fexp(zeta * (-1)) if MUTATE != "noezeta" else Fld.const(1))
    G0 = [[Ginvhalf[i][j].map(lambda v: subs(v, _IT, TO)) for j in range(3)] for i in range(3)]
    n = [Fld.const(c) for c in NV]
    e0 = [sum((G0[i][j] * (-NV[j]) for j in range(3)), Fld()) + NV[i] for i in range(3)]
    e0 = [e0[i].without_const() for i in range(3)]
    dX = [e0[i] * (T - TO) for i in range(3)]
    for it in range(Ctx.maxdeg):
        e = [n[i] * (-1) + dX[i].dt() for i in range(3)]
        Sr = [[taylor_on_ray(S[i][j], dX) for j in range(3)] for i in range(3)]
        Gr = [[[taylor_on_ray(Gam[a][i][j], dX) for j in range(3)] for i in range(3)] for a in range(3)]
        sr = [[taylor_on_ray(sig[i][j], dX) for j in range(3)] for i in range(3)]
        see = sum((sr[i][j] * e[i] * e[j] for i in range(3) for j in range(3)), Fld())
        A = []
        for i in range(3):
            Ai = sum((Sr[i][j] * e[j] for j in range(3)), Fld()) * (-2)
            Ai = Ai - sum((Gr[i][j][k] * e[j] * e[k] for j in range(3) for k in range(3)), Fld())
            Ai = Ai + e[i] * see
            A.append(Ai)
        if MUTATE == "noray":
            A = [Fld() for _ in range(3)]
        dX = []
        for i in range(3):
            d = {}
            for Sx in set(e0[i].d) | set(A[i].d):
                val = norm(e0[i].d.get(Sx, ZERO) * (T - TO))
                if Sx in A[i].d:
                    Atp = subs(A[i].d[Sx], _IT, TP)
                    val = val + integ(norm((TP - T) * Atp), _ITP, T, TO)
                d[Sx] = norm(val)
            dX.append(Fld(d))
        LOG(f"ray iteration {it + 1}")
    e = [n[i] * (-1) + dX[i].dt() for i in range(3)]
    sr = [[taylor_on_ray(sig[i][j], dX) for j in range(3)] for i in range(3)]
    F = sum((sr[i][j] * e[i] * e[j] for i in range(3) for j in range(3)), Fld())
    phi = taylor_on_ray(nu, dX)
    delta = solve_root(phi, TS, ITS)
    Iexp = Fld({Sx: integ(v, _IT, TS, TO) for Sx, v in F.d.items()})
    corr = Fld()
    dn = Fld.const(1)
    G = F
    for m in range(1, Ctx.maxdeg + 1):
        dn = dn * delta
        corr = corr + dn * G.map(lambda v: subs(v, _IT, TS)) * rat(1, sp.factorial(m))
        G = G.dt()
    if MUTATE != "noend":
        Iexp = Iexp - corr
    LOG("exponent I built (ray, emission root, endpoint)")
    Tn = fexp(Iexp * (-1))
    m1 = (Tn.map(sky_mean) - 1).without_const()
    inv, term = Fld.const(1), Fld.const(1)
    for _ in range(Ctx.maxdeg):
        term = term * m1 * (-1)
        inv = inv + term
    Thc = (Tn * inv - 1).without_const() if MUTATE != "nonorm" else (Iexp * (-1)).without_const()
    omega = Fld(dict(nu.d))  # nu(t, 0): no phase at the origin
    tI_shift = solve_root(omega, TI, ITI)
    eom = fexp(omega)
    Wint = Fld({Sx: integ(norm(T * v), _IT, ZERO, TP) for Sx, v in (eom - 1).without_const().d.items()})
    Wf = Wint.map(lambda v: subs(v, _ITP, T))
    D = Fld()
    for _ in range(Ctx.maxdeg + 1):
        Rr = tI_shift * (TI * 2) + tI_shift * tI_shift - shift_eval(Wf, D, TOB) * 2 + shift_eval(Wf, tI_shift, TI) * 2
        D = (Rr - D * D) * (ITOB * rat(1, 2))
    Th = Thc.map(lambda v: subs(v, _ITO, TOB))
    dn = Fld.const(1)
    G = Thc
    for m in range(1, Ctx.maxdeg + 1):
        dn = dn * D
        G = G.map(lambda v: v.diff(TO))
        Th = Th + dn * G.map(lambda v: subs(v, _ITO, TOB)) * rat(1, sp.factorial(m))
    LOG("normalised, clock composed")
    return Th


# ------------------------------------------------------------------ fields from the fast engine
def build_fields(labels, zero_q=False):
    """Fields for labels at k = (p,0,0), q = (b1,b2,0), -q, from fast_series (exact)."""
    second = {}
    pairs = [pr for pr in ((1, 2), (1, 3), (2, 3)) if set(pr) <= set(labels)]
    for pr in pairs:
        second[pr] = fs.solve(list(pr), {})
        LOG(f"fast engine: pair {pr} solved")
    top = None
    if len(labels) == 3:
        top = fs.solve(list(labels), second)
        LOG("fast engine: third order solved")
    nu_s, zeta_s, H_s = fs.build(labels, second, top)
    base = {fs.P: P, fs.B1: B1, fs.B2: B2, fs.IP: IP, fs.U: U, fs.V: V}
    del base

    def conv(series):
        out = ZERO
        for k, c in series.d.items():
            assert k >= 0, "negative power of t in a field kernel"
            if k > KMAX:
                continue
            out = out + c.set_ring(R) * KAP**k * T**k
        if zero_q:
            out = out.compose([(B1, ZERO), (B2, ZERO)])
            out = out.compose([(U, IP**2), (V, IP**2)])
        return norm(out)

    tofld = lambda Fs: Fld({k: conv(v) for k, v in Fs.d.items()})
    nu, zeta = tofld(nu_s), tofld(zeta_s)
    H = [[tofld(H_s[i][j]) for j in range(3)] for i in range(3)]
    q = (ZERO, ZERO, ZERO) if zero_q else (KAP * B1, KAP * B2, ZERO)
    Ctx.momenta = {1: (KAP * P, ZERO, ZERO), 2: q, 3: tuple(-c for c in q)}
    Ctx.maxdeg = len(labels)
    return nu, zeta, H


def to_sympy(poly):
    """Ring element -> SymPy expression with the formal inverses restored."""
    p, q, b1, b2 = sp.symbols("p q", positive=True) + sp.symbols("b1 b2", real=True)
    kap, t_, to, ts, tI, tob = sp.symbols("kappa t t_o tbar_star tbar_I tbar_o", positive=True)
    n1, n2, n3 = sp.symbols("n1 n2 n3", real=True)
    sym = dict(
        kap=kap, t=t_, to=to, ts=ts, tI=tI, tob=tob, its=1 / ts, itI=1 / tI, itob=1 / tob,
        n1=n1, n2=n2, n3=n3, p=p, b1=b1, b2=b2, ip=1 / p,
        u=1 / ((p + b1) ** 2 + b2**2), v=1 / ((p - b1) ** 2 + b2**2), tp=sp.Symbol("tp"),
    )  # fmt: skip
    e = poly.as_expr()
    return sp.expand(e.xreplace({s_: sym[s_.name] for s_ in e.free_symbols})), sym


# ------------------------------------------------------------------ checks
def superhorizon_prediction(labels, sym):
    """Multilinear coefficient of -tau^2 e^{-2 zeta} TF[zeta_ij - zeta_i zeta_j] n n at the origin."""
    kap, p, b1, b2 = sym["kap"], sym["p"], sym["b1"], sym["b2"]
    ks = {1: sp.Matrix([kap * p, 0, 0]), 2: sp.Matrix([kap * b1, kap * b2, 0]), 3: sp.Matrix([-kap * b1, -kap * b2, 0])}
    zs = sp.symbols("z1:4")
    Z = sum(zs[a - 1] for a in labels)
    Zi = sum((zs[a - 1] * sp.I * ks[a] for a in labels), sp.zeros(3, 1))
    Zij = sum((zs[a - 1] * (-ks[a] * ks[a].T) for a in labels), sp.zeros(3, 3))
    Mm = sp.exp(-2 * Z) * (Zij - Zi * Zi.T)
    nn = sp.Matrix([sym["n1"], sym["n2"], sym["n3"]])
    val = (nn.T * Mm * nn)[0, 0] - Mm.trace() * (nn.T * nn)[0, 0] / 3
    tau2 = (sym["tob"] ** 2 - sym["ts"] ** 2) / 6
    return sp.expand(sp.diff(-tau2 * val, *[zs[a - 1] for a in labels]).subs({z_: 0 for z_ in zs}))


def reduce_q(e, sym):
    """Exact zero test helper: clear the formal denominators and impose nothing else (q^2 = b1^2+b2^2 is implicit)."""
    return sp.simplify(sp.together(e))


def check_superhorizon(labels):
    nu, zeta, H = build_fields(labels)
    LOG("fields ready")
    Th = observable(nu, zeta, H)
    got, sym = to_sympy(Th.coeff(labels))
    kap = sym["kap"]
    pred = superhorizon_prediction(labels, sym)
    low = [sp.simplify(got.coeff(kap, m)) for m in (0, 1)]
    r = reduce_q(sp.expand(got - pred).coeff(kap, 2), sym)
    assert low == [0, 0] and r == 0, (low, r)
    if len(labels) == 2:
        got1, _ = to_sympy(Th.coeff([1]))
        pred1 = superhorizon_prediction([1], sym)
        assert sp.simplify(sp.expand(got1 - pred1).coeff(kap, 2)) == 0
    return Th


def check_ward(order):
    labels = [1, 2] if order == 2 else [1, 2, 3]
    nu, zeta, H = build_fields(labels, zero_q=True)
    LOG("fields ready")
    Th = observable(nu, zeta, H)
    base, sym = to_sympy(Th.coeff([1]))
    kap = sym["kap"]
    Dk = lambda e: sp.expand(kap * sp.diff(e, kap))
    pred = -Dk(base) if order == 2 else Dk(Dk(base))
    got, _ = to_sympy(Th.coeff(labels))
    assert base != 0
    assert sp.simplify(sp.expand(got - pred)) == 0
    return Th


if __name__ == "__main__":
    if MODE == "SH":
        check_superhorizon([1, 2])
        print("[ok] SH (fast): superhorizon closed form at O(kappa^2), orders 1 and 2, general (k, q)")
    elif MODE == "SH3":
        check_superhorizon([1, 2, 3])
        print("[ok] SH3 (fast): superhorizon closed form at O(kappa^2), third order at (k, q, -q)")
    elif MODE in ("W2", "W3"):
        check_ward(int(MODE[1]))
        print(f"[ok] {MODE} (fast): observable dilation Ward identity through O(kappa^{KMAX})")
    else:
        raise SystemExit(f"unknown mode {MODE}")
