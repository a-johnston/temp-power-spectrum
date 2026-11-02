"""Closed-form time integrals in ENTIRE special functions (resonance-safe).

Every second-order field kernel is a finite sum of integrals
    P(t) = int_0^t B(s) ds,   B(s) = sum_w e^{i w s} R_w(s),
with R_w Laurent polynomials in s and B analytic at s = 0.  The frequencies w are integer
combinations of p/sqrt3, q/sqrt3, K/sqrt3 and K, so they vanish on resonance surfaces (collinear
momenta, the sound cone).  Term-by-term integration is done with
    s^m, m >= 0:   int_0^t s^m e^{i w s} ds = t^{m+1} Ephi(m, i w t),
                   Ephi(m, z) = int_0^1 u^m e^{z u} du                 (entire in z)
    s^{-M}, M >= 1: int_0^t s^{-M} [e^{i w s} - sum_{j<M} (i w s)^j / j!] ds = J_M(w, t),
                   J_1 = -Cin(w t) + i Si(w t),   Cin(x) = int_0^x (1 - cos u)/u du   (entire, even)
                   J_M = t^{1-M} E_M /(1-M) + (i w/(M-1)) J_{M-1},  E_M = e^{i w t} - sum_{j<M} (i w t)^j/j!
The subtracted Taylor pieces are proportional to s^{j-M} (j < M); their total over all frequencies
is the principal part of B at s = 0, which must vanish (asserted).  No division by w ever occurs,
so the result is manifestly regular on every resonance surface: there are no ln|w| terms to
recombine.
"""

import sympy as sp

I = sp.I
s3 = sp.sqrt(3)


class Cin(sp.Function):
    """Cin(x) = int_0^x (1 - cos u)/u du  (entire, even)."""

    @classmethod
    def eval(cls, x):
        if x.is_zero:
            return sp.Integer(0)
        if x.could_extract_minus_sign():
            return cls(-x)

    def fdiff(self, argindex=1):
        x = self.args[0]
        return (1 - sp.cos(x)) / x


class Ephi(sp.Function):
    """Ephi(m, z) = int_0^1 u^m e^{z u} du  (entire in z)."""

    nargs = 2

    @classmethod
    def eval(cls, m, z):
        if z.is_zero:
            return 1 / (m + 1)

    def fdiff(self, argindex=2):
        m, z = self.args
        return Ephi(m + 1, z)

    def elementary(self):
        m, z = self.args
        m = int(sp.Integer(m))
        return (
            sp.factorial(m)
            * (-1) ** (m + 1)
            / z ** (m + 1)
            * (1 - sp.exp(z) * sum((-z) ** j / sp.factorial(j) for j in range(m + 1)))
        )


def J(M, w, t):
    if M == 1:
        return -Cin(w * t) + I * sp.Si(w * t)
    EM = sp.exp(I * w * t) - sum((I * w * t) ** j / sp.factorial(j) for j in range(M))
    return t ** (1 - M) * EM / (1 - M) + I * w / (M - 1) * J(M - 1, w, t)


class Freq:
    """Map exp(i * t * (a p + b q + c K)/sqrt3 + i d K t) <-> integer tuple (a, b, c, d)."""

    def __init__(self, syms, t):
        self.syms = syms  # (p, q, K)
        self.t = t

    def key(self, arg):
        p, q, K = self.syms
        lin = sp.expand(arg / (I * self.t))
        assert not lin.has(self.t), arg
        ca = sp.nsimplify(sp.expand(lin.coeff(p) * s3))
        cb = sp.nsimplify(sp.expand(lin.coeff(q) * s3))
        ck = sp.expand(lin.coeff(K) * s3)  # = c + d sqrt3
        d = sp.nsimplify(ck.coeff(s3))
        c = sp.nsimplify(sp.expand(ck - d * s3))
        rest = sp.expand(lin - (ca * p + cb * q) / s3 - (c / s3 + d) * K)
        assert rest == 0, (arg, rest)
        for v in (ca, cb, c, d):
            assert v.is_integer, (arg, v)
        return (int(ca), int(cb), int(c), int(d))

    def omega(self, key):
        p, q, K = self.syms
        a, b, c, d = key
        return (a * p + b * q + c * K) / s3 + d * K


def split_frequencies(expr, fr, s):
    """expr(s) (trig in s, rational in s) -> {key: Laurent coefficient dict {m: coeff}}."""
    e = sp.expand(sp.expand(expr).rewrite(sp.exp))
    out = {}
    for term in sp.Add.make_args(e):
        num, den = sp.fraction(term)
        if den.has(sp.exp):
            den = sp.factor(den)
            dfac = []
            for f in sp.Mul.make_args(den):
                base, ex = f.as_base_exp()
                if isinstance(base, sp.exp):
                    num = num * sp.exp(-base.args[0] * ex)
                elif base == sp.E:
                    num = num * sp.exp(-ex)
                else:
                    dfac.append(f)
            den = sp.Mul(*dfac)
            assert not den.has(sp.exp), term
            term = sp.expand(num) / den
            if isinstance(term, sp.Add) or sp.expand(num).is_Add:
                # re-split after normalisation
                sub = split_frequencies(sp.expand(sp.expand(num) / den), fr, s)
                for k_, lau in sub.items():
                    dd = out.setdefault(k_, {})
                    for m_, c_ in lau.items():
                        dd[m_] = dd.get(m_, 0) + c_
                continue
        fac = sp.Mul.make_args(term)
        coef = 1
        kk = [0, 0, 0, 0]
        for f in fac:
            base, ex = f.as_base_exp()
            if isinstance(base, sp.exp) or (base == sp.E):
                arg = base.args[0] * ex if isinstance(base, sp.exp) else ex
                if arg.has(s):
                    k_ = fr.key(arg.subs(s, fr.t))
                    kk = [kk[i] + k_[i] for i in range(4)]
                    continue
            coef *= f
        key = tuple(kk)
        # coef = (polynomial in s) / (s^m * s-independent factor)
        n_, d_ = sp.fraction(sp.together(coef))
        d_ = sp.factor(d_)
        ind, dep = d_.as_independent(s, as_Add=False)
        mden = 0
        if dep != 1:
            b_, e_ = dep.as_base_exp()
            assert b_ == s and e_.is_integer, ("non-monomial s-denominator", dep)
            mden = int(e_)
        cpow = sp.Poly(sp.expand(n_), s)
        for (m,), c in cpow.terms():
            dd = out.setdefault(key, {})
            dd[m - mden] = dd.get(m - mden, 0) + c / ind
    return out


def integrate0(expr, s, t, syms, reduce=lambda e: sp.cancel(e)):
    """Exact P(t) = int_0^t expr(s) ds in entire functions; asserts regularity at s = 0."""
    fr = Freq(syms, s)
    parts = split_frequencies(expr, fr, s)
    res = 0
    principal = {}
    for key, lau in parts.items():
        w = fr.omega(key)
        for m, c in lau.items():
            c = reduce(c)
            if c == 0:
                continue
            if m >= 0:
                res += c * t ** (m + 1) * Ephi(m, I * w * t)
            else:
                M = -m
                res += c * J(M, w, t)
                for j in range(M):
                    n = j - M
                    principal[n] = principal.get(n, 0) + c * (I * w) ** j / sp.factorial(j)
    for n, c in principal.items():
        assert reduce(sp.expand(c)) == 0, f"principal part s^{n} of the integrand does not vanish"
    return res


def to_elementary(e):
    return e.replace(lambda x: isinstance(x, Ephi), lambda x: x.elementary())


def check_antiderivative(P, expr, s, t, syms, reduce):
    """Exact check  dP/dt = expr(t): Si, Cin as independent transcendentals, exp by frequency."""
    d = sp.diff(P, t)
    d = to_elementary(d) - expr.subs(s, t)
    d = sp.expand(sp.expand(d).rewrite(sp.exp))
    fr = Freq(syms, t)
    parts = split_frequencies(d, fr, t)
    for key, lau in parts.items():
        for m, c in lau.items():
            assert reduce(c) == 0, ("antiderivative check failed", key, m)
    return True
