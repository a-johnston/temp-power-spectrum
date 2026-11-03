"""Independent forward cubic temperature map at finite wavelength (VALIDATION ONLY, not evidence).

Everything is solved from definitions, for an arbitrary momentum triple, with no adjoint,
stationary, second-ray-elimination or constraint-row formula:

  fields    regular solutions of (N) and the Cartesian rows (T) through third order (cubic_fields)
  frame     Fermi-Walker triad dE/dt = -S(t,0) E, E(tbar_I) = gamma^(-1/2), e(t_o) = -E(t_o) n
  ray       e' = -2 S e - Gamma(e,e) + e sigma(e,e),  X' = e,  X(t_o) = 0  (Picard, multilinear)
  emission  t_e exp(nu(t_e, X(t_e))) = tbar_*
  exponent  I = int_{t_e}^{t_o} sigma(e,e) dt,  Theta = e^{-I}/<e^{-I}> - 1
  clock     t_I exp(nu(t_I,0)) = tbar_I,  int_{t_I}^{t_o} t e^{nu(t,0)} dt = (tbar_o^2 - tbar_I^2)/2

Multilinear convention: a quantity is a dict {label set: array}; products only between disjoint
label sets, so the {1,2,3} coefficient is the polarised cubic kernel with no factorials.
Time integrals along the ray are spectral (Chebyshev-Lobatto), sky averages Gauss-Legendre x
uniform azimuth, clock composition a five-point stencil in the observation time.
"""

from __future__ import annotations

import numpy as np

EMPTY = frozenset()


class F(dict):
    def __add__(self, o):
        r = F(self)
        for k, v in o.items():
            r[k] = r[k] + v if k in r else v
        return r

    def __sub__(self, o):
        return self + o.scale(-1.0)

    def scale(self, c):
        return F({k: c * v for k, v in self.items()})

    def drop0(self):
        return F({k: v for k, v in self.items() if k != EMPTY})

    def map(self, fn):
        return F({k: fn(v) for k, v in self.items()})


def mul(A, B, op=lambda a, b: a * b, maxlab=3):
    r = F()
    for ka, va in A.items():
        for kb, vb in B.items():
            if ka & kb:
                continue
            k = ka | kb
            if len(k) > maxlab:
                continue
            v = op(va, vb)
            r[k] = r[k] + v if k in r else v
    return r


def const(x):
    return F({EMPTY: np.asarray(x, dtype=complex)})


def fexp(X, order=3):
    out, term = const(1.0), const(1.0)
    X = X.drop0()
    for n in range(1, order + 1):
        term = mul(term, X).scale(1.0 / n)
        out = out + term
    return out


mm = lambda a, b: a @ b  # noqa: E731
mv = lambda a, b: np.einsum("...ij,...j->...i", a, b)  # noqa: E731
vv = lambda a, b: np.einsum("...i,...i->...", a, b)  # noqa: E731
smat = lambda a, b: a[..., None, None] * b  # noqa: E731
svec = lambda a, b: a[..., None] * b  # noqa: E731


def mexp(H, sign=1.0):
    Hs = H.drop0().scale(sign)
    eye = const(np.eye(3))
    H2 = mul(Hs, Hs, mm)
    return eye + Hs + H2.scale(0.5) + mul(H2, Hs, mm).scale(1 / 6)


def mexp_dot(H, Hd):
    """d/dt exp(H) along Hd, through third order."""
    H, Hd = H.drop0(), Hd.drop0()
    return (
        Hd
        + (mul(Hd, H, mm) + mul(H, Hd, mm)).scale(0.5)
        + (mul(mul(Hd, H, mm), H, mm) + mul(mul(H, Hd, mm), H, mm) + mul(mul(H, H, mm), Hd, mm)).scale(1 / 6)
    )


# ----------------------------------------------------------------------------- spectral tools
class Spectral:
    """Chebyshev-Lobatto nodes on [lo, hi] with node 0 at hi; cumulative integral from hi."""

    def __init__(self, lo, hi, n):
        x = np.cos(np.pi * np.arange(n) / (n - 1))
        self.x, self.t = x, lo + (hi - lo) * (1 + x) / 2
        V = np.polynomial.chebyshev.chebvander(x, n - 1)
        Vi = np.linalg.inv(V)
        half = (hi - lo) / 2
        I = np.eye(n)
        Cint = np.stack([np.polynomial.chebyshev.chebint(I[:, j], lbnd=1.0) for j in range(n)], axis=1)  # (n+1, n)
        Vint = np.polynomial.chebyshev.chebvander(x, n)
        self.Q = half * Vint @ Cint @ Vi  # (Q f)_j = int_hi^{t_j} f
        Cd = np.stack([np.concatenate([np.polynomial.chebyshev.chebder(I[:, j]), [0.0]]) for j in range(n)], axis=1)
        self.D = V @ Cd @ Vi / half
        self.full = -(self.Q[-1])  # weights for int_lo^hi f

    def cumint(self, f, axis):
        f = np.moveaxis(f, axis, -1)
        return np.moveaxis(f @ self.Q.T, -1, axis)

    def deriv(self, f, axis):
        f = np.moveaxis(f, axis, -1)
        return np.moveaxis(f @ self.D.T, -1, axis)


def sphere(nth, nph):
    mu, wmu = np.polynomial.legendre.leggauss(nth)
    ph = 2 * np.pi * (np.arange(nph) + 0.5) / nph
    MU, PH = np.meshgrid(mu, ph, indexing="ij")
    st = np.sqrt(1 - MU**2)
    n = np.stack([st * np.cos(PH), st * np.sin(PH), MU], axis=-1).reshape(-1, 3)
    w = (wmu[:, None] * np.full(nph, 1.0 / (2 * nph))[None, :]).reshape(-1)  # weights of the average
    return n, w


# ----------------------------------------------------------------------------- the forward map
class ForwardMap:
    def __init__(self, momenta, profiles, e, o, iota, labels_on=None, nth=18, nph=36, nt=44, rotate=True):
        """profiles: {label set: cubic_fields.Profile}; labels_on: which label sets carry fields.

        rotate=False replaces the Fermi-Walker triad by the symmetric one (negative control)."""
        self.rotate = rotate
        from itertools import combinations

        labs = sorted(momenta)
        self.K = {
            frozenset(c): sum((np.asarray(momenta[a], dtype=float) for a in c), np.zeros(3))
            for r in range(1, len(labs) + 1)
            for c in combinations(labs, r)
        }
        self.prof = {S: p for S, p in profiles.items() if labels_on is None or S in labels_on}
        self.e, self.o, self.iota = e, o, iota
        self.n, self.w = sphere(nth, nph)
        self.nt = nt
        self.single = sorted(S for S in profiles if len(S) == 1)

    # time profiles (no spatial phase) at times tt: nu, H, H', zeta
    def profiles_at(self, tt):
        nu, H, Hp, z = F(), F(), F(), F()
        for S, p in self.prof.items():
            if p.kind == "first":
                from cubic_fields import first_profiles

                f = first_profiles(p.kmag, tt)
                nu[S] = f["a"][0]
                H[S] = f["b"][0][..., None, None] * p.Qhat
                Hp[S] = f["b"][1][..., None, None] * p.Qhat
                z[S] = np.ones_like(tt)
            else:
                nu[S] = p.nu(tt)
                H[S] = p.A(tt)
                Hp[S] = p.Ap(tt)
        return nu, H, Hp, z

    def on_ray(self, P, dX, lam, deriv=False):
        """Fields with time profiles P (arrays (Nt,...)) at X = n lam + dX; optional gradient index."""
        out = F()
        y_cache = {}
        for S, v in P.items():
            KS = self.K[S] if S else np.zeros(3)
            if S == EMPTY:
                if not deriv:  # time-independent constant (identity, unit scalar)
                    out = out + F({EMPTY: np.broadcast_to(v, (len(self.n), len(lam)) + v.shape)})
                continue
            ph = np.exp(1j * (self.n @ KS)[:, None] * lam[None, :])
            base = v[None] * ph.reshape(ph.shape + (1,) * (v.ndim - 1))
            if deriv:
                base = base[..., None] * (1j * KS)
            key = tuple(KS)
            if key not in y_cache:
                y = F({kk: 1j * np.einsum("...i,i->...", d, KS) for kk, d in dX.items() if kk != EMPTY})
                y_cache[key] = fexp(y)
            out = out + mul(F({S: base}), y_cache[key], lambda a, b: a * b.reshape(b.shape + (1,) * (a.ndim - b.ndim)))
        return out

    def frame(self, t_o, nf=40):
        """Fermi-Walker triad at the observer, E(t_o), multilinear through second order."""
        sp_ = Spectral(self.iota, t_o, nf)  # node 0 at t_o, last node at iota
        nu, H, Hp, z = self.profiles_at(sp_.t)
        Mi = mexp(H, -1.0)
        S = mul(Mi, mexp_dot(H, Hp), mm).scale(0.5)
        # E(iota) = e^{-zeta} e^{-H(iota)/2}
        HI = F({k: v[-1] for k, v in H.items()})
        zI = F({k: v[-1] for k, v in z.items()})
        EI = mul(fexp(zI.scale(-1.0)), mexp(HI, -0.5), smat)
        E = F({k: np.broadcast_to(v, (nf, 3, 3)).copy() for k, v in EI.items()})
        for _ in range(3):
            SE = mul(S, E, mm, maxlab=2).drop0()  # third-order frame acts on the background sky only
            # E(t) = E(iota) - int_iota^t S E = E(iota) - [int_{t_o}^t - int_{t_o}^iota] S E
            newE = F({k: np.broadcast_to(v, (nf, 3, 3)).copy() for k, v in EI.items()})
            for k, v in SE.items():
                cum = sp_.cumint(v, 0)
                newE[k] = newE.get(k, 0) - (cum - cum[-1])
            E = newE
        return F({k: v[0] for k, v in E.items()})

    def exponent(self, t_o):
        """I(n) multilinear at coordinate observation time t_o (before normalisation)."""
        sp_ = Spectral(self.e, t_o, self.nt)
        tr = sp_.t
        lam = t_o - tr
        Ns = len(self.n)
        nu, H, Hp, z = self.profiles_at(tr)
        H0, Hp0 = H.drop0(), Hp.drop0()
        M, Mi = mexp(H0, 1.0), mexp(H0, -1.0)
        Md = mexp_dot(H0, Hp0)
        e2z = fexp(z.scale(2.0))
        grad = lambda P: F({k: v[..., None] * (1j * self.K[k]) for k, v in P.items() if k != EMPTY})  # noqa: E731
        dM, dz = grad(M), grad(z)
        Sprof = mul(Mi, Md, mm).scale(0.5)
        sig = mul(e2z, Md, smat).scale(0.5)
        P1 = mul(Mi, dM, lambda a, b: np.einsum("...ia,...ajl->...ijl", a, b))
        P2 = mul(Mi, dM, lambda a, b: np.einsum("...il,...jkl->...ijk", a, b))
        Q3 = mul(Mi, dz, mv)
        if self.rotate:
            E = self.frame(t_o)
        else:
            Ho = F({k: v[0] for k, v in H.items()})
            zo = F({k: v[0] for k, v in z.items()})
            E = mul(fexp(zo.scale(-1.0)), mexp(Ho, -0.5), smat)
        e_o = F({k: -np.einsum("ij,sj->si", v, self.n) for k, v in E.items()})
        nvec = np.broadcast_to(self.n[:, None, :], (Ns, len(tr), 3))
        de_o = e_o.drop0()
        dX, de = F(), F()
        for _ in range(3):
            ray = lambda P, dX=dX: self.on_ray(P, dX, lam)  # noqa: E731
            e = const(-nvec) + de
            Se = mul(ray(Sprof), e, mv)
            see = mul(mul(ray(sig), e, mv), e, vv)
            t1a = mul(mul(ray(P1), e, lambda a, b: np.einsum("...ijl,...l->...ij", a, b)), e, mv)
            t1b = mul(mul(ray(P2), e, lambda a, b: np.einsum("...ijk,...k->...ij", a, b)), e, mv)
            edz = mul(e, ray(dz), vv)
            eMe = mul(mul(e, ray(M), lambda a, b: np.einsum("...i,...ij->...j", a, b)), e, vv)
            Gam = (
                t1a
                - t1b.scale(0.5)
                + mul(e, edz, lambda a, b: a * b[..., None]).scale(2.0)
                - mul(eMe, ray(Q3), lambda a, b: a[..., None] * b)
            )
            acc = (Se.scale(-2.0) - Gam + mul(e, see, lambda a, b: a * b[..., None])).drop0()
            new_de = F()
            for k, a in acc.items():
                new_de[k] = sp_.cumint(a, 1)
            for k, v in de_o.items():
                new_de[k] = new_de.get(k, 0) + np.broadcast_to(v[:, None, :], (Ns, len(tr), 3))
            new_dX = F({k: sp_.cumint(v, 1) for k, v in new_de.items()})
            de, dX = new_de, new_dX
        Fint = see.drop0()
        Itot = F({k: np.einsum("st,t->s", v, sp_.full) for k, v in Fint.items()})
        phi = self.on_ray(nu.drop0(), dX, lam).drop0()
        phi0 = F({k: v[:, -1] for k, v in phi.items()})
        phid = F({k: sp_.deriv(v, 1)[:, -1] for k, v in phi.items()})
        delta = F()
        for _ in range(3):
            delta = (fexp((phi0 + mul(phid, delta)).scale(-1.0)) - const(1.0)).drop0().scale(self.e)
        F0 = F({k: v[:, -1] for k, v in Fint.items()})
        F1 = F({k: sp_.deriv(v, 1)[:, -1] for k, v in Fint.items()})
        Iend = mul(F0, delta) + mul(mul(F1, delta), delta).scale(0.5)
        return Itot - Iend

    def _nu_and_derivative(self, T_):
        from cubic_fields import first_profiles

        nu, nup = F(), F()
        for k, p in self.prof.items():
            if p.kind == "first":
                f = first_profiles(p.kmag, np.array([T_]))
                nu[k], nup[k] = f["a"][0][0], f["a"][1][0]
            else:
                nu[k], nup[k] = p.nu(np.array([T_]))[0], p.nup(np.array([T_]))[0]
        return nu, nup

    def clock(self):
        """D = t_o - tbar_o (orders 1 and 2), from the physical initial event and proper time."""
        sp_ = Spectral(self.iota, self.o, 60)
        nu, *_ = self.profiles_at(sp_.t)
        g = (fexp(nu) - const(1.0)).drop0()
        Gdiff = F({k: np.dot(sp_.full, sp_.t * v) for k, v in g.items()})  # int_iota^o t g
        nuI, nuIp = self._nu_and_derivative(self.iota)
        nuo, _ = self._nu_and_derivative(self.o)
        dI = F()
        for _ in range(3):
            dI = (fexp((nuI + mul(nuIp, dI)).scale(-1.0)) - const(1.0)).drop0().scale(self.iota)
        gpo = (fexp(nuo) - const(1.0)).drop0().scale(self.o)
        gpI = (fexp(nuI) - const(1.0)).drop0().scale(self.iota)
        D = F()
        for _ in range(3):
            rhs = mul(D, D).scale(-0.5) - Gdiff - mul(gpo, D) + dI.scale(self.iota) + mul(dI, dI).scale(0.5) + mul(gpI, dI)
            D = rhs.drop0().scale(1.0 / self.o)
        return D

    def map(self, h=0.02, clock=True):
        """Normalised Theta(n), multilinear, with the physical clock composition."""
        Ls = {}
        for j in (-2, -1, 0, 1, 2):
            I = self.exponent(self.o + j * h)
            Ls[j] = F({k: -(v - np.sum(self.w * v)) for k, v in I.items()})
        keys = set().union(*Ls.values())
        L0 = F({k: Ls[0][k] for k in keys})
        g = lambda k, j: Ls[j].get(k, 0)  # noqa: E731
        d1 = F({k: (g(k, -2) - 8 * g(k, -1) + 8 * g(k, 1) - g(k, 2)) / (12 * h) for k in keys})
        d2 = F({k: (-g(k, -2) + 16 * g(k, -1) - 30 * g(k, 0) + 16 * g(k, 1) - g(k, 2)) / (12 * h * h) for k in keys})
        L = L0
        if clock:
            D = self.clock()
            L = L0 + mul(D, d1) + mul(mul(D, D), d2).scale(0.5)
        num = fexp(L)
        den = F({k: np.sum(self.w * v) if np.ndim(v) else v for k, v in num.items()})
        d = den.drop0()
        inv = const(1.0) - d + mul(d, d) - mul(mul(d, d), d)
        return (mul(num, inv) - const(1.0)).drop0()

    def quadrupole(self, Theta, C, labels):
        w = 7.5 * np.einsum("si,ij,sj->s", self.n, C, self.n)
        return np.sum(self.w * w * Theta.get(frozenset(labels), 0))
