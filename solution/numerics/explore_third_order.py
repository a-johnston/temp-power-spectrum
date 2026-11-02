"""EXPLORATORY numerical third-order map at (k, q, -q) and the missing (1/3) Re C[theta, L3].

NOT EVIDENCE (see explore_fourth_order.py).  A numerical multilinear engine: every quantity is a
dict {label set -> array}, labels 1 -> k, 2 -> q, 3 -> -q, products only between disjoint label
sets (the multilinear rule of the report, full-series coefficients).  The map is solved from its
DEFINITIONS, not from hand-expanded formulas:
  ray      e' = -2 S e - Gamma(e,e) + e sigma(e,e),  X' = e,  X(t_o) = 0,
           e(t_o) = -e^{-zeta} exp(-H/2) n  (symmetric frame; FW rotation vanishes on the shell),
  fields   evaluated at X = Xbar + dX as sum_S f_S(t) e^{i K_S.Xbar} exp(i K_S.dX),
  emission ln t_e + nu(t_e, X(t_e)) = ln tbar_*,
  exponent I = int_{t_e}^{t_o} sigma(e,e) dt,  L_c = -(I - <I>_sky),
  clock    ln t_I + nu(t_I,0) = ln tbar_I,  int_{t_I}^{t_o} t e^{nu(t,0)} dt = (tbar_o^2 - tbar_I^2)/2,
  compose  L = L_c(tbar_o + D)  (finite differences in t_o),  Theta = e^L/<e^L>_sky - 1.
Second-order pair fields come from the Green solver of explore_fourth_order.py (pair (q,-q) is the
exact K = 0 kernel); the third-order field A_3 from the same Green formulas with the exact
third-order sources (third_order.py, including the Hamiltonian source).

Modes:
  python explore_third_order.py check2   independent second-order map vs explore_fourth_order
  python explore_third_order.py ward3    K3(k,q,-q) -> (k d_k)^2 K1 as q -> 0 (O(q^2) residual)
  python explore_third_order.py c13      (1/3) Re C_2[theta, L3] / (sigma^2 C_2^(11)) vs x_o
"""

import os
import pickle
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))  # solution/numerics
SOL = os.path.dirname(HERE)  # solution/: exact modules and results/
sys.path.insert(0, HERE)
sys.path.insert(0, SOL)
import explore_fourth_order as ef  # noqa: E402

s3 = np.sqrt(3.0)
FULL = frozenset((1, 2, 3))
EMPTY = frozenset()


# ------------------------------------------------------------------ multilinear numeric fields
class F(dict):
    """{frozenset(labels): ndarray}.  Arrays broadcast; tensor indices are trailing."""

    def __add__(self, o):
        r = F(self)
        for k, v in o.items():
            r[k] = r[k] + v if k in r else v
        return r

    def __sub__(self, o):
        return self + o.scale(-1)

    def scale(self, c):
        return F({k: c * v for k, v in self.items()})

    def map(self, fn):
        return F({k: fn(v) for k, v in self.items()})

    def drop0(self):
        return F({k: v for k, v in self.items() if k != EMPTY})


def mul(A, B, op=lambda a, b: a * b):
    r = F()
    for ka, va in A.items():
        for kb, vb in B.items():
            if ka & kb:
                continue
            k = ka | kb
            v = op(va, vb)
            r[k] = r[k] + v if k in r else v
    return r


def const(x):
    return F({EMPTY: np.asarray(x, dtype=complex)})


def fexp(X, order=3):
    """exp of a field with no EMPTY component, truncated (multilinear)."""
    out = const(1.0)
    term = const(1.0)
    for n in range(1, order + 1):
        term = mul(term, X).scale(1.0 / n)
        out = out + term
    return out


mm = lambda a, b: a @ b  # matrix-matrix
mv = lambda a, b: np.einsum("...ij,...j->...i", a, b)  # matrix-vector
vv = lambda a, b: np.einsum("...i,...i->...", a, b)  # vector-vector
sv = lambda a, b: a[..., None] * b  # scalar-vector
smat = lambda a, b: a[..., None, None] * b  # scalar-matrix


def mexp_and_derivative(H, Hd, sign=1.0):
    """exp(sign H) and its derivative along Hd (d/dt or d/dx_l), truncated at third order."""
    Hs = H.scale(sign)
    Hds = Hd.scale(sign)
    E = const(np.eye(3)) + Hs + mul(Hs, Hs, mm).scale(0.5) + mul(mul(Hs, Hs, mm), Hs, mm).scale(1 / 6)
    dE = (
        Hds
        + (mul(Hds, Hs, mm) + mul(Hs, Hds, mm)).scale(0.5)
        + (mul(mul(Hds, Hs, mm), Hs, mm) + mul(mul(Hs, Hds, mm), Hs, mm) + mul(mul(Hs, Hs, mm), Hds, mm)).scale(1 / 6)
    )
    return E, dE


# ------------------------------------------------------------------ third-order sources
_S3 = None


def third_sources():
    global _S3
    if _S3 is not None:
        return _S3
    path = os.path.join(SOL, "results", "third_order_sources_numeric.pkl")
    if not os.path.exists(path):
        # exact third-order sources at (k, q, -q), including the Hamiltonian source (about 20 s)
        import third_order as to

        E, unk, _ = to.build()
        lab = frozenset((1, 2, 3))
        fnl = [unk["n"]] + [e for e in set(unk["A"]) if isinstance(e, sp.Function)]
        src = lambda e: e.subs({d: 0 for d in e.atoms(sp.Derivative) if d.expr in fnl}).subs({f: 0 for f in fnl})
        raw = dict(
            N=src(E["N"].coeff(lab)),
            M=[src(E["M"][i].coeff(lab)) for i in range(3)],
            T=sp.Matrix(3, 3, lambda i, j: src(E["T"][i][j].coeff(lab))),
            H=src(E["Ham"].coeff(lab)),
        )
        with open(path, "wb") as fh:
            pickle.dump(raw, fh)
    with open(path, "rb") as fh:
        raw = pickle.load(fh)
    exprs = [raw["N"], raw["H"]] + list(raw["M"]) + [raw["T"][i, j] for i in range(3) for j in range(3)]
    funcs = set()
    for e in exprs:
        funcs |= {f for f in e.atoms(sp.Function) if f.func.__name__.startswith(("n1", "n2", "A1", "A2"))}
    tsym = [s for s in raw["N"].free_symbols if str(s) == "t"][0]
    names = []
    rep = {}
    for f in sorted(funcs, key=lambda f: f.func.__name__):
        nm = f.func.__name__
        for d in (2, 1):
            sym = sp.Symbol(f"{nm}_d{d}")
            rep[sp.Derivative(f, (tsym, d))] = sym
            names.append(sym)
        sym = sp.Symbol(f"{nm}_d0")
        rep[f] = sym
        names.append(sym)
    by = {s.name: s for e in exprs for s in e.free_symbols}
    base = [by.get(nm, sp.Symbol(nm)) for nm in ("t", "p", "q", "b1", "b2")]
    lam = [sp.lambdify(base + names, e.xreplace(rep), "numpy") for e in exprs]
    _S3 = dict(lam=lam, names=[s.name for s in names])
    return _S3


# ------------------------------------------------------------------ second-order pair fields
def rot_z(phi):
    c, s = np.cos(phi), np.sin(phi)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def pair_fields(ka, kb, tg):
    """nu, nu', nu'', A, A', A'' of the second-order kernel for the pair (ka, kb) on tg (global frame).

    Uses the frame k1 = (|ka|,0,0) of second_order.py after rotating ka onto the x axis.  For a
    zero-total-momentum pair the exact K = 0 Green functions are used.
    """
    if ef.SRC is None:
        ef.SRC = ef.load_sources()
    pa = np.linalg.norm(ka)
    phi = np.arctan2(ka[1], ka[0])
    R = rot_z(phi)
    kb_loc = R.T @ kb
    p, b1, b2 = pa, kb_loc[0], kb_loc[1]
    q = np.linalg.norm(kb)
    KV = np.array([p + b1, b2, 0.0])
    K = np.linalg.norm(KV)
    ev = lambda fn: np.broadcast_to(np.asarray(fn(tg, p, q, b1, b2), dtype=complex), tg.shape)
    cN = ev(ef.SRC["N"])
    cT = np.array([[ev(ef.SRC["T"][i][j]) for j in range(3)] for i in range(3)])
    c2 = (2 * (p**2 + q**2) - p * b1) / 9
    if K < 1e-12:
        # exact K = 0 kernel: nu = c2 t^2 - (1/3) int (t^3 - s^3)/(3 t s) c_N, A = -2 int (s/t)(t-s) c_T
        I1 = ef.cumint(cN / tg, tg)  # int c_N / s
        I2 = ef.cumint(tg**2 * cN, tg)  # int s^2 c_N
        nu = c2 * tg**2 - (tg**2 * I1 / 3 - I2 / (3 * tg)) / 3
        nup = 2 * c2 * tg - (2 * tg * I1 / 3 + I2 / (3 * tg**2)) / 3
        nupp = (6 * nu / tg**2 - cN) / 3
        J1 = ef.cumint(tg * cT, tg)
        J2 = ef.cumint(tg**2 * cT, tg)
        A = -2 * (J1 - J2 / tg)
        Ap = -2 * J2 / tg**2
        App = -2 * Ap / tg - 2 * cT
    else:
        tg_, A, Ap = ef.second_order_A2(p, q, b1, b2, tg[-1], ng=len(tg), lo=tg[0] / tg[-1])
        assert np.allclose(tg_, tg)
        x = K / s3 * tg
        u1 = np.sin(x) / x - np.cos(x)
        u2 = -np.cos(x) / x - np.sin(x)
        P1 = ef.cumint(u1 * cN, tg)
        P2 = ef.cumint(u2 * cN, tg)
        c = K / s3
        nu = 9 * c2 / K**2 * u1 - (u2 * P1 - u1 * P2) / (3 * c)
        u1p = np.cos(x) / x - np.sin(x) / x**2 + np.sin(x)
        u2p = np.sin(x) / x + np.cos(x) / x**2 - np.cos(x)
        nup = 9 * c2 / K**2 * c * u1p - (u2p * P1 - u1p * P2) / 3
        nupp = (-(K**2) * nu + 6 * nu / tg**2 - cN) / 3
        R1 = 0.5 * (K**2 * A - np.einsum("a,i,ajn->ijn", KV, KV, A) - np.einsum("a,j,ain->ijn", KV, KV, A))
        R1 = R1 - np.eye(3)[:, :, None] * np.einsum("iin->n", R1)[None, None, :] / 3
        KK = np.outer(KV, KV) - np.eye(3) * K**2 / 3
        App = -2 * Ap / tg - 2 * (R1 + 2 * KK[:, :, None] * nu + cT)
    rot = lambda M: np.einsum("ia,abn,jb->ijn", R, M, R)
    tomat = lambda M: np.moveaxis(M, -1, 0)
    return dict(nu=nu, nup=nup, nupp=nupp, A=tomat(rot(A)), Ap=tomat(rot(Ap)), App=tomat(rot(App)))


# ------------------------------------------------------------------ third-order field
def third_field(k, q, tg, pairs, scalar_only=False):
    """A_3 and A_3' at (k, q, -q) with k = (p,0,0), q = (b1,b2,0); pairs = pair_fields dicts."""
    S3 = third_sources()
    p = k[0]
    b1, b2 = q[0], q[1]
    qm = np.linalg.norm(q)
    vals = {}
    keymap = {"12": (1, 2), "13": (1, 3), "23": (2, 3)}
    for nm in S3["names"]:
        base, d = nm.rsplit("_d", 1)
        d = int(d)
        if base.startswith("n"):
            pr = pairs[keymap[base[1:]]]
            arr = [pr["nu"], pr["nup"], pr["nupp"]][d]
        else:
            pair, ij = base[1:].split("_")
            pr = pairs[keymap[pair]]
            arr = [pr["A"], pr["Ap"], pr["App"]][d][:, int(ij[0]), int(ij[1])]
        vals[nm] = arr
    args = [tg, p, qm, b1, b2] + [vals[nm] for nm in S3["names"]]
    out = [np.broadcast_to(np.asarray(fn(*args), dtype=complex), tg.shape) for fn in S3["lam"]]
    cN, cH = out[0], out[1]
    cM = np.array(out[2:5])
    cT = np.array(out[5:14]).reshape(3, 3, -1)
    KV = np.array(k, dtype=float)
    K = np.linalg.norm(KV)
    c3 = 2 / 9 * (-2 * (p**2 + 2 * qm**2) + (np.dot(k, q) - np.dot(k, q) - qm**2))
    x = K / s3 * tg
    u1 = np.sin(x) / x - np.cos(x)
    u2 = -np.cos(x) / x - np.sin(x)
    P1 = ef.cumint(u1 * cN, tg)
    P2 = ef.cumint(u2 * cN, tg)
    c = K / s3
    nu = 9 * c3 / K**2 * u1 - (u2 * P1 - u1 * P2) / (3 * c)
    u1p = np.cos(x) / x - np.sin(x) / x**2 + np.sin(x)
    u2p = np.sin(x) / x + np.cos(x) / x**2 - np.cos(x)
    nup = 9 * c3 / K**2 * c * u1p - (u2p * P1 - u1p * P2) / 3
    nupp = (-(K**2) * nu + 6 * nu / tg**2 - cN) / 3
    hS = 3 / (2 * K**2) * (4 * K**2 * nu + 12 * nup / tg + 12 * nu / tg**2 + cH)
    hSp = 3 / (2 * K**2) * (4 * K**2 * nup + 12 * nupp / tg - 24 * nu / tg**3 + np.gradient(cH, tg))
    Kh = KV / K
    proj = (np.outer(Kh, Kh) - np.eye(3) / 3)[:, :, None]
    A = proj * hS
    Ap = proj * hSp
    if not scalar_only:
        Pp = np.eye(3) - np.outer(Kh, Kh)
        w = 2j / K**2 * ef.cumint(Pp @ cM, tg)
        wp = 2j / K**2 * (Pp @ cM)
        PMP = np.einsum("ia,abn,bj->ijn", Pp, cT, Pp)
        TT = PMP - Pp[:, :, None] * np.einsum("aan->n", PMP)[None, None, :] / 2
        Ic = ef.cumint(tg * np.cos(K * tg) * TT, tg)
        Is = ef.cumint(tg * np.sin(K * tg) * TT, tg)
        h = -2 / (K * tg) * (np.sin(K * tg) * Ic - np.cos(K * tg) * Is)
        hp = 2 / (K * tg**2) * (np.sin(K * tg) * Ic - np.cos(K * tg) * Is) - 2 / tg * (np.cos(K * tg) * Ic + np.sin(K * tg) * Is)
        A = A + KV[:, None, None] * w[None] + w[:, None] * KV[None, :, None] + h
        Ap = Ap + KV[:, None, None] * wp[None] + wp[:, None] * KV[None, :, None] + hp
    return dict(A=np.moveaxis(A, -1, 0), Ap=np.moveaxis(Ap, -1, 0), nu=nu)


# ------------------------------------------------------------------ the map
class Map:
    def __init__(self, x_o, ks, R=3.0, RI=0.5, nt=300, nth=16, nph=32, ng=12000, third=True, scalar_only=False, xstar=None):
        self.ks = {a: np.array(v, dtype=float) for a, v in ks.items()}
        self.t_o = s3 * x_o
        Rx = R if xstar is None else x_o / xstar
        self.ts = self.t_o / Rx
        self.tI = self.ts * RI
        self.nt = nt
        self.n, self.w = ef.sphere(nth, nph)
        self.tg = np.linspace(self.t_o * 2e-4, self.t_o * 1.02, ng)  # extends past t_o for the clock
        KS = lambda S: sum((self.ks[a] for a in S), np.zeros(3))
        self.KS = KS
        tg = self.tg
        # time profiles of the primordial and metric fields, per label set
        self.zeta = F({frozenset([a]): np.ones_like(tg) for a in (1, 2, 3)})
        self.H, self.Hp, self.nu = F(), F(), F()
        for a in (1, 2, 3):
            k = self.ks[a]
            Q, Qp, nu = ef.mode(np.linalg.norm(k), tg)
            M = np.outer(k, k) - np.dot(k, k) * np.eye(3) / 3
            self.H[frozenset([a])] = Q[:, None, None] * M
            self.Hp[frozenset([a])] = Qp[:, None, None] * M
            self.nu[frozenset([a])] = nu
        self.pairs = {}
        for pr in ((1, 2), (1, 3), (2, 3)):
            f2 = pair_fields(self.ks[pr[0]], self.ks[pr[1]], tg)
            self.pairs[pr] = f2
            S = frozenset(pr)
            self.H[S] = f2["A"]
            self.Hp[S] = f2["Ap"]
            self.nu[S] = f2["nu"]
        if third:
            f3 = third_field(self.ks[1], self.ks[2], tg, self.pairs, scalar_only=scalar_only)
            self.H[FULL] = f3["A"]
            self.Hp[FULL] = f3["Ap"]

    # interpolate a profile F onto times tt
    def at_times(self, Fp, tt):
        def interp(arr):
            flat = arr.reshape(len(self.tg), -1)
            out = np.stack(
                [
                    np.interp(tt, self.tg, flat[:, j].real) + 1j * np.interp(tt, self.tg, flat[:, j].imag)
                    for j in range(flat.shape[1])
                ],
                axis=-1,
            )
            return out.reshape((len(tt),) + arr.shape[1:])

        return F({k: interp(v) for k, v in Fp.items()})

    def on_ray(self, prof_t, dX, tr, deriv=False):
        """Field with time profiles prof_t (already on tr) at X = Xbar + dX.  If deriv, return the
        spatial gradient (extra trailing index l) instead."""
        lam = self.t_o_cur - tr
        out = F()
        for S, v in prof_t.items():
            if S == EMPTY:  # constant (no time axis, zero wavevector)
                v = np.asarray(v)
                if deriv:
                    continue
                out = out + F(
                    {
                        EMPTY: np.broadcast_to(
                            v, (len(self.n), len(tr)) + v.shape[-2:] if v.ndim >= 2 else (len(self.n), len(tr)) + v.shape
                        )
                    }
                )
                continue
            KS = self.KS(S)
            ph = np.exp(1j * (self.n @ KS)[:, None] * lam[None, :])  # (Ns, Nt)
            base = v[None] * ph.reshape(ph.shape + (1,) * (v.ndim - 1))
            if deriv:
                base = base[..., None] * (1j * KS)
            y = F({kk: 1j * np.einsum("...i,i->...", vv_, KS) for kk, vv_ in dX.items() if kk != EMPTY})
            ex = fexp(y)
            piece = mul(F({S: base}), ex, lambda a, b: a * b.reshape(b.shape + (1,) * (a.ndim - b.ndim)))
            out = out + piece
        return out

    def build(self, t_o, drop=()):
        self.t_o_cur = t_o
        tr = np.linspace(self.ts, t_o, self.nt)
        n = self.n
        Ns, Nt = len(n), len(tr)
        Hr, Hpr, zr, nur = (self.at_times(P, tr) for P in (self.H, self.Hp, self.zeta, self.nu))
        # observer values (X = 0, t = t_o): profiles at t_o
        H_o = F({k: v[0] for k, v in self.at_times(self.H, np.array([t_o])).items()})
        z_o = F({k: v[0] for k, v in self.at_times(self.zeta, np.array([t_o])).items()})
        nvec = n[:, None, :] * np.ones((1, Nt, 1))
        e_bg = const(-nvec)
        # initial direction e(t_o) = -e^{-zeta} exp(-H/2) n at the observer
        Eh, _ = mexp_and_derivative(H_o, H_o, sign=-0.5)
        ez = fexp(z_o.scale(-1))
        Mo = mul(ez, Eh, smat)  # matrices (3, 3) per label set
        e_o = F({kk: -np.einsum("ij,sj->si", v, n) for kk, v in Mo.items()})
        de_o = F({k: v for k, v in e_o.items() if k != EMPTY})  # perturbation at t_o, (Ns, 3)
        dX = F()
        dX[EMPTY] = np.zeros((Ns, Nt, 3))
        de = F()
        de[EMPTY] = np.zeros((Ns, Nt, 3))
        # --- algebra on time-only profiles (products at a point commute with evaluation at X) ---
        Hr0, Hpr0 = Hr.drop0(), Hpr.drop0()
        Mprof, Mpprof = mexp_and_derivative(Hr0, Hpr0, 1.0)
        Minvprof, _ = mexp_and_derivative(Hr0, Hpr0, -1.0)
        e2zprof = fexp(zr.drop0().scale(2.0))
        grad = lambda P: F({kk: v[..., None] * (1j * self.KS(kk)) for kk, v in P.items() if kk != EMPTY})
        dMprof = grad(Mprof)  # (Nt,3,3,l)
        dzprof = grad(zr)  # (Nt,l)
        Sprof = mul(Minvprof, Mpprof, mm).scale(0.5)
        sigprof = mul(e2zprof, Mpprof, smat).scale(0.5)
        P1prof = mul(Minvprof, dMprof, lambda a, b: np.einsum("...ia,...ajl->...ijl", a, b))  # (M^-1 d_l M)_ij
        P2prof = mul(Minvprof, dMprof, lambda a, b: np.einsum("...il,...jkl->...ijk", a, b))  # M^{il} d_l M_jk
        Q3prof = mul(Minvprof, dzprof, mv)  # M^{-1} d zeta
        for _ in range(3):
            ray = lambda P, dX=dX: self.on_ray(P, dX, tr)
            e = e_bg + de.drop0()
            Se = mul(ray(Sprof), e, mv)
            see = mul(mul(ray(sigprof), e, mv), e, vv)  # sigma(e,e)
            t1a = mul(mul(ray(P1prof), e, lambda a, b: np.einsum("...ijl,...l->...ij", a, b)), e, mv)
            t1b = mul(mul(ray(P2prof), e, lambda a, b: np.einsum("...ijk,...k->...ij", a, b)), e, mv)
            term1 = t1a - t1b.scale(0.5)
            edz = mul(e, ray(dzprof), vv)
            term2 = mul(e, edz, lambda a, b: a * b[..., None]).scale(2.0)
            eMe = mul(mul(e, ray(Mprof), lambda a, b: np.einsum("...i,...ij->...j", a, b)), e, vv)
            term3 = mul(eMe, ray(Q3prof), lambda a, b: a[..., None] * b)
            Gam = term1 + term2 - term3
            acc = Se.scale(-2.0) - Gam + mul(e, see, lambda a, b: a * b[..., None])
            acc = acc.drop0()
            # integrate from t_o backwards
            new_de, new_dX = F(), F()
            for kk, a in acc.items():
                cum = ef.cumint(np.moveaxis(a, 1, -1), tr)
                R0 = np.moveaxis(cum[..., -1:] - cum, -1, 1)
                v = -R0
                if kk in de_o:
                    v = v + de_o[kk][:, None, :]
                new_de[kk] = v
            for kk, v in de_o.items():
                if kk not in new_de:
                    new_de[kk] = np.broadcast_to(v[:, None, :], (Ns, Nt, 3)).copy()
            for kk, v in new_de.items():
                cum = ef.cumint(np.moveaxis(v, 1, -1), tr)
                new_dX[kk] = -np.moveaxis(cum[..., -1:] - cum, -1, 1)
            new_dX[EMPTY] = np.zeros((Ns, Nt, 3))
            new_de[EMPTY] = np.zeros((Ns, Nt, 3))
            if "ray" in drop:
                new_dX = F({EMPTY: np.zeros((Ns, Nt, 3))})
                new_de = F({EMPTY: np.zeros((Ns, Nt, 3))})
            de, dX = new_de, new_dX
        Fint = see.drop0()
        Itot = F({kk: np.trapezoid(v, tr, axis=1) for kk, v in Fint.items()})
        # emission root: delta = t*(exp(-phi(t*+delta)) - 1), phi = nu along the ray
        phi = self.on_ray(nur, dX, tr).drop0()
        phi0 = F({kk: v[:, 0] for kk, v in phi.items()})
        phid = F({kk: np.gradient(v, tr, axis=1)[:, 0] for kk, v in phi.items()})
        delta = F()
        for _ in range(3):
            arg = phi0 + mul(phid, delta)
            delta = (fexp(arg.scale(-1.0)) - const(1.0)).drop0().scale(self.ts)
        if "endpoint" in drop:
            delta = F()
        F0 = F({kk: v[:, 0] for kk, v in Fint.items()})
        F1 = F({kk: np.gradient(v, tr, axis=1)[:, 0] for kk, v in Fint.items()})
        Iend = mul(F0, delta) + mul(mul(F1, delta), delta).scale(0.5)
        I = Itot - Iend
        Lc = F({kk: -(v - ef.smean(v, self.w)) for kk, v in I.items()})
        return Lc

    def clock(self):
        """D = t_o - tbar_o as a multilinear field (orders 1, 2)."""
        tg = self.tg
        om = self.nu
        g = (fexp(om) - const(1.0)).drop0()
        gT = F({kk: ef.cumint(tg * v, tg) for kk, v in g.items()})
        val = lambda Fp, T: F({kk: np.interp(T, tg, v.real) + 1j * np.interp(T, tg, v.imag) for kk, v in Fp.items()})
        gp = F({kk: tg * v for kk, v in g.items()})
        omI, omIp = val(om, self.tI), F({kk: np.interp(self.tI, tg, np.gradient(v, tg)) for kk, v in om.items()})
        dI = F()
        for _ in range(3):
            dI = (fexp((omI + mul(omIp, dI)).scale(-1.0)) - const(1.0)).drop0().scale(self.tI)
        to = self.t_o
        D = F()
        for _ in range(3):
            rhs = (
                mul(D, D).scale(-0.5)
                - val(gT, to)
                - mul(val(gp, to), D)
                + dI.scale(self.tI)
                + mul(dI, dI).scale(0.5)
                + val(gT, self.tI)
                + mul(val(gp, self.tI), dI)
            )
            D = rhs.drop0().scale(1.0 / to)
        return D

    def L(self, h=2e-3, drop=()):
        L0 = self.build(self.t_o, drop)
        Lp = self.build(self.t_o + h, drop)
        Lm = self.build(self.t_o - h, drop)
        d1 = F({kk: (Lp[kk] - Lm[kk]) / (2 * h) for kk in L0})
        d2 = F({kk: (Lp[kk] - 2 * L0[kk] + Lm[kk]) / h**2 for kk in L0})
        D = F() if "clock" in drop else self.clock()
        return L0 + mul(D, d1) + mul(mul(D, D), d2).scale(0.5)

    def Theta(self, L):
        num = fexp(L.drop0())
        den = F({kk: ef.smean(v, self.w) if np.ndim(v) else v for kk, v in num.items()})
        # 1/den with den = 1 + d
        d = den.drop0()
        inv = const(1.0) - d + mul(d, d) - mul(mul(d, d), d)
        return (mul(num, inv) - const(1.0)).drop0()


def symmetric_cfg(q, mu):
    return {1: (1.0, 0, 0), 2: (q * mu, q * np.sqrt(1 - mu**2), 0), 3: (-q * mu, -q * np.sqrt(1 - mu**2), 0)}


# ------------------------------------------------------------------ checks and spectrum
def check2(x_o=1.5):
    """Independent second-order map: engine L at {1,2} vs explore_fourth_order K^L."""
    mu = 0.3
    cfg = symmetric_cfg(1.0, mu)
    Mp = Map(x_o, cfg, third=False, nt=320, nth=24, nph=48)
    L = Mp.L()
    S = ef.Setup(x_o, nt=320, nth=24, nph=48)
    a, b = S.first(np.array(cfg[1])), S.first(np.array(cfg[2]))
    k2 = S.second(a, b, ef.second_order_A2(1.0, 1.0, mu, np.sqrt(1 - mu**2), S.t_o))
    eng = L[frozenset((1, 2))]
    ref = k2["KL"]
    err = np.sqrt(ef.smean(np.abs(eng - ref) ** 2, S.w) / ef.smean(np.abs(ref) ** 2, S.w))
    e1 = np.sqrt(ef.smean(np.abs(L[frozenset([1])] - a["theta"]) ** 2, S.w) / ef.smean(np.abs(a["theta"]) ** 2, S.w))
    print(f"  check2 x_o={x_o}: first order rel err {e1:.1e};  second-order L2 kernel rel err {err:.1e}")


def ward3(x_o=1.5, qs=(0.08, 0.04)):
    k = 1.0
    h = 2e-3

    def K1(kk):
        Mp = Map(x_o, {1: (kk, 0, 0), 2: (0, 0.01, 0), 3: (0, -0.01, 0)}, third=False, nth=12, nph=24)
        return Mp.Theta(Mp.L())[frozenset([1])], Mp.w

    km, w = K1(k * np.exp(-h))
    k0, _ = K1(k)
    kp, _ = K1(k * np.exp(h))
    target = (kp - 2 * k0 + km) / h**2  # (k d_k)^2 K1
    for q in qs:
        Mp = Map(x_o, symmetric_cfg(q, 0.6), nth=12, nph=24)
        K3 = Mp.Theta(Mp.L())[FULL]
        err = np.sqrt(ef.smean(np.abs(K3 - target) ** 2, w) / ef.smean(np.abs(target) ** 2, w))
        print(f"  ward3 x_o={x_o} q={q}: rel L2 error of K3(k,q,-q) vs (k d_k)^2 K1 = {err:.2e}", flush=True)


def c13(xos=(0.15, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0), nmu=16, l=2, xstar=None, scalar_only=True):
    mus, wm = np.polynomial.legendre.leggauss(nmu)
    rows = []
    for xo in xos:
        nth = max(16, int(8 + 4 * xo))
        acc = 0
        th = None
        for mu, w in zip(mus, wm, strict=True):
            Mp = Map(
                xo, symmetric_cfg(1.0, mu), nth=nth, nph=2 * nth, nt=max(240, int(100 * xo)), scalar_only=scalar_only, xstar=xstar
            )
            L = Mp.L()
            th = L[frozenset([1])]
            acc += 0.5 * w * ef.S_l(th, L[FULL], Mp.n, Mp.w, l).real
        SS = ef.S_l(th, th, Mp.n, Mp.w, l).real
        rows.append((xo, acc / SS, SS / (2 * l + 1)))
        print(
            f"  x_o={xo}: (1/3) Re C_{l}[theta,L3] / (sigma^2 C11) = {acc / SS:.6g};  (1/3)Re C[theta,L3]/sigma^4 = {acc / (2 * l + 1):.5g};  C11/sigma^2 = {SS / (2 * l + 1):.4g}",
            flush=True,
        )
    return rows


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else "check2"
    {"check2": check2, "ward3": ward3, "c13": c13}[m]()
