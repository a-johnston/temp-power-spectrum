"""Assembled cubic quadrupole from the report's reduced formulas, evaluated numerically (VALIDATION).

  J_abc = J_111 - sum_{a|bc} S_meas,a[nu_bc, A_bc] + J_lin[y_abc]          (Prop. combined)

  J_111      the second-ray-free direct first-field measurement, built by the repository's own
             direct_cubic_measurement.first_only_exponent_triple and normalised as in
             normalized_cubic_raw; its nested definite integrals are evaluated by Gauss quadrature.
  S_meas     the no-probe-ray measurement Hessian of the first field-adjoint source (report
             eq. measurement-source), implemented here from the formula, applied to the solved
             second-order pair fields.
  J_lin      the constraint-row linear measurement (report eq. combined-linear): acoustic lapse,
             LOCAL momentum source and TT wave response of the complete third-order source.

All quantities use averaged sky weights: J = <w Theta>, w = (15/2) n.C.n.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)
import cubic_fields as cf  # noqa: E402
from forward_cubic import Spectral, sphere  # noqa: E402

# ----------------------------------------------------------------------------- numeric evaluator
_FUN = {sp.exp: np.exp, sp.sin: np.sin, sp.cos: np.cos, sp.sqrt: np.sqrt}


class Evaluator:
    """Vectorised evaluation of a SymPy expression with nested definite Integrals."""

    def __init__(self, ng=20):
        self.x, self.w = np.polynomial.legendre.leggauss(ng)

    def __call__(self, expr, env):
        return self._ev(expr, env, {})

    def _ev(self, e, env, cache):
        if e in cache:
            return cache[e]
        if e.is_Number:
            out = complex(e)
        elif e is sp.I:
            out = 1j
        elif e.is_Symbol:
            out = env[e]
        elif e.is_Add:
            out = 0
            for a in e.args:
                out = out + self._ev(a, env, cache)
        elif e.is_Mul:
            out = 1
            for a in e.args:
                out = out * self._ev(a, env, cache)
        elif e.is_Pow:
            b = self._ev(e.args[0], env, cache)
            ex = e.args[1]
            out = b ** (int(ex) if ex.is_Integer else self._ev(ex, env, cache))
        elif isinstance(e, sp.Integral):
            f, lim = e.args[0], e.args[1:]
            assert all(len(li) == 3 for li in lim), "definite limits expected"
            if len(lim) > 1:  # innermost limit first: peel off the outermost
                f = sp.Integral(f, *lim[:-1])
            x, lo, hi = lim[-1]
            lo, hi = np.asarray(self._ev(lo, env, cache)), np.asarray(self._ev(hi, env, cache))
            new = {k: np.asarray(v)[..., None] for k, v in env.items()}
            half = (hi - lo) / 2
            new[x] = lo[..., None] + half[..., None] * (self.x + 1)
            val = np.asarray(self._ev(f, new, {}))
            val = np.broadcast_to(val, np.broadcast_shapes(val.shape, new[x].shape))
            out = half * np.tensordot(val, self.w, axes=([-1], [0]))
        elif e.func in _FUN:
            out = _FUN[e.func](self._ev(e.args[0], env, cache))
        else:
            raise TypeError(f"unsupported node {e.func}")
        cache[e] = out
        return out


# ----------------------------------------------------------------------------- J_111
def direct_measurement(momenta, e, o, iota, C, nth=18, nph=36, ng=20, chunk=36):
    from direct_cubic_measurement import FourierFirst, first_only_exponent_triple

    t = sp.Symbol("t", positive=True)
    nsym = sp.Matrix(sp.symbols("n1 n2 n3", real=True))

    def seed(k):
        k = sp.Matrix([sp.nsimplify(c) for c in k])
        km2 = k.dot(k)
        x = sp.sqrt(km2) * t / sp.sqrt(3)
        a = sp.cos(x) - sp.sin(x) / x
        b = -6 * (2 * sp.sin(x) / x - sp.cos(x) - 1)
        return FourierFirst(k, 1, a, b * (k * k.T / km2 - sp.eye(3) / 3), t)

    R = sp.nsimplify
    data = first_only_exponent_triple([seed(momenta[a]) for a in (1, 2, 3)], nsym, R(e), R(o), R(iota))
    n, w = sphere(nth, nph)
    ev = Evaluator(ng)
    first = [np.zeros(len(n), complex) for _ in range(3)]
    second = {k: np.zeros(len(n), complex) for k in data.second}
    third = np.zeros(len(n), complex)
    for lo in range(0, len(n), chunk):
        sl = slice(lo, lo + chunk)
        env = {nsym[i]: n[sl, i] for i in range(3)}
        for a in range(3):
            first[a][sl] = ev(data.first[a], env)
        for k, x in data.second.items():
            second[k][sl] = ev(x, env)
        third[sl] = ev(data.third, env)
    avg = lambda f: np.sum(w * f)  # noqa: E731
    u = [f - avg(f) for f in first]
    raw = -third - u[0] * u[1] * u[2]
    for a in range(3):
        b, c = (i for i in range(3) if i != a)
        pair = second[(b, c)]
        raw = raw + u[a] * (pair - avg(pair)) + avg(u[b] * u[c]) * u[a]
    wq = 7.5 * np.einsum("si,ij,sj->s", n, C, n)
    return avg(wq * raw)


# ----------------------------------------------------------------------------- S_meas
def stf_sym(u, v):
    """TF sym(u v^T) for arrays (..., 3)."""
    M = 0.5 * (u[..., :, None] * v[..., None, :] + v[..., :, None] * u[..., None, :])
    return M - np.eye(3) * np.trace(M, axis1=-2, axis2=-1)[..., None, None] / 3


def measurement_hessian(p, q, first_kmag, Qhat, probe, e, o, iota, C, nth=18, nph=36, nt=44):
    """S_meas,a[m, j] (averaged weights) for first direction (z=1, a, b Qhat) at p and probe at q."""
    n, wsky = sphere(nth, nph)
    avg = lambda f: np.tensordot(wsky, f, axes=([0], [0]))  # noqa: E731
    w = 7.5 * np.einsum("si,ij,sj->s", n, C, n)
    N = n[:, :, None] * n[:, None, :] - np.eye(3) / 3
    sp_ = Spectral(e, o, nt)  # node 0 at o, last node at e
    tt = sp_.t
    f1 = cf.first_profiles(first_kmag, tt)
    a, h, hp = f1["a"][0], f1["b"][0][:, None, None] * Qhat, f1["b"][1][:, None, None] * Qhat
    pn, qn = n @ p, n @ q
    php = np.exp(1j * pn[:, None] * (o - tt)[None, :])  # (s, t)
    phq = np.exp(1j * qn[:, None] * (o - tt)[None, :])
    G = 2 * np.eye(3) + h  # (t,3,3)
    Pn = np.eye(3) - n[:, :, None] * n[:, None, :]  # (s,3,3)
    hpn = np.einsum("tij,sj->sti", hp, n)
    hpnn = np.einsum("sti,si->st", hpn, n)
    Gn = np.einsum("tij,sj->sti", G, n)
    Gnn = np.einsum("sti,si->st", Gn, n)
    # first ray, fixed coordinate observer: X(o) = v(o) = 0 (cumint from o is exactly that)
    vdot = php[..., None] * (
        np.einsum("sij,stj->sti", Pn, hpn / 2)
        + 1j * np.einsum("sij,j->si", Pn, p)[:, None, :]
        + 0.5j * (p[None, None, :] * (Gnn - 2)[..., None] - pn[:, None, None] * (Gn - 2 * n[:, None, :]))
    )
    v = sp_.cumint(vdot, 1)
    Xdot = v + php[..., None] * (n[:, None, :] + np.einsum("tij,sj->sti", h, n) / 2)
    X = sp_.cumint(Xdot, 1)
    r = 0.5j * p[None, None, :] * (php * hpnn)[..., None]
    sray = -php[..., None] * np.einsum("sij,stj->sti", Pn, hpn)
    cum_e = lambda f: sp_.cumint(f, 1) - sp_.cumint(f, 1)[:, -1:, ...]  # noqa: E731  int_e^t
    Rr = cum_e(r)
    Kk = cum_e(sray - Rr)
    Fa = 0.5 * php * hpnn
    Ua = np.einsum("st,t->s", Fa, sp_.full)
    Th = -(Ua - avg(Ua))
    lam = w * Th - avg(w * Th)
    ti_sp = Spectral(iota, o, 44)
    at = cf.first_profiles(first_kmag, ti_sp.t)["a"][0]
    aI = cf.first_profiles(first_kmag, np.array([iota]))["a"][0][0]
    Da = -(iota**2 * aI + np.dot(ti_sp.full, ti_sp.t * at)) / o
    hpo = cf.first_profiles(first_kmag, np.array([o]))["b"][1][0] * Qhat
    dUo = 0.5 * np.einsum("si,ij,sj->s", n, hpo, n) + np.einsum("st,t->s", 0.5j * pn[:, None] * php * hpnn, sp_.full)
    Aa = avg(w * dUo)
    PK = np.einsum("sij,stj->sti", Pn, Kk)
    E1 = (
        0.5j * np.einsum("j,stj->st", q, X)[..., None, None] * N[:, None]
        - stf_sym(v, n[:, None, :])
        - 0.5 * stf_sym(PK, n[:, None, :])
    )
    E0 = (
        -0.5 * stf_sym(Rr, n[:, None, :])
        - 0.5j * np.einsum("j,stj->st", q, Kk)[..., None, None] * N[:, None]
        + 0.5j * qn[:, None, None, None] * stf_sym(Kk, n[:, None, :])
    )
    MH = phq[..., None, None] * w[:, None, None, None] * E0
    MHp = phq[..., None, None] * (
        w[:, None, None, None] * E1
        + 0.5 * lam[:, None, None, None] * N[:, None]
        + 0.5j * (w * Da * qn)[:, None, None, None] * N[:, None]
    )
    MHbar, MHpbar = avg(MH), avg(MHp)  # (t,3,3)
    nu_p, A_p, Ap_p = probe  # callables: lapse, spatial field and its time derivative
    j, jp = A_p(tt), Ap_p(tt)
    bulk = np.dot(sp_.full, np.einsum("tij,tij->t", MHbar, j) + np.einsum("tij,tij->t", MHpbar, jp))
    ae = php[:, -1] * a[-1]
    me = nu_p(np.array([e]))[0]
    jpe = Ap_p(np.array([e]))[0]
    jpo = Ap_p(np.array([o]))[0]
    emission = avg(phq[:, -1] * (e * w * Fa[:, -1] * me + (e / 2) * w * ae * np.einsum("sij,ij->s", N, jpe)))
    observer = Da / 2 * np.einsum("ij,ij->", avg(w[:, None, None] * N), jpo)
    mI = nu_p(np.array([iota]))[0]
    clock = -Aa / o * (iota**2 * mI + np.dot(ti_sp.full, ti_sp.t * nu_p(ti_sp.t)))
    return bulk + emission + observer + clock


# ----------------------------------------------------------------------------- J_lin (constraint rows)
def projector_T(K, W):
    K2 = K @ K
    return 2 * np.einsum("...ij,j->...i", W, K) / K2 - np.einsum("i,...ij,j->...", K, W, K)[..., None] * K / (2 * K2**2)


def linear_constraint_rows(momenta, profiles, triple, e, o, C, nth=18, nph=36, nt=44):
    """J_lin[y_abc] from the acoustic lapse, the local momentum source and the TT wave response."""
    P = sum((np.asarray(momenta[a], dtype=float) for a in triple), np.zeros(3))
    n, wsky = sphere(nth, nph)
    w = 7.5 * np.einsum("si,ij,sj->s", n, C, n)
    sp_ = Spectral(e, o, nt)
    tt = sp_.t
    ph = np.exp(1j * (n @ P)[:, None] * (o - tt)[None, :])
    N = n[:, :, None] * n[:, None, :] - np.eye(3) / 3
    CP = np.einsum("s,st,sij->tij", wsky * w, ph, N)  # <w TF(nn) e^{iP.n(o-t)}>
    P2 = P @ P
    ell = np.einsum("i,tij,j->t", P, CP, P) / P2
    pr = profiles[frozenset(triple)]
    nu, nup = pr.nu(tt), pr.nup(tt)
    members = {S: v for S, v in profiles.items() if frozenset(triple) != S}
    _, _, cM, cT = cf.evaluate_sources(momenta, triple, members, tt)
    lapse = -3 * np.dot(sp_.full, ell * (nup - nu / tt))
    momentum = -1j * np.dot(sp_.full, np.einsum("ti,ti->t", projector_T(P, CP), cM))
    # TT wave response h(t) = -2 int_0^t (s/t) sin(|P|(t-s))/|P| TT[c_T](s) ds, and h'(t)
    Pm = np.sqrt(P2)
    Pp = np.eye(3) - np.outer(P, P) / P2

    def TT(M):
        PMP = np.einsum("ia,...ab,bj->...ij", Pp, M, Pp)
        return PMP - Pp * np.trace(PMP, axis1=-2, axis2=-1)[..., None, None] / 2

    s = 0.5 * tt[:, None] * (cf.GL_X[None, :] + 1)
    ws = 0.5 * tt[:, None] * cf.GL_W[None, :]
    cTs = cf.evaluate_sources(momenta, triple, members, s)[3]
    d = tt[:, None] - s
    kern_d = -2 * (s / tt[:, None]) * (np.cos(Pm * d) - np.sin(Pm * d) / (Pm * tt[:, None]))
    hp = np.einsum("tg,tgij,tg->tij", kern_d, TT(cTs), ws)
    tensor = -0.5 * np.dot(sp_.full, np.einsum("tij,tij->t", TT(CP), hp))
    return lapse + momentum + tensor, dict(lapse=lapse, momentum=momentum, tensor=tensor)


def linear_cartesian(momenta, profiles, triple, e, o, C, nth=18, nph=36, nt=44):
    """-1/2 int C_P : A_3' with the Cartesian (hyperbolic-sector) third-order field, for comparison."""
    P = sum((np.asarray(momenta[a], dtype=float) for a in triple), np.zeros(3))
    n, wsky = sphere(nth, nph)
    w = 7.5 * np.einsum("si,ij,sj->s", n, C, n)
    sp_ = Spectral(e, o, nt)
    tt = sp_.t
    ph = np.exp(1j * (n @ P)[:, None] * (o - tt)[None, :])
    N = n[:, :, None] * n[:, None, :] - np.eye(3) / 3
    CP = np.einsum("s,st,sij->tij", wsky * w, ph, N)
    Ap = profiles[frozenset(triple)].Ap(tt)
    return -0.5 * np.dot(sp_.full, np.einsum("tij,tij->t", CP, Ap))
