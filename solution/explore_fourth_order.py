"""EXPLORATORY numerical evaluation of the second-order part of the fourth-order spectrum.

NOT EVIDENCE.  The problem forbids numerics in place of derivation; this script is a guide for
finding structure, to be proved analytically afterwards.  It uses the exact second-order sources
(second_order.py) and the triangular Green-function formulas of the report, integrated by
quadrature, and the map formulas of the report's section "The observable through third order".

Computed on the monochromatic shell (k0 = 1), for l = 2:
  C11            first-order spectrum (checked against the entire-function T_2)
  N_l            lognormal term, from T_l (closed form, Thm. log-map)
  (1/4) C[L2,L2] split into the local part and the A_2 (resonant) part
  M_l            collapsed log-map bispectrum
  (1/4) C^(22)   full, and the internal check (1/4) C[theta^2, theta^2] = first part of N_l
Missing: (1/3) Re C[theta, L3] (needs third order).

Validation modes:
  python explore_fourth_order.py ward    dilation Ward identity K2(k, q->0) = -k d_k K1 at finite k
  python explore_fourth_order.py scan    spectrum pieces versus x_o
"""
import sys
import pickle
import os
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'results', 'second_order_sources_numeric.pkl')
s3 = np.sqrt(3.0)


# ------------------------------------------------------------------ exact sources -> numpy
def load_sources():
    if os.path.exists(CACHE):
        out = pickle.load(open(CACHE, 'rb'))
    else:
        sys.path.insert(0, HERE)
        import second_order as so
        src = so.sources()
        unk = src['unk']
        fnl = [unk['n']] + [e for e in set(unk['A']) if isinstance(e, sp.Function)]
        full = src['Ham_full']
        cH = full.subs({d: 0 for d in full.atoms(sp.Derivative)}).subs({f: 0 for f in fnl})
        out = dict(N=src['N'], M=list(src['M']), T=src['T'], H=sp.expand(cH))
        pickle.dump(out, open(CACHE, 'wb'))
    def f(e):
        e = sp.sympify(e)
        by_name = {s.name: s for s in e.free_symbols}
        syms = [by_name.get(nm, sp.Symbol(nm)) for nm in ('t', 'p', 'q', 'b1', 'b2')]
        return sp.lambdify(syms, e, 'numpy')
    return dict(N=f(out['N']), H=f(out['H']), M=[f(e) for e in out['M']],
                T=[[f(out['T'][i, j]) for j in range(3)] for i in range(3)])


SRC = None


def cumint(y, tg):
    """Cumulative trapezoid from tg[0]; y shape (..., len(tg))."""
    dt = np.diff(tg)
    inc = 0.5 * (y[..., 1:] + y[..., :-1]) * dt
    return np.concatenate([np.zeros(y.shape[:-1] + (1,), dtype=inc.dtype), np.cumsum(inc, axis=-1)], axis=-1)


def second_order_A2(p, q, b1, b2, t_o, ng=12000):
    """A_2(t) (3x3) and its derivative on a dense grid (0, t_o]; frame k1=(p,0,0), k2=(b1,b2,0)."""
    global SRC
    if SRC is None:
        SRC = load_sources()
    tg = np.linspace(t_o * 2e-4, t_o, ng)
    KV = np.array([p + b1, b2, 0.0])
    K = np.linalg.norm(KV)
    c = K / s3
    x = c * tg
    u1 = np.sin(x) / x - np.cos(x)
    u2 = -np.cos(x) / x - np.sin(x)
    ev = lambda fn: np.broadcast_to(np.asarray(fn(tg, p, q, b1, b2), dtype=complex), tg.shape)
    cN = ev(SRC['N'])
    cH = ev(SRC['H'])
    P1 = cumint(u1 * cN, tg)
    P2 = cumint(u2 * cN, tg)
    c2 = (2 * (p**2 + q**2) - p * b1) / 9
    nu = 9 * c2 / K**2 * u1 - (u2 * P1 - u1 * P2) / (3 * c)
    u1p = np.cos(x) / x - np.sin(x) / x**2 + np.sin(x)
    u2p = np.sin(x) / x + np.cos(x) / x**2 - np.cos(x)
    nup = 9 * c2 / K**2 * c * u1p - (u2p * P1 - u1p * P2) / 3
    hS = 3 / (2 * K**2) * (4 * K**2 * nu + 12 * nup / tg + 12 * nu / tg**2 + cH)
    Kh = KV / K
    Pp = np.eye(3) - np.outer(Kh, Kh)
    cM = np.array([ev(fn) for fn in SRC['M']])                      # (3, ng)
    w = 2j / K**2 * cumint(Pp @ cM, tg)                             # (3, ng)
    cT = np.array([[ev(SRC['T'][i][j]) for j in range(3)] for i in range(3)])   # (3,3,ng)
    PMP = np.einsum('ia,abn,bj->ijn', Pp, cT, Pp)
    TT = PMP - Pp[:, :, None] * np.einsum('aan->n', PMP)[None, None, :] / 2
    Ic = cumint(tg * np.cos(K * tg) * TT, tg)
    Is = cumint(tg * np.sin(K * tg) * TT, tg)
    h = -2 / (K * tg) * (np.sin(K * tg) * Ic - np.cos(K * tg) * Is)
    A = (np.outer(Kh, Kh) - np.eye(3) / 3)[:, :, None] * hS + KV[:, None, None] * w[None, :, :] \
        + w[:, None, :] * KV[None, :, None] + h
    Ap = np.gradient(A, tg, axis=-1)
    return tg, A, Ap


# ------------------------------------------------------------------ first-order mode functions
def mode(k, t):
    x = k * t / s3
    s = 2 * np.sin(x) / x - np.cos(x)
    sp_ = np.sin(x) + 2 * np.cos(x) / x - 2 * np.sin(x) / x**2
    Q = -6 * (s - 1) / k**2
    Qp = -6 * sp_ * (k / s3) / k**2
    nu = -(np.sin(x) / x - np.cos(x))
    return Q, Qp, nu


def sphere(nth=40, nph=80):
    mu, wmu = np.polynomial.legendre.leggauss(nth)
    ph = 2 * np.pi * np.arange(nph) / nph
    MU, PH = np.meshgrid(mu, ph, indexing='ij')
    st = np.sqrt(1 - MU**2)
    n = np.stack([st * np.cos(PH), st * np.sin(PH), MU], axis=-1).reshape(-1, 3)
    w = (wmu[:, None] * np.full(nph, 2 * np.pi / nph)[None, :]).reshape(-1)
    return n, w


def smean(f, w):
    return np.tensordot(w, f, axes=(0, 0)) / (4 * np.pi)


class Setup:
    def __init__(self, x_o, R=3.0, RI=0.5, nt=320, nth=40, nph=80):
        self.t_o = s3 * x_o
        self.ts = self.t_o / R
        self.tI = self.ts * RI
        self.tr = np.linspace(self.ts, self.t_o, nt)
        self.lam = self.t_o - self.tr
        self.n, self.w = sphere(nth, nph)

    def first(self, kv):
        """First-order quantities for mode with wavevector kv along the ray grid."""
        k = np.linalg.norm(kv)
        n, tr, lam = self.n, self.tr, self.lam
        Q, Qp, nu = mode(k, tr)
        kn = n @ kv                                                   # (Ns,)
        ph = np.exp(1j * kn[:, None] * lam[None, :])                  # (Ns, Nt)
        Mn = np.outer(kn, kv) - k**2 * n / 3                          # (Ns, 3)  M n, M = kk - delta k^2/3
        nMn = kn**2 - k**2 / 3
        snn = 0.5 * Qp[None, :] * nMn[:, None]                        # sigma_nn (no phase)
        # ray acceleration A = H' n - Gamma(n,n) - n sigma_nn
        Gam = 1j * (2 * kn[:, None, None] * n[:, None, :] - kv[None, None, :]) * np.ones_like(Q)[None, :, None] \
            + 0.5j * (2 * kn[:, None, None] * Q[None, :, None] * Mn[:, None, :]
                      - kv[None, None, :] * (Q[None, :] * nMn[:, None])[:, :, None])
        Acc = (Qp[None, :, None] * Mn[:, None, :] - Gam - n[:, None, :] * snn[:, :, None]) * ph[:, :, None]
        Q_o, _, _ = mode(k, np.array([self.t_o]))
        e_o = n + 0.5 * Q_o[0] * Mn                                   # (Ns, 3)
        # reverse cumulative integrals from t to t_o
        tot0 = cumint(np.moveaxis(Acc, 1, -1), tr)                    # (Ns,3,Nt) from ts
        tot1 = cumint(np.moveaxis(Acc * tr[None, :, None], 1, -1), tr)
        R0 = (tot0[..., -1:] - tot0)                                  # int_t^{t_o} A
        R1 = (tot1[..., -1:] - tot1)
        R0 = np.moveaxis(R0, -1, 1)
        R1 = np.moveaxis(R1, -1, 1)
        X = e_o[:, None, :] * (tr - self.t_o)[None, :, None] + R1 - tr[None, :, None] * R0
        e = e_o[:, None, :] - R0
        sig_n = 0.5 * Qp[None, :, None] * Mn[:, None, :]              # sigma n (no phase)
        I1 = np.trapezoid(snn * ph, tr, axis=1)
        dI1 = snn[:, -1] * ph[:, -1] + np.trapezoid(snn * ph * (1j * kn[:, None]), tr, axis=1)
        theta = -(I1 - smean(I1, self.w))
        dtheta = -(dI1 - smean(dI1, self.w))
        # clock
        tg = np.linspace(1e-6, self.t_o, 20000)
        _, _, om = mode(k, tg)
        W = cumint(tg * om, tg)
        Wo, WI = W[-1], np.interp(self.tI, tg, W)
        omI = mode(k, np.array([self.tI]))[2][0]
        dI = -self.tI * omI
        D = (self.tI * dI - Wo + WI) / self.t_o
        _, _, nus = mode(k, np.array([self.ts]))
        d_star = -self.ts * nus[0] * ph[:, 0]
        return dict(k=k, kv=kv, kn=kn, ph=ph, Q=Q, Qp=Qp, Mn=Mn, nMn=nMn, snn=snn, sig_n=sig_n,
                    X=X, e=e, theta=theta, dtheta=dtheta, D=D, d_star=d_star, I1=I1)

    def second(self, a, b, A2data, drop=()):
        """Second-order map kernels for modes a, b (dicts from first())."""
        tg, A, Ap = A2data
        tr = self.tr
        Apr = np.array([[np.interp(tr, tg, Ap[i, j].real) + 1j * np.interp(tr, tg, Ap[i, j].imag)
                         for j in range(3)] for i in range(3)])     # (3,3,Nt)
        n = self.n
        phab = a['ph'] * b['ph']
        A2nn = np.einsum('si,ijt,sj->st', n, Apr, n)
        Qa, Qpa, Qb, Qpb = a['Q'], a['Qp'], b['Q'], b['Qp']
        MaMb = np.einsum('si,si->s', a['Mn'], b['Mn'])
        sig_loc = 0.5 * (Qpa * Qb + Qa * Qpb)[None, :] * MaMb[:, None] \
            + Qpa[None, :] * a['nMn'][:, None] + Qpb[None, :] * b['nMn'][:, None]
        integ_loc = sig_loc * phab \
            + 1j * np.einsum('i,sti->st', b['kv'], a['X']) * b['snn'] * b['ph'] \
            + 1j * np.einsum('i,sti->st', a['kv'], b['X']) * a['snn'] * a['ph'] \
            - 2 * np.einsum('sti,sti->st', b['sig_n'] * np.ones_like(tr)[None, :, None], a['e']) * b['ph'] \
            - 2 * np.einsum('sti,sti->st', a['sig_n'] * np.ones_like(tr)[None, :, None], b['e']) * a['ph']
        endpoint = a['snn'][:, 0] * a['ph'][:, 0] * b['d_star'] + b['snn'][:, 0] * b['ph'][:, 0] * a['d_star']
        if 'endpoint' in drop:
            endpoint = 0 * endpoint
        if 'ray' in drop:
            integ_loc = sig_loc * phab
        I_loc = np.trapezoid(integ_loc, tr, axis=1) - endpoint
        I_A = np.trapezoid(0.5 * A2nn * phab, tr, axis=1)
        w = self.w
        KA = -(I_A - smean(I_A, w))
        clock = 0 if 'clock' in drop else a['D'] * b['dtheta'] + b['D'] * a['dtheta']
        Kloc = -(I_loc - smean(I_loc, w)) + clock
        if 'A2' in drop:
            KA = 0 * KA
        KL = Kloc + KA
        K11 = a['theta'] * b['theta']
        return dict(KL=KL, KA=KA, Kloc=Kloc, K11=K11, Theta2=KL + K11)


_PCACHE = {}


def S_l(f, g, n, w, l=2):
    """sum_m f_lm g_lm^* = (2l+1)/(4pi) int int f(n) g(n')^* P_l(n.n')."""
    key = (id(n), l)
    if key not in _PCACHE:
        _PCACHE[key] = (2 * l + 1) / (4 * np.pi) * np.polynomial.legendre.legval(n @ n.T, [0] * l + [1])
    P = _PCACHE[key]
    return (w * f) @ P @ (w * np.conj(g))


# ------------------------------------------------------------------ modes
def ward(x_o=1.5, qs=(0.04, 0.02), drop=()):
    """Dilation Ward identity of the observable at finite k: [K2(k,q)+K2(k,-q)]/2 -> -k d_k K1."""
    S = Setup(x_o, nt=400, nth=24, nph=48)
    k = 1.0
    h = 1e-4
    th = lambda kk: S.first(np.array([kk, 0, 0]))['theta']
    target = -k * (th(k + h) - th(k - h)) / (2 * h)
    for q in qs:
        res = []
        for sgn in (+1, -1):
            qv = sgn * q * np.array([0.6, 0.8, 0.0])
            a, b = S.first(np.array([k, 0, 0])), S.first(qv)
            A2 = second_order_A2(k, q, qv[0], qv[1], S.t_o)
            res.append(S.second(a, b, A2, drop)['Theta2'])
        avg = 0.5 * (res[0] + res[1])
        err = np.sqrt(smean(np.abs(avg - target)**2, S.w) / smean(np.abs(target)**2, S.w))
        print(f'  ward x_o={x_o} q={q} drop={drop}: relative L2 error of [K2(k,q)+K2(k,-q)]/2 vs -k dK1/dk = {err:.2e}')


def T_l_closed(l, xs, xo):
    sys.path.insert(0, HERE)
    import transfer_entire as te
    V, co = te.coefficients(l)
    return te.numeric(l, co, xs, xo)


def scan(xos=(0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0), nmu=24, l=2, lmax=None):
    from sympy.physics.wigner import wigner_3j
    mus, wmus = np.polynomial.legendre.leggauss(nmu)
    rows = []
    for xo in xos:
        S = Setup(xo, nth=max(24, int(8 + 4 * xo)), nph=max(48, int(16 + 8 * xo)))
        n, w = S.n, S.w
        a = S.first(np.array([1.0, 0, 0]))
        acc = dict(LL=0, LL_loc=0, LL_A=0, LL_x=0, LT=0, TLT=0, TT=0, Th=0)
        for mu, wm in zip(mus, wmus):
            qv = np.array([mu, np.sqrt(1 - mu**2), 0.0])
            b = S.first(qv)
            A2 = second_order_A2(1.0, 1.0, mu, np.sqrt(1 - mu**2), S.t_o)
            k2 = S.second(a, b, A2)
            acc['LL'] += wm * S_l(k2['KL'], k2['KL'], n, w, l).real
            acc['LL_loc'] += wm * S_l(k2['Kloc'], k2['Kloc'], n, w, l).real
            acc['LL_A'] += wm * S_l(k2['KA'], k2['KA'], n, w, l).real
            acc['LL_x'] += wm * 2 * S_l(k2['Kloc'], k2['KA'], n, w, l).real
            acc['LT'] += wm * S_l(k2['KL'], k2['K11'], n, w, l).real
            acc['TLT'] += wm * S_l(a['theta'], np.conj(b['theta']) * k2['KL'], n, w, l).real
            acc['TT'] += wm * S_l(k2['K11'], k2['K11'], n, w, l).real
            acc['Th'] += wm * S_l(k2['Theta2'], k2['Theta2'], n, w, l).real
        pre = 1.0 / (2 * l + 1)         # sigma^4 = 1
        C11 = S_l(a['theta'], a['theta'], n, w, l).real / (2 * l + 1)   # sigma^2 = 1
        T2_pipeline = np.sqrt(C11 / (4 * np.pi))
        T2_closed = abs(float(T_l_closed(l, xo / 3.0, xo)))
        # lognormal term from T_l of the pipeline (first part) and the closed form (-4 T_l^4 part)
        Ts = {}
        lmax_ = lmax or int(6 + 3 * xo)
        for L in range(1, lmax_ + 1):
            Ts[L] = S_l(a['theta'], a['theta'], n, w, L).real / (4 * np.pi * (2 * L + 1))
        conv = sum((2 * l1 + 1) * (2 * l2 + 1) * float(wigner_3j(l1, l2, l, 0, 0, 0))**2 * Ts[l1] * Ts[l2]
                   for l1 in Ts for l2 in Ts)
        N = 2 * np.pi * (conv - 4 * Ts[l]**2)
        quarterTT = 0.25 * pre * acc['TT']
        LL = 0.25 * pre * acc['LL']
        M = 0.5 * pre * acc['LT'] + pre * acc['TLT']
        row = dict(xo=xo, C11=C11, T2_pipe=T2_pipeline, T2_closed=T2_closed,
                   N=N / C11, quarterTT_check=quarterTT / (2 * np.pi * conv),
                   LL=LL / C11, LL_loc=0.25 * pre * acc['LL_loc'] / C11, LL_A=0.25 * pre * acc['LL_A'] / C11,
                   LL_x=0.25 * pre * acc['LL_x'] / C11, M=M / C11, C22q=0.25 * pre * acc['Th'] / C11)
        rows.append(row)
        print('  ' + '  '.join(f'{k}={v:.5g}' for k, v in row.items()), flush=True)
    return rows


if __name__ == '__main__':
    mode_ = sys.argv[1] if len(sys.argv) > 1 else 'scan'
    if mode_ == 'ward':
        ward()
        for d in ('clock', 'endpoint', 'ray', 'A2'):
            ward(qs=(0.02,), drop=(d,))
    else:
        scan()
