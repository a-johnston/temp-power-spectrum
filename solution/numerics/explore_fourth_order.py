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

Modes:
  python explore_fourth_order.py ward     dilation Ward identity K2(k, q->0) = -k d_k K1 at finite k,
                                          with negative controls (drop clock/endpoint/ray/A2)
  python explore_fourth_order.py harden   analytic vs finite-difference A2', early-time cutoff,
                                          moving-endpoint term vs direct root of the emission
                                          condition, ray-grid convergence of T_2 at x_o = 12
  python explore_fourth_order.py scan     spectrum pieces versus x_o (fixed ratio x* = x_o/3);
                                          scan(..., xstar=1.0) follows fixed emission, l=3 allowed.
At a zero of T_l(k0) the C^(13) term vanishes exactly, so the scan's quarter C22 is then the
complete fourth-order spectrum (report Prop. "nodes").
"""
import sys
import pickle
import os
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))           # solution/numerics
SOL = os.path.dirname(HERE)                                    # solution/: exact modules and results/
CACHE = os.path.join(SOL, 'results', 'second_order_sources_numeric.pkl')
s3 = np.sqrt(3.0)


# ------------------------------------------------------------------ exact sources -> numpy
def load_sources():
    if os.path.exists(CACHE):
        out = pickle.load(open(CACHE, 'rb'))
    else:
        sys.path.insert(0, SOL)
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
    tsym = [x for x in sp.sympify(out['H']).free_symbols if x.name == 't'][0]
    return dict(N=f(out['N']), H=f(out['H']), Hp=f(sp.diff(out['H'], tsym)), M=[f(e) for e in out['M']],
                T=[[f(out['T'][i, j]) for j in range(3)] for i in range(3)])


SRC = None


def cumint(y, tg):
    """Cumulative trapezoid from tg[0]; y shape (..., len(tg))."""
    dt = np.diff(tg)
    inc = 0.5 * (y[..., 1:] + y[..., :-1]) * dt
    return np.concatenate([np.zeros(y.shape[:-1] + (1,), dtype=inc.dtype), np.cumsum(inc, axis=-1)], axis=-1)


CFG = dict(ng=12000, lo=2e-4, analytic=True)


def second_order_A2(p, q, b1, b2, t_o, ng=None, lo=None, analytic=None):
    """A_2(t) (3x3) and its derivative on a dense grid [lo*t_o, t_o]; frame k1=(p,0,0), k2=(b1,b2,0).

    The derivative is analytic in the Green representation (analytic=True): the integrals'
    boundary terms cancel, nu'' comes from the lapse equation, and c_Ham' is differentiated
    symbolically.  analytic=False uses np.gradient (kept as a cross-check)."""
    ng = ng or CFG['ng']
    lo = CFG['lo'] if lo is None else lo
    analytic = CFG['analytic'] if analytic is None else analytic
    global SRC
    if SRC is None:
        SRC = load_sources()
    tg = np.linspace(t_o * lo, t_o, ng)
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
    if not analytic:
        return tg, A, np.gradient(A, tg, axis=-1)
    cHp = ev(SRC['Hp'])
    nupp = (-K**2 * nu + 6 * nu / tg**2 - cN) / 3                    # lapse equation
    hSp = 3 / (2 * K**2) * (4 * K**2 * nup + 12 * nupp / tg - 24 * nu / tg**3 + cHp)
    wp = 2j / K**2 * (Pp @ cM)
    hp = 2 / (K * tg**2) * (np.sin(K * tg) * Ic - np.cos(K * tg) * Is) - 2 / tg * (np.cos(K * tg) * Ic + np.sin(K * tg) * Is)
    Ap = (np.outer(Kh, Kh) - np.eye(3) / 3)[:, :, None] * hSp + KV[:, None, None] * wp[None, :, :] \
        + wp[:, None, :] * KV[None, :, None] + hp
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
    key = (n.shape, hash(n.tobytes()), l)      # content key: id() can be reused after garbage collection
    if key not in _PCACHE:
        if len(_PCACHE) >= 2:                    # each entry is a dense Ns x Ns matrix: keep memory bounded
            _PCACHE.clear()
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
    sys.path.insert(0, SOL)
    import transfer_entire as te
    V, co = te.coefficients(l)
    return te.numeric(l, co, xs, xo)


def scan(xos=(0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0), nmu=24, l=2, lmax=None, verbose=True, R=3.0, xstar=None, RI=0.5, nt=None):
    from sympy.physics.wigner import wigner_3j
    mus, wmus = np.polynomial.legendre.leggauss(nmu)
    rows = []
    for xo in xos:
        Rx = R if xstar is None else xo / xstar
        S = Setup(xo, R=Rx, RI=RI, nt=nt or max(400, int(160 * xo)), nth=max(24, int(8 + 4 * xo)), nph=max(48, int(16 + 8 * xo)))
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
        T2_closed = abs(float(T_l_closed(l, xo / Rx, xo)))
        # lognormal term from T_l of the pipeline (first part) and the closed form (-4 T_l^4 part)
        Ts = {}
        lmax_ = lmax or int(6 + 3 * xo)
        # theta for a single mode is axisymmetric about k = x: T_L = (1/(4 pi i^L)) int dOmega theta P_L(n_x)
        for L in range(1, lmax_ + 1):
            PL = np.polynomial.legendre.legval(n[:, 0], [0] * L + [1])
            TL = np.sum(w * a['theta'] * PL) / (4 * np.pi * 1j**L)
            Ts[L] = abs(TL)**2
        conv = sum((2 * l1 + 1) * (2 * l2 + 1) * float(wigner_3j(l1, l2, l, 0, 0, 0))**2 * Ts[l1] * Ts[l2]
                   for l1 in Ts for l2 in Ts)
        N = 2 * np.pi * (conv - 4 * Ts[l]**2)
        quarterTT = 0.25 * pre * acc['TT']
        LL = 0.25 * pre * acc['LL']
        M = 0.5 * pre * acc['LT'] + pre * acc['TLT']
        row = dict(xo=xo, C11=C11, T2_pipe=T2_pipeline, T2_closed=T2_closed,
                   N=N / C11, quarterTT_check=quarterTT / (2 * np.pi * conv),
                   LL=LL / C11, LL_loc=0.25 * pre * acc['LL_loc'] / C11, LL_A=0.25 * pre * acc['LL_A'] / C11,
                   LL_x=0.25 * pre * acc['LL_x'] / C11, M=M / C11, C22q=0.25 * pre * acc['Th'] / C11,
                   C22q_abs=0.25 * pre * acc['Th'], LL_abs=LL, M_abs=M, N_abs=N,
                   LT_abs=0.5 * pre * acc['LT'], TLT_abs=pre * acc['TLT'])
        rows.append(row)
        if verbose:
            print('  ' + '  '.join(f'{k}={v:.5g}' for k, v in row.items()), flush=True)
    return rows


def harden():
    """Pipeline diagnostics: analytic vs finite-difference A2', cutoff, moving endpoint, ray-grid convergence."""
    # (a) analytic A2' versus finite differences
    for mu in (-0.9, -1 / 3, 0.5):
        b2 = np.sqrt(1 - mu**2)
        tg, A, Ap = second_order_A2(1.0, 1.0, mu, b2, s3 * 3.0, analytic=True)
        _, _, Ag = second_order_A2(1.0, 1.0, mu, b2, s3 * 3.0, analytic=False)
        sel = tg > tg[-1] / 3
        err = np.max(np.abs(Ap - Ag)[..., sel]) / np.max(np.abs(Ap)[..., sel])
        print(f"  (a) mu={mu:+.3f}: max |A2'(analytic) - A2'(gradient)| / max|A2'| on [t*, t_o] = {err:.1e}")
    # (b) early-time cutoff of the Green integrals, independently varied
    for lo in (1e-3, 2e-4, 2e-5):
        CFG['lo'] = lo
        r = scan(xos=(3.0,), nmu=16, verbose=False)[0]
        print(f"  (b) cutoff lo={lo:.0e}: quarter C22/C11 = {r['C22q']:.7g}, LL_A = {r['LL_A']:.7g}")
    CFG['lo'] = 2e-4
    # (c) moving emission endpoint against a direct root of ln t + nu = ln t*
    S = Setup(2.0, nt=200, nth=6, nph=12)
    ka, kb = np.array([1.0, 0, 0]), np.array([0.3, 0.8, 0.0])
    code = 0
    for sa in (1, -1):
        for sb in (1, -1):
            a, b = S.first(sa * ka), S.first(sb * kb)
            code = code + 0.25 * (a['snn'][:, 0] * a['ph'][:, 0] * b['d_star'] + b['snn'][:, 0] * b['ph'][:, 0] * a['d_star'])
    for i in (3, 17, 40):
        n = S.n[i]
        def fields(k, t):
            Q, Qp, nu = mode(np.linalg.norm(k), t)
            kn = n @ k
            M = np.outer(k, k) - np.dot(k, k) * np.eye(3) / 3
            return 0.5 * Qp * (n @ M @ n) * np.cos(kn * (S.t_o - t)), nu * np.cos(kn * (S.t_o - t))
        def I(za, zb):
            phi = lambda t: za * fields(ka, t)[1] + zb * fields(kb, t)[1]
            te = S.ts
            for _ in range(80):
                te = S.ts * np.exp(-phi(np.array([te]))[0])
            tt = np.linspace(te, S.t_o, 40001)
            return np.trapezoid(za * fields(ka, tt)[0] + zb * fields(kb, tt)[0], tt)
        h = 1e-4
        mixed = (I(h, h) - I(h, -h) - I(-h, h) + I(-h, -h)) / (4 * h * h)
        print(f"  (c) n#{i}: mixed z_a z_b coefficient from root-solving = {mixed:.8f}; code endpoint term = {-code[i].real:.8f}")
    # (d) the x_o = 12 first-order mismatch versus ray and sphere resolution
    T_cl = abs(float(T_l_closed(2, 4.0, 12.0)))
    for nt, nth in ((320, 56), (960, 56), (320, 96), (960, 96)):
        S = Setup(12.0, nt=nt, nth=nth, nph=2 * nth)
        a = S.first(np.array([1.0, 0, 0]))
        T = np.sqrt(S_l(a['theta'], a['theta'], S.n, S.w, 2).real / (4 * np.pi * 5))
        print(f"  (d) x_o=12 nt={nt} nth={nth}: T2 pipeline {T:.7f}  closed {T_cl:.7f}  rel {abs(T - T_cl) / T_cl:.1e}")


if __name__ == '__main__':
    mode_ = sys.argv[1] if len(sys.argv) > 1 else 'scan'
    if mode_ == 'harden':
        harden()
    elif mode_ == 'ward':
        ward()
        for d in ('clock', 'endpoint', 'ray', 'A2'):
            ward(qs=(0.02,), drop=(d,))
    else:
        scan()


# ------------------------------------------------------------------ amplitude-level recombination
def kernel_pieces(S, a, b, A2data):
    """Pieces of the second-order log-map kernel K^L(n) (l >= 1; sky-mean subtraction omitted).

    Base split:  loc (product terms of sigma^(2)), A (A_2' term), ray, end (moving emission),
    clk (observation clock).  Metric variable gamma~ = A_2 + (H_a H_b + H_b H_a)/2 + 2(H_a + H_b),
    with sigma^(2) = (1/2) d_t gamma~, and its integration by parts
        -int (1/2) gamma~'_nn e^{iK.n(t_o-t)} = g_obs + g_emit + g_bulk,
        g_obs = -(1/2) gamma~_nn(t_o),  g_emit = +(1/2) gamma~_nn(t*) e^{iK.n lam*},
        g_bulk = -(i/2) K.n int gamma~_nn e^{iK.n(t_o-t)} dt.
    """
    tg, A, Ap = A2data
    tr, n = S.tr, S.n
    ip = lambda M: np.array([[np.interp(tr, tg, M[i, j].real) + 1j * np.interp(tr, tg, M[i, j].imag)
                              for j in range(3)] for i in range(3)])
    Ar, Apr = ip(A), ip(Ap)
    phab = a['ph'] * b['ph']
    Qa, Qpa, Qb, Qpb = a['Q'], a['Qp'], b['Q'], b['Qp']
    MaMb = np.einsum('si,si->s', a['Mn'], b['Mn'])
    nn = lambda M: np.einsum('si,ijt,sj->st', n, M, n)
    sig_loc = 0.5 * (Qpa * Qb + Qa * Qpb)[None, :] * MaMb[:, None] + Qpa[None, :] * a['nMn'][:, None] + Qpb[None, :] * b['nMn'][:, None]
    gam_loc = (Qa * Qb)[None, :] * MaMb[:, None] + 2 * (Qa[None, :] * a['nMn'][:, None] + Qb[None, :] * b['nMn'][:, None])
    gam = nn(Ar) + gam_loc                                        # gamma~_nn (no phase)
    ray = (1j * np.einsum('i,sti->st', b['kv'], a['X']) * b['snn'] * b['ph']
           + 1j * np.einsum('i,sti->st', a['kv'], b['X']) * a['snn'] * a['ph']
           - 2 * np.einsum('sti,sti->st', b['sig_n'], a['e']) * b['ph']
           - 2 * np.einsum('sti,sti->st', a['sig_n'], b['e']) * a['ph'])
    endpoint = a['snn'][:, 0] * a['ph'][:, 0] * b['d_star'] + b['snn'][:, 0] * b['ph'][:, 0] * a['d_star']
    Kn = n @ (a['kv'] + b['kv'])
    sig_zeta = Qpa[None, :] * a['nMn'][:, None] + Qpb[None, :] * b['nMn'][:, None]     # from e^{2 zeta}
    rayX = (1j * np.einsum('i,sti->st', b['kv'], a['X']) * b['snn'] * b['ph']
            + 1j * np.einsum('i,sti->st', a['kv'], b['X']) * a['snn'] * a['ph'])
    # first-order fluid displacement removing the scalar part of H_1:  Psi_a = -(i/2) Q_a k_a e^{i k_a.Xbar}
    Psi_a = -0.5j * Qa[None, :, None] * a['kv'][None, None, :] * a['ph'][:, :, None]
    Psi_b = -0.5j * Qb[None, :, None] * b['kv'][None, None, :] * b['ph'][:, :, None]
    shiftX = (1j * np.einsum('i,sti->st', b['kv'], Psi_a) * b['snn'] * b['ph']
              + 1j * np.einsum('i,sti->st', a['kv'], Psi_b) * a['snn'] * a['ph'])
    # full first-order relabelling Y = X + Psi:  e = (1 - dPsi)(e^Y - v),  v = d_t Psi
    v_a = -0.5j * Qpa[None, :, None] * a['kv'][None, None, :] * a['ph'][:, :, None]
    v_b = -0.5j * Qpb[None, :, None] * b['kv'][None, None, :] * b['ph'][:, :, None]
    kna, knb = (n @ a['kv'])[:, None], (n @ b['kv'])[:, None]
    sn = lambda sig_n, vec: np.einsum('sti,sti->st', sig_n, vec)            # sigma(n, vec), no phase
    jac = (-2j * kna * sn(b['sig_n'], Psi_a) * b['ph'] - 2j * knb * sn(a['sig_n'], Psi_b) * a['ph'])
    vel = (2 * sn(b['sig_n'], v_a) * b['ph'] + 2 * sn(a['sig_n'], v_b) * a['ph'])
    dirE = (ray - rayX) - vel - jac
    P = dict(
        jac=-np.trapezoid(jac, tr, axis=1),
        vel=-np.trapezoid(vel, tr, axis=1),
        dirE=-np.trapezoid(dirE, tr, axis=1),
        dispE=-np.trapezoid(rayX + shiftX, tr, axis=1),        # photon displacement relative to the fluid
        fluidshift=+np.trapezoid(shiftX, tr, axis=1),          # minus the fluid's own displacement
        loc_zeta=-np.trapezoid(sig_zeta * phab, tr, axis=1),
        loc_HH=-np.trapezoid((sig_loc - sig_zeta) * phab, tr, axis=1),
        rayX=-np.trapezoid(rayX, tr, axis=1),
        raydir=-np.trapezoid(ray - rayX, tr, axis=1),
        loc=-np.trapezoid(sig_loc * phab, tr, axis=1),
        A=-np.trapezoid(0.5 * nn(Apr) * phab, tr, axis=1),
        ray=-np.trapezoid(ray, tr, axis=1),
        end=+endpoint,
        clk=a['D'] * b['dtheta'] + b['D'] * a['dtheta'],
        g_obs=-0.5 * gam[:, -1] * phab[:, -1],
        g_emit=+0.5 * gam[:, 0] * phab[:, 0],
        g_bulk=-0.5j * Kn * np.trapezoid(gam * phab, tr, axis=1),
    )
    return P


def recombination(xos=(1.0, 3.0, 6.0, 12.0), nmu=16, l=2):
    names = ['loc', 'A', 'ray', 'end', 'clk', 'g_obs', 'g_emit', 'g_bulk', 'loc_zeta', 'loc_HH', 'rayX', 'raydir',
             'dispE', 'fluidshift', 'jac', 'vel', 'dirE']
    splits = {
        'local | A2 | ray | end | clock': [['loc'], ['A'], ['ray'], ['end'], ['clk']],
        'gamma~ (=loc+A) | ray | end | clock': [['loc', 'A'], ['ray'], ['end'], ['clk']],
        'IBP: obs | emit | bulk | ray | end | clock': [['g_obs'], ['g_emit'], ['g_bulk'], ['ray'], ['end'], ['clk']],
        'obs | emit+end | bulk+ray | clock': [['g_obs'], ['g_emit', 'end'], ['g_bulk', 'ray'], ['clk']],
        'obs+clock | emit+end | bulk+ray': [['g_obs', 'clk'], ['g_emit', 'end'], ['g_bulk', 'ray']],
        'observer (obs+clock) | rest': [['g_obs', 'clk'], ['g_emit', 'end', 'g_bulk', 'ray']],
        'A+HH | zeta+ray | end | clock': [['A', 'loc_HH'], ['loc_zeta', 'ray'], ['end'], ['clk']],
        'A+HH | zeta+dir | displacement | end | clock': [['A', 'loc_HH'], ['loc_zeta', 'raydir'], ['rayX'], ['end'], ['clk']],
        'A | HH | zeta | dir | displacement | end | clock': [['A'], ['loc_HH'], ['loc_zeta'], ['raydir'], ['rayX'], ['end'], ['clk']],
        'A+HH+zeta+dir | displacement | end | clock': [['A', 'loc_HH', 'loc_zeta', 'raydir'], ['rayX'], ['end'], ['clk']],
        'EULER: A+HH+zeta+dir-Psi.dsigma | (X+Psi).dsigma | end | clock':
            [['A', 'loc_HH', 'loc_zeta', 'raydir', 'fluidshift'], ['dispE'], ['end'], ['clk']],
        'FULL EULER: fieldE | velocity | dispE | dirE+zeta | end | clock':
            [['A', 'loc_HH', 'fluidshift', 'jac'], ['vel'], ['dispE'], ['dirE', 'loc_zeta'], ['end'], ['clk']],
        'FULL EULER: fieldE+velocity | dispE | dirE+zeta | end | clock':
            [['A', 'loc_HH', 'fluidshift', 'jac', 'vel'], ['dispE'], ['dirE', 'loc_zeta'], ['end'], ['clk']],
        'FULL EULER: fieldE+velocity+dirE+zeta | dispE | end | clock':
            [['A', 'loc_HH', 'fluidshift', 'jac', 'vel', 'dirE', 'loc_zeta'], ['dispE'], ['end'], ['clk']],
        'EULER: gamma~-Psi.dsigma | dir | (X+Psi).dsigma | end | clock':
            [['loc', 'A', 'fluidshift'], ['raydir'], ['dispE'], ['end'], ['clk']],
    }
    mus, wm = np.polynomial.legendre.leggauss(nmu)
    for xo in xos:
        S = Setup(xo, nt=max(400, int(160 * xo)), nth=max(24, int(8 + 4 * xo)), nph=max(48, int(16 + 8 * xo)))
        a = S.first(np.array([1.0, 0, 0]))
        G = np.zeros((len(names), len(names)))
        idcheck = 0
        for mu, w in zip(mus, wm):
            qv = np.array([mu, np.sqrt(1 - mu**2), 0.0])
            b = S.first(qv)
            P = kernel_pieces(S, a, b, second_order_A2(1.0, 1.0, mu, np.sqrt(1 - mu**2), S.t_o))
            lhs, rhs = P['loc'] + P['A'], P['g_obs'] + P['g_emit'] + P['g_bulk']
            idcheck = max(idcheck, np.max(np.abs(lhs - rhs)) / np.max(np.abs(lhs)))
            for i, ni in enumerate(names):
                for j, nj in enumerate(names):
                    if j >= i:
                        G[i, j] += w * S_l(P[ni], P[nj], S.n, S.w, l).real
        G = G + np.triu(G, 1).T
        idx = {nm: i for i, nm in enumerate(names)}
        def norm(group):
            v = np.zeros(len(names)); v[[idx[g] for g in group]] = 1
            return v @ G @ v
        total = norm(['loc', 'A', 'ray', 'end', 'clk'])
        print(f'  x_o={xo}: IBP identity residual {idcheck:.1e};  |K^L|^2 = {total:.4g}')
        for label, groups in splits.items():
            parts = [norm(g) for g in groups]
            print(f'      kappa = {sum(parts) / total:7.3g}   [{label}]   parts: ' + ', '.join(f'{p / total:.3g}' for p in parts))


# ------------------------------------------------------------------ observer-dilation decomposition
def X_centred(xo, nmu=16, l=2, h=1e-4):
    """X_l = (1/(2l+1)) Re sum_m int int P P A^c_2 A^(1)* on the shell (sigma = 1), fixed-ratio times.
    K2^c = Theta2 + (k d_k) theta(k) + (q d_q) theta(q), derivatives at fixed n and times."""
    S = Setup(xo, nt=max(400, int(160 * xo)), nth=max(24, int(8 + 4 * xo)), nph=max(48, int(16 + 8 * xo)))
    kv = np.array([1.0, 0, 0])
    a = S.first(kv)
    Dth = lambda v: (S.first(v * np.exp(h))['theta'] - S.first(v * np.exp(-h))['theta']) / (2 * h)
    Dk = Dth(kv)
    mus, wm = np.polynomial.legendre.leggauss(nmu)
    acc = 0
    for mu, w in zip(mus, wm):
        qv = np.array([mu, np.sqrt(1 - mu**2), 0.0])
        b = S.first(qv)
        k2 = S.second(a, b, second_order_A2(1.0, 1.0, mu, np.sqrt(1 - mu**2), S.t_o))
        K2c = k2['Theta2'] + Dk + Dth(qv)
        acc += 0.5 * w * S_l(K2c, a['theta'], S.n, S.w, l).real
    return acc / (2 * l + 1)


def cross_term(xo, dl=0.01, **kw):
    """-2 D_{k0} X at fixed times: along the fixed-ratio trajectory D_{k0} = d/d ln x_o."""
    return -2 * (X_centred(xo * np.exp(dl), **kw) - X_centred(xo * np.exp(-dl), **kw)) / (2 * dl)
