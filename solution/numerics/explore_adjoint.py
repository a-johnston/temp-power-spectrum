"""EXPLORATORY: the top-order cubic response through observable-weighted (adjoint) Green kernels.

NOT EVIDENCE (see explore_fourth_order.py).  At (k, q, -q) the third-order field kernels enter the
map only through the line-of-sight term -int (1/2) A_3'_nn e^{i K.n (t_o - t)}: the emission and clock
shifts of order 3 multiply background quantities that vanish.  On the shell only the scalar sector
h_S of A_3 survives the q-average.  Its l-projection is the transfer integral
    T_l[h] = int_{t*}^{t_o} h'(t) w_l(K (t_o - t)) dt,      w_l = j_l / 3 + j_l'',
and h_S = Hop[nu_3] + 3 c_Ham / (2K^2) with Hop[nu] = (3/2K^2)(4K^2 nu + 12 nu'/t + 12 nu/t^2), and
nu_3 = 9 c_3 u_1 / K^2 + Green[c_N].  Exchanging the time integrals gives
    T_l[h_S3] = (c_3/c_1) T_l[h_S1] + int_0^{t_o} W_N(s) c_N(s) ds + (3/2K^2) int_{t*}^{t_o} w_l c_Ham' dt,
    W_N(s) = -(1/3c) [u_1(cs) U_2(s) - u_2(cs) U_1(s)],   U_i(s) = int_{max(s,t*)}^{t_o} w_l d_t Hop[u_i(c.)] dt,
so the adjoint kernel is the Wronskian of the PARTIAL transfer integrals of the growing (u_1) and
decaying (u_2) lapse modes.  Checked here against
  (a) the direct third-order Green solution (kernel level, every mu), and
  (b) the full numerical third-order map of explore_third_order.py (map level, l = 2 projection).
"""

import sys

import explore_fourth_order as ef
import explore_third_order as et
import numpy as np

s3 = np.sqrt(3.0)


def jl_w(l, y):
    """w_l(y) = j_l(y)/3 + j_l''(y) for l = 2, from the closed form."""
    assert l == 2
    y = np.asarray(y, dtype=float)
    out = np.empty_like(y)
    small = y < 1e-2
    ys = y[small]
    # series: j_2 = y^2/15 - y^4/210 + ...,  j_2'' = 2/15 - 12 y^2/210 + ...
    out[small] = (ys**2 / 15 - ys**4 / 210) / 3 + (2 / 15 - 12 * ys**2 / 210 + 30 * ys**4 / 7560)
    yb = y[~small]
    s, c = np.sin(yb), np.cos(yb)
    j2 = (3 / yb**3 - 1 / yb) * s - 3 * c / yb**2
    # j_2'' from the spherical Bessel equation: j'' = -2 j'/y - (1 - 6/y^2) j
    j2p = (-9 / yb**4 + 1 / yb**2) * s + (3 / yb**3 - 1 / yb) * c + 6 * c / yb**3 + 3 * s / yb**2
    j2pp = -2 * j2p / yb - (1 - 6 / yb**2) * j2
    out[~small] = j2 / 3 + j2pp
    return out


def hop(u, up, t, K):
    return 3 / (2 * K**2) * (4 * K**2 * u + 12 * up / t + 12 * u / t**2)


def adjoint_vs_direct(x_o, mu, R=3.0, ng=12000, l=2):
    t_o = s3 * x_o
    ts = t_o / R
    tg = np.linspace(t_o * 2e-4, t_o * 1.02, ng)
    cfg = et.symmetric_cfg(1.0, mu)
    ks = {a: np.array(v, dtype=float) for a, v in cfg.items()}
    pairs = {pr: et.pair_fields(ks[pr[0]], ks[pr[1]], tg) for pr in ((1, 2), (1, 3), (2, 3))}
    f3 = et.third_field(ks[1], ks[2], tg, pairs, scalar_only=True)
    K = 1.0
    c = K / s3
    win = (tg >= ts) & (tg <= t_o)
    w = np.where(win, jl_w(l, K * np.clip(t_o - tg, 0, None)), 0.0)
    T = lambda hp: np.trapezoid(hp * w, tg)
    # first order: h_S1 = Q k^2
    _, Qp, _ = ef.mode(K, tg)
    T1 = T(Qp * K**2)
    direct = T(f3["hSp"])
    # adjoint
    x = c * tg
    u1, u2 = np.sin(x) / x - np.cos(x), -np.cos(x) / x - np.sin(x)
    u1p = c * (np.cos(x) / x - np.sin(x) / x**2 + np.sin(x))
    u2p = c * (np.sin(x) / x + np.cos(x) / x**2 - np.cos(x))
    dH1 = np.gradient(hop(u1, u1p, tg, K), tg)
    dH2 = np.gradient(hop(u2, u2p, tg, K), tg)
    rev = lambda f: (lambda C: C[-1] - C)(ef.cumint(f * w, tg))  # int_s^{end} f w
    U1, U2 = rev(dH1), rev(dH2)
    # Wronskian part plus the local part from the jump d_t G(t,t) = -1/3 (Hop contains nu'')
    WN = -(u1 * U2 - u2 * U1) / (3 * c) - 6 / (K**2 * tg) * w
    tsel = tg <= t_o
    c1 = -(K**2) / 9
    hom = f3["c3"] / c1 * T1
    srcN = np.trapezoid((WN * f3["cN"])[tsel], tg[tsel])
    srcH = 3 / (2 * K**2) * T(np.gradient(f3["cH"], tg))
    return dict(direct=direct, adjoint=hom + srcN + srcH, hom=hom, srcN=srcN, srcH=srcH, T1=T1, WN=(tg, WN))


def kernel_check(xos=(1.0, 3.0), mus=(0.9, 0.3, -0.5)):
    for xo in xos:
        for mu in mus:
            r = adjoint_vs_direct(xo, mu)
            rel = abs(r["adjoint"] - r["direct"]) / abs(r["direct"])
            print(
                f"  x_o={xo} mu={mu:+.2f}: T[h_S3]/T[h_S1] direct {r['direct'] / r['T1']:.6g}  adjoint {r['adjoint'] / r['T1']:.6g}"
                f"  rel {rel:.1e}   [hom {r['hom'] / r['T1']:.4g}, c_N {r['srcN'] / r['T1']:.4g}, c_Ham {r['srcH'] / r['T1']:.4g}]",
                flush=True,
            )


def map_check(x_o=1.5, mus=(0.7, -0.4), l=2):
    """(1/2) Re S_l[theta, L3_top] / S_l[theta, theta] from the full map vs Re T[h_S3]/T[h_S1]."""
    for mu in mus:
        nth = max(16, int(8 + 4 * x_o))
        args = dict(nth=nth, nph=2 * nth, nt=max(240, int(100 * x_o)), scalar_only=True)
        M3 = et.Map(x_o, et.symmetric_cfg(1.0, mu), third=True, **args)
        M0 = et.Map(x_o, et.symmetric_cfg(1.0, mu), third=False, **args)
        L3, L0 = M3.L(), M0.L()
        th = L3[frozenset([1])]
        top = L3[et.FULL] - L0[et.FULL]
        ratio_map = ef.S_l(th, top, M3.n, M3.w, l).real / ef.S_l(th, th, M3.n, M3.w, l).real
        r = adjoint_vs_direct(x_o, mu)
        print(
            f"  x_o={x_o} mu={mu:+.2f}: map top-order share {ratio_map:.6g}   adjoint Re T[h_S3]/T[h_S1] {(r['adjoint'] / r['T1']).real:.6g}",
            flush=True,
        )


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else "kernel"
    {"kernel": kernel_check, "map": map_check}[m]()
