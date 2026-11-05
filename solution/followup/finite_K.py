"""Finite-K comparison: the repository's general-momentum second-order Green solution vs the K->0 closed forms.

Uses the exact sources in solution/results (via numerics/explore_fourth_order.py); independent of k0_fields.py.
Part 1: K_R(q1,q2;t) at decreasing K with generic angle (including |q1| != |q2|).
Part 2: P22(K) by direct 2D quadrature at finite K against the closed K->0 coefficient A0.
"""

import os
import sys

import numpy as np
from k0_closed import W_R, iso_R

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "numerics"))
sys.path.insert(0, os.path.join(HERE, ".."))
import explore_fourth_order as ef  # noqa: E402

s3 = np.sqrt(3.0)


def TF(M):
    return M - np.eye(3) * np.trace(M) / 3


def pair_fields(p, q, b1, b2, t_o, ng=20000, lo=1e-3):
    """nu_2 and A_2 at (k1=(p,0,0), k2=(b1,b2,0)) at time t_o by a 4th-order cumulative rule."""
    if ef.SRC is None:
        ef.SRC = ef.load_sources()
    SRC = ef.SRC
    tg = np.linspace(t_o * lo, t_o, ng)
    dt = tg[1] - tg[0]

    def cum(y):
        out = np.zeros_like(y)
        inc = (y[..., 1:] + y[..., :-1]) / 2 * dt
        out[..., 1:] = np.cumsum(inc, axis=-1)
        g = np.gradient(y, dt, axis=-1, edge_order=2)
        return out - dt**2 / 12 * (g - g[..., :1])

    KV = np.array([p + b1, b2, 0.0])
    K = np.linalg.norm(KV)
    c = K / s3
    x = c * tg
    # Series-safe u1, u2-combination: use the Green form at small x
    u1 = np.sin(x) / x - np.cos(x)
    u2 = -np.cos(x) / x - np.sin(x)
    ev = lambda fn: np.broadcast_to(np.asarray(fn(tg, p, q, b1, b2), dtype=complex), tg.shape)
    cN, cH = ev(SRC["N"]), ev(SRC["H"])
    # first part of the integrals from 0 to t_lo is negligible: sources are O(1) x O(s^2), lo small
    P1, P2 = cum(u1 * cN), cum(u2 * cN)
    c2 = (2 * (p**2 + q**2) - p * b1) / 9
    # leading regular behaviour below t_lo: nu = c2 t^2 handled by the homogeneous piece
    nu = 9 * c2 / K**2 * u1 - (u2 * P1 - u1 * P2) / (3 * c)
    u1p = np.cos(x) / x - np.sin(x) / x**2 + np.sin(x)
    u2p = np.sin(x) / x + np.cos(x) / x**2 - np.cos(x)
    nup = 9 * c2 / K**2 * c * u1p - (u2p * P1 - u1p * P2) / 3
    hS = 3 / (2 * K**2) * (4 * K**2 * nu + 12 * nup / tg + 12 * nu / tg**2 + cH)
    Kh = KV / K
    Pp = np.eye(3) - np.outer(Kh, Kh)
    cM = np.array([ev(fn) for fn in SRC["M"]])
    w = 2j / K**2 * cum(Pp @ cM)
    cT = np.array([[ev(SRC["T"][i][j]) for j in range(3)] for i in range(3)])
    PMP = np.einsum("ia,abn,bj->ijn", Pp, cT, Pp)
    TT = PMP - Pp[:, :, None] * np.einsum("aan->n", PMP)[None, None, :] / 2
    Ic, Is = cum(tg * np.cos(K * tg) * TT), cum(tg * np.sin(K * tg) * TT)
    h = -2 / (K * tg) * (np.sin(K * tg) * Ic - np.cos(K * tg) * Is)
    A = (np.outer(Kh, Kh) - np.eye(3) / 3) * hS[-1] + np.outer(KV, w[:, -1]) + np.outer(w[:, -1], KV) + h[:, :, -1]
    return nu[-1], A


def Q_nu(k, t):
    Q, _, nu = ef.mode(k, t)
    return Q, nu


def KR_finite(p, q, b1, b2, t):
    """Full-series z1 z2 coefficient of R at K (the formula proved off-shell in check_chart_R.py)."""
    k1, k2 = np.array([p, 0.0, 0.0]), np.array([b1, b2, 0.0])
    K = k1 + k2
    nu2, A = pair_fields(p, q, b1, b2, t)
    Q1, n1 = Q_nu(p, t)
    Q2, n2 = Q_nu(q, t)
    b_1, b_2 = Q1 / 2, Q2 / 2
    psi1, psi2 = 1 + n1, 1 + n2
    kk = k1 @ k2
    xi_dpsi = -kk * (b_1 * psi2 + b_2 * psi1)
    trm2 = 2 * b_1 * b_2 * kk**2
    xi_dH = -kk * (b_1 * Q2 * TF(np.outer(k2, k2)) + b_2 * Q1 * TF(np.outer(k1, k1)))
    m2 = b_1 * b_2 * kk * (np.outer(k1, k2) + np.outer(k2, k1))
    S = A + xi_dH - TF(m2)
    Kh = K / np.linalg.norm(K)
    return nu2 + xi_dpsi - trm2 / 6 - Kh @ S @ Kh / 4


def KR_limit(q, mu, t):
    x = q * t / s3
    return iso_R(x) - (mu**2 - 1 / 3) * W_R(x)[0] / 4


def part1():
    print("K_R(q, K-q; t) -> closed K->0 limit. q=1, x = q t/sqrt3, mu = Khat.qhat")
    for x_ in (0.5, 2.0, 6.0):
        t = s3 * x_
        for theta in (np.pi / 2, np.pi / 3, 0.3):
            row = []
            for K in (0.4, 0.2, 0.1, 0.05):
                KV = K * np.array([np.cos(theta), np.sin(theta)])
                b1, b2 = KV[0] - 1.0, KV[1]
                q2 = np.hypot(b1, b2)
                row.append(KR_finite(1.0, q2, b1, b2, t).real)
            lim = KR_limit(1.0, np.cos(theta), t)
            print(
                f"x={x_:4} mu={np.cos(theta):+.3f} K_R: "
                + " ".join(f"{v:+.6f}" for v in row)
                + f"  limit {lim:+.6f}  err(K=.05) {row[-1] - lim:+.1e}",
                flush=True,
            )


def bump(q, q0=1.0, w=0.3):
    u = (q - q0) / w
    out = np.zeros_like(q)
    m = np.abs(u) < 1
    out[m] = np.exp(1 - 1 / (1 - u[m] ** 2))
    return out


def A0_closed(t, q0=1.0, w=0.3):
    """A0 = (1/pi^2) int q^2 dq P^2 [alpha^2 + beta^2/5], F2 -> alpha + beta P2(mu), alpha = I/2, beta = -W/12."""
    qs, ws = np.polynomial.legendre.leggauss(200)
    q = q0 + w * qs
    x = q * t / s3
    al, be = iso_R(x) / 2, -W_R(x) / 12
    return np.sum(w * ws * q**2 * bump(q) ** 2 * (al**2 + be**2 / 5)) / np.pi**2


def P22_direct(K, t, nq=14, nmu=14, q0=1.0, w=0.3):
    qs, wq = np.polynomial.legendre.leggauss(nq)
    mus, wm = np.polynomial.legendre.leggauss(nmu)
    tot = 0.0
    for qv, wv in zip(q0 + w * qs, w * wq, strict=True):
        for mu, wmu in zip(mus, wm, strict=True):
            # q = qv at polar angle mu to K. Frame: k1 = (qv,0,0), K = K (mu, s, 0), k2 = K - k1
            s = np.sqrt(1 - mu**2)
            b1, b2 = K * mu - qv, K * s
            q2 = np.hypot(b1, b2)
            Pq2 = bump(np.array([q2]))[0]
            if Pq2 == 0:
                continue
            F2 = KR_finite(qv, q2, b1, b2, t).real / 2
            tot += wv * wmu * qv**2 * bump(np.array([qv]))[0] * Pq2 * F2**2
    return 2 * tot * 2 * np.pi / (2 * np.pi) ** 3


if __name__ == "__main__":
    part1()
    print("P22(K)/A0 at x0 = q0 t/sqrt3:")
    for x_ in (1.0, 4.0):
        t = s3 * x_
        A0 = A0_closed(t)
        vals = [P22_direct(K, t, nq=28, nmu=10) / A0 for K in (0.2, 0.1, 0.05)]
        print(f"x0={x_}: A0={A0:.6e}; P22/A0 at K=0.2,0.1,0.05: " + " ".join(f"{v:.5f}" for v in vals), flush=True)
