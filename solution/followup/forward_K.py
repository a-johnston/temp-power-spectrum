"""Forward test at K != 0: exact plane-symmetric PDE and exact rays for a collinear near-cancelling pair.

Units K = 1. Seed zeta = 2 eps [cos(n y) + cos((n-1) y)], so q = n, q - K = n - 1, and qhat || Khat (mu_K = 1).
Observers are labelled by y0. At O(eps^2) the observer-label harmonic y0 -> cos(y0) contains only the pair
(q, K-q), with coefficient 2 eps^2 Theta_K; the map's first-order part has no such harmonic for n >= 3, and
for n >= 4 neither does its O(eps^3) part. Harmonic 0 contains the two K=0 pairs, Q(n) + Q(n-1), which
is an internal control of the new seed and initial data.

Initial data: exact first-order closed forms per mode, plus the exact nonlinear gradient data
    nu = (t^2/9) e^{-2 zeta} (zeta_yy + zeta_y^2 / 2),   f = (t^2/3) e^{-2 zeta} (zeta_yy - zeta_y^2),
minus their linear parts, started at x <= 0.1. The K channel also gets the exact series through t^8 (k_series.py):
at small K, an initial-data error in that channel excites the constraint-violating mode nu ~ t^2, which grows and
enters f' through (M). With gradient data alone this produced a spurious K->0 offset, about -0.44 at (x_*,x_o) = (3,10).
"""

import argparse
from functools import cache

import numpy as np
from forward_plane import Evol, observe
from k0_closed import fields
from k_series import pair_series as _pair_series
from map_k0 import theta_pair
from numpy.polynomial.legendre import leggauss

s3 = np.sqrt(3.0)
pair_series = cache(_pair_series)


class EvolK(Evol):
    def __init__(self, eps, n, tstart, tend, dt):
        self.n = n
        super().__init__(eps, tstart, tend, N=max(32, 8 * n), dt=dt, mmax=2 * n + 2)

    def seed(self):
        self.zeta, self.zy = self.zeta_at(self.y)
        e, n = self.eps, self.n
        self.zyy = -2 * e * (n**2 * np.cos(n * self.y) + (n - 1) ** 2 * np.cos((n - 1) * self.y))

    def zeta_at(self, y):
        e, n = self.eps, self.n
        return (
            2 * e * (np.cos(n * y) + np.cos((n - 1) * y)),
            -2 * e * (n * np.sin(n * y) + (n - 1) * np.sin((n - 1) * y)),
        )

    def initial(self, t):
        f = np.zeros(self.N)
        nu = np.zeros(self.N)
        nup = np.zeros(self.N)
        for q in (self.n, self.n - 1):
            F = fields(np.array([q * t]))
            g = 2 * self.eps * np.cos(q * self.y)
            f += F["h"][0] * g
            nu += F["n1"][0] * g
            nup += q * F["n1p"][0] * g
        ex = np.exp(-2 * self.zeta)
        nu2 = t**2 / 9 * (ex * (self.zyy + self.zy**2 / 2) - self.zyy)
        f2 = t**2 / 3 * (ex * (self.zyy - self.zy**2) - self.zyy)
        # K channel (harmonic y, amplitude 2 eps^2): exact series beyond the t^2 gradient term
        c, d = pair_series(self.n, -(self.n - 1), order=8)
        cK = 2 * self.eps**2 * np.cos(self.y)
        dnu = sum(float(c[m]) * t**m for m in c if m > 2)
        dnup = sum(float(c[m]) * m * t ** (m - 1) for m in c if m > 2)
        df = sum(float(d[m]) * t**m for m in d if m > 2)
        return np.array([f + f2 + df * cK, nu + nu2 + dnu * cK, nup + 2 * nu2 / t + dnup * cK]), 2 * np.mean(f2) / t


def channels(eps, n, xe, xo):
    """P2 coefficient of the K and K=0 observer-label harmonics; times from x = (n - 1/2) t / sqrt3."""
    qm = n - 0.5
    te, to, tI = s3 * xe / qm, s3 * xo / qm, s3 * xe / qm / 2
    tstart = min(0.4 * tI, 0.1 * s3 / n)
    ev = EvolK(eps, n, tstart, to * (1 + 200 * eps) + 0.5 / n, 0.01 / n)
    mus, wts = leggauss(int(max(32, 5 * s3 * (xo - xe) + 16)))
    wts = wts / 2
    m = 4 * n + 6
    y0s = 2 * np.pi * np.arange(m) / m
    Iv, _ = observe(ev, y0s, mus, te, tI, to, step=0.04 / n)
    th = np.exp(-Iv) / np.sum(wts * np.exp(-Iv), axis=1, keepdims=True) - 1
    p2 = (3 * mus**2 - 1) / 2
    proj = lambda y: 5 * np.sum(wts * y * p2)
    return proj((th * np.cos(y0s)[:, None]).mean(axis=0)) / eps**2, proj(th.mean(axis=0)) / eps**2, te, to, tI


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", default="1,3:8,16,32,64;3,10:16,32,64;5,20:32,64")
    args = ap.parse_args()
    print("mu_K=1 pair (q, K-q), K/q_m = 1/(n-1/2). Q_K: forward K-channel; Q_lim: K->0 kernel at q_m;")
    print("ctrl: forward K=0 channel minus Q(n)+Q(n-1), which should vanish.")
    for case in args.cases.split(";"):
        xs, ns = case.split(":")
        xe, xo = map(float, xs.split(","))
        ql = theta_pair(s3 * xe, s3 * xo, s3 * xe / 2, nmu=160)["Q2"]
        for n in map(int, ns.split(",")):
            a = channels(2e-3, n, xe, xo)
            b = channels(1e-3, n, xe, xo)
            QK = (4 * b[0] - a[0]) / 3
            Q0 = (4 * b[1] - a[1]) / 3
            te, to, tI = b[2:]
            q0 = sum(theta_pair(q * te, q * to, q * tI, nmu=160)["Q2"] for q in (n, n - 1))
            print(
                f"x*={xe} xo={xo} n={n:3} K t_o={to:.3f}: Q_K {QK:+.5f}  Q_lim {ql:+.5f}  diff {QK - ql:+.2e}  ctrl {Q0 - q0:+.1e}",
                flush=True,
            )
