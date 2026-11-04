"""Does the leading-structure projection of K^c_2 obey -(1/2) D deltaK_1 at all orders? (NUMERICAL)

Exact through O(k0^4) (check_observer_dilation_k6.py). With the fixed template T = TF(khat qhat + qhat khat),
compare f = <T . K^c_2(k,q)>_mu and g = <T . (-(1/2) D deltaK_1)>_mu, deltaK_1 = K_1(k+q) - K_1(k) - K_1(q).
(f - g)/f ~ kappa^4 would mean the relation fails at O(k0^6); ~1e-9 would mean it is exact.
Run: .venv/bin/python -B solution/numerics/half_dilation_projection.py
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from centred22_ratio import BASIS, D_single, quad  # noqa: E402

E, O, IOTA = 1.2, 3.0, 0.5


def projections(kappa, n_mu=8):
    khat = np.array([0.0, 0.0, 1.0])
    k = kappa * khat
    Dk = D_single(k, E, O, IOTA)
    mus, wmu = np.polynomial.legendre.leggauss(n_mu)
    f = g = 0.0
    for mu, w in zip(mus, wmu, strict=True):
        qhat = np.array([np.sqrt(1 - mu**2), 0.0, mu])
        q = kappa * qhat
        M = np.outer(khat, qhat) + np.outer(qhat, khat)
        T = np.array([np.sum(C * (M - np.eye(3) * np.trace(M) / 3)) for C in BASIS])
        Dq = D_single(q, E, O, IOTA)
        Kc = quad({1: tuple(k), 2: tuple(q)}, (1, 2), E, O, IOTA) + Dk + Dq
        DdK = D_single(k + q, E, O, IOTA) - Dk - Dq  # D(deltaK_1): dilation of each single kernel
        f += 0.5 * w * (T @ Kc)
        g += 0.5 * w * (T @ (-0.5 * DdK))
    return f, g


if __name__ == "__main__":
    for kappa in (0.1, 0.2, 0.3, 0.6):
        f, g = projections(kappa)
        print(f"kappa={kappa}: f={f.real:+.10e} g={g.real:+.10e} (f-g)/f={((f - g) / f).real:+.3e}", flush=True)
