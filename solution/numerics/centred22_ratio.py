"""Is the centred 22 spectrum proportional to C11 order by order? (NUMERICAL test, not evidence)

Exact results: C^{c,22}/C11 = 5/3 at k0^4 and the k0^6 coefficients have ratio 28/3 with identical
time dependence. This script fits the k0^4, k0^6, k0^8 coefficients of
  R_c(kappa) = (1/2) < |K^c_2(k,q)|^2_{l=2} >_mu   and   R_1(kappa) = |K_1(k)|^2_{l=2},
K^c_2 = K_2 + D K_1(k) + D K_1(q), from the validated forward map at small kappa, for given times.
If the k0^8 ratio is the same for different (e, o), the proportionality persists at that order.
Run: .venv/bin/python -B solution/numerics/centred22_ratio.py e o [iota]
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import cubic_fields as cf  # noqa: E402
import forward_cubic as fc  # noqa: E402

r2 = np.sqrt(2.0)
BASIS = [np.diag([1, -1, 0]) / r2, np.diag([1, 1, -2]) / np.sqrt(6.0)]
for i, j in ((0, 1), (0, 2), (1, 2)):
    B = np.zeros((3, 3))
    B[i, j] = B[j, i] = 1 / r2
    BASIS.append(B)
H = 0.01


def quad(ks, labels, e, o, iota):
    """Five orthonormal quadrupole components of the multilinear map coefficient."""
    prof = cf.build_profiles(ks, 4.5 * o / 3, third=False)
    fm = fc.ForwardMap(ks, prof, e, o, iota)
    Th = fm.map()
    return np.array([fm.quadrupole(Th, C, list(labels)) for C in BASIS])


def D_single(k, e, o, iota):
    v = {j: quad({1: tuple(np.exp(j * H) * np.asarray(k))}, (1,), e, o, iota) for j in (-2, -1, 1, 2)}
    return (v[-2] - 8 * v[-1] + 8 * v[1] - v[2]) / (12 * H)


def ratios(kappa, e, o, iota, n_mu=10):
    k = kappa * np.array([0.0, 0.0, 1.0])
    R1 = np.sum(np.abs(quad({1: tuple(k)}, (1,), e, o, iota)) ** 2)
    Dk = D_single(k, e, o, iota)
    mus, wmu = np.polynomial.legendre.leggauss(n_mu)
    Rc = 0.0
    for mu, w in zip(mus, wmu, strict=True):
        q = kappa * np.array([np.sqrt(1 - mu**2), 0.0, mu])
        K2 = quad({1: tuple(k), 2: tuple(q)}, (1, 2), e, o, iota)
        Kc = K2 + Dk + D_single(q, e, o, iota)
        Rc += 0.5 * w * 0.5 * np.sum(np.abs(Kc) ** 2)
    return Rc, R1


if __name__ == "__main__":
    e, o = float(sys.argv[1]), float(sys.argv[2])
    iota = float(sys.argv[3]) if len(sys.argv) > 3 else 0.5
    ks = tuple(float(x) for x in os.environ.get("KAPPAS", "0.08 0.12 0.16 0.2 0.24").split())
    data = [ratios(kp, e, o, iota) for kp in ks]
    M = np.array([[kp**4, kp**6, kp**8, kp**10, kp**12] for kp in ks])
    c = np.linalg.solve(M, [d[0] for d in data])
    l = np.linalg.solve(M, [d[1] for d in data])  # noqa: E741
    print(
        f"e={e} o={o}: ratio k^4 {c[0] / l[0]:.6f} (5/3={5 / 3:.6f}), k^6 {c[1] / l[1]:.6f} (28/3={28 / 3:.6f}), k^8 {c[2] / l[2]:.5f}"
    )
