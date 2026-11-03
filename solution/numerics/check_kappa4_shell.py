"""Numerical check of the analytic kappa^4 shell cubic quadrupole (report eq. stationary-four-gradient).

  Jbar_abc = kappa^2 4(o^2-e^2)/9 + kappa^4 [3984e^4-6256e^3o+4692e^2o^2-1853o^4-567iota^4]/3402 + O(kappa^6)

Independent of all analytic routes: the validated forward map at small kappa, the shell q average
by Gauss-Legendre in mu (exact for the low-degree polynomial kernel; q -> -q symmetry halves the
nodes), and Richardson extrapolation in kappa^2. Two iota values isolate the observer-local term.
Run: .venv/bin/python -B solution/numerics/check_kappa4_shell.py iota
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

E, O = 1.2, 3.0
CQ = np.diag([2 / 3, -1 / 3, -1 / 3])  # C = Q_k for k along x: w = 5 P_2(n_x)


def shell_average(kappa, iota, n_mu=6):
    mus, wmu = np.polynomial.legendre.leggauss(n_mu)
    total = 0.0
    for mu, w in zip(mus, wmu, strict=True):
        if mu < 0:
            continue  # K3(k,q,-q) is even in mu
        q = kappa * np.array([mu, np.sqrt(1 - mu**2), 0.0])
        ks = {1: (kappa, 0.0, 0.0), 2: tuple(q), 3: tuple(-q)}
        prof = cf.build_profiles(ks, 4.0)
        fm = fc.ForwardMap(ks, prof, E, O, iota)
        total += 2 * 0.5 * w * fm.quadrupole(fm.map(), CQ, [1, 2, 3])
    return total


if __name__ == "__main__":
    iota = float(sys.argv[1])
    a2 = 4 * (O**2 - E**2) / 9
    a4 = (3984 * E**4 - 6256 * E**3 * O + 4692 * E**2 * O**2 - 1853 * O**4 - 567 * iota**4) / 3402
    ks = (0.04, 0.06, 0.08, 0.1)
    vals = [shell_average(k, iota).real / k**2 for k in ks]
    # fit vals = A + B k^2 + C k^4 + D k^6 exactly through four points (truncation ~ kappa^8)
    M = np.array([[1, k**2, k**4, k**6] for k in ks])
    A, B, Cc, _ = np.linalg.solve(M, vals)
    print(f"iota={iota}: kappa^2 coeff forward {A:+.9f} analytic {a2:+.9f} rel {abs(A - a2) / abs(a2):.1e}")
    print(
        f"iota={iota}: kappa^4 coeff forward {B:+.9f} analytic {a4:+.9f} rel {abs(B - a4) / abs(a4):.1e}  (kappa^6 fit {Cc:+.4f})"
    )
