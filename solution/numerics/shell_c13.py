"""Monochromatic-shell cubic quadrupole at finite k0 with the validated forward map (NUMERICAL).

For an isotropic shell, (1/3) Re C_2^(13) / (sigma^2 C_2^(11)) = <K_3(k,q,-q)>_qhat / K_1(k) for the
axisymmetric weight w = 5 P_2(n.khat); by rotational symmetry the q average reduces to mu = khat.qhat.
The observer-centred split K_3 = K^c_3 - sum D K_2 - sum D^2 K_1 is evaluated at each mu.
Fixed-ratio trajectory of the exploratory appendix: k0 = 1, x_o = k0 o / sqrt 3, x_* = x_o/3, x_I = x_o/6.

Run: .venv/bin/python -B solution/numerics/shell_c13.py x_o [n_mu]
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import cubic_fields as cf  # noqa: E402
import forward_cubic as fc  # noqa: E402
import observer_centred as oc  # noqa: E402

CQ = np.diag([-1 / 3, -1 / 3, 2 / 3])  # w = (15/2) n.C.n = 5 P_2(n_z)


def run(x_o, n_mu=16, dilations=True):
    o = np.sqrt(3) * x_o
    e, iota = o / 3, o / 6
    T = 1.15 * o
    sky = dict(nth=18, nph=36) if x_o <= 1.5 else dict(nth=26, nph=52)
    nt = 44 if x_o <= 1.5 else 64
    k = np.array([0.0, 0.0, 1.0])

    def fwd(ks, labels):
        prof = cf.build_profiles(ks, T, third=len(labels) == 3)
        fm = fc.ForwardMap(ks, prof, e, o, iota, nt=nt, **sky)
        return fm.quadrupole(fm.map(), CQ, list(labels))

    K1 = fwd({1: tuple(k)}, (1,))
    mus, wmu = np.polynomial.legendre.leggauss(n_mu)
    K3s, Kcs = [], []
    start = time.time()
    for mu in mus:
        q = np.array([np.sqrt(1 - mu**2), 0.0, mu])
        ks = {1: tuple(k), 2: tuple(q), 3: tuple(-q)}
        K3 = fwd(ks, (1, 2, 3))
        K3s.append(K3)
        if dilations:

            def kern(ks_, labels, e_, o_, iota_, C_, scale=1.0):
                m = {a: tuple(scale * np.asarray(ks_[a], dtype=float)) for a in labels}
                return fwd(m, labels)

            Kcs.append(oc.split(ks, e, o, iota, CQ, K3=K3, kernel_fn=kern)["centred"])
        print(f"  mu={mu:+.4f}  K3={K3.real:+.6e}  {time.time() - start:.0f}s", flush=True)
    avg = lambda v: 0.5 * np.dot(wmu, np.asarray(v))  # noqa: E731
    ratio = avg(K3s) / K1
    print(f"x_o={x_o}: K1={K1.real:+.6e}  <K3>/K1 = (1/3)Re C13_2/(sigma^2 C11_2) = {ratio.real:+.6f}", flush=True)
    if dilations:
        rc = avg(Kcs) / K1
        print(f"x_o={x_o}: observer-centred part {rc.real:+.6f}, dilation part {(ratio - rc).real:+.6f}", flush=True)
    return ratio


if __name__ == "__main__":
    x = float(sys.argv[1])
    nm = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    run(x, nm, dilations="--no-dilation" not in sys.argv)
