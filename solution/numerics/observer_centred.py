"""Observer-centred split of the cubic quadrupole at finite wavelength (NUMERICAL, not evidence).

Exact identity (report, Sec. observer-centred dilation): with D the simultaneous momentum dilation
of a kernel at fixed sky direction,

    K_3(k1,k2,k3) = K^c_3 - sum_{pairs} D K_2(k_b,k_c) - sum_a D^2 K_1(k_a),

and K^c_3 vanishes whenever any momentum is zero.  The dilation pieces need only the linear and
quadratic maps.  This script measures how much of the finite-wavelength cubic quadrupole they carry.
"""

from __future__ import annotations

import os
import sys
from itertools import combinations

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import cubic_fields as cf  # noqa: E402
import forward_cubic as fc  # noqa: E402

H = 0.01
STENCIL = (-2, -1, 0, 1, 2)


def kernel(ks, labels, e, o, iota, C, scale=1.0):
    m = {a: tuple(scale * np.asarray(ks[a], dtype=float)) for a in labels}
    prof = cf.build_profiles(m, 4.0, third=len(labels) == 3)
    fm = fc.ForwardMap(m, prof, e, o, iota)
    return fm.quadrupole(fm.map(), C, list(labels))


def dilate(ks, labels, e, o, iota, C, order, kernel_fn=kernel):
    v = {j: kernel_fn(ks, labels, e, o, iota, C, np.exp(j * H)) for j in STENCIL}
    if order == 1:
        return (v[-2] - 8 * v[-1] + 8 * v[1] - v[2]) / (12 * H)
    return (-v[-2] + 16 * v[-1] - 30 * v[0] + 16 * v[1] - v[2]) / (12 * H * H)


def split(ks, e, o, iota, C, K3=None, kernel_fn=kernel):
    labs = sorted(ks)
    if K3 is None:
        K3 = kernel_fn(ks, labs, e, o, iota, C)
    pairs = sum(dilate(ks, pr, e, o, iota, C, 1, kernel_fn) for pr in combinations(labs, 2))
    singles = sum(dilate(ks, (a,), e, o, iota, C, 2, kernel_fn) for a in labs)
    centred = K3 + pairs + singles
    return dict(K3=K3, pair_dilation=-pairs, single_dilation=-singles, centred=centred)


if __name__ == "__main__":
    from compare_cubic_quadrupole import IOTA, TRIPLES, C, E, O

    for name, ks in TRIPLES.items():
        r = split(ks, E, O, IOTA, C)
        print(
            f"{name}: K3 {r['K3'].real:+.8f} = centred {r['centred'].real:+.8f}"
            f" + pair dilation {r['pair_dilation'].real:+.8f} + single dilation {r['single_dilation'].real:+.8f}"
            f"   |centred/K3| = {abs(r['centred'] / r['K3']):.3f}",
            flush=True,
        )
