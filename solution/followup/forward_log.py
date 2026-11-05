"""Forward test in the logarithmic regime: exact PDE and exact rays at K = 0 out to large x_o.

Compares forward Q (Richardson in eps) with the map kernel. Smaller amplitudes than forward_plane.py, because the
second-order coefficient needs eps * x_o << 1 at the latest time.
"""

import numpy as np
from forward_plane import theta_avg
from map_k0 import s3, theta_pair


def P2(mus, wts, y):
    return 5 * np.sum(wts * y * (3 * mus**2 - 1) / 2)


if __name__ == "__main__":
    for xe, xo in [(5, 40), (2, 40), (5, 80), (2, 80), (5, 160)]:
        q2 = theta_pair(s3 * xe, s3 * xo, s3 * xe / 2, nmu=int(5 * s3 * (xo - xe) + 60), nt=int(100 * s3 * (xo - xe)) + 2000)[
            "Q2"
        ]
        e0 = 2e-4 * 40 / xo * 5
        vals = []
        for eps in (2 * e0, e0):
            mus, wts, th, _ = theta_avg(eps, xe, xo, nmu=int(4 * s3 * (xo - xe) + 16))
            vals.append(P2(mus, wts, th) / eps**2)
        rich = (4 * vals[1] - vals[0]) / 3
        print(
            f"xe={xe} xo={xo} eps={e0:.1e}: forward {vals[0]:+.5f}, {vals[1]:+.5f}, extrap {rich:+.5f}; kernel {q2:+.5f}; diff {rich - q2:+.2e}",
            flush=True,
        )
