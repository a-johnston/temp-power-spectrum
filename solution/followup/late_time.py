"""Late-time behaviour of the K->0 temperature quadrupole and of the field balance.

Numerical regression (not a proof). Q2(x_o) at fixed x_* is fitted to c0 + c1 ln x_o + c2 cos 2x_o + c3 sin 2x_o + (decaying terms).
"""

import numpy as np
from k0_closed import W_R
from map_k0 import s3, theta_pair


def fit(xe, xos):
    rows, ys, res = [], [], []
    for xo in xos:
        r = theta_pair(s3 * xe, s3 * xo, s3 * xe / 2, nmu=int(4 * s3 * (xo - xe) + 60), nt=int(70 * s3 * (xo - xe)))
        y = xo - xe
        rows.append(
            [
                1,
                np.log(xo),
                np.cos(2 * xo),
                np.sin(2 * xo),
                np.cos(s3 * y) / (s3 * y),
                np.sin(s3 * y) / (s3 * y),
                np.cos(xo) / xo,
                np.sin(xo) / xo,
            ]
        )
        ys.append(r["Q2"])
        res.append(r["Q2"] + (W_R(np.array([xo]))[0] - W_R(np.array([xe]))[0]) / 3)
    A = np.array(rows)
    coef, *_ = np.linalg.lstsq(A, np.array(ys), rcond=None)
    resid = np.max(np.abs(A @ coef - np.array(ys)))
    return coef, resid, np.array(res)


if __name__ == "__main__":
    rng = np.random.default_rng(1)
    for xe in (2.0, 5.0, 8.0):
        xos = np.sort(rng.uniform(xe + 15, xe + 110, 36))
        coef, resid, res = fit(xe, xos)
        C = coef[0] + 8 * np.log(xe) - 2
        print(
            f"x*={xe}: ln coeff {coef[1]:+.4f}, cos2x_o {coef[2]:+.4f}, sin2x_o {coef[3]:+.4f}; C(x*)={coef[0] + coef[1] * (-np.log(xe)) * -1:+.4f} [Q2 = c0 + c1 ln x_o]; c0={coef[0]:+.4f}; max resid {resid:.3f}",
            flush=True,
        )
