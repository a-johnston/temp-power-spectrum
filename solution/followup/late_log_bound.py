"""Bound on a logarithm in the non-generated pieces of Q, over a long lever arm.

Nongen(x_o) = Q - (-[W]/3) - em is fitted at fixed x_* over log-spaced x_o, to
    c0 + c1 ln x_o + c2 cos 2x_o + c3 sin 2x_o + decaying terms.
c1 = 0 would mean that the whole logarithm of Q is -8 ln x_o. Spread: repeated fits on subsets of the points.
Chunks of 2 sky directions keep the peak memory under about 300 MB at the largest x_o.
"""

import numpy as np
from k0_closed import W_R
from map_k0 import s3, theta_pair


def nongen(xe, xo):
    lam = s3 * (xo - xe)
    r = theta_pair(s3 * xe, s3 * xo, s3 * xe / 2, nmu=int(1.5 * lam + 80), nt=int(24 * lam) + 2000, chunk=2)
    W = lambda x: W_R(np.array([float(x)]))[0]
    return r["Q2"] - r["Q2pieces"]["em"] + (W(xo) - W(xe)) / 3


def design(xos, xe):
    y = xos - xe
    return np.column_stack(
        [
            np.ones_like(xos),
            np.log(xos),
            np.cos(2 * xos),
            np.sin(2 * xos),
            np.cos(xos) / xos,
            np.sin(xos) / xos,
            np.cos(s3 * y) / (s3 * y),
            np.sin(s3 * y) / (s3 * y),
            1 / xos,
        ]
    )


if __name__ == "__main__":
    print("accuracy check (5,160): low-resolution Nongen vs values from the Table 1 settings")
    lo = nongen(5.0, 160.0)
    hi = theta_pair(s3 * 5, s3 * 160, s3 * 5 / 2, nmu=int(5 * s3 * 155 + 60), nt=int(100 * s3 * 155) + 2000)
    W = lambda x: W_R(np.array([float(x)]))[0]
    print(f"  {lo:+.6f}  {hi['Q2'] - hi['Q2pieces']['em'] + (W(160.0) - W(5.0)) / 3:+.6f}", flush=True)
    rng = np.random.default_rng(3)
    for xe in (2.0, 5.0, 8.0):
        xos = np.exp(rng.uniform(np.log(30.0), np.log(1000.0), 48))
        ys = np.array([nongen(xe, x) for x in xos])
        A = design(xos, xe)
        coef, *_ = np.linalg.lstsq(A, ys, rcond=None)
        res = ys - A @ coef
        subs = []
        for _ in range(200):
            idx = rng.integers(0, len(xos), len(xos))
            subs.append(np.linalg.lstsq(A[idx], ys[idx], rcond=None)[0][1])
        print(
            f"x*={xe}: ln coeff of non-generated part {coef[1]:+.4f} (bootstrap sd {np.std(subs):.4f}); "
            f"cos2x_o {coef[2]:+.4f} sin2x_o {coef[3]:+.4f}; rms resid {np.sqrt(np.mean(res**2)):.3f}",
            flush=True,
        )
