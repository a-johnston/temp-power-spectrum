"""Numerical evaluation of the closed K=0 forms proved in k0_fields.py (q = 1, x = t/sqrt3)."""

import mpmath as mpm
import numpy as np

s3 = np.sqrt(3.0)


def cin(z):
    """Cin(z) = int_0^z (1-cos u)/u du, for arrays of z >= 0."""
    z = np.atleast_1d(np.asarray(z, float))
    out = np.empty_like(z)
    small = z < 4.0
    zs = z[small]
    s, p = np.zeros_like(zs), np.ones_like(zs)
    for m in range(1, 40):
        p = p * zs * zs / ((2 * m - 1) * (2 * m))
        s += (-1) ** (m + 1) * p / (2 * m)
    out[small] = s
    if (~small).any():
        out[~small] = [float(mpm.euler + mpm.log(v) - mpm.ci(v)) for v in z[~small]]
    return out


def j0(x):
    return np.sin(x) / x


def fields(t):
    x = np.asarray(t, float) / s3
    sx, cx = np.sin(x), np.cos(x)
    j = sx / x
    jp = (cx - j) / x
    n1 = -(j - cx)
    n1p = -(jp + sx) / s3
    h = -6 * (2 * j - cx - 1)
    hp = -6 * (2 * jp + sx) / s3
    return dict(n1=n1, n1p=n1p, h=h, hp=hp)


def a_exact(t):
    x = np.atleast_1d(np.asarray(t, float)) / s3
    r = (
        -78
        - 36 * np.cos(x)
        - 6 * np.cos(2 * x)
        + 72 * np.sin(x) / x
        + 36 * np.sin(2 * x) / x
        - 24 * (np.sin(x) / x) ** 2
        + 24 * cin(2 * x)
    )
    return r if np.ndim(t) else r[0]


def ap_exact(t):
    """da/dt from x^2 a_x = 36x^2 sin x + 72x cos x - 72 sin x + 12x^2 sin 2x + 48x cos 2x - 60 sin 2x + 24x + 48 sin^2 x / x."""
    x = np.asarray(t, float) / s3
    sx, cx, s2, c2 = np.sin(x), np.cos(x), np.sin(2 * x), np.cos(2 * x)
    ax = 36 * sx + 72 * cx / x - 72 * sx / x**2 + 12 * s2 + 48 * c2 / x - 60 * s2 / x**2 + 24 / x + 48 * sx**2 / x**3
    return ax / s3


def n2_exact(t):
    x = np.asarray(t, float) / s3
    return -np.cos(2 * x) - 6 * np.cos(x) + 6 * np.sin(x) / x + 3 * np.sin(2 * x) / x - 5 * (np.sin(x) / x) ** 2


def n2p_exact(t):
    x = np.asarray(t, float) / s3
    sx, cx, s2, c2 = np.sin(x), np.cos(x), np.sin(2 * x), np.cos(2 * x)
    d = 2 * s2 + 6 * sx + 6 * (cx / x - sx / x**2) + 3 * (2 * c2 / x - s2 / x**2) - 5 * (s2 / x**2 - 2 * sx**2 / x**3)
    return d / s3


def iso_R(x):
    """Isotropic part of K_R = 2 F_2 at K->0: I(x) = 4 + cos^2 x - 5 j0^2 (proved in k0_fields.py)."""
    return 4 + np.cos(x) ** 2 - 5 * j0(x) ** 2


def W_R(x):
    """Coefficient W in K_R = I - (mu^2-1/3) W/4: W(x) = 24 Cin(2x) + 6 cos^2 x - 54 + 48 j0^2."""
    return 24 * cin(2 * x) + 6 * np.cos(x) ** 2 - 54 + 48 * j0(x) ** 2
