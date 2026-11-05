"""Analytic structure of the K->0 temperature quadrupole: the identities behind the boundedness theorem.

q = 1, c = 1/sqrt3, x = c t, lam = t_o - t, phase E_s = exp(i s mu lam), w2 = 1 - mu^2, v = h'/2.
Phi = -4 j0'(x)/x = 4 j1(x)/x, twice the Newtonian potential: Phi = O(x^-2) and Phi' = O(x^-2).

Part 1 (exact, SymPy). The constant-amplitude sound terms are exact t-derivatives along the unperturbed ray:
  sigma_s := v (mu^2 - 1/3) E_s             = d/dt[(h''/2 + i s mu v) E_s] + Phi' E_s
  dmu_s integrand := [i s (1 - h/6) + mu v] E_s = d/dt[i s v E_s] + i s Phi E_s
  F_s := [-i s w2 v - mu (1 + h/3)] E_s      = d/dt[(-i s h/2 - mu v) E_s] - mu Phi E_s
and W = a + h^2/2.

Part 2 (numerical check of the algebra). The reduced ray terms, computed from these forms, reproduce the complete
second-order exponent of map_k0.py.
"""

import numpy as np
import sympy as sp
from k0_closed import W_R, a_exact, fields
from map_k0 import cum_from_right, s3, theta_pair
from numpy.polynomial.legendre import leggauss, legval


def part1():
    t, mu, to = sp.symbols("t mu t_o", positive=True)
    c = 1 / sp.sqrt(3)
    x = c * t
    j0 = sp.sin(x) / x
    h = 6 * (1 + sp.cos(x) - 2 * j0)
    v = sp.diff(h, t) / 2
    j0x = (sp.cos(x) - j0) / x  # j0'(x)
    Phi = -4 * j0x / x
    res = []
    for s in (1, -1):
        E = sp.exp(sp.I * s * mu * (to - t))
        sigma = v * (mu**2 - sp.Rational(1, 3)) * E
        res.append(sigma - sp.diff((sp.diff(h, t, 2) / 2 + sp.I * s * mu * v) * E, t) - sp.diff(Phi, t) * E)
        dmu = (sp.I * s * (1 - h / 6) + mu * v) * E
        res.append(dmu - sp.diff(sp.I * s * v * E, t) - sp.I * s * Phi * E)
        F = (-sp.I * s * (1 - mu**2) * v - mu * (1 + h / 3)) * E
        res.append(F - sp.diff((-sp.I * s * h / 2 - mu * v) * E, t) + mu * Phi * E)
    for r in res:
        assert sp.simplify(sp.expand(sp.simplify(r / sp.exp(sp.I * mu * (to - t))))) == 0 or sp.simplify(r) == 0
    # W = a + h^2/2 (Cin(2x) is common to both, so it is omitted)
    C = sp.Symbol("C")
    X = sp.Symbol("X", positive=True)
    jj, cc = sp.sin(X) / X, sp.cos(X)
    a = -78 - 36 * cc - 6 * sp.cos(2 * X) + 72 * jj + 36 * sp.sin(2 * X) / X - 24 * jj**2 + 24 * C
    W = 24 * C + 6 * cc**2 - 54 + 48 * jj**2
    hh = 6 * (1 + cc - 2 * jj)
    assert sp.simplify(sp.expand(sp.expand_trig(W - a - hh**2 / 2))) == 0
    print("PASS three total-derivative identities and W = a + h^2/2 (exact)")


def phi_fields(t):
    x = t / s3
    j = np.sin(x) / x
    jx = (np.cos(x) - j) / x
    jxx = -2 * jx / x - j
    Phi = -4 * jx / x
    Phip = -4 * (jxx / x - jx / x**2) / s3
    return Phi, Phip


def reduced_exponent(mus, t, to):
    """Second-order ray exponent without the generated a-term, from the exact plane-symmetric reduction."""
    dt = t[1] - t[0]
    F = fields(t)
    h, v = F["h"], F["hp"] / 2
    lam = to - t
    w2 = 1 - mus**2
    out_direct = np.zeros(len(mus), complex)
    out_reduced = np.zeros(len(mus), complex)
    out_print = np.zeros(len(mus), complex)
    Phi, _ = phi_fields(t)
    for i, mu in enumerate(mus):
        Tdirect = 0
        Tred = 0
        Tprint = 0
        for s in (1, -1):
            Es, Em = np.exp(1j * s * mu * lam), np.exp(-1j * s * mu * lam)
            sig = v * (mu**2 - 1 / 3) * Es
            # Direct form, opposite mode m = -s: dmu_m = w2 int_t^to [i m (1-h/6) + mu v] E_m, dy_m = int_t^to (dmu_m - mu A_m)
            dmu_m = w2[i] * cum_from_right((-1j * s * (1 - h / 6) + mu * v) * Em, dt)
            dy_m = cum_from_right(dmu_m - mu * (1 + h / 3) * Em, dt)
            Tdirect += np.trapezoid(v * Es * 2 * mu * dmu_m + 1j * s * sig * dy_m, dx=dt)
            # Reduced form: the same quantities from the total-derivative identities
            W_m = cum_from_right(Phi * Em, dt)
            Wt_m = cum_from_right(W_m, dt)
            dmu_r = w2[i] * (-1j * s) * (v[-1] - v * Em + W_m)
            Y_m = 1j * s * h / 2 - mu * v
            dy_r = w2[i] * (-1j * s) * (v[-1] * lam + Wt_m) + Y_m[-1] - Y_m * Em - mu * W_m
            Tred += np.trapezoid(v * Es * 2 * mu * dmu_r + 1j * s * sig * dy_r, dx=dt)
            # The six nonlocal terms printed in the report as R(mu); the local h^2 term is added back below
            J_s = np.trapezoid(sig, dx=dt)
            Tprint += (
                -2j * s * mu * w2[i] * np.trapezoid(v * Es * (v[-1] + W_m), dx=dt)
                + w2[i] * np.trapezoid(sig * (v[-1] * lam + Wt_m), dx=dt)
                + 1j * s * Y_m[-1] * J_s
                - 1j * s * mu * np.trapezoid(sig * W_m, dx=dt)
            )
        out_direct[i], out_reduced[i] = Tdirect, Tred
        out_print[i] = Tprint + (h[-1] ** 2 - h[0] ** 2) / 4 * (mu**2 - 1 / 3)
    return out_direct, out_reduced, out_print


def part2():
    print("(x*, xo): Q2 gen+ray from map_k0 | direct reduction | total-derivative form | printed R(mu) | B = prod+clk+rem")
    for xe, xo in [(0.5, 1.0), (1.0, 3.0), (3.0, 10.0), (5.0, 20.0), (5.0, 40.0)]:
        te, to = s3 * xe, s3 * xo
        lam = to - te
        r = theta_pair(te, to, te / 2, nmu=int(5 * lam + 60), nt=int(100 * lam) + 2000)
        pm = r["Q2pieces"]["gen"] + r["Q2pieces"]["ray"]
        mus, wts = leggauss(int(3 * lam + 40))
        wts = wts / 2
        t = np.linspace(te, to, int(120 * lam) + 4001)
        d, red, prt = reduced_exponent(mus, t, to)
        p2 = legval(mus, [0, 0, 1])
        agen = a_exact(to) - a_exact(te)
        qd, qr, qp = (-agen / 3 - 5 * np.sum(wts * y * p2).real for y in (d, red, prt))
        P = r["Q2pieces"]
        B = pm + (W_R(np.array([xo]))[0] - W_R(np.array([xe]))[0]) / 3 + P["prod"] + P["clk"]
        print(f"({xe},{xo}): {pm:+.7f} | {qd:+.7f} | {qr:+.7f} | printed R {qp:+.7f} | B = {B:+.4f}", flush=True)
        assert abs(qd - pm) < 2e-5 and abs(qr - pm) < 2e-5 and abs(qp - pm) < 2e-5
    print("PASS reduced ray terms reproduce the map")


if __name__ == "__main__":
    part1()
    part2()
