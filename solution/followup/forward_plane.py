"""Forward test (no perturbative map formula): exact plane-symmetric Einstein-radiation evolution + exact rays.

Field solve, zeta = 2 eps cos y (q=1):
  nu    : hyperbolic equation (N), pseudo-spectral RK4;
  f'    : momentum constraint (M) at nonzero wavenumber (evolving (T) instead would excite the -K^2/3 scalar
          instability that the constraints exclude);
  <f'>  : the zero mode of (T), which has no constraint.
Initial data: exact first-order closed forms, exact K=0 pair kernels, gradient data for the 2q channel
(whose error reaches the K=0 channel only at O(eps^4)).

Rays are integrated backward from the observer at label y0, with e = -gamma^{-1/2} n, to the root of
ln t + nu = ln tbar_*. Observation occurs at fixed proper time after T = T_I on the worldline.
Theta is exactly normalised. Averaging over the observer label y0 isolates the zero-momentum channel.
"""

import numpy as np
from k0_closed import a_exact, ap_exact, fields, n2_exact, n2p_exact
from numpy.polynomial.legendre import leggauss

s3 = np.sqrt(3.0)


class Evol:
    def __init__(self, eps, tstart, tend, N=32, dt=0.01, mmax=10):
        self.eps, self.N = eps, N
        self.y = 2 * np.pi * np.arange(N) / N
        self.kk = np.fft.rfftfreq(N, 1 / N)
        self.zeta = 2 * eps * np.cos(self.y)
        self.zy = -2 * eps * np.sin(self.y)
        self.zyy = -2 * eps * np.cos(self.y)
        n = int(round((tend - tstart) / dt))
        self.dt = (tend - tstart) / n
        self.ts = tstart + self.dt * np.arange(n + 1)
        self.mmax = mmax
        coef = np.zeros((n + 1, 4, mmax + 1), complex)  # f, f', nu, nu'
        s = self.initial(tstart)
        coef[0] = self.coefs(s, tstart)
        self.ham = np.zeros(n + 1)
        for i in range(n):
            s = self.rk4(s, self.ts[i], self.dt)
            coef[i + 1] = self.coefs(s, self.ts[i + 1])
        self.c = coef

    def dy(self, u):
        return np.fft.irfft(1j * self.kk * np.fft.rfft(u), self.N)

    def dyy(self, u):
        return np.fft.irfft(-(self.kk**2) * np.fft.rfft(u), self.N)

    def initial(self, t):
        eps = self.eps
        F = fields(np.array([t]))
        h, n1, n1p = F["h"][0], F["n1"][0], F["n1p"][0]
        g = 2 * np.cos(self.y)
        c2 = 2 * np.cos(2 * self.y)
        f = eps * h * g + eps**2 * (a_exact(t) + t**2 * c2)
        nu = eps * n1 * g + eps**2 * (n2_exact(t) + t**2 / 6 * c2)
        nup = eps * n1p * g + eps**2 * (n2p_exact(t) + t / 3 * c2)
        return np.array([f, nu, nup]), eps**2 * ap_exact(t) + 0 * eps * F["hp"][0]

    def fprime(self, f, nu, nup, mean_fp, t):
        """Solve (M): (E w)_y = 3 E (2 nup_y - 2 nu_y/t - 2 nup nu_y), E = exp(3B + 2nu), with <w> = mean_fp."""
        ny = self.dy(nu)
        E = np.exp(3 * self.zeta - f / 2 + 2 * nu)
        rhs = 3 * E * (2 * self.dy(nup) - 2 * ny / t - 2 * nup * ny)
        R = np.fft.rfft(rhs)
        R[0] = 0  # compatibility: the zero mode of (M) is identically satisfied by the exact solution
        Phi = np.zeros_like(R)
        Phi[1:] = R[1:] / (1j * self.kk[1:])
        phi = np.fft.irfft(Phi, self.N)
        c = (mean_fp - np.mean(phi / E)) / np.mean(1 / E)
        return (c + phi) / E

    def rhs_all(self, s, m, t):
        f, nu, nup = s
        w = self.fprime(f, nu, nup, m, t)
        fy, fyy, ny, nyy = self.dy(f), self.dyy(f), self.dy(nu), self.dyy(nu)
        zy, zyy = self.zy, self.zyy
        ex = np.exp(-2 * f / 3 - 2 * self.zeta)
        Ff = ex * (-(fy**2) / 18 + fy * ny / 3 + fy * zy / 6 + fyy / 6 + 2 * ny**2 + 4 * ny * zy - 2 * nyy + zy**2 - zyy)
        mp = np.mean(-2 * (1 / t + nup) * w - 2 * Ff)
        Ay, By = zy + fy / 3, zy - fy / 6
        lap = ex * (nyy - Ay * ny + 2 * By * ny)
        grad = ex * ny**2
        nupp = (lap + 2 * grad - w**2 / 6 - 3 / t**2 * (np.exp(-2 * nu) - 1)) / 3
        return np.array([w, nup, nupp]), mp

    def rk4(self, sm, t, h):
        s, m = sm if isinstance(sm, tuple) else (sm, None)
        k1, l1 = self.rhs_all(s, m, t)
        k2, l2 = self.rhs_all(s + h / 2 * k1, m + h / 2 * l1, t + h / 2)
        k3, l3 = self.rhs_all(s + h / 2 * k2, m + h / 2 * l2, t + h / 2)
        k4, l4 = self.rhs_all(s + h * k3, m + h * l3, t + h)
        return s + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4), m + h / 6 * (l1 + 2 * l2 + 2 * l3 + l4)

    def coefs(self, sm, t):
        s, m = sm
        f, nu, nup = s
        w = self.fprime(f, nu, nup, m, t)
        arr = np.array([f, w, nu, nup])
        c = np.fft.rfft(arr, axis=-1) / self.N
        c[:, 1:] *= 2
        return c[:, : self.mmax + 1]

    # --- local-in-time cubic Lagrange interpolation of the Fourier coefficients
    def coef_at(self, t):
        u = (t - self.ts[0]) / self.dt
        i = int(np.clip(np.floor(u), 1, len(self.ts) - 3))
        s = u - i
        w = np.array(
            [-s * (s - 1) * (s - 2) / 6, (s + 1) * (s - 1) * (s - 2) / 2, -(s + 1) * s * (s - 2) / 2, (s + 1) * s * (s - 1) / 6]
        )
        return w[0] * self.c[i - 1] + w[1] * self.c[i] + w[2] * self.c[i + 1] + w[3] * self.c[i + 2]

    def eval_fields(self, t, y):
        c = self.coef_at(t)
        m = np.arange(self.mmax + 1)
        E = np.exp(1j * np.outer(y, m))
        vals = (E @ c.T).real  # (nray, 4)
        dv = (E * (1j * m) @ c.T).real
        return vals, dv

    def worldline(self, y0):
        m = np.arange(self.mmax + 1)
        return (self.c[:, 2] @ np.exp(1j * m * y0)).real


def geo_rhs(ev, t, st):
    y, ey, ep, intI = st
    vals, dv = ev.eval_fields(t, y)
    f, fp, nu, nup = vals.T
    fy, ny = dv[:, 0], dv[:, 2]
    z, zy = 2 * ev.eps * np.cos(y), -2 * ev.eps * np.sin(y)
    A, B = z + f / 3, z - f / 6
    At, Bt = fp / 3, -fp / 6
    Ay, By = zy + fy / 3, zy - fy / 6
    sig = At * np.exp(2 * A) * ey**2 + Bt * np.exp(2 * B) * ep**2
    dey = -2 * At * ey - Ay * ey**2 + np.exp(2 * B - 2 * A) * By * ep**2 + ey * sig
    dep = -2 * Bt * ep - 2 * By * ey * ep + ep * sig
    g = np.log(t) + nu
    gd = 1 / t + nup + ny * ey
    return np.array([ey, dey, dep, sig]), g, gd


def observer_time(ev, y0, tbar_I, tbar_o):
    ts = ev.ts
    nw = ev.worldline(y0)
    gI = np.log(ts) + nw - np.log(tbar_I)
    iI = np.searchsorted(gI, 0.0)
    # high-accuracy root and clock with cubic-interpolated nu
    lo, hi = ts[iI - 1], ts[iI]
    for _ in range(60):
        mid = (lo + hi) / 2
        c = ev.coef_at(mid)
        val = np.log(mid) + (c[2] @ np.exp(1j * np.arange(ev.mmax + 1) * y0)).real - np.log(tbar_I)
        lo, hi = (mid, hi) if val < 0 else (lo, mid)
    tI = (lo + hi) / 2

    # tau(t) = int t e^nu dt (a1 = 1): Simpson on the grid, with the interpolant on the partial cells
    def integ(a, b, n):
        tt = np.linspace(a, b, n + 1)
        vals = np.array([tv * np.exp((ev.coef_at(tv)[2] @ np.exp(1j * np.arange(ev.mmax + 1) * y0)).real) for tv in tt])
        return (b - a) / n / 3 * (vals[0] + vals[-1] + 4 * vals[1:-1:2].sum() + 2 * vals[2:-1:2].sum())

    target = (tbar_o**2 - tbar_I**2) / 2
    # Newton iteration on the proper-time equation, starting at the background value
    to = tbar_o
    n = max(400, int((tbar_o - tI) / (ev.dt * 2)) // 2 * 2)
    for _ in range(4):
        F = integ(tI, to, n) - target
        c = ev.coef_at(to)
        to = to - F / (to * np.exp((c[2] @ np.exp(1j * np.arange(ev.mmax + 1) * y0)).real))
    return to


def observe(ev, y0s, mus, tbar_e, tbar_I, tbar_o, step=0.04):
    ny, nm = len(y0s), len(mus)
    Y0, MU = np.meshgrid(y0s, mus, indexing="ij")
    Y0, MU = Y0.ravel(), MU.ravel()
    tos = {y0: observer_time(ev, y0, tbar_I, tbar_o) for y0 in y0s}
    out = np.zeros(ny * nm)
    # each y0 has its own t_o: integrate per observer
    for y0 in y0s:
        to = tos[y0]
        sel = y0 == Y0
        n_rays = sel.sum()
        y = np.full(n_rays, y0)
        mu = MU[sel]
        w = np.sqrt(1 - mu**2)
        vals, _ = ev.eval_fields(to, y)
        f = vals[:, 0]
        z = 2 * ev.eps * np.cos(y)
        ey, ep = -np.exp(-(z + f / 3)) * mu, -np.exp(-(z - f / 6)) * w
        st = np.array([y, ey, ep, np.zeros(n_rays)])
        t = to
        nsteps = int(np.ceil((to - tbar_e * 0.6) / step)) + 2
        h = -(to - tbar_e * 0.7) / nsteps
        done = np.zeros(n_rays, bool)
        Ival = np.zeros(n_rays)
        k1, g0, gd0 = geo_rhs(ev, t, st)
        lg = np.log(tbar_e)
        for _ in range(nsteps):
            k2, _, _ = geo_rhs(ev, t + h / 2, st + h / 2 * k1)
            k3, _, _ = geo_rhs(ev, t + h / 2, st + h / 2 * k2)
            k4, _, _ = geo_rhs(ev, t + h, st + h * k3)
            new = st + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            kn, g1, gd1 = geo_rhs(ev, t + h, new)
            crossing = (g1 <= lg) & ~done
            if crossing.any():
                # cubic Hermite root for g, then Hermite for I
                lo, hi = np.zeros(crossing.sum()), np.ones(crossing.sum())
                a0, a1, d0, d1 = g0[crossing] - lg, g1[crossing] - lg, gd0[crossing] * h, gd1[crossing] * h

                def H(s, v0, v1, m0, m1):
                    return (
                        (2 * s**3 - 3 * s**2 + 1) * v0
                        + (s**3 - 2 * s**2 + s) * m0
                        + (-2 * s**3 + 3 * s**2) * v1
                        + (s**3 - s**2) * m1
                    )

                for _b in range(60):
                    mid = (lo + hi) / 2
                    v = H(mid, a0, a1, d0, d1)
                    hi = np.where(v <= 0, mid, hi)
                    lo = np.where(v <= 0, lo, mid)
                s = (lo + hi) / 2
                Ival[crossing] = H(s, st[3][crossing], new[3][crossing], k1[3][crossing] * h, kn[3][crossing] * h)
                done |= crossing
            st, t, k1, g0, gd0 = new, t + h, kn, g1, gd1
            if done.all():
                break
        if not done.all():
            raise RuntimeError("some rays did not reach the emission surface")
        out[sel] = -Ival  # the march runs backward, so the accumulator holds int_{t_o}^{t_e} = -I
    return out.reshape(ny, nm), tos


def theta_avg(eps, xe, xo, xI_ratio=0.5, ny0=8, nmu=32, N=32, dt=0.01, step=0.04):
    te, to, tI = s3 * xe, s3 * xo, s3 * xe * xI_ratio
    tstart = min(tI, te) * 0.4
    ev = Evol(eps, tstart, to * (1 + 200 * eps) + 0.5, N=N, dt=dt)
    mus, wts = leggauss(nmu)
    wts = wts / 2
    y0s = 2 * np.pi * np.arange(ny0) / ny0
    Iv, _ = observe(ev, y0s, mus, te, tI, to, step=step)
    th = np.exp(-Iv) / np.sum(wts * np.exp(-Iv), axis=1, keepdims=True) - 1
    return mus, wts, th.mean(axis=0), ev


if __name__ == "__main__":
    from map_k0 import theta_pair

    def P2(mus, wts, y):
        return 5 * np.sum(wts * y * (3 * mus**2 - 1) / 2)

    cases = [(0.1, 0.2), (0.5, 1.0), (1, 2), (1, 3), (1, 5), (3, 10), (5, 20)]
    for xe, xo in cases:
        q2 = theta_pair(s3 * xe, s3 * xo, s3 * xe / 2, nmu=160)["Q2"]
        vals = []
        for eps in (2e-3, 1e-3):
            mus, wts, th, ev = theta_avg(eps, xe, xo, nmu=int(max(32, 5 * s3 * (xo - xe) + 16)))
            vals.append(P2(mus, wts, th) / eps**2)
        rich = (4 * vals[1] - vals[0]) / 3
        print(
            f"xe={xe} xo={xo}: forward 2e-3 {vals[0]:+.6f}, 1e-3 {vals[1]:+.6f}, extrap {rich:+.6f}; kernel {q2:+.6f}; diff {rich - q2:+.2e}",
            flush=True,
        )
