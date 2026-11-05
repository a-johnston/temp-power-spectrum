"""Complete second-order temperature kernel in the K->0 channel, from the report's map formulas.

q = 1, x = t/sqrt3. Seed zeta = z1 e^{iy} + z2 e^{-iy}; Theta_ab(mu) is the full-series z1 z2 coefficient
for an observer at the origin, mu = qhat.nhat. Fields: closed forms of k0_fields.py (via k0_closed.py).
Exposes the additive pieces of the l=2 coefficient so that the generated part can be separated.

The sky is processed in chunks of `chunk` directions: only 1D per-direction results are kept, so peak
memory is about 30 complex arrays of shape (chunk, nt+1) and does not grow with the number of directions.
"""

import numpy as np
from k0_closed import a_exact, fields
from numpy.polynomial.legendre import leggauss, legval

s3 = np.sqrt(3.0)


def cum_from_right(y, dt):
    """int_t^{t_o} y dt' along the last axis on a uniform grid: trapezoid with the Gregory end correction."""
    c = np.zeros_like(y)
    c[..., :-1] = np.cumsum(((y[..., :-1] + y[..., 1:]) / 2)[..., ::-1], axis=-1)[..., ::-1] * dt
    g = np.gradient(y, dt, axis=-1, edge_order=2)
    return c - (dt**2 / 12) * (g[..., -1:] - g)


def _directions(mus, t, dt, to, h, hp):
    """Per-direction first-order exponents I_a, I_b and the ray integral, for one chunk of mu."""
    lam = to - t
    mu = mus[:, None]
    w = np.sqrt(1 - mu**2)

    def mode(s):
        """first-order objects for mode k = s y; arrays have shape (nmu, nt+1)"""
        ph = np.exp(1j * s * mu * lam)
        n = (mu, w)
        Hn = (2 / 3 * mu, -1 / 3 * w)  # TF(qq) contracted with n, in-plane components
        Hnn = mu**2 - 1 / 3
        kn = s * mu
        sig_nn = 0.5 * hp * Hnn * ph
        # A_1 = 2 S n - Gamma(n,n) - n sig_nn, S = H'/2
        ry = hp * Hn[0] * ph - (1j * (2 * kn * n[0] - s) + 0.5j * h * (2 * kn * Hn[0] - s * Hnn)) * ph - n[0] * sig_nn
        rp = hp * Hn[1] * ph - (1j * (2 * kn * n[1]) + 0.5j * h * (2 * kn * Hn[1])) * ph - n[1] * sig_nn
        eoy = n[0] + 0.5 * h[-1] * Hn[0]  # e_a(t_o) = (1 + H_a/2) n, zeta_a = 1, unit phase at the observer
        eop = n[1] + 0.5 * h[-1] * Hn[1]
        By, Bp = cum_from_right(ry, dt), cum_from_right(rp, dt)
        Cy = cum_from_right(ry * t, dt)
        Xy = eoy * (t - to) + Cy - t * By  # X = e_o (t-t_o) + int_t^{t_o} (t'-t) A dt'
        vy = eoy - By - (n[0] + 0.5 * h * Hn[0]) * ph  # orthonormal v_a = e_a - (1 + H_a/2) n e^{i lam k.n}
        vp = eop - Bp - (n[1] + 0.5 * h * Hn[1]) * ph
        return dict(ph=ph, sig_nn=sig_nn, Xy=Xy, vy=vy, vp=vp, I=cum_from_right(sig_nn, dt)[:, 0], Hn=Hn)

    a, b = mode(+1), mode(-1)

    def H1p_nv(mo, other):  # H_b'(n, v_a)
        return hp * (mo["Hn"][0] * other["vy"] + mo["Hn"][1] * other["vp"]) * mo["ph"]

    ray_integrand = 1j * (-1) * a["Xy"] * b["sig_nn"] + 1j * (+1) * b["Xy"] * a["sig_nn"] - H1p_nv(b, a) - H1p_nv(a, b)
    Iray = cum_from_right(ray_integrand, dt)[:, 0]
    return a["I"], b["I"], Iray, a["sig_nn"][:, 0] * b["ph"][:, 0], b["sig_nn"][:, 0] * a["ph"][:, 0]


def theta_pair(te, to, tI, nmu=None, nt=None, chunk=8):
    lam_max = to - te
    nmu = nmu or int(max(60, 3 * lam_max + 40))
    nt = nt or int(max(2000, 60 * lam_max))
    mus, wts = leggauss(nmu)
    wts = wts / 2  # sky mean of a function of mu
    t = np.linspace(te, to, nt + 1)
    dt = t[1] - t[0]
    F = fields(t)
    h, hp, n1 = F["h"], F["hp"], F["n1"]
    res = [_directions(mus[i : i + chunk], t, dt, to, h, hp) for i in range(0, nmu, chunk)]
    Ia, Ib, Iray, sa_pb, sb_pa = (np.concatenate(c) for c in zip(*res, strict=True))
    sigo = 0.5 * hp[-1] * (mus**2 - 1 / 3)
    dIa = sigo + 1j * mus * Ia  # d/dt_o at fixed t_e
    dIb = sigo - 1j * mus * Ib

    def sky(fn):
        return np.sum(wts * fn)

    # emission root: -[F_a d_b + F_b d_a] at t_e, with d = -t_e nu_1(t_e) x phase
    Iem = te * n1[0] * (sa_pb + sb_pa)
    # generated field term
    agen = a_exact(to) - a_exact(te)
    Igen = 0.5 * agen * (mus**2 - 1 / 3)
    # clock: D_a = D_b = D for this pair
    n1I = fields(np.array([tI]))["n1"][0]
    dI = -tI * n1I
    tt = np.linspace(tI, to, 20001)
    FF = fields(tt)["n1"] * tt
    W_int = (tt[1] - tt[0]) / 3 * (FF[0] + FF[-1] + 4 * FF[1:-1:2].sum() + 2 * FF[2:-1:2].sum())
    D = (tI * dI - W_int) / to

    pieces = dict(
        gen=-(Igen - sky(Igen)),
        ray=-(Iray - sky(Iray)),
        em=-(Iem - sky(Iem)),
        prod=Ia * Ib - sky(Ia) * Ib - sky(Ib) * Ia + 2 * sky(Ia) * sky(Ib) - sky(Ia * Ib),
        clk=-D * (dIa - sky(dIa) + dIb - sky(dIb)),
    )
    total = sum(pieces.values())

    def proj(y, ell=2):
        c = np.zeros(ell + 1)
        c[ell] = 1
        return (2 * ell + 1) * np.sum(wts * y * legval(mus, c)).real

    return dict(
        mu=mus,
        w=wts,
        total=total,
        pieces=pieces,
        Q2=proj(total),
        Q2pieces={k: proj(v) for k, v in pieces.items()},
        proj=proj,
        imag=float(np.max(np.abs(total.imag))),
        D=D,
    )


if __name__ == "__main__":
    from k0_closed import W_R

    W = lambda x: W_R(np.array([float(x)]))[0]
    print("xe xo Q2 | -dW/3 rem em prod clk | l=4 coefficient")
    for xe, xo in [(0.1, 0.2), (0.5, 1), (1, 2), (1, 3), (1, 5), (3, 10), (5, 20), (5, 40), (5, 80), (5, 160)]:
        te, to = s3 * xe, s3 * xo
        r = theta_pair(te, to, te / 2, nmu=int(5 * s3 * (xo - xe) + 60), nt=int(100 * s3 * (xo - xe)) + 2000)
        p = r["Q2pieces"]
        gW = -(W(xo) - W(xe)) / 3
        rem = p["gen"] + p["ray"] - gW
        print(
            f"{xe:5} {xo:5} {r['Q2']:+10.6f} | {gW:+10.6f} {rem:+10.6f} {p['em']:+10.6f} {p['prod']:+10.6f} {p['clk']:+10.6f} | {r['proj'](r['total'], 4):+10.6f} | imag {r['imag']:.0e}",
            flush=True,
        )
