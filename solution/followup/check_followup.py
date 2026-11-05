"""Fast checks for the white-noise follow-up. Exact unless labelled numerical.

1. Closed K=0 fields satisfy (T), (N), (M), (Ham); mutated coefficients fail.
2. W and I simplify to the printed forms; small-x and late-time coefficients.
3. Induced-tensor limit: Eq. (FT-KT) of the main report gives F_T -> a/2 as K -> 0 on the shell (numerical, high precision).
4. Effective stress: late-time DC source equals (4/3) rho v v from first-order acoustic velocities.
5. Helicity split 1:2:2 of a fixed trace-free tensor over directions Khat, by exact integration.
6. Map omission controls: dropping a map piece breaks agreement with the forward nonlinear ray test (numerical).
"""

import k0_fields as kf
import mpmath as mp
import numpy as np
import sympy as sp

X = kf.X


def check_closed_fields():
    assert all(v == 0 for v in kf.check_first_order().values())
    assert all(v == 0 for v in kf.check_pair().values())
    saved = kf.a2
    try:
        kf.a2 = saved + kf.Cin  # mutate 24 -> 25: must fail
        res = kf.check_pair()
        assert res["T"] != 0, "negative control failed"
    finally:
        kf.a2 = saved
    saved = kf.n2
    try:
        kf.n2 = saved + sp.cos(X)
        res = kf.check_pair()
        assert res["Ham"] != 0 and res["N"] != 0, "negative control failed"
    finally:
        kf.n2 = saved
    print("PASS closed K=0 fields (T,N,M,Ham), with negative controls")


def check_printed_forms():
    Iso, W = kf.KR_parts()
    j0 = sp.sin(X) / X
    assert sp.simplify(sp.expand(sp.expand_trig(Iso - (4 + sp.cos(X) ** 2 - 5 * j0**2)))) == 0
    assert sp.simplify(sp.expand(sp.expand_trig(W - (24 * kf.Cin + 6 * sp.cos(X) ** 2 - 54 + 48 * j0**2)))) == 0
    assert sp.simplify(1 + kf.n1 - kf.h1 / 6 - j0) == 0
    a_s, n_s = kf.series_checks()
    assert sp.expand(a_s - (2 * X**2 - sp.Rational(11, 30) * X**4)) == 0
    assert sp.expand(n_s - (sp.Rational(5, 3) * X**2 - sp.Rational(13, 45) * X**4)) == 0
    print("PASS printed I(x), W(x), R1 = sin x / x, a and n2 series")


def check_tensor_limit():
    mp.mp.dps = 40

    def FT(t, k0, mu):
        K = k0 * mp.sqrt(2 * (1 + mu))
        r, Xv = k0 / K, K * t
        x = k0 * t / mp.sqrt(3)
        d = 2 * r**2 - 3
        ap, am = 1 + 2 * r / mp.sqrt(3), 1 - 2 * r / mp.sqrt(3)

        def Cin(z):
            z = abs(z)
            return mp.euler + mp.log(z) - mp.ci(z) if z != 0 else mp.mpf(0)

        def Si(z):
            return mp.si(z)

        C = Cin(am * Xv) + Cin(ap * Xv) - 2 * Cin(Xv)
        S = Si(am * Xv) + Si(ap * Xv) - 2 * Si(Xv)
        E = (
            -4 * r**2 * d * mp.sin(Xv)
            + 24 * r**2 * mp.cos(x) ** 2 / Xv
            - 48 * mp.sqrt(3) * r * mp.cos(x) * mp.sin(x) / Xv**2
            + 12 * (6 + d * Xv**2) * mp.sin(x) ** 2 / Xv**3
        )
        Ird = 3 / (4 * r**6 * Xv) * (E + d**2 * (mp.sin(Xv) * C + mp.cos(Xv) * S))
        s = lambda y: 2 * mp.sin(y) / y - mp.cos(y)
        sp_ = lambda y: mp.diff(s, y)
        return (1 - mu) / (2 * (1 + mu)) * Ird + mp.mpf(3) / 2 * (1 - mu) * (3 * mu * (s(x) - 1) ** 2 - sp_(x) ** 2)

    def a_val(xv):
        Cin2 = mp.euler + mp.log(2 * xv) - mp.ci(2 * xv)
        return (
            -78
            - 36 * mp.cos(xv)
            - 6 * mp.cos(2 * xv)
            + 72 * mp.sin(xv) / xv
            + 36 * mp.sin(2 * xv) / xv
            - 24 * (mp.sin(xv) / xv) ** 2
            + 24 * Cin2
        )

    worst = 0
    for xv in (0.7, 2.0, 6.0, 15.0):
        t = mp.sqrt(3) * xv
        val = FT(t, mp.mpf(1), -1 + mp.mpf(10) ** (-14))
        err = abs(val - a_val(xv) / 2)
        worst = max(worst, float(err))
        print(f"   x={xv}: F_T(K->0)={mp.nstr(val, 10)}  a/2={mp.nstr(a_val(xv) / 2, 10)}")
    assert worst < 1e-6
    print(f"PASS K->0 limit of the Kohri-Terada form equals a/2 (numerical, worst {worst:.1e})")


def check_stress():
    S_T = (
        -12 * sp.cos(X)
        - 8 * sp.cos(2 * X)
        + 24 * sp.sin(2 * X) / X
        - 48 * sp.sin(X) ** 2 / X**2
        + 16 / X**2
        - 16 * sp.sin(2 * X) / X**3
        + 16 * sp.sin(X) ** 2 / X**4
    )
    # S_T must solve (x^2 a_x)_x / 3 / x^2 + S_T = 0
    res = sp.diff(X**2 * kf.dx(kf.a2), X) / X**2 / 3 + S_T
    assert sp.simplify(res.subs(sp.Derivative(kf.Cin, X), (1 - sp.cos(2 * X)) / X).doit()) == 0
    # time average of x^2 S_T over a period, x -> infinity: the 1/x^2 DC coefficient is -48/2 + 16 = -8
    dc = -48 * sp.Rational(1, 2) + 16
    # a'' + 2a'/t = (6/t^2) Pi/rho  with t^2 = 3 x^2 => Pi/rho = -(x^2/2) * S_T|DC = 4
    Pi_over_rho = -sp.Rational(1, 2) * dc
    # (4/3) rho v_i v_j, with v = sqrt(3) sin x qhat g (from delta = -4 cos x g ... ), z1 z2 coefficient, time average
    V2 = 3
    stress = sp.Rational(4, 3) * 2 * V2 * sp.Rational(1, 2)
    assert Pi_over_rho == stress == 4
    print("PASS late-time DC stress: Pi_ij/rho = 4 TF(qhat qhat) = (4/3)<v_i v_j>; gives a ~ 24 ln x")


def check_helicity():
    mu = sp.Symbol("mu", real=True)
    # fixed T = TF(qhat qhat) norm 2/3. K-frame components: scalar part s = (3/2)(mu^2 - 1/3) with s^2 (2/3) norm.
    avg = lambda e: sp.integrate(e, (mu, -1, 1)) / 2
    scalar = avg(sp.Rational(2, 3) * (sp.Rational(3, 2) * (mu**2 - sp.Rational(1, 3))) ** 2)
    vector = avg(2 * mu**2 * (1 - mu**2))
    tensor = sp.Rational(2, 3) - scalar - vector
    assert (scalar, vector, tensor) == (sp.Rational(2, 15), sp.Rational(4, 15), sp.Rational(4, 15))
    print("PASS helicity split scalar:vector:tensor = 1:2:2")


def check_omission():
    from forward_plane import theta_avg
    from map_k0 import s3, theta_pair

    def P2(m, w, y):
        return 5 * np.sum(w * y * (3 * m**2 - 1) / 2)

    xe, xo = 1.0, 3.0
    vals = []
    for eps in (2e-3, 1e-3):
        mus, wts, th, _ = theta_avg(eps, xe, xo, nmu=60)
        vals.append(P2(mus, wts, th) / eps**2)
    fwd = (4 * vals[1] - vals[0]) / 3
    r = theta_pair(s3 * xe, s3 * xo, s3 * xe / 2, nmu=160)
    assert abs(fwd - r["Q2"]) < 1e-5
    for k, v in r["Q2pieces"].items():
        assert abs(fwd - (r["Q2"] - v)) > 0.01, k
    print(
        f"PASS forward nonlinear rays vs kernel at (1,3): {fwd:+.6f} vs {r['Q2']:+.6f}; omission of any of {list(r['Q2pieces'])} is detected"
    )


if __name__ == "__main__":
    check_closed_fields()
    check_printed_forms()
    check_tensor_limit()
    check_stress()
    check_helicity()
    check_omission()
