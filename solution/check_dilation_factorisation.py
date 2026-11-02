"""Dilation factorisation of the leading-gradient quadrupole (exact symbolic checks).

1. Theta_lead[zeta] = tau^2 f (d_i d_j f)_TF n^i n^j, f = e^{-zeta}, obeys
     Theta_lead[zeta + c] = e^{-2c} Theta_lead[zeta]   and   Theta_lead[zeta(lambda X)] = lambda^2 Theta_lead[zeta](lambda X),
   so a constant shift of zeta is exactly the momentum dilation k -> e^{-c} k.
2. With zeta = mu + eps g (g Gaussian, zero mean), the exact Gaussian average factorises as
     C_2,lead = <e^{-4 zeta_o}> * (16 pi/225) tau^4 eps^2 (m_4 + (5/3) eps^2 m_2^2),
     <e^{-4 zeta_o}> = e^{-4 mu + 8 eps^2 m_0},
   so the mean (the residual dilation convention) multiplies the whole spectrum and the ratio
   (5/3) eps^2 m_2^2 / m_4 is convention independent.
3. The dilation factor is the Gaussian average of the dilation operator: on functions of
   homogeneous degree 4 in the momenta, <exp(-eps g_o D)> = exp(eps^2 m_0 D^2 / 2) = e^{8 eps^2 m_0},
   D = k d/dk.
"""

import sympy as sp

X = sp.symbols("X1:4", real=True)
n = sp.symbols("n1:4", real=True)
c, lam, tau, eps, mu, m0, m2, m4 = sp.symbols("c lambda tau epsilon mu m0 m2 m4", positive=True)
Z = sp.Function("Z")


def theta_lead(zeta):
    f = sp.exp(-zeta)
    H = sp.Matrix(3, 3, lambda i, j: sp.diff(f, X[i], X[j]))
    TF = H - sp.eye(3) * H.trace() / 3
    nv = sp.Matrix(n)
    return tau**2 * f * (nv.T * TF * nv)[0]


def check_dilation():
    zeta = Z(*X)
    lhs = theta_lead(zeta + c)
    assert sp.simplify(lhs - sp.exp(-2 * c) * theta_lead(zeta)) == 0
    # scaling: generic cubic polynomial seed (all coefficients independent symbols)
    import itertools

    monos = [m for d in range(4) for m in itertools.combinations_with_replacement(range(3), d)]
    coef = sp.symbols(f"a0:{len(monos)}")
    seed = lambda Y: sum(a * sp.Mul(*[Y[i] for i in m]) for a, m in zip(coef, monos, strict=True))
    scaled = theta_lead(seed([lam * xi for xi in X]))
    target = lam**2 * theta_lead(seed(X)).subs({xi: lam * xi for xi in X}, simultaneous=True)
    assert sp.expand(sp.powsimp(sp.expand(scaled - target))) == 0
    print("[ok] constant shift of zeta = momentum dilation (exact in zeta, generic seed)")


def check_mean_factorisation():
    # g, its Hessian and gradient are independent at a point (isotropy); see check_report_identities.
    AA, BB = sp.Rational(2, 3) * m4, sp.Rational(10, 9) * m2**2  # <A:A>, <B:B>; <A:B> = 0
    gauss = sp.exp(-4 * mu) * sp.exp(8 * eps**2 * m0)  # <e^{-4(mu + eps g)}>
    C2 = sp.Rational(8, 1) * sp.pi / 75 * tau**4 * gauss * (eps**2 * AA + eps**4 * BB)
    target = gauss * sp.Rational(16, 225) * sp.pi * tau**4 * eps**2 * (m4 + sp.Rational(5, 3) * eps**2 * m2**2)
    assert sp.simplify(C2 - target) == 0
    # Moving the mean is the same as dilating the spectrum: m_j -> e^{-j mu} m_j.
    dilated = target.subs({mu: 0, m4: sp.exp(-4 * mu) * m4, m2: sp.exp(-2 * mu) * m2}, simultaneous=True)
    assert sp.simplify(dilated - target) == 0
    ratio = sp.simplify((target / (gauss * sp.Rational(16, 225) * sp.pi * tau**4 * eps**2 * m4)) - 1)
    assert sp.simplify(ratio - sp.Rational(5, 3) * eps**2 * m2**2 / m4) == 0
    print("[ok] mean of zeta factors out as <e^{-4 zeta_o}>; (5/3) eps^2 m2^2/m4 is convention independent")


def check_operator():
    k, s = sp.symbols("k s", positive=True)
    F = k**4
    # <exp(-s D)> for Gaussian s of variance v: average of F(e^{-s} k) over s
    v = eps**2 * m0
    avg = sp.integrate(F.subs(k, sp.exp(-s) * k) * sp.exp(-(s**2) / (2 * v)) / sp.sqrt(2 * sp.pi * v), (s, -sp.oo, sp.oo))
    assert sp.simplify(avg - sp.exp(8 * v) * k**4) == 0
    # exp(v D^2/2) on k^4: D^2 k^4 = 16 k^4
    assert sp.simplify(sp.exp(v / 2 * 16) - sp.exp(8 * v)) == 0
    print("[ok] Gaussian average of the dilation operator: exp(eps^2 m0 D^2/2) = e^{8 eps^2 m0} on degree 4")


def check_centred_kernels():
    """Superhorizon kernels of the observer-centred field (Prop. obs-dilation):
    K2^c = K2 + D K1(k) + D K1(q) = -tau^2 TF(kq+qk)_nn  and  K3^c(k, q, -q) = 0."""
    k = sp.Matrix(sp.symbols("k1:4", real=True))
    q = sp.Matrix(sp.symbols("q1:4", real=True))
    nv = sp.Matrix(n)
    z = sp.symbols("z1:4")
    lam = sp.Symbol("lam", positive=True)

    def kernels(kv):
        # leading map at X = 0 for zeta = sum_a z_a e^{i k_a.X}: zeta = sum z, zeta_i = i sum z k, zeta_ij = -sum z k k
        zeta = sum(z)
        zi = sp.I * sum((za * ka for za, ka in zip(z, kv, strict=True)), sp.zeros(3, 1))
        zij = -sum((za * ka * ka.T for za, ka in zip(z, kv, strict=True)), sp.zeros(3))
        TF = lambda M: M - sp.eye(3) * M.trace() / 3
        th = -(tau**2) * sp.exp(-2 * eps * zeta) * (nv.T * TF(eps * zij - eps**2 * zi * zi.T) * nv)[0]
        th = sp.expand(sp.series(th, eps, 0, 4).removeO())
        return th

    th = kernels([k, q, -q])
    coef = lambda mons: sp.expand(sp.Poly(th, *z, eps).coeff_monomial(mons))
    z1, z2, z3 = z
    K1k = coef(z1 * eps)
    K1q = coef(z2 * eps)
    K2kq = coef(z1 * z2 * eps**2)
    K2kmq = coef(z1 * z3 * eps**2)
    K2qmq = coef(z2 * z3 * eps**2)
    K3 = coef(z1 * z2 * z3 * eps**3)
    Dil = lambda e: sp.expand(
        sp.diff(e.subs({**{c: lam * c for c in k}, **{c: lam * c for c in q}}, simultaneous=True), lam).subs(lam, 1)
    )
    K1mq = K1q  # even in q
    K2c = lambda K2, a, b: sp.expand(K2 + Dil(a) + Dil(b))
    K2c_kq = K2c(K2kq, K1k, K1q)
    target = -(tau**2) * sp.expand((nv.T * ((k * q.T + q * k.T) - sp.eye(3) * 2 * (k.T * q)[0] / 3) * nv)[0])
    assert sp.simplify(K2c_kq - target) == 0
    K3c = sp.expand(
        K3
        + Dil(K2c_kq)
        + Dil(K2c(K2kmq, K1k, K1mq))
        + Dil(K2c(K2qmq, K1q, K1mq))
        - Dil(Dil(K1k))
        - Dil(Dil(K1q))
        - Dil(Dil(K1mq))
    )
    assert sp.simplify(K3c) == 0
    print("[ok] superhorizon centred kernels: K2^c = -tau^2 TF(kq+qk)_nn (pure gradient^2), K3^c(k,q,-q) = 0")


if __name__ == "__main__":
    check_dilation()
    check_mean_factorisation()
    check_operator()
    check_centred_kernels()
    print("DILATION FACTORISATION CHECKS PASSED (exact symbolic algebra)")
