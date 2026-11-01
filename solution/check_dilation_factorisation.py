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

X = sp.symbols('X1:4', real=True)
n = sp.symbols('n1:4', real=True)
c, lam, tau, eps, mu, m0, m2, m4 = sp.symbols('c lambda tau epsilon mu m0 m2 m4', positive=True)
Z = sp.Function('Z')


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
    coef = sp.symbols(f'a0:{len(monos)}')
    seed = lambda Y: sum(a * sp.Mul(*[Y[i] for i in m]) for a, m in zip(coef, monos))
    scaled = theta_lead(seed([lam * xi for xi in X]))
    target = lam**2 * theta_lead(seed(X)).subs({xi: lam * xi for xi in X}, simultaneous=True)
    assert sp.expand(sp.powsimp(sp.expand(scaled - target))) == 0
    print('[ok] constant shift of zeta = momentum dilation (exact in zeta, generic seed)')


def check_mean_factorisation():
    # g, its Hessian and gradient are independent at a point (isotropy); see check_report_identities.
    g = sp.Symbol('g')
    AA, BB = sp.Rational(2, 3) * m4, sp.Rational(10, 9) * m2**2     # <A:A>, <B:B>; <A:B> = 0
    gauss = sp.exp(-4 * mu) * sp.exp(8 * eps**2 * m0)                # <e^{-4(mu + eps g)}>
    C2 = sp.Rational(8, 1) * sp.pi / 75 * tau**4 * gauss * (eps**2 * AA + eps**4 * BB)
    target = gauss * sp.Rational(16, 225) * sp.pi * tau**4 * eps**2 * (m4 + sp.Rational(5, 3) * eps**2 * m2**2)
    assert sp.simplify(C2 - target) == 0
    # Moving the mean is the same as dilating the spectrum: m_j -> e^{-j mu} m_j.
    dilated = target.subs({mu: 0, m4: sp.exp(-4 * mu) * m4, m2: sp.exp(-2 * mu) * m2}, simultaneous=True)
    assert sp.simplify(dilated - target) == 0
    ratio = sp.simplify((target / (gauss * sp.Rational(16, 225) * sp.pi * tau**4 * eps**2 * m4)) - 1)
    assert sp.simplify(ratio - sp.Rational(5, 3) * eps**2 * m2**2 / m4) == 0
    print('[ok] mean of zeta factors out as <e^{-4 zeta_o}>; (5/3) eps^2 m2^2/m4 is convention independent')


def check_operator():
    k, s = sp.symbols('k s', positive=True)
    F = k**4
    # <exp(-s D)> for Gaussian s of variance v: average of F(e^{-s} k) over s
    v = eps**2 * m0
    avg = sp.integrate(F.subs(k, sp.exp(-s) * k) * sp.exp(-s**2 / (2 * v)) / sp.sqrt(2 * sp.pi * v),
                       (s, -sp.oo, sp.oo))
    assert sp.simplify(avg - sp.exp(8 * v) * k**4) == 0
    # exp(v D^2/2) on k^4: D^2 k^4 = 16 k^4
    assert sp.simplify(sp.exp(v / 2 * 16) - sp.exp(8 * v)) == 0
    print('[ok] Gaussian average of the dilation operator: exp(eps^2 m0 D^2/2) = e^{8 eps^2 m0} on degree 4')


if __name__ == '__main__':
    check_dilation()
    check_mean_factorisation()
    check_operator()
    print('DILATION FACTORISATION CHECKS PASSED (exact symbolic algebra)')
