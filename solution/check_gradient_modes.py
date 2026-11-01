"""Exact gradient-mode checks motivated by Mirbabayi & Zaldarriaga, arXiv:1409.4777v2.

The physical calculation remains radiation-only. No matter-era transfer coefficient
is imported. These tests check the special-conformal null solution exactly, its
quadratic initial completion, the plane-wave soft generator, and the leading-gradient
quadratic observable. They do NOT substitute the arbitrary-time quadratic kernel into
the Ward identity; that identity is derived by covariance in the report.

Run: python -B solution/check_gradient_modes.py
"""
import sympy as sp


def zero(expr):
    assert sp.expand(expr) == 0, sp.expand(expr)


def tf(matrix):
    return matrix-sp.eye(3)*sp.trace(matrix)/3


def check_exact_null_profile():
    x = sp.Matrix(sp.symbols('x1:4', real=True))
    b = sp.Matrix(sp.symbols('b1:4', real=True))
    x2, b2, bx = x.dot(x), b.dot(b), b.dot(x)
    denominator = 1-2*bx+b2*x2
    numerator = x-b*x2
    dD = sp.Matrix([sp.diff(denominator, xi) for xi in x])
    # Jacobian of Y=(X-b X^2)/D, with a common D^2 denominator.
    jac_num = denominator*numerator.jacobian(x)-numerator*dD.T
    for entry in jac_num.T*jac_num-denominator**2*sp.eye(3):
        zero(entry)
    origin = dict.fromkeys(x, 0)
    assert numerator.subs(origin) == sp.zeros(3, 1)
    assert jac_num.subs(origin) == sp.eye(3)
    assert denominator.subs(origin) == 1

    # Exact Ricci of gamma=f^{-2} delta, multiplied by f^2:
    # f Hess(f) + delta [f Lap(f) - 2 |grad f|^2].
    Hess_f = sp.hessian(denominator, x)
    Ric_num = denominator*Hess_f+sp.eye(3)*(denominator*sp.trace(Hess_f)-2*dD.dot(dD))
    for entry in Ric_num:
        zero(entry)
    assert tf(Hess_f) == sp.zeros(3)
    t = sp.Symbol('t', positive=True)
    gamma = sp.eye(3)/denominator**2
    assert gamma.diff(t) == sp.zeros(3)
    # Thus S=sigma=0 and Ric=0. With nu=0 every exact field residual
    # vanishes, I=0 on every ray, and the normalised map is identically zero.
    print('[ok] exact special-conformal Jacobian, flat seed, fixed observer/frame, and zero shear')


def check_completion_and_physical_gradient():
    x = sp.Matrix(sp.symbols('X1:4', real=True))
    p = sp.Matrix(sp.symbols('p1:4', real=True))
    q = sp.Matrix(sp.symbols('q1:4', real=True))
    n = sp.Matrix(sp.symbols('n1:4', real=True))
    eps, tau = sp.symbols('epsilon tau', real=True)
    denominator = 1-eps*p.dot(x)+eps**2*p.dot(p)*x.dot(x)/4
    zeta_ad = sp.series(-sp.log(denominator), eps, 0, 3).removeO().expand()
    completion = p.dot(x)**2/2-p.dot(p)*x.dot(x)/4
    zero(zeta_ad.coeff(eps, 1)-p.dot(x))
    zero(zeta_ad.coeff(eps, 2)-completion)
    # Our second Taylor coefficient is twice the displayed quadratic bracket.
    assert sp.diff(zeta_ad, eps, 2).subs(eps, 0).expand() == sp.expand(2*completion)
    for entry in tf(sp.hessian(completion, x)-p*p.T):
        zero(entry)
    # zeta=eps p.X alone does not have the cancelling Hessian.
    quad = tau**2*eps**2*(n.T*tf(p*p.T)*n)[0]
    expected = tau**2*eps**2*(p.dot(n)**2-p.dot(p)*n.dot(n)/3)
    zero(quad-expected)
    assert quad != 0
    # Only the INITIAL PROFILE is converted to the paper's matter-era notation.
    substitutions = {p[i]: -sp.Rational(5, 3)*q[i] for i in range(3)}
    paper_phi2 = -sp.Rational(3, 5)*completion.subs(substitutions)
    zero(paper_phi2-(sp.Rational(5, 12)*q.dot(q)*x.dot(x)-sp.Rational(5, 6)*q.dot(x)**2))
    print('[ok] nonlinear completion, nonzero uncompensated gradient, and the paper\'s primordial phi2')


def check_soft_generator():
    x = sp.Matrix(sp.symbols('x1:4', real=True))
    p = sp.Matrix(sp.symbols('p1:4', real=True))
    k = sp.Matrix(sp.symbols('k1:4', real=True))
    q = sp.Matrix(sp.symbols('q1:4', real=True))
    n = sp.Matrix(sp.symbols('n1:4', real=True))
    tau = sp.Symbol('tau', real=True)

    def D(i, expression):
        return sum(k[j]*sp.diff(expression, k[i], k[j]) for j in range(3)) \
            - k[i]*sum(sp.diff(expression, kj, 2) for kj in k)/2

    phase = sp.exp(sp.I*k.dot(x))
    xi = p.dot(x)*x-x.dot(x)*p/2
    deformation = sum(xi[i]*sp.diff(phase, x[i]) for i in range(3))
    zero(deformation+sp.I*sum(p[i]*D(i, phase) for i in range(3)))
    soft_phase = sp.exp(sp.I*q.dot(x))
    at_q0 = dict.fromkeys(q, 0)
    zero(-sp.I*sum(p[i]*sp.diff(soft_phase, q[i]) for i in range(3)).subs(at_q0)-p.dot(x))

    # Extract the two plane-wave coefficients DIRECTLY from the local
    # curvature observable, rather than assuming the proposed Ward identity.
    a, b = sp.symbols('a b')
    zeta_at_origin = a+b
    grad_at_origin = sp.I*(a*k+b*q)
    Hess_at_origin = -a*k*k.T-b*q*q.T
    contracted = (n.T*tf(Hess_at_origin-grad_at_origin*grad_at_origin.T)*n)[0]
    theta = -tau**2*sp.exp(-2*zeta_at_origin)*contracted
    K1 = sp.expand(sp.diff(theta, a).subs({a: 0, b: 0}))
    K2 = sp.expand(sp.diff(theta, a, b).subs({a: 0, b: 0}))
    dilation = sum(k[i]*sp.diff(K1, k[i]) for i in range(3))
    zero(K2.subs(at_q0)+dilation)
    for i in range(3):
        zero(sp.diff(K2, q[i]).subs(at_q0)+D(i, K1))
    remainder = sp.expand(K2+dilation+sum(q[i]*D(i, K1) for i in range(3)))
    zero(remainder.subs(at_q0))
    for qi in q:
        zero(sp.diff(remainder, qi).subs(at_q0))
    # Negative control: incorrectly dropping grad(zeta)^2 breaks the test.
    wrong = -tau**2*sp.exp(-2*zeta_at_origin)*(n.T*tf(Hess_at_origin)*n)[0]
    wrong_K2 = sp.expand(sp.diff(wrong, a, b).subs({a: 0, b: 0}))
    assert sp.expand(sp.diff(wrong_K2, q[0]).subs(at_q0)+D(0, K1)) != 0
    print('[ok] plane-wave generator, leading-gradient Ward identity, O(q^2) remainder, and negative control')


if __name__ == '__main__':
    check_exact_null_profile()
    check_completion_and_physical_gradient()
    check_soft_generator()
    print('GRADIENT-MODE CHECKS PASSED (exact symbolic algebra; no matter-era transfer imported)')
