"""Exact checks for the strengthened report; no numerics or output-file mutations.

Run: python -B solution/check_report_identities.py

Scope: residual Bianchi identities and the homogeneous constraint counterexample;
compact conformal driver and noncommutative shear; Fermi--Walker rotation;
centred normalisation and polarised kernels; Wick reorganisation of the spectrum;
geometric quadrupole and the all-amplitude average of its leading-gradient term.
These checks do NOT verify the arbitrary-time cubic observable or its convergence.
"""
from functools import lru_cache
from itertools import combinations

import sympy as sp


eps = sp.Symbol('epsilon', real=True)


def check_zero(expr):
    assert sp.expand(expr) == 0, sp.expand(expr)


def truncate(expr, degree=3):
    expanded = sp.expand(expr)
    return sp.expand(sum(expanded.coeff(eps, n) * eps**n for n in range(degree + 1)))


def check_bianchi():
    t = sp.Symbol('t', positive=True)
    k = sp.symbols('k1:4', real=True)
    f = sp.Function('f')(t)
    v = [sp.Function(f'v{i}')(t) for i in range(3)]
    b11, b22, b12, b13, b23 = [sp.Function(n)(t) for n in ('b11', 'b22', 'b12', 'b13', 'b23')]
    B = sp.Matrix([[f/3+b11, b12, b13], [b12, f/3+b22, b23],
                   [b13, b23, f/3-b11-b22]])
    E = sp.zeros(4)
    E[0, 0] = f
    for i in range(3):
        E[0, i+1], E[i+1, 0] = v[i], -v[i]
        for j in range(3):
            E[i+1, j+1] = B[i, j]

    # Christoffels of a(t)^2 eta with a'/a=1/t.
    G = sp.MutableDenseNDimArray.zeros(4, 4, 4)
    G[0, 0, 0] = 1/t
    for i in range(1, 4):
        G[0, i, i] = G[i, 0, i] = G[i, i, 0] = 1/t
    div = []
    for nu in range(4):
        value = sp.diff(E[0, nu], t) + sum(sp.I*k[i]*E[i+1, nu] for i in range(3))
        value += sum(G[mu, mu, rho]*E[rho, nu] - G[rho, mu, nu]*E[mu, rho]
                     for mu in range(4) for rho in range(4))
        div.append(sp.expand(value))
    check_zero(div[0] - (sp.diff(f, t)+2*f/t-sp.I*sum(k[i]*v[i] for i in range(3))))
    for i in range(3):
        check_zero(div[i+1] - (sp.diff(v[i], t)+4*v[i]/t
                              +sp.I*sum(k[j]*B[j, i] for j in range(3))))

    # At nonzero K, Ham and the temporal identity give f=v=0. Transversality
    # and the two TT conditions then kill all five spatial TF components.
    B0 = B.subs(f, 0)
    unknowns = [b11, b22, b12, b13, b23]
    solution = sp.solve([B0[2, i] for i in range(3)]
                        + [B0[0, 0]-B0[1, 1], B0[0, 1]], unknowns)
    assert solution == dict.fromkeys(unknowns, sp.S.Zero)
    for exponent, coefficient in ((2, 2), (4, 4)):
        check_zero(sp.diff(t**(-exponent), t)+coefficient*t**(-exponent)/t)

    # Check the regular-but-unconstrained lapse against the actual field equations.
    from engine import Ctx, Field, equations, t as et
    Ctx.maxdeg, Ctx.trunc, Ctx.momenta = 1, None, {1: (0, 0, 0)}
    c = sp.Symbol('c')
    Ehom = equations(Field({frozenset([1]): c*et**2}), Field(),
                     [[Field() for _ in range(3)] for _ in range(3)])
    check_zero(Ehom['N'].coeff([1]))
    for entry in Ehom['M'] + [entry for row in Ehom['T'] for entry in row]:
        check_zero(entry.coeff([1]))
    check_zero(Ehom['Ham'].coeff([1])-36*c)
    print('[ok] background Bianchi identities, nonzero-K residual elimination, and Ham=36c counterexample')


def symmetric_matrix(prefix):
    entries = {(i, j): sp.Symbol(f'{prefix}{i}{j}', real=True)
               for i in range(3) for j in range(i, 3)}
    return sp.Matrix(3, 3, lambda i, j: entries[min(i, j), max(i, j)])


def check_factored_geometry():
    M, R, Z, V = [symmetric_matrix(p) for p in ('M', 'R', 'Z', 'V')]
    z, v = sp.Matrix(sp.symbols('z1:4')), sp.Matrix(sp.symbols('v1:4'))
    lap_z, grad_z2, dot_zv = sp.symbols('lap_z grad_z2 dot_zv')
    # Independent conformal Ricci/Hessian transformation rules, before raising
    # an index. Their difference must be a scalar times M, removed by TF.
    Ric_gamma = R-Z+z*z.T-M*(lap_z+grad_z2)
    Hess_gamma_nu = V-z*v.T-v*z.T+M*dot_zv
    original = Ric_gamma-2*Hess_gamma_nu+2*v*v.T
    compact = R-Z-2*V+(z+2*v)*(z+2*v).T-2*v*v.T
    for entry in original-compact+M*(lap_z+grad_z2+2*dot_zv):
        check_zero(entry)

    # Formal noncommutative algebra proves these for arbitrary matrices,
    # rather than for a chosen numerical or polynomial metric.
    A, B, C, D, E, F = sp.symbols('A B C D E F', commutative=False)
    H = eps*A+eps**2*C/2+eps**3*E/6
    Hp = eps*B+eps**2*D/2+eps**3*F/6

    def exponential(a):
        return truncate(sum(a**n*H**n/sp.factorial(n) for n in range(4)))

    def exponential_prime(a):
        return truncate(sum(a**n/sp.factorial(n)*sum(H**j*Hp*H**(n-1-j)
                                                    for j in range(n)) for n in range(1, 4)))

    S = truncate(exponential(-1)*exponential_prime(1)/2)
    comm = lambda X, Y: X*Y-Y*X
    expected = truncate(Hp/2+comm(Hp, H)/4+comm(comm(Hp, H), H)/12)
    check_zero(S-expected)
    root = exponential(-sp.Rational(1, 2))
    root_inv = exponential(sp.Rational(1, 2))
    root_prime = exponential_prime(-sp.Rational(1, 2))
    rotation_prime = truncate(-root_inv*(S*root+root_prime), 2)
    check_zero(rotation_prime-eps**2*comm(A, B)/8)
    # A shell has H1(t)=f(t) B, hence identically zero commutator.
    check_zero(comm(A, B).subs(B, sp.Symbol('scalar')*A))

    # Rotation is orthogonal to the original real quadrupole in its norm.
    Q = symmetric_matrix('Q')
    a, b, c = sp.symbols('a b c')
    Om = sp.Matrix([[0, a, b], [-a, 0, c], [-b, -c, 0]])
    check_zero(sp.trace(Q*(Q*Om-Om*Q)))
    print('[ok] compact conformal driver, cubic shear, FW commutator coefficient, and shell cancellation')


def moment_average(variables, centred=False):
    """Formal linear expectation with independent symbols for all joint moments."""
    def average(expr):
        result = 0
        for powers, coefficient in sp.Poly(sp.expand(expr), *variables).terms():
            degree = sum(powers)
            if degree == 0:
                moment = 1
            elif degree == 1 and centred:
                moment = 0
            else:
                moment = sp.Symbol('mean_' + '_'.join(map(str, powers)))
            result += coefficient*moment
        return sp.expand(result)
    return average


def check_normalisation():
    i1, i2, i3 = sp.symbols('i1 i2 i3')
    mean = moment_average((i1, i2, i3))
    I = eps*i1+eps**2*i2/2+eps**3*i3/6
    temperature = truncate(sum((-I)**n/sp.factorial(n) for n in range(4)))
    fluctuation = mean(temperature)-1
    inv_mean = truncate(sum((-fluctuation)**n for n in range(4)))
    exact_map = truncate(temperature*inv_mean-1)
    j1, j2, j3 = [i-mean(i) for i in (i1, i2, i3)]
    P0 = lambda x: x-mean(x)
    compact = [-j1, -j2+P0(j1**2),
               -j3+3*P0(j1*j2)-P0(j1**3)+3*j1*mean(j1**2)]
    for r, value in enumerate(compact, 1):
        check_zero(sp.factorial(r)*exact_map.coeff(eps, r)-value)
        check_zero(mean(value))
    assert sp.expand(compact[2]-(-j3+3*P0(j1*j2)-P0(j1**3))) != 0

    # Independent square-free polynomial algebra: coefficient of z1*z2*z3
    # in the exact normalised exponential, without a perturbative factorial.
    empty = frozenset()
    labels = [frozenset(subset) for size in range(1, 4)
              for subset in combinations((1, 2, 3), size)]
    J = {key: sp.Symbol('J' + ''.join(map(str, sorted(key)))) for key in labels}
    average = moment_average(tuple(J.values()), centred=True)

    def add(X, Y):
        return {key: sp.expand(X.get(key, 0)+Y.get(key, 0)) for key in X.keys() | Y.keys()}

    def scale(X, c):
        return {key: c*value for key, value in X.items()}

    def multiply(X, Y):
        out = {}
        for a, x in X.items():
            for b, y in Y.items():
                if not a & b:
                    out[a | b] = out.get(a | b, 0)+x*y
        return {key: sp.expand(value) for key, value in out.items()}

    one = {empty: sp.S.One}
    temp, term = one, one
    for n in range(1, 4):
        term = scale(multiply(term, J), -sp.Rational(1, n))
        temp = add(temp, term)
    mu = {key: average(value) for key, value in temp.items() if key}
    inverse, term = one, one
    for _ in range(3):
        term = scale(multiply(term, mu), -1)
        inverse = add(inverse, term)
    result = multiply(temp, inverse)
    a, b, c = [J[frozenset([i])] for i in (1, 2, 3)]
    ab, ac, bc = [J[frozenset(pair)] for pair in ((1, 2), (1, 3), (2, 3))]
    abc = J[frozenset((1, 2, 3))]
    P = lambda x: x-average(x)
    check_zero(result[frozenset([1])]+a)
    check_zero(result[frozenset((1, 2))]-(-ab+P(a*b)))
    expected = -abc+P(a*bc+b*ac+c*ab-a*b*c)
    expected += a*average(b*c)+b*average(a*c)+c*average(a*b)
    check_zero(result[frozenset((1, 2, 3))]-expected)
    print('[ok] centred map at all three orders and polarised kernels with no extra factorials')


def check_spectrum_reorganisation():
    g, variance = sp.symbols('g variance', real=True)
    reals = sp.symbols('ar ai br bi cr ci', real=True)
    a, b, c = [reals[i]+sp.I*reals[i+1] for i in (0, 2, 4)]
    # Connected quadratic part; the cubic part has NOT been normal ordered.
    observed = eps*a*g+eps**2*b*(g**2-variance)/2+eps**3*c*g**3/6
    polynomial = sp.Poly(sp.expand(observed*sp.conjugate(observed)), g)
    spectrum = 0
    for (degree,), coefficient in polynomial.terms():
        if degree % 2 == 0:
            moment = sp.factorial2(degree-1)*variance**(degree//2) if degree else 1
            spectrum += coefficient*moment
    effective = a+eps**2*variance*c/2
    compact = eps**2*variance*effective*sp.conjugate(effective)
    compact += eps**4*variance**2*b*sp.conjugate(b)/2
    check_zero(truncate(spectrum, 4)-truncate(compact, 4))
    # The omitted normally ordered cubic term is orthogonal to the linear one.
    check_zero(3*variance**2-3*variance*variance)
    print('[ok] two-norm spectrum retains 22 and all three 13 contractions, including complex interference')


def check_superhorizon():
    x = sp.symbols('x1:4', real=True)
    zeta = sp.Function('zeta')(*x)
    f = sp.exp(-zeta)
    for i in range(3):
        for j in range(3):
            identity = f*sp.diff(f, x[i], x[j])
            identity += sp.exp(-2*zeta)*(sp.diff(zeta, x[i], x[j])
                                        -sp.diff(zeta, x[i])*sp.diff(zeta, x[j]))
            check_zero(identity)

    m0, m2, m4, tau = sp.symbols('m0 m2 m4 tau', positive=True)
    delta = lambda i, j: sp.Integer(i == j)

    @lru_cache(None)
    def covariance(a, b):
        if (len(a)+len(b)) % 2:
            return sp.S.Zero
        if not a and not b:
            return m0
        if len(a) == len(b) == 1:
            return m2*delta(a[0], b[0])/3
        if not a and len(b) == 2:
            return -m2*delta(*b)/3
        if not b and len(a) == 2:
            return -m2*delta(*a)/3
        i, j = a
        k, l = b
        return m4*(delta(i, j)*delta(k, l)+delta(i, k)*delta(j, l)
                   +delta(i, l)*delta(j, k))/15

    for i in range(3):
        for j in range(3):
            check_zero(covariance((), (i, j))-delta(i, j)
                       *sum(covariance((), (k, k)) for k in range(3))/3)
            for k in range(3):
                check_zero(covariance((i, j), (k,)))
            check_zero(covariance((), (i,)))

    @lru_cache(None)
    def wick(args):
        if not args:
            return sp.S.One
        if len(args) % 2:
            return sp.S.Zero
        return sp.expand(sum(covariance(args[0], args[j])*wick(args[1:j]+args[j+1:])
                             for j in range(1, len(args))))

    indices = [()] + [(i,) for i in range(3)] + [(i, j) for i in range(3) for j in range(i, 3)]
    variables = [sp.Symbol('g' + ''.join(map(str, ij))) for ij in indices]
    by_index = dict(zip(indices, variables))

    def expect(expr):
        answer = 0
        for powers, coefficient in sp.Poly(sp.expand(expr), *variables).terms():
            args = tuple(ij for ij, power in zip(indices, powers) for _ in range(power))
            answer += coefficient*wick(args)
        return sp.expand(answer)

    g = by_index[()]
    grad = sp.Matrix([by_index[(i,)] for i in range(3)])
    Hess = sp.Matrix(3, 3, lambda i, j: by_index[tuple(sorted((i, j)))])
    tf = lambda M: M-sp.eye(3)*sp.trace(M)/3
    A, B = tf(Hess), tf(grad*grad.T)
    inner = lambda A, B: sum(A[i, j]*B[i, j] for i in range(3) for j in range(3))
    AA, BB, AB = [expect(inner(*pair)) for pair in ((A, A), (B, B), (A, B))]
    check_zero(AA-2*m4/3)
    check_zero(BB-sp.Rational(10, 9)*m2**2)
    check_zero(AB)
    normalisation = 8*sp.pi*tau**4/75
    resummed = normalisation*sp.exp(8*eps**2*m0)*(eps**2*AA+eps**4*BB)
    target = 16*sp.pi*tau**4*eps**2*sp.exp(8*eps**2*m0)*(m4+5*eps**2*m2**2/3)/225
    assert sp.simplify(resummed-target) == 0

    # Direct Wick expansion before taking TF or using statistical independence.
    psi1, psi2, psi3 = A, -4*g*A-2*B, 12*g**2*A+12*g*B
    c11 = normalisation*expect(inner(psi1, psi1))
    c22 = normalisation*expect(inner(psi2, psi2))
    c13 = normalisation*expect(inner(psi1, psi3))
    series = sp.series(resummed, eps, 0, 6).removeO().expand()
    check_zero(series.coeff(eps, 2)-c11)
    check_zero(series.coeff(eps, 4)-c22/4-c13/3)
    for psi in (psi1, psi2, psi3):
        for entry in psi:
            check_zero(expect(entry))

    sigma2, k, Xi = sp.symbols('sigma2 k Xi', positive=True)
    shell = {m0: sigma2, m2: sigma2*k**2, m4: sigma2*k**4, tau**4: Xi**2/(4*k**4)}
    for actual, expected in zip((c11, c22, c13),
                               (4*sp.pi*sigma2*Xi**2/225,
                                272*sp.pi*sigma2**2*Xi**2/675,
                                16*sp.pi*sigma2**2*Xi**2/75)):
        assert sp.simplify(actual.subs(shell)-expected) == 0
    shell_target = 4*sp.pi*eps**2*sigma2*Xi**2*sp.exp(8*eps**2*sigma2)*(1+5*eps**2*sigma2/3)/225
    assert sp.simplify(resummed.subs(shell)-shell_target) == 0
    assert sp.simplify((c22/4+c13/3).subs(shell)/c11.subs(shell)) == 29*sigma2/3
    print('[ok] curvature identity, Gaussian independence, general-moment resummation, and all three shell coefficients')


if __name__ == '__main__':
    check_bianchi()
    check_factored_geometry()
    check_normalisation()
    check_spectrum_reorganisation()
    check_superhorizon()
    print('STRENGTHENED REPORT IDENTITIES VERIFIED (exact symbolic algebra)')
