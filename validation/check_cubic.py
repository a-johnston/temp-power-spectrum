"""Exact checks of the new cubic *dynamical* construction, not of the sky map."""
import unittest
import sympy as s

from amplitude_jets import (
    AmplitudeJet, add, coefficient, connection_linear, connection_product,
    derivative, from_modes, identity, product, ricci_derivative, ricci_product,
    scale, trace,
)
from cubic_source import cubic_regular_source, cubic_source
from initial_data import (I3, linear_boundary_jet, linear_mode,
                          quadratic_boundary_coefficient, spatial_operator, time_operator)


def exact_rhs_polarization(time, momenta, linear, quadratic):
    """Expand the full nonlinear metric RHS; no use of the cubic source formula."""
    h, j = from_modes(linear, quadratic)
    I = identity()
    perturbation = add(h, scale(s.Rational(1, 2), j))
    squared = product(perturbation, perturbation)
    cubed = product(squared, perturbation)
    q = add(I, perturbation)
    inverse = add(I, scale(-1, perturbation), squared, scale(-1, cubed))
    zeta = AmplitudeJet({1 << a: 1 for a in range(3)})
    chi = trace(add(perturbation, scale(-s.Rational(1, 2), squared),
                    scale(s.Rational(1, 3), cubed))) / 6 - zeta
    cp = chi.dt(time)
    ep = 1 + 2 * chi + 2 * chi * chi + s.Rational(4, 3) * chi * chi * chi
    em = 1 - 2 * chi + 2 * chi * chi - s.Rational(4, 3) * chi * chi * chi
    qp = derivative(q, time)
    gamma = connection_product(inverse, connection_linear(q, momenta))
    ricci = add(ricci_derivative(gamma, momenta), ricci_product(gamma, gamma))
    spatial = [[
        chi.dx(l, momenta).dx(i, momenta)
        - sum(gamma[a][i][l] * chi.dx(a, momenta) for a in range(3))
        + chi.dx(i, momenta) * chi.dx(l, momenta) - ricci[i][l]
        for l in range(3)] for i in range(3)]
    rhs = add(product(product(qp, inverse), qp), scale(-2 * (1 / time + cp), qp),
              scale(2 * ((em - 1) / time**2 - 2 * cp / time), q), scale(2 * ep, spatial))
    return coefficient(rhs, 7)


class CubicTests(unittest.TestCase):
    def zero(self, expr):
        self.assertEqual(s.simplify(s.expand(expr)), 0)

    def matrix_zero(self, mat):
        for item in mat:
            self.zero(item)

    def test_cubic_determinant_and_fluid_composition(self):
        # Noncommuting symbolic symmetric matrix coefficients, no spatial assumptions.
        eps = s.Symbol('epsilon')
        def matrix(prefix):
            return s.Matrix(3, 3, lambda i, j: s.Symbol(f'{prefix}{min(i,j)}{max(i,j)}'))
        h, j, l = matrix('h'), matrix('j'), matrix('l')
        q = s.eye(3) + eps * h + eps**2 * j / 2 + eps**3 * l / 6
        det = s.Poly(q.det(), eps)
        d1, d2, d3 = det.nth(1), det.nth(2), det.nth(3)
        # Cubic coefficient of log(1+d1 eps+d2 eps²+d3 eps³).
        log3 = 6 * (d3 - d1 * d2 + d1**3 / 3)
        self.zero(log3 - s.trace(l) + 3 * s.trace(h * j) - 2 * s.trace(h**3))
        c, d, e = s.symbols('c d e')
        chi = eps * c + eps**2 * d / 2 + eps**3 * e / 6
        rho = s.series(s.exp(-4 * chi), eps, 0, 4).removeO().expand()
        lapse = s.series(s.exp(chi), eps, 0, 4).removeO().expand()
        self.zero(6 * rho.coeff(eps, 3) + 4 * e - 48 * c * d + 64 * c**3)
        self.zero(6 * lapse.coeff(eps, 3) - e - 3 * c * d - c**3)

    def test_cubic_source_against_exact_metric_rhs(self):
        eta = s.Symbol('eta', positive=True)
        momenta = [s.Matrix([s.Symbol(f'k{a}{i}') for i in range(3)]) for a in range(3)]
        a = s.Matrix([[1, 2, -1], [2, -1, 3], [-1, 3, 2]])
        b = s.Matrix([[2, -1, 1], [-1, 3, 2], [1, 2, -2]])
        # Noncommuting time-dependent probes; all nine momentum components free.
        linear = [2 * I3 + eta * a, 2 * I3 + eta**2 * b,
                  2 * I3 + eta * b + eta**2 * a]
        quadratic = {(0, 1): 4 * I3 + eta * b,
                     (0, 2): 4 * I3 + eta**2 * a,
                     (1, 2): 4 * I3 + eta * (a + b)}
        self.matrix_zero(cubic_source(eta, momenta, linear, quadratic)
                         - exact_rhs_polarization(eta, momenta, linear, quadratic))

    def test_two_constant_dilations_exact_finite_time(self):
        eta, p = s.symbols('eta p', positive=True)
        zero, momentum = s.zeros(3, 1), s.Matrix([0, 0, p])
        mode, lapse = linear_mode(eta, momentum)
        pair = 2 * mode - p * mode.diff(p)
        triple = 4 * mode - 3 * p * mode.diff(p) + p**2 * mode.diff(p, 2)
        momenta = [zero, zero, momentum]
        linear = [2 * I3, 2 * I3, mode]
        quadratic = {(0, 1): 4 * I3, (0, 2): pair, (1, 2): pair}
        source = cubic_source(eta, momenta, linear, quadratic)
        self.matrix_zero(time_operator(eta, triple) + spatial_operator(momentum, triple) - source)
        h, j = from_modes(linear, quadratic)
        w = trace(product(h, j)) / 2 - trace(product(product(h, h), h)) / 3
        chi3 = s.trace(triple) / 6 - w.coefficient(7) / 6
        self.zero(chi3 - p * s.diff(lapse, p) - p**2 * s.diff(lapse, p, 2))
        self.matrix_zero(triple.applyfunc(lambda x: s.limit(x, eta, 0, dir='+')) - 8 * I3)

    def test_cubic_regular_boundary_general_momenta(self):
        eta = s.Symbol('eta', positive=True)
        momenta = [s.Matrix([s.Symbol(f'k{a}{i}') for i in range(3)]) for a in range(3)]
        linear = [linear_boundary_jet(eta, k)[0] for k in momenta]
        quadratic = {(a, b): 4 * I3 + eta**2 * quadratic_boundary_coefficient(momenta[a], momenta[b])
                     for a in range(3) for b in range(a + 1, 3)}
        source = cubic_regular_source(eta, momenta, linear, quadratic)
        for value in source:
            expanded = s.expand(value)
            self.zero(expanded.coeff(eta, -2))
            self.zero(expanded.coeff(eta, 0))
            # The finite boundary jet must contain no other negative powers.
            self.assertTrue(s.Poly(expanded, eta).is_univariate)

    def test_cubic_lapse_boundary_and_internal_zero_pair(self):
        eta = s.Symbol('eta', positive=True)
        momenta = [s.Matrix([s.Symbol(f'k{a}{i}') for i in range(3)]) for a in range(3)]
        linear = [linear_boundary_jet(eta, k)[0] for k in momenta]
        quadratic = {(a, b): 4 * I3 + eta**2 * quadratic_boundary_coefficient(momenta[a], momenta[b])
                     for a in range(3) for b in range(a + 1, 3)}
        h, j = from_modes(linear, quadratic)
        w = trace(product(h, j)) / 2 - trace(product(product(h, h), h)) / 3
        # l has seed kernel 8 I and no eta² term; its tr(l)/6 at mask 7 is 24.
        chi3 = 4 - w.coefficient(7) / 6
        leading = s.expand(chi3).coeff(eta, 2)
        expected = (-s.Rational(4, 9) * sum(k.dot(k) for k in momenta)
                    + s.Rational(2, 9) * sum(momenta[a].dot(momenta[b])
                                             for a in range(3) for b in range(a + 1, 3)))
        self.zero(leading - expected)
        k, q, r = momenta
        restricted = leading.subs(dict(zip(r, -q)), simultaneous=True)
        self.zero(restricted + s.Rational(4, 9) * k.dot(k) + s.Rational(10, 9) * q.dot(q))
        self.zero(s.expand(chi3).coeff(eta, 0))


if __name__ == '__main__':
    unittest.main(verbosity=2)
