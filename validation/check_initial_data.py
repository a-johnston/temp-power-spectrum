"""Exact checks for the quadratic initial-data construction.

No quadrature or numerical evolution. The direct source test uses a separate
nilpotent-amplitude expansion of the exact metric equation, not the source
formula under test. See README for its probe coverage and remaining gaps.
"""
import unittest

import sympy as s

from initial_data import (
    I3, gamma_linear, green_zero, indicial_inverse, linear_boundary_jet,
    linear_mode, quadratic_boundary_coefficient, quadratic_source,
    regular_source, ricci_linear, spatial_operator, time_operator,
)


class Jet:
    """Two commuting nilpotent amplitudes z,w: coefficient order 1,z,w,zw."""
    def __init__(self, coefficients):
        self.c = tuple(map(s.sympify, coefficients))

    @staticmethod
    def lift(value):
        return value if isinstance(value, Jet) else Jet((value, 0, 0, 0))

    def __add__(self, other):
        other = self.lift(other)
        return Jet(tuple(a + b for a, b in zip(self.c, other.c)))

    __radd__ = __add__

    def __neg__(self):
        return Jet(tuple(-a for a in self.c))

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) - self

    def __mul__(self, other):
        a, b, c, d = self.c
        e, f, g, h = self.lift(other).c
        return Jet((a * e, a * f + b * e, a * g + c * e,
                    a * h + b * g + c * f + d * e))

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self * (s.S.One / other)

    def diff_time(self, eta):
        return Jet(tuple(s.diff(c, eta) for c in self.c))

    def diff_space(self, i, k, p):
        return Jet((0, s.I * k[i] * self.c[1], s.I * p[i] * self.c[2],
                    s.I * (k[i] + p[i]) * self.c[3]))


def direct_metric_rhs_cross(eta, k, p, h, f):
    """Independently expand exact equation (3) with q=I+z h+w f.

Uses nilpotent ring arithmetic and the full spatial connection, not R_Q or S.
Seed zeta=z+w has zero quadratic coefficient.
"""
    q = [[Jet((I3[i, j], h[i, j], f[i, j], 0)) for j in range(3)]
         for i in range(3)]
    inverse_cross = h * f + f * h
    qi = [[Jet((I3[i, j], -h[i, j], -f[i, j], inverse_cross[i, j]))
           for j in range(3)] for i in range(3)]
    chi = Jet((0, s.trace(h) / 6 - 1, s.trace(f) / 6 - 1, -s.trace(h * f) / 6))
    cp = chi.diff_time(eta)
    ep = 1 + 2 * chi + 2 * chi * chi
    em = 1 - 2 * chi + 2 * chi * chi
    dx = lambda value, i: value.diff_space(i, k, p)
    qp = [[q[i][j].diff_time(eta) for j in range(3)] for i in range(3)]
    gamma = [[[sum(qi[a][b] * (dx(q[b][j], i) + dx(q[b][i], j) - dx(q[i][j], b))
                      for b in range(3)) * s.Rational(1, 2)
               for j in range(3)] for i in range(3)] for a in range(3)]
    out = s.zeros(3)
    for i in range(3):
        for j in range(3):
            ricci = sum(dx(gamma[a][i][j], a) - dx(gamma[a][i][a], j)
                        for a in range(3))
            ricci += sum(gamma[a][a][b] * gamma[b][i][j]
                         - gamma[a][j][b] * gamma[b][i][a]
                         for a in range(3) for b in range(3))
            cov_hessian = dx(dx(chi, j), i) - sum(
                gamma[a][i][j] * dx(chi, a) for a in range(3))
            kinetic = sum(qp[i][a] * qi[a][b] * qp[b][j]
                          for a in range(3) for b in range(3))
            rhs = (kinetic - 2 * (1 / eta + cp) * qp[i][j]
                   + 2 * ((em - 1) * (1 / eta**2) - 2 * cp * (1 / eta)) * q[i][j]
                   + 2 * ep * (cov_hessian + dx(chi, i) * dx(chi, j) - ricci))
            out[i, j] = rhs.c[3]
    return out


def direct_constraints_cross(eta, k, p, h, f, j):
    """Constraints (4),(5) assembled from full matrix jets, independently of S."""
    q = [[Jet((I3[i, l], h[i, l], f[i, l], j[i, l])) for l in range(3)]
         for i in range(3)]
    inverse_cross = h * f + f * h - j
    qi = [[Jet((I3[i, l], -h[i, l], -f[i, l], inverse_cross[i, l]))
           for l in range(3)] for i in range(3)]
    chi = Jet((0, s.trace(h) / 6 - 1, s.trace(f) / 6 - 1,
               (s.trace(j) - s.trace(h * f)) / 6))
    dx = lambda value, i: value.diff_space(i, k, p)
    gamma = [[[sum(qi[a][b] * (dx(q[b][l], i) + dx(q[b][i], l) - dx(q[i][l], b))
                      for b in range(3)) * s.Rational(1, 2)
               for l in range(3)] for i in range(3)] for a in range(3)]
    ricci = [[sum(dx(gamma[a][i][l], a) - dx(gamma[a][i][a], l)
                   for a in range(3))
              + sum(gamma[a][a][b] * gamma[b][i][l]
                    - gamma[a][l][b] * gamma[b][i][a]
                    for a in range(3) for b in range(3))
              for l in range(3)] for i in range(3)]
    scalar = sum(qi[i][l] * ricci[i][l] for i in range(3) for l in range(3))
    M = [[sum(qi[i][a] * q[a][l].diff_time(eta) for a in range(3)) / 2
          for l in range(3)] for i in range(3)]
    m = sum(M[i][i] for i in range(3))
    ep, em = 1 + 2 * chi + 2 * chi * chi, 1 - 2 * chi + 2 * chi * chi
    hamiltonian = (ep * scalar + 6 / eta**2 + 4 * m * (1 / eta) + m * m
                   - sum(M[i][l] * M[l][i] for i in range(3) for l in range(3))
                   - 6 * em * (1 / eta**2))
    momentum = []
    for i in range(3):
        div = sum(dx(M[l][i], l) for l in range(3)) + sum(
            gamma[l][l][a] * M[a][i] - gamma[a][l][i] * M[l][a]
            for l in range(3) for a in range(3))
        momentum.append(div - dx(m, i) - sum(dx(chi, l) * M[l][i] for l in range(3))
                        + (2 / eta + m) * dx(chi, i))
    return hamiltonian.c[3], [value.c[3] for value in momentum]


class InitialDataTests(unittest.TestCase):
    def zero(self, expr):
        self.assertEqual(s.simplify(s.expand(expr)), 0)

    def matrix_zero(self, matrix):
        for item in matrix:
            self.zero(item)

    def test_exact_fluid_reduction(self):
        eta = s.Symbol('eta', positive=True)
        chi = s.Function('chi')(eta)
        # sqrt(gamma) is proportional to a^3 exp(3 chi), up to the seed.
        rho = eta**-4 * s.exp(-4 * chi)
        volume = eta**3 * s.exp(3 * chi)
        self.zero(s.diff(rho, eta) + s.Rational(4, 3) * rho
                  * s.diff(volume, eta) / volume)
        self.zero(s.diff(s.log(eta * s.exp(chi)) + s.log(rho) / 4, eta))
        # Spatial Euler uses the same algebra with eta held constant.
        x = s.Symbol('x')
        f = s.Function('f')(x)
        self.zero(s.diff(f, x) + s.diff(-4 * f, x) / 4)

    def test_seed_ricci_from_connection(self):
        x = s.symbols('x0:3')
        z = s.Function('zeta')(*x)
        # Direct connection of exp(2 zeta) delta, followed by its definition of R.
        gamma = [[[s.KroneckerDelta(a, j) * s.diff(z, x[i])
                   + s.KroneckerDelta(a, i) * s.diff(z, x[j])
                   - s.KroneckerDelta(i, j) * s.diff(z, x[a])
                   for j in range(3)] for i in range(3)] for a in range(3)]
        lap = sum(s.diff(z, xi, 2) for xi in x)
        grad2 = sum(s.diff(z, xi)**2 for xi in x)
        ricci = s.zeros(3)
        for i in range(3):
            for j in range(3):
                ricci[i, j] = sum(s.diff(gamma[a][i][j], x[a])
                                  - s.diff(gamma[a][i][a], x[j]) for a in range(3))
                ricci[i, j] += sum(gamma[a][a][b] * gamma[b][i][j]
                                  - gamma[a][j][b] * gamma[b][i][a]
                                  for a in range(3) for b in range(3))
                expected = (-s.diff(z, x[i], x[j]) - I3[i, j] * lap
                            + s.diff(z, x[i]) * s.diff(z, x[j]) - I3[i, j] * grad2)
                self.zero(ricci[i, j] - expected)
        self.zero(s.exp(-2 * z) * s.trace(ricci)
                  - s.exp(-2 * z) * (-4 * lap - 2 * grad2))

    def test_leading_constraints_and_evolution(self):
        # General seed Ricci tensor: equations determine C without projection.
        r = s.Matrix([[s.Symbol(f'R{min(i,j)}{max(i,j)}') for j in range(3)]
                      for i in range(3)])
        scalar = s.trace(r)
        c = -r / 3 + I3 * scalar / 18
        t = s.trace(c)
        self.matrix_zero(3 * c + t * I3 + r)
        self.zero(scalar + 6 * t)
        # Contracted Bianchi: div R^j_i = d_i R/2, for each arbitrary component.
        dR = s.Symbol('dR')
        div_c = -dR / 6 + dR / 18
        dt = -dR / 6
        self.zero(div_c - s.Rational(2, 3) * dt)

    def test_linear_full_equations(self):
        eta, k = s.symbols('eta k', positive=True)
        momentum = s.Matrix([0, 0, k])
        h, c = linear_mode(eta, momentum)
        # Rotational covariance lifts this arbitrary k>0 scalar check to any direction.
        residual = (h.diff(eta, 2) + 2 * h.diff(eta) / eta
                    + 4 * (c / eta**2 + s.diff(c, eta) / eta) * I3
                    + 2 * (momentum * momentum.T) * c + 2 * ricci_linear(momentum, h))
        self.matrix_zero(residual)
        self.zero(s.trace(h) / 6 - 1 - c)
        self.zero(s.trace(ricci_linear(momentum, h))
                  + 12 * s.diff(c, eta) / eta + 12 * c / eta**2)
        self.matrix_zero(h.diff(eta) * momentum / 2
                         - momentum * s.trace(h.diff(eta)) / 2 + 2 * momentum * c / eta)

    def test_linear_zero_mode_and_boundary(self):
        eta, k = s.symbols('eta k', positive=True)
        h, c = linear_mode(eta, [0, 0, k])
        h0, c0 = linear_mode(eta, [0, 0, 0])
        self.matrix_zero(h.applyfunc(lambda v: s.limit(v, k, 0, dir='+')) - h0)
        self.zero(s.limit(c, k, 0, dir='+') - c0)
        hj, cj = linear_boundary_jet(eta, [0, 0, k])
        self.matrix_zero(h.applyfunc(lambda v: s.series(v, eta, 0, 4).removeO()) - hj)
        self.zero(s.series(c, eta, 0, 4).removeO() - cj)

    def test_quadratic_source_against_direct_metric_expansion(self):
        eta = s.Symbol('eta', positive=True)
        k, p = s.Matrix(s.symbols('k0:3')), s.Matrix(s.symbols('p0:3'))
        # Noncommuting symmetric probes with independent time profiles; momenta
        # remain completely symbolic. This is a regression, not a tensor theorem.
        a = s.Matrix([[1, 2, -1], [2, 0, 3], [-1, 3, -2]])
        b = s.Matrix([[2, -1, 0], [-1, 3, 2], [0, 2, 1]])
        h = 2 * I3 + eta * a + eta**2 * b
        f = 2 * I3 - eta * b + eta**3 * a
        c, d = s.trace(h) / 6 - 1, s.trace(f) / 6 - 1
        expected = direct_metric_rhs_cross(eta, k, p, h, f)
        actual = quadratic_source(eta, k, p, h, f, c, d)
        self.matrix_zero(expected - actual)

    def test_regular_source_boundary_general_momenta(self):
        eta = s.Symbol('eta', positive=True)
        k, p = s.Matrix(s.symbols('k0:3')), s.Matrix(s.symbols('p0:3'))
        h, c = linear_boundary_jet(eta, k)
        f, d = linear_boundary_jet(eta, p)
        source = regular_source(eta, k, p, h, f, c, d)
        at_zero = source.applyfunc(lambda x: s.expand(x).subs(eta, 0))
        expected = 2 * (k * p.T + p * k.T) - 4 * I3 * k.dot(p)
        self.matrix_zero(at_zero - expected)
        self.matrix_zero(indicial_inverse(2, at_zero) - quadratic_boundary_coefficient(k, p))

    def test_source_symmetry(self):
        eta = s.Symbol('eta', positive=True)
        k, p = s.Matrix(s.symbols('k0:3')), s.Matrix(s.symbols('p0:3'))
        h, c = linear_boundary_jet(eta, k)
        f, d = linear_boundary_jet(eta, p)
        source = quadratic_source(eta, k, p, h, f, c, d)
        self.matrix_zero(source - source.T)
        self.matrix_zero(source - quadratic_source(eta, p, k, f, h, d, c))

    def test_zero_total_momentum_forced_tracefree_response(self):
        k = s.Matrix(s.symbols('k0:3'))
        b = quadratic_boundary_coefficient(k, -k)
        self.zero(s.trace(b) - 2 * k.dot(k) / 3)
        self.matrix_zero(b - I3 * s.trace(b) / 3
                         + s.Rational(2, 3) * (k * k.T - I3 * k.dot(k) / 3))
        x = s.Matrix([[s.Symbol(f'X{min(i,j)}{max(i,j)}') for j in range(3)]
                      for i in range(3)])
        self.matrix_zero(spatial_operator([0, 0, 0], x))

    def test_lapse_density_boundary_kernel(self):
        eta = s.Symbol('eta', positive=True)
        k, p = s.Matrix(s.symbols('k0:3')), s.Matrix(s.symbols('p0:3'))
        h, c = linear_boundary_jet(eta, k)
        f, d = linear_boundary_jet(eta, p)
        j = 4 * I3 + eta**2 * quadratic_boundary_coefficient(k, p)
        chi2 = s.expand(s.trace(j) / 6 - s.trace(h * f) / 6).coeff(eta, 2)
        self.zero(chi2 - (2 * (k.dot(k) + p.dot(p)) - k.dot(p)) / 9)
        antipodal = chi2.subs(dict(zip(p, -k)), simultaneous=True)
        self.zero(antipodal - 5 * k.dot(k) / 9)
        self.zero(-4 * antipodal + 20 * k.dot(k) / 9)

    def test_matrix_green_boundary_and_equation(self):
        eta, time = s.symbols('eta time', positive=True)
        x = s.Matrix([[s.Symbol(f'X{min(i,j)}{max(i,j)}') for j in range(3)]
                      for i in range(3)])
        green = green_zero(eta, time, x)
        self.matrix_zero(green.subs(eta, time))
        self.matrix_zero(green.diff(eta).subs(eta, time) - x)
        self.matrix_zero(time_operator(eta, green))
        self.matrix_zero(green.subs(time, 0))

    def test_indicial_inverse_and_exact_monomial_integral(self):
        eta, time = s.symbols('eta time', positive=True)
        x = s.Matrix([[s.Symbol(f'X{min(i,j)}{max(i,j)}') for j in range(3)]
                      for i in range(3)])
        for m in (0, 2, 4):
            result = eta**(m + 2) * indicial_inverse(m + 2, x)
            self.matrix_zero(time_operator(eta, result) - eta**m * x)
            integral = green_zero(eta, time, time**m * x).applyfunc(
                lambda v: s.integrate(v, (time, 0, eta)))
            self.matrix_zero(integral - result)

    def test_quadratic_constraints_through_next_time_coefficient(self):
        eta = s.Symbol('eta', positive=True)
        a, b, c = s.symbols('a b c', real=True)
        # Any two nonzero vectors can be rotated into this plane; no assumption
        # about their relative angle or magnitudes is made, and no division by a,b,c.
        k, p = s.Matrix([a, 0, 0]), s.Matrix([b, c, 0])

        def fourth_jet(momentum):
            metric, lapse = linear_boundary_jet(eta, momentum)
            k2 = momentum.dot(momentum)
            return (metric + eta**4 * (I3 * k2**2 / 540
                                       + (momentum * momentum.T) * k2 / 60),
                    lapse + eta**4 * k2**2 / 270)

        h, lapse_h = fourth_jet(k)
        f, lapse_f = fourth_jet(p)
        source = regular_source(eta, k, p, h, f, lapse_h, lapse_f)
        f2 = source.applyfunc(lambda value: s.expand(value).coeff(eta, 2))
        j2 = quadratic_boundary_coefficient(k, p)
        j4 = indicial_inverse(4, f2 - spatial_operator(k + p, j2))
        j = 4 * I3 + eta**2 * j2 + eta**4 * j4
        hamiltonian, momentum = direct_constraints_cross(eta, k, p, h, f, j)
        hamiltonian = s.expand(hamiltonian)
        for power in (-2, 0, 2):
            self.zero(hamiltonian.coeff(eta, power))
        for component in momentum:
            component = s.expand(component)
            for power in (-1, 1, 3):
                self.zero(component.coeff(eta, power))

    def test_constraint_propagation_energy_identity(self):
        hr, hi = s.symbols('hr hi', real=True)
        mr, mi = s.symbols('mr0:3', real=True), s.symbols('mi0:3', real=True)
        k = s.symbols('k0:3', real=True)
        eta = s.Symbol('eta', positive=True)
        h = hr + s.I * hi
        m = [mr[i] + s.I * mi[i] for i in range(3)]
        hp = -4 * h / eta - 2 * s.I * sum(k[i] * m[i] for i in range(3))
        mp = [-2 * m[i] / eta - s.I * k[i] * h / 2 for i in range(3)]
        derivative = 2 * s.re(s.conjugate(h) * hp) + 8 * sum(
            s.re(s.conjugate(m[i]) * mp[i]) for i in range(3))
        self.zero(derivative + 8 * (hr**2 + hi**2) / eta
                  + 16 * sum(mr[i]**2 + mi[i]**2 for i in range(3)) / eta)


if __name__ == '__main__':
    unittest.main(verbosity=2)
