#!/usr/bin/env python3
"""Exact checks of selected identities in APPROACH.md, NOT a spectrum solver.

No random sampling, numerical evolution, quadrature, or floating-point evidence.
Run: python -m unittest discover -s validation -p 'check_*.py' -v
"""

import unittest
from math import factorial

import sympy as s


def pairings(items):
    """Enumerate labelled perfect matchings independently of the kernel formulas."""
    if not items:
        yield ()
        return
    first, *rest = items
    for index, partner in enumerate(rest):
        remaining = rest[:index] + rest[index + 1 :]
        for tail in pairings(remaining):
            yield ((first, partner),) + tail


class ExactAlgebraTests(unittest.TestCase):
    def zero(self, expr):
        self.assertEqual(s.simplify(s.expand(expr)), 0)

    def test_radiation_background(self):
        eta, A, G, alpha = s.symbols('eta A G alpha', positive=True)
        a = A * eta
        H = s.diff(a, eta) / a
        rho = 3 / (8 * s.pi * G * A**2 * eta**4)
        temperature = (rho / alpha) ** s.Rational(1, 4)
        self.zero(3 * H**2 - 8 * s.pi * G * a**2 * rho)
        self.zero(s.diff(rho, eta) + 4 * H * rho)
        self.zero(s.diff(temperature, eta) + H * temperature)
        self.zero(s.diff(A * eta**2 / 2, eta) - a)

    def test_sky_division_and_zero_mean(self):
        e = s.Symbol('epsilon')
        t1, t2, t3, b1, b2, b3 = s.symbols('t1 t2 t3 b1 b2 b3')
        raw = 1 + e * t1 + e**2 * t2 / 2 + e**3 * t3 / 6
        mean = 1 + e * b1 + e**2 * b2 / 2 + e**3 * b3 / 6
        derived = s.series(raw / mean - 1, e, 0, 4).removeO().expand()
        expected = (
            t1 - b1,
            t2 - b2 - 2 * b1 * (t1 - b1),
            t3 - b3 - 3 * b1 * (t2 - b2)
            + (6 * b1**2 - 3 * b2) * (t1 - b1),
        )
        for order, value in enumerate(expected, 1):
            self.zero(factorial(order) * derived.coeff(e, order) - value)
            # Each value is linear in the raw t_r; b_r are sky constants.
            self.zero(value.subs({t1: b1, t2: b2, t3: b3}))

    def test_log_to_temperature(self):
        e = s.Symbol('epsilon')
        q1, q2, q3 = s.symbols('q1 q2 q3')
        log_map = e * q1 + e**2 * q2 / 2 + e**3 * q3 / 6
        derived = s.series(s.exp(log_map), e, 0, 4).removeO().expand()
        for order, value in enumerate(
            (q1, q2 + q1**2, q3 + 3 * q1 * q2 + q1**3), 1
        ):
            self.zero(factorial(order) * derived.coeff(e, order) - value)

    def test_common_temperature_factor_cancels(self):
        e = s.Symbol('epsilon')
        t1, t2, t3, b1, b2, b3, h1, h2, h3 = s.symbols(
            't1 t2 t3 b1 b2 b3 h1 h2 h3'
        )
        raw = 1 + e * t1 + e**2 * t2 / 2 + e**3 * t3 / 6
        mean = 1 + e * b1 + e**2 * b2 / 2 + e**3 * b3 / 6
        common = 1 + e * h1 + e**2 * h2 / 2 + e**3 * h3 / 6
        multiplied_raw = s.Poly(s.expand(common * raw), e)
        multiplied_mean = s.Poly(s.expand(common * mean), e)
        trunc_raw = sum(multiplied_raw.nth(i) * e**i for i in range(4))
        trunc_mean = sum(multiplied_mean.nth(i) * e**i for i in range(4))
        # Cross-multiply instead of doing a more expensive rational expansion.
        residual = s.expand(trunc_raw * mean - trunc_mean * raw)
        for order in range(4):
            self.zero(residual.coeff(e, order))

    def test_cubic_emission_root(self):
        e, delta = s.symbols('epsilon delta')
        # F[r,j] = partial_epsilon^r partial_lambda^j F at (0, lambda_0).
        F = {(r, j): s.Symbol(f'F{r}{j}')
             for r in range(4) for j in range(4 - r)}
        F[0, 0] = s.S.Zero
        d = F[0, 1]  # Assumed nonzero: transverse intersection.
        l1 = -F[1, 0] / d
        l2 = -(F[2, 0] + 2 * l1 * F[1, 1] + l1**2 * F[0, 2]) / d
        l3 = -(
            F[3, 0] + 3 * l1 * F[2, 1] + 3 * l2 * F[1, 1]
            + 3 * l1**2 * F[1, 2] + 3 * l1 * l2 * F[0, 2]
            + l1**3 * F[0, 3]
        ) / d
        taylor = sum(value * e**r * delta**j / (factorial(r) * factorial(j))
                     for (r, j), value in F.items())
        # Substitute independent jet coefficients before inserting root formulas.
        a1, a2, a3 = s.symbols('lambda1 lambda2 lambda3')
        composed = s.Poly(s.expand(taylor.subs(
            delta, e * a1 + e**2 * a2 / 2 + e**3 * a3 / 6
        )), e)
        for order in range(4):
            self.zero(composed.nth(order).subs({a1: l1, a2: l2, a3: l3}))

    def test_taylor_power_factors(self):
        e = s.Symbol('epsilon', real=True)
        a = s.symbols('a1:6')
        b = s.symbols('b1:6')  # Independent symbols for conjugates.
        left = sum(e**r * a[r - 1] / factorial(r) for r in range(1, 6))
        right = sum(e**r * b[r - 1] / factorial(r) for r in range(1, 6))
        power = s.expand(left * right)
        self.zero(power.coeff(e, 2) - a[0] * b[0])
        self.zero(power.coeff(e, 4)
                  - (a[1] * b[1] / 4 + (a[0] * b[2] + a[2] * b[0]) / 6))
        # Any degree 3 or 5 contribution has an odd number of Gaussian seeds.
        for total in (3, 5):
            for r in range(1, total):
                self.assertEqual(len(list(pairings(list(range(total))))), 0)
                self.assertEqual((r + (total - r)) % 2, 1)

    def test_wick_pairings(self):
        pairs = list(pairings(list(range(4))))
        self.assertEqual(len(pairs), 3)
        # For 22, slots 0,1 belong to one factor; slots 2,3 to the other.
        crosses = [sum((a < 2) != (b < 2) for a, b in p) for p in pairs]
        self.assertEqual(crosses.count(0), 1)  # Product of means: subtract.
        self.assertEqual(crosses.count(2), 2)  # Connected quadratic covariance.
        # For 13, slot 0 belongs to the linear factor: all 3 pairings remain.
        for p in pairs:
            self.assertEqual(sum(0 in pair for pair in p), 1)
            self.assertEqual(sum(0 not in pair for pair in p), 1)
        self.assertEqual(len(list(pairings([0, 1]))), 1)

    def test_scalar_fundamental_solutions(self):
        y = s.Symbol('y', positive=True)
        f = 3 * (s.sin(y) - y * s.cos(y)) / y**3
        g = (s.cos(y) + y * s.sin(y)) / y**3
        for solution in (f, g):
            self.zero(s.diff(solution, y, 2) + 4 * s.diff(solution, y) / y
                      + solution)
        self.assertEqual(s.limit(f, y, 0, dir='+'), 1)
        self.zero(f * s.diff(g, y) - s.diff(f, y) * g + 3 / y**4)

    def test_scalar_retarded_kernel(self):
        eta, time, w = s.symbols('eta time w', positive=True)

        def f(x):
            return 3 * (s.sin(w * x) - w * x * s.cos(w * x)) / (w * x)**3

        def g(x):
            return (s.cos(w * x) + w * x * s.sin(w * x)) / (w * x)**3

        W = -3 / (w**3 * time**4)
        green = (f(time) * g(eta) - g(time) * f(eta)) / W
        self.zero(green.subs(eta, time))
        self.zero(s.diff(green, eta).subs(eta, time) - 1)
        self.zero(s.diff(green, eta, 2) + 4 * s.diff(green, eta) / eta
                  + w**2 * green)
        # The Heaviside extension is retarded. These two boundary conditions
        # establish a unit delta source, with no delta-prime term.
        zero_mode = time / 3 * (1 - (time / eta)**3)
        self.zero(s.limit(green, w, 0, dir='+') - zero_mode)
        self.zero(s.diff(zero_mode, eta, 2) + 4 * s.diff(zero_mode, eta) / eta)

    def test_tensor_fundamental_solutions(self):
        y, eta = s.symbols('y eta', positive=True)
        for solution in (s.sin(y) / y, s.cos(y) / y):
            self.zero(s.diff(solution, y, 2) + 2 * s.diff(solution, y) / y
                      + solution)
        for solution in (s.S.One, 1 / eta):
            self.zero(s.diff(solution, eta, 2) + 2 * s.diff(solution, eta) / eta)

    def test_shell_measure_and_linear_angular_factor(self):
        k0, sigma2 = s.symbols('k0 sigma2', positive=True)
        # delta(log(k/k0)) = k0 delta(k-k0), for k>0 and k0>0.
        k = s.Symbol('k', positive=True)
        self.zero(1 / s.diff(s.log(k / k0), k).subs(k, k0) - k0)
        shell_weight = k0**2 / (2 * s.pi)**3 * (2 * s.pi**2 * sigma2 / k0**2)
        self.zero(shell_weight - sigma2 / (4 * s.pi))
        self.zero(4 * s.pi * shell_weight - sigma2)
        ell = s.Symbol('ell', integer=True, nonnegative=True)
        addition_theorem = (2 * ell + 1) / (4 * s.pi)
        radial_factor = (4 * s.pi)**2 * addition_theorem / (2 * ell + 1)
        radial_factor *= 4 * s.pi / (2 * s.pi)**3
        self.zero(radial_factor - 2 / s.pi)
        self.zero(radial_factor * 2 * s.pi**2 - 4 * s.pi)


if __name__ == '__main__':
    unittest.main(verbosity=2)
