"""Independent finite-time Newtonian/temperature-frame tensor cross-check.

Literature: Kohri & Terada, arXiv:1804.08577v2, Eqs. (1), (7)-(9),
(15), (18)-(20), (22). See report/tensor_newtonian_crosscheck.md.

No project modules, cached sources, numerical evolution, or quadrature are used.
The temperature-frame TT source is rebuilt from elementary Christoffel/Ricci
formulas. The comparison is an exact finite-time differential identity plus
regular initial data, not a truncated-series comparison. Series are used only
for boundary conditions and negative controls.

Run: .venv/bin/python -B solution/check_tensor_newtonian.py
"""

import sympy as sp


class Cin(sp.Function):
    """Entire Cin(z) = integral_0^z (1-cos u)/u du."""

    def fdiff(self, argindex=1):
        z = self.args[0]
        return (1 - sp.cos(z)) / z


class SinIntegral(sp.Function):
    """Entire Si; avoid branch-dependent automatic rewrites in symbolic checks."""

    def fdiff(self, argindex=1):
        z = self.args[0]
        return sp.sin(z) / z


def zero(expr):
    result = sp.simplify(sp.trigsimp(sp.expand(expr)))
    assert result == 0, result


def polarization(matrix):
    """Coefficient of vv-ww with K along z, v along x, w along y."""
    return sp.expand((matrix[0, 0] - matrix[1, 1]) / 2)


def temperature_source():
    """Rebuild the z1*z2 source with the unknown second-order H set to zero."""
    b, c, h, nu, nudot, hdot = sp.symbols("b c h nu nudot hdot", real=True)
    p, q = sp.Matrix([b, 0, c]), sp.Matrix([-b, 0, c])
    total = p + q
    identity = sp.eye(3)
    k2 = b**2 + c**2
    P, Q = p * p.T - identity * k2 / 3, q * q.T - identity * k2 / 3
    gp, gq = 2 * identity + h * P, 2 * identity + h * Q
    gpq = 4 * identity + 2 * h * (P + Q) + h**2 * (P * Q + Q * P) / 2

    def connection(metric, momentum):
        return [
            [
                [
                    sp.I * (momentum[i] * metric[a, j] + momentum[j] * metric[a, i] - momentum[a] * metric[i, j]) / 2
                    for j in range(3)
                ]
                for i in range(3)
            ]
            for a in range(3)
        ]

    cp, cq = connection(gp, p), connection(gq, q)
    linear_pair = connection(gpq, total)
    cpair = [
        [
            [
                sp.expand(linear_pair[a][i][j] - sum(gp[a, l] * cq[l][i][j] + gq[a, l] * cp[l][i][j] for l in range(3)))
                for j in range(3)
            ]
            for i in range(3)
        ]
        for a in range(3)
    ]

    def linear_ricci(christoffel, momentum):
        return sp.Matrix(
            3,
            3,
            lambda i, j: sum(
                sp.I * momentum[a] * christoffel[a][i][j] - sp.I * momentum[j] * christoffel[a][i][a] for a in range(3)
            ),
        )

    rp, rq = linear_ricci(cp, p), linear_ricci(cq, q)
    rpair = linear_ricci(cpair, total) + sp.Matrix(
        3,
        3,
        lambda i, j: sum(
            cp[a][a][l] * cq[l][i][j] + cq[a][a][l] * cp[l][i][j] - cp[a][j][l] * cq[l][i][a] - cq[a][j][l] * cp[l][i][a]
            for a in range(3)
            for l in range(3)
        ),
    )
    mixed_ricci = rpair - gp * rq - gq * rp
    mixed_hessian = sp.Matrix(
        3, 3, lambda i, j: -sp.I * nu * sum(cp[a][i][j] * q[a] + cq[a][i][j] * p[a] for a in range(3))
    ) + nu * (gp * q * q.T + gq * p * p.T)
    # On the shell H_1=h(t)*fixed_matrix, so the quadratic commutator shear is zero.
    source = nudot * hdot * (P + Q) + mixed_ricci - 2 * mixed_hessian - 2 * nu**2 * (p * q.T + q * p.T)
    forcing = sp.factor(-2 * polarization(source))

    k, mu = sp.symbols("k mu", real=True)
    shell = {b**2: k**2 * (1 - mu) / 2, c**2: k**2 * (1 + mu) / 2}
    forcing = sp.expand(forcing).subs(shell, simultaneous=True)
    expected = (
        k**2
        * (1 - mu)
        / 36
        * (
            k**4 * h**2 * (9 * mu**2 + 12 * mu + 1)
            - 6 * k**2 * h * (6 * mu * nu + 3 * mu + 2)
            - 36 * hdot * nudot
            - 72 * nu**2
            + 36
        )
    )
    zero(forcing - expected)
    print("[ok] independent Christoffel/Ricci reconstruction of the full shell TT source")

    # Directly check the quadratic pullback, conformal normalization, and log.
    B, A = sp.symbols("B A", real=True)
    C = 1 + k2 * B / 3  # A/t - Phi - nu = g - Delta B/3, per unit seed.
    product = p * q.T + q * p.T
    pullback = (p.dot(q) * B**2 + A**2) * product - 4 * C * B * (p * p.T + q * q.T)
    log_subtraction = gpq.subs(h, -2 * B)
    gauge = polarization(pullback - log_subtraction)
    zero(gauge - b**2 * ((c**2 - b**2) * B**2 - A**2))
    print("[ok] quadratic coordinate pullback + conformal factor + matrix log give the elementary gauge correction")
    return expected, (k, mu, h, nu, hdot, nudot)


def finite_time_comparison(forcing, symbols):
    k, mu, h, nu_symbol, hdot, nudot = symbols
    x = sp.Symbol("x", positive=True)
    s = 2 * sp.sin(x) / x - sp.cos(x)
    B = 3 * (s - 1) / k**2
    dt = lambda expr: sp.diff(expr, x) * k / sp.sqrt(3)
    A = dt(B)
    phi = -2 * (sp.sin(x) / x - sp.cos(x)) / x**2
    nu = sp.cos(x) - sp.sin(x) / x
    time = sp.sqrt(3) * x / k
    zero(nu - phi - dt(A) - A / time)
    zero(A / time - phi - nu - k**2 * B / 3 - 1)
    zero(A + time * (phi + time * dt(phi)) / 2)
    print("[ok] first-order lapse, determinant, zero-shift and comoving gauge relations")

    f = 2 * phi**2 + (phi + x * sp.diff(phi, x)) ** 2
    # KT Eq. (20), specialized to equal scalar momenta; acoustic variable x.
    f_kt = (
        4
        / x**6
        * (6 * x**2 * sp.cos(x) ** 2 + (6 - 4 * x**2 + x**4) * sp.sin(x) ** 2 + 4 * x * (x**2 - 3) * sp.cos(x) * sp.sin(x))
    )
    zero(f - f_kt)
    print("[ok] KT radiation source, including the primordial 4/9 factor and velocity term")

    total_k2 = 2 * k**2 * (1 + mu)
    wave = lambda expr: dt(dt(expr)) + 2 * dt(expr) / time + total_k2 * expr
    gauge_space = k**4 * mu * (1 - mu) * B**2 / 2
    gauge_time = -(k**2) * (1 - mu) * A**2 / 2
    predicted = k**2 * (1 - mu) * f + wave(gauge_space + gauge_time)
    actual = forcing.subs({h: -2 * B, nu_symbol: nu, hdot: -2 * A, nudot: dt(nu)}, simultaneous=True)
    zero((actual - predicted) * x**6 / k**2)
    print("[ok] exact arbitrary-time, arbitrary-angle identity: L F_T = source_T (not a time series)")

    # Boundary data: I/K^2 = 2 t^2/9 + O(t^4), A=t/3+O(t^3), B=t^2/6+O(t^4).
    leading = k**2 * (1 - mu) * 2 * time**2 / 9 + sp.series(gauge_space + gauge_time, x, 0, 3).removeO()
    zero(leading - (1 - mu) * x**2 / 2)
    print("[ok] F_T=(1-mu) k^2 t^2/6+O(t^4), no constant or decaying tensor seed")

    # Negative controls: omitting either gauge term changes the equation.
    for name, term in [("time-shift gradient", gauge_time), ("spatial/log correction", gauge_space)]:
        residual = sp.simplify(sp.series(wave(term), x, 0, 3).removeO())
        assert residual != 0
        print(f"[ok] negative control: dropping {name} fails")
    assert sp.limit(f, x, 0) == sp.Rational(4, 3)
    print("[ok] negative control: confusing h_KT with h_KT/2 changes the nonzero tensor source")


def closed_kt_kernel():
    """Check KT Eq. (22) on u=v, with all removable logarithms eliminated."""
    X, r = sp.symbols("X r", positive=True)
    x = r * X / sp.sqrt(3)
    d = 2 * r**2 - 3
    a, b = 1 - 2 * r / sp.sqrt(3), 1 + 2 * r / sp.sqrt(3)
    C = Cin(a * X) + Cin(b * X) - 2 * Cin(X)
    S = SinIntegral(a * X) + SinIntegral(b * X) - 2 * SinIntegral(X)
    E = (
        -4 * r**2 * d * sp.sin(X)
        + 24 * r**2 * sp.cos(x) ** 2 / X
        - 48 * sp.sqrt(3) * r * sp.cos(x) * sp.sin(x) / X**2
        + 12 * (6 + d * X**2) * sp.sin(x) ** 2 / X**3
    )
    integral = 3 / (4 * r**6 * X) * (E + d**2 * (sp.sin(X) * C + sp.cos(X) * S))
    phi = -2 * (sp.sin(x) / x - sp.cos(x)) / x**2
    f = 2 * phi**2 + (phi + X * sp.diff(phi, X)) ** 2
    residual = sp.expand(sp.diff(integral, X, 2) + 2 * sp.diff(integral, X) / X + integral - f)
    assert not residual.has(Cin, SinIntegral)
    residual = residual.replace(lambda e: e.func in (sp.sin, sp.cos), lambda e: e.func(sp.expand(e.args[0])))
    zero(sp.expand_trig(residual))
    series_integral = integral.replace(
        lambda e: e.func == Cin,
        lambda e: sum((-1) ** (j + 1) * e.args[0] ** (2 * j) / (2 * j * sp.factorial(2 * j)) for j in range(1, 6)),
    ).replace(
        lambda e: e.func == SinIntegral,
        lambda e: sum((-1) ** j * e.args[0] ** (2 * j + 1) / ((2 * j + 1) * sp.factorial(2 * j + 1)) for j in range(6)),
    )
    zero(sp.series(series_integral, X, 0, 4).removeO() - 2 * X**2 / 9)
    print("[ok] entire-function KT shell closed form: exact wave equation and regular initial data")
    resonant = integral.subs(r, sp.sqrt(3) / 2).subs({Cin(0): 0, SinIntegral(0): 0})
    resonant_residual = sp.diff(resonant, X, 2) + 2 * sp.diff(resonant, X) / X + resonant - f.subs(r, sp.sqrt(3) / 2)
    zero(sp.expand_trig(resonant_residual))
    print("[ok] exact resonant wave equation at r=sqrt(3)/2 with Cin(0)=Si(0)=0")


if __name__ == "__main__":
    source, variables = temperature_source()
    finite_time_comparison(source, variables)
    closed_kt_kernel()
    print("INDEPENDENT NEWTONIAN/TEMPERATURE TENSOR MATCH PASSED (exact symbolic algebra)")
