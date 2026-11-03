"""Exact constant-mode source identities and an augmented-adjoint Ward fixture.

The field residual test is exact Einstein-radiation algebra at arbitrary time.
The finite matrix experiment is explicitly NOT a radiation solution. The full
physical Ward identity follows from spatial dilation covariance, not that fixture.
"""

import sympy as s
from entire_sphere import RadialBessel
from stationary_cubic import direct_field_residuals


def zero(expr):
    if isinstance(expr, s.MatrixBase):
        for value in expr:
            zero(value)
    else:
        assert s.cancel(s.expand(expr)) == 0, expr


def check_field_sources():
    import engine as eng

    p = s.Matrix(s.symbols("p0:3", real=True))
    saved = eng.Ctx.maxdeg, eng.Ctx.momenta, eng.Ctx.trunc
    try:
        eng.Ctx.maxdeg = 3
        eng.Ctx.momenta = {1: [p[i] for i in range(3)], 2: [0, 0, 0], 3: [0, 0, 0]}
        eng.Ctx.trunc = s.expand
        a = s.Function("a")(eng.t)
        values = [s.Function("h" + str(i))(eng.t) for i in range(5)]
        H = s.Matrix(
            [[values[0], values[2], values[3]], [values[2], values[1], values[4]], [values[3], values[4], -values[0] - values[1]]]
        )
        one = lambda value: eng.Field({frozenset([1]): value})
        zeta = eng.Field({frozenset([i]): s.S.One for i in range(1, 4)})
        N, T, c = direct_field_residuals(one(a), zeta, [[one(H[i, j]) for j in range(3)] for i in range(3)])

        def D(expr):
            return sum(p[i] * s.diff(expr, p[i]) for i in range(3))

        for field in (N, c, *(T[i][j] for i in range(3) for j in range(3))):
            linear = field.coeff([1])
            pair = field.coeff([1, 2])
            third = field.coeff([1, 2, 3])
            zero(pair + D(linear))
            zero(third + 2 * pair)
            zero(field.coeff([2, 3]))
        zero(c.coeff([1, 2]) - 2 * p.dot(p) / 9)
        zero(c.coeff([1, 2, 3]) + 4 * p.dot(p) / 9)
        assert N.coeff([1, 2]) != 0
    finally:
        eng.Ctx.maxdeg, eng.Ctx.momenta, eng.Ctx.trunc = saved
    print(
        "[ok] exact C_(p00)=-2 C_(p0), c_(p00)=-2 c_(p0), and zero constant-constant sources at symbolic momentum and arbitrary field profiles"
    )


def check_augmented_adjoint():
    k = s.Symbol("k", positive=True)
    D = lambda M: k * M.diff(k)
    A0 = s.Matrix([[1, 1, 0], [0, 2, 1], [0, 0, 3]])
    A2 = s.Matrix([[1, 0, 1], [0, 2, 0], [0, 0, 1]])
    assert A0 * A2 != A2 * A0
    A = (A0 + k * k * A2).row_join(s.Matrix([1, -1, 2])).col_join(s.Matrix([[0, 0, 0, 1]]))
    B = k * k * s.Matrix([1, 2, -1, s.Rational(1, 9)])
    j = s.Matrix([RadialBessel(0, k * k), k * k * RadialBessel(1, k * k), RadialBessel(2, k * k), 0])
    y = -A.inv() * B
    lam = -A.T.inv() * j
    # Independent source-equation solve, with the fourth initial-data row kept.
    first = A.T.inv() * (D(A).T * lam + D(j))
    zero(first + D(lam))
    pair = -D(A) * y - D(B)
    third = D(D(A)) * y + D(D(B))
    zero(third + 2 * pair)
    direct = (D(D(j)).T * y)[0]
    assembled = direct + (lam.T * third)[0] + 2 * (first.T * pair)[0]
    measured = (j.T * y)[0]
    zero(assembled - D(D(measured)))
    initial = lam[3] * third[3] + 2 * first[3] * pair[3]
    assert s.simplify(initial.subs(k, 1)) != 0
    print(
        "[ok] lambda_constant=-D lambda0 and complete stationary cubic=D^2 linear measurement, including the nonzero initial row (finite algebraic fixture)"
    )


if __name__ == "__main__":
    check_field_sources()
    check_augmented_adjoint()
    print("DILATION SOURCE/ADJOINT IDENTITIES VERIFIED; generic finite-wavelength cubic spectrum remains open")
