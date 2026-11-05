"""Exact plane-symmetric reduction of the temperature-frame field equations.

gamma = e^{2 zeta} exp(H), H = f diag(2/3,-1/3,-1/3), so gamma = diag(e^{2A}, e^{2B}, e^{2B}) with
A = zeta + f/3, B = zeta - f/6, functions of (t, y). S^i_j = (1/2) f' diag(2/3,-1/3,-1/3).
Equation (T) is reported as f'' + 2(1/t + nu') f' + 2 F_f = 0, where TF[...] = F_f diag(2/3,-1/3,-1/3).
Every expression is obtained from the generic Ricci formula, not typed in by hand.
"""

import sympy as sp

t, y = sp.symbols("t y", positive=True)
zeta = sp.Function("zeta")(y)
f = sp.Function("f")(t, y)
nu = sp.Function("nu")(t, y)


def christoffel(g, coords):
    ginv = g.inv()
    n = len(coords)
    return [
        [
            [
                sum(
                    ginv[i, l] * (sp.diff(g[l, j], coords[k]) + sp.diff(g[l, k], coords[j]) - sp.diff(g[j, k], coords[l]))
                    for l in range(n)
                )
                / 2
                for k in range(n)
            ]
            for j in range(n)
        ]
        for i in range(n)
    ]


def ricci(g, coords):
    G = christoffel(g, coords)
    n = len(coords)
    R = sp.zeros(n, n)
    for j in range(n):
        for k in range(n):
            R[j, k] = sum(
                sp.diff(G[i][j][k], coords[i])
                - sp.diff(G[i][j][i], coords[k])
                + sum(G[i][i][l] * G[l][j][k] - G[i][k][l] * G[l][j][i] for l in range(n))
                for i in range(n)
            )
    return R, G


def reduced_equations():
    x2, x3 = sp.symbols("x2 x3")
    coords = (y, x2, x3)
    A = zeta + f / 3
    B = zeta - f / 6
    g = sp.diag(sp.exp(2 * A), sp.exp(2 * B), sp.exp(2 * B))
    Ric, G = ricci(g, coords)
    ginv = g.inv()
    Ric_mixed = (ginv * Ric).applyfunc(sp.simplify)
    dn = [sp.diff(nu, c) for c in coords]
    # D^i D_j nu and D^i nu D_j nu
    hess = sp.Matrix(3, 3, lambda i, j: sp.diff(nu, coords[i], coords[j]) - sum(G[k][i][j] * dn[k] for k in range(3)))
    hess_mixed = ginv * hess
    grad_mixed = sp.Matrix(3, 3, lambda i, j: ginv[i, i] * dn[i] * dn[j])
    M = Ric_mixed - 2 * hess_mixed + 2 * grad_mixed
    tr = M.trace()
    F_tf = M - sp.eye(3) * tr / 3
    # TF part is F_f * diag(2/3,-1/3,-1/3)
    F_f = sp.simplify(sp.Rational(3, 2) * F_tf[0, 0])
    assert sp.simplify(F_tf[1, 1] + F_f / 3) == 0
    assert sp.simplify(F_tf[0, 1]) == 0
    lap = hess_mixed.trace()
    gradsq = grad_mixed.trace()
    fp = sp.diff(f, t)
    trS2 = fp**2 / 6
    N_eq = 3 * sp.diff(nu, t, 2) - lap - 2 * gradsq + trS2 + 3 / t**2 * (sp.exp(-2 * nu) - 1)
    T_eq = sp.diff(f, t, 2) + 2 * (1 / t + sp.diff(nu, t)) * fp + 2 * F_f
    # (M): D_j S^j_y - 2 d_y nu' + (2/t) d_y nu + 2 S^j_y d_j nu + 2 nu' d_y nu
    S = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3)) * fp / 2
    divS = sum(sp.diff(S[j, 0], coords[j]) + sum(G[j][j][k] * S[k, 0] - G[k][j][0] * S[j, k] for k in range(3)) for j in range(3))
    M_eq = (
        divS - 2 * sp.diff(nu, t, y) + 2 / t * sp.diff(nu, y) + 2 * S[0, 0] * sp.diff(nu, y) + 2 * sp.diff(nu, t) * sp.diff(nu, y)
    )
    Rs = Ric_mixed.trace()
    Ham = Rs - 4 * lap - 2 * gradsq - trS2 + 6 * (1 / t + sp.diff(nu, t)) ** 2 - 6 * sp.exp(-2 * nu) / t**2
    return dict(N=N_eq, T=T_eq, M=M_eq, Ham=Ham, F_f=F_f, Rs=Rs)


if __name__ == "__main__":
    e = reduced_equations()
    for k in ("F_f", "Rs"):
        print(k, "=", sp.simplify(e[k]))
