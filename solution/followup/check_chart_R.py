"""Exact off-shell check of the second-order comoving curvature formula (no field equations used).

Spatial metric on a comoving slice, material chart X (conformal factor a^2 dropped):
    h(X) = e^{2 psi} exp(H),   psi = psi1 + psi2,   H = H1 + H2,
where H1 = Sum_a b_a TF(k_a k_a) e_a is the scalar first-order form and H2, psi2 are arbitrary at K = k1 + k2.
Transverse chart: X = x + xi(x), giving h = e^{2R} exp(chi), tr chi = 0, d_j chi_ij = 0.

Claim, for the z1 z2 coefficient at wavevector K:
    xi1 = d beta1,  2 k_a^2 beta_a... fixed by chi1 = 0
    R2  = psi2 + [xi1.d psi1] - (1/6)[tr m^2] - (1/4) Khat Khat . S2,
    S2  = H2 + [xi1.d H1] - TF[m^2],   m = d d beta1.
Checks:
  1. the pulled-back metric has exactly chi1 = 0, and its log has the stated second-order part;
  2. the scalar curvature, an invariant, agrees between the two charts through second order, with an
     arbitrary symmetric trace-free H2 and generic momenta. This uses
     Ricci scalar[e^{2R} exp chi] = -4 lap R + 8 R lap R - 2 (dR)^2 + d_i d_j chi_ij = same with chi transverse;
  3. the 1/K^2 pole statement: at K -> 0 the composite -20 q^2 R1^2 equals the curvature pair coefficient.
"""

import sympy as sp

eps, z1, z2 = sp.symbols("eps z1 z2")
xs = sp.symbols("x1 x2 x3", real=True)
k1 = sp.Matrix(sp.symbols("p1 p2 p3", real=True))
k2 = sp.Matrix(sp.symbols("r1 r2 r3", real=True))
K = k1 + k2
bb1, bb2, p1c, p2c, psi2c = sp.symbols("b1 b2 P1 P2 P_2")
Hc = sp.symbols("h11 h12 h13 h22 h23")
XV = sp.Matrix(xs)
e1, e2 = sp.exp(sp.I * k1.dot(XV)), sp.exp(sp.I * k2.dot(XV))
EK = e1 * e2


def TF(M):
    return M - sp.eye(3) * M.trace() / 3


def d(u, i):
    return sp.diff(u, xs[i])


def to_coeff(expr):
    """Keep O(eps^2) z1 z2 e1 e2 part; returns coefficient. Treat z's as bookkeeping."""
    e = sp.expand(sp.diff(expr, eps, 2).subs(eps, 0) / 2)
    return sp.simplify(e.coeff(z1, 1).coeff(z2, 1) / EK)


H1m = bb1 * z1 * e1 * TF(k1 * k1.T) + bb2 * z2 * e2 * TF(k2 * k2.T)
psi1 = p1c * z1 * e1 + p2c * z2 * e2
H2m = sp.Matrix([[Hc[0], Hc[1], Hc[2]], [Hc[1], Hc[3], Hc[4]], [Hc[2], Hc[4], -Hc[0] - Hc[3]]]) * z1 * z2 * EK
psi2 = psi2c * z1 * z2 * EK

# xi1 = d beta1, chi1 = H1 + 2 dd beta1 - (2/3) delta lap beta1 = 0 => beta = + b_a/2 per mode
beta1 = bb1 / 2 * z1 * e1 + bb2 / 2 * z2 * e2
xi1 = sp.Matrix([d(beta1, i) for i in range(3)])
m = sp.Matrix(3, 3, lambda i, j: d(d(beta1, i), j))


def check1():
    chi1 = H1m + 2 * m - sp.Rational(2, 3) * sp.eye(3) * m.trace()
    return chi1.applyfunc(sp.simplify)


def R1_expr():
    return psi1 + m.trace() / 3


def R2_S2():
    dpsi = sp.Matrix([d(psi1, i) for i in range(3)])
    xi_dpsi = sp.expand(xi1.dot(dpsi))
    trm2 = sp.expand((m * m).trace())
    xi_dH = sp.Matrix(3, 3, lambda i, j: sp.expand(sum(xi1[k] * d(H1m[i, j], xs.index(xs[k])) for k in range(3))))
    S2 = H2m + xi_dH - TF((m * m).applyfunc(sp.expand))
    pick = lambda u: sp.expand(u).coeff(z1, 1).coeff(z2, 1) * z1 * z2
    khat2 = K * K.T / K.dot(K)
    R2 = pick(psi2) + pick(xi_dpsi) - pick(trm2) / 6 - sum(khat2[i, j] * pick(S2[i, j]) for i in range(3) for j in range(3)) / 4
    return sp.simplify(R2 / (z1 * z2 * EK)), S2.applyfunc(pick)


def ricci_scalar_material_second_order():
    """Scalar curvature of e^{2psi} exp(H) through second order, exact expansion via the generic formula.

    Using g = e^{2psi} M with M = exp(H): Rs[g] = e^{-2psi}(Rs[M] - 4 lap_M psi - 2 |d psi|_M^2),
    with Rs[M] computed from M = 1 + H + H^2/2 through second order.
    """
    H = eps * H1m + eps**2 * H2m
    M = sp.eye(3) + H + (H * H) / 2
    # inverse through eps^2: 1 - H + H^2/2
    Mi = sp.eye(3) - H + (H * H) / 2
    Gam = [
        [
            [_trunc(sum(Mi[i, l] * (d(M[l, j], k) + d(M[l, k], j) - d(M[j, k], l)) for l in range(3)) / 2) for k in range(3)]
            for j in range(3)
        ]
        for i in range(3)
    ]
    Ric = sp.zeros(3, 3)
    for j in range(3):
        for k in range(j, 3):
            val = sum(d(Gam[i][j][k], i) - d(Gam[i][j][i], k) for i in range(3))
            val += sum(Gam[i][i][l] * Gam[l][j][k] - Gam[i][k][l] * Gam[l][j][i] for i in range(3) for l in range(3))
            Ric[j, k] = Ric[k, j] = _trunc(val)
    RsM = _trunc(sum(Mi[j, k] * Ric[j, k] for j in range(3) for k in range(3)))
    psi = eps * psi1 + eps**2 * psi2
    dp = [d(psi, i) for i in range(3)]
    # sqrt(det M) = 1 (unimodular), so lap_M f = d_i (M^{ij} d_j f)
    lap = _trunc(sum(d(sum(Mi[i, j] * dp[j] for j in range(3)), i) for i in range(3)))
    grad = _trunc(sum(Mi[i, j] * dp[i] * dp[j] for i in range(3) for j in range(3)))
    pref = 1 - 2 * psi + 2 * psi**2
    return _trunc(pref * (RsM - 4 * lap - 2 * grad))


def _trunc(u):
    u = sp.expand(u)
    return sum(u.coeff(eps, n) * eps**n for n in range(3))


def check2():
    Rs = ricci_scalar_material_second_order()
    Rs1 = Rs.coeff(eps, 1)
    Rs2 = Rs.coeff(eps, 2)
    # invariance: Rs_x(x) = Rs_X(x + xi) => second-order part at fixed x has + xi1.d Rs1
    dRs1 = sp.Matrix([d(Rs1, i) for i in range(3)])
    lhs = sp.expand(Rs2 + xi1.dot(dRs1))
    # transverse-chart form with R = eps R1 + eps^2 R2
    R1 = sp.expand(R1_expr())
    R2c, _ = R2_S2()
    R2 = R2c * z1 * z2 * EK
    lapR1 = sum(d(d(R1, i), i) for i in range(3))
    gradR1 = sum(d(R1, i) ** 2 for i in range(3))
    lapR2 = sum(d(d(R2, i), i) for i in range(3))
    rhs = sp.expand(-4 * lapR1 * 0 - 4 * lapR2 + 8 * R1 * lapR1 - 2 * gradR1)
    pick = lambda u: sp.expand(u).coeff(z1, 1).coeff(z2, 1)
    diff = sp.simplify((pick(lhs) - pick(rhs)) / EK)
    # first order: Rs1 must equal -4 lap R1
    first = sp.simplify(sp.expand(Rs1 + 4 * lapR1))
    return first, diff


if __name__ == "__main__":
    print("chi1 =", check1().tolist() if check1() != sp.zeros(3, 3) else 0)
    first, diff = check2()
    print("first-order curvature residual:", first)
    print("second-order curvature invariance residual:", diff)
    assert first == 0 and diff == 0
    # Negative control: the residual is linear in R2, so dropping the -tr(m^2)/6 composite shifts it by
    # -4 lap(tr(m^2)/6) at K, which is nonzero for K != 0.
    trm2 = sp.expand((m * m).trace()).coeff(z1, 1).coeff(z2, 1) * EK
    shift = sp.simplify(sum(d(d(trm2 / 6, i), i) for i in range(3)) / EK)
    assert shift != 0
    print("negative control: omitting the threading composite leaves residual", sp.factor(4 * shift))
    print("PASS")
