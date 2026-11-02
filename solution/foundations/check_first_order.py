"""First-order (epsilon^1) solution and transfer function, verified exactly.

Units 8 pi G = 1; background a(t) = a1 t (t = conformal time = temperature-frame time
at zeroth order), rho = 3 a1^2 / N^4.  One Fourier mode zeta_k e^{i k z}.
x = k t / sqrt(3).

(A) temperature frame  g = N^2(-dt^2 + gamma_ij dX dX):
      N = a1 t (1 + eps nu e^{ikz}),
      gamma_ij = delta_ij + 2 eps (C delta_ij + d_i d_j E) e^{ikz}
    regular adiabatic growing mode (no decaying mode, gamma -> e^{2 zeta} delta, nu -> 0):
      nu = -zeta (sin x / x - cos x),  C = zeta (2 sin x / x - cos x),
      k^2 E = 3 zeta (2 sin x/x - cos x - 1)   [so 3C - k^2 E = 3 zeta: expansion-free]
(B) Newtonian gauge  ds^2 = a^2[-(1+2Phi)dt^2 + (1-2Psi)dx^2], Phi = Psi = -2 zeta (sin x - x cos x)/x^3,
    with delta_N and velocity potential v_N fixed by the constraints.

Checks: all 10 Einstein equations at O(eps) for (A) and (B); regularity and t->0 limits;
the exact line-of-sight transfer functions of both routes coincide for every l >= 1
(proved for symbolic l via d/dt_* of the difference = 0 plus the coincidence limit).
"""

import sympy as sp
from gr import assert_zero, christoffel, ricci

t, z, xx, yy = sp.symbols("t z x y", real=True)
k, a1, zeta, eps = sp.symbols("k a1 zeta epsilon", positive=True)
X = [t, xx, yy, z]
s3 = sp.sqrt(3)
xv = k * t / s3
ph = sp.exp(sp.I * k * z)


def einstein_residual(g, Tdn, label):
    Gam, ginv = christoffel(g, X)
    Ric = ricci(Gam, X)
    R = sum(ginv[m, n] * Ric[m, n] for m in range(4) for n in range(4))
    for m in range(4):
        for n in range(m, 4):
            res = Ric[m, n] - g[m, n] * R / 2 - Tdn[m, n]
            res1 = sp.diff(res, eps).subs(eps, 0)
            res0 = res.subs(eps, 0)
            assert sp.simplify(res0) == 0, (label, m, n, "background")
            r = sp.simplify(sp.expand((res1 / ph).rewrite(sp.exp)))
            assert r == 0, (label, m, n, r)
    print(f"[ok] {label}: G_mu nu = T_mu nu at O(eps^0) and O(eps^1), all 10 components")


# ------------------------------------------------------------- (A) temperature frame
nu = -zeta * (sp.sin(xv) / xv - sp.cos(xv))
C = zeta * (2 * sp.sin(xv) / xv - sp.cos(xv))
k2E = 3 * zeta * (2 * sp.sin(xv) / xv - sp.cos(xv) - 1)
N = a1 * t * (1 + eps * nu * ph)
gam = sp.diag(1 + 2 * eps * C * ph, 1 + 2 * eps * C * ph, 1 + 2 * eps * (C - k2E) * ph)
g = sp.zeros(4, 4)
g[0, 0] = -(N**2)
for i in range(3):
    g[i + 1, i + 1] = N**2 * gam[i, i]
rho = 3 * a1**2 / N**4
p = rho / 3
udn = [-N, 0, 0, 0]
Tdn = sp.Matrix(4, 4, lambda m, n: (rho + p) * udn[m] * udn[n] + p * g[m, n])
einstein_residual(g, Tdn, "(A) temperature-frame linear solution")
assert_zero(3 * C - k2E - 3 * zeta, "(A) expansion-free: 3C - k^2E = 3 zeta for all t")
for f, lim, nm in [(nu, 0, "nu"), (C, zeta, "C"), (k2E, 0, "k^2 E")]:
    ser = sp.expand(sp.series(f, t, 0, 10).removeO())
    assert sp.simplify(ser.subs(t, 0) - lim) == 0
    assert all(sp.simplify(ser.coeff(t, p_)) == 0 for p_ in (1, 3, 5, 7, 9)), nm
    assert not any(sp.simplify(ser.coeff(t, -p_)) != 0 for p_ in range(1, 6)), nm
    print(f"[ok] (A) {nm} regular: no negative and no odd powers of t through t^9, limit {lim} as t->0")

# ------------------------------------------------------------- (B) Newtonian gauge
a = a1 * t
Phi = -2 * zeta * (sp.sin(xv) - xv * sp.cos(xv)) / xv**3
H = 1 / t
lapPsi = -(k**2) * Phi
vN = -(t**2 / 2) * (sp.diff(Phi, t) + H * Phi)  # v^i = d_i vN  (0i constraint)
dN = (2 * t**2 / 3) * (lapPsi - 3 * H * (sp.diff(Phi, t) + H * Phi))  # 00 constraint
gN = sp.diag(-(a**2) * (1 + 2 * eps * Phi * ph), *([a**2 * (1 - 2 * eps * Phi * ph)] * 3))
rhoN = 3 / (a1**2 * t**4) * (1 + eps * dN * ph)
pN = rhoN / 3
# u^mu = a^-1 (1 - eps Phi, 0, 0, eps i k vN)  -> lower with gN
uup = [(1 - eps * Phi * ph) / a, 0, 0, eps * sp.I * k * vN * ph / a]
udnN = [sum(gN[m, n] * uup[n] for n in range(4)) for m in range(4)]
TdnN = sp.Matrix(4, 4, lambda m, n: (rhoN + pN) * udnN[m] * udnN[n] + pN * gN[m, n])
einstein_residual(gN, TdnN, "(B) Newtonian-gauge linear solution (Phi=Psi, delta_N, v_N)")

# ------------------------------------------------------------- transfer functions
# Theta^(1)(n) = int d^3k/(2pi)^3 zeta_k e^{ik.x_o} sum_l (2l+1) i^l T_l(k) P_l(khat.n)
# (A): T^A_l = - int_{t*}^{to} dt k^2 E'(t) [ j_l/3 + j_l'' ](k (to - t))
# (B): T^B_l = [delta_N/4 + Phi]_* j_l(k r*) - k vN(t*) j_l'(k r*) + int 2 Phi' j_l + k vN(to) j_l'(0)
ts, to, L = sp.symbols("t_star t_o L", positive=True)  # L = l(l+1)
J0, J1, J2, J3 = sp.symbols("J0 J1 J2 J3")  # j_l, j_l', j_l'', j_l''' at y = k(to - t*)
yv = k * (to - ts)
bessel_J2 = -(2 / yv) * J1 - (1 - L / yv**2) * J0  # spherical Bessel ODE
# d/dt* of T^A: + k^2 E'(t*) [J0/3 + J2]
dTA = (sp.diff(k2E, t) * (J0 / 3 + J2)).subs(t, ts)
# d/dt* of T^B: d/dt*[F j(y)] = F' j - k F j' ; d/dt*[-k vN j'(y)] = -k vN' j' + k^2 vN j''; ISW: -2Phi'(t*) j
Fb = dN / 4 + Phi
dTB = (sp.diff(Fb, t) * J0 - k * Fb * J1 - k * sp.diff(vN, t) * J1 + k**2 * vN * J2 - 2 * sp.diff(Phi, t) * J0).subs(t, ts)
# structural identity (found by the independent auditor): the temperature-frame shear rate
# equals the Newtonian-gauge velocity gradient, (k^2 E)' = k^2 v_N
assert sp.simplify(sp.diff(k2E, t) - k**2 * vN) == 0
print("[ok] (k^2 E)' = k^2 v_N  (temperature-frame shear rate = Newtonian fluid velocity gradient)")
Draw = sp.expand(dTA - dTB)
# The j'' (J2) terms cancel identically by the identity above, so the result holds for ANY
# kernel j (stronger than needed); we still substitute the Bessel ODE to be explicit.
assert sp.simplify(Draw.coeff(J2)) == 0
D = sp.expand(Draw.subs(J2, bessel_J2))
cJ0 = sp.simplify(D.coeff(J0))
cJ1 = sp.simplify(D.coeff(J1))
assert cJ0 == 0 and cJ1 == 0, (cJ0, cJ1)
print("[ok] d/dt*(T^A_l - T^B_l) = 0 identically, for symbolic l (indeed for any radial kernel j)")
# coincidence limit t* -> to (y -> 0): j_l(0)=delta_l0, j_l'(0)=delta_l1/3. T^A -> 0.
# T^B(t*=to) = [dN/4+Phi](to) delta_l0 - k vN(to) delta_l1/3 + k vN(to) delta_l1/3
yy0 = sp.Symbol("y", positive=True)
tt_ = sp.Symbol("tt", positive=True)


def TB_full(l, tstar):
    """Complete Newtonian-gauge transfer function as ONE expression (no retyping)."""
    j = lambda arg: sp.expand_func(sp.jn(l, arg))
    jp = lambda arg: sp.diff(sp.expand_func(sp.jn(l, yy0)), yy0).subs(yy0, arg)
    jp0 = sp.limit(sp.diff(sp.expand_func(sp.jn(l, yy0)), yy0), yy0, 0)
    em = (Fb * j(k * (to - t)) - k * vN * jp(k * (to - t))).subs(t, tstar)
    isw = sp.Integral((2 * sp.diff(Phi, t)).subs(t, tt_) * j(k * (to - tt_)), (tt_, tstar, to))
    obs = k * vN.subs(t, to) * jp0
    return em, isw, obs


def TA_full(l, tstar):
    w = sp.expand_func(sp.jn(l, yy0)) / 3 + sp.diff(sp.expand_func(sp.jn(l, yy0)), yy0, 2)
    return -sp.Integral((sp.diff(k2E, t) / zeta).subs(t, tt_) * zeta * w.subs(yy0, k * (to - tt_)), (tt_, tstar, to))


for l in range(1, 5):
    em, isw, obs = TB_full(l, ts)
    lim_em = sp.limit(em, ts, to)  # emission-endpoint terms as t* -> to
    tb = sp.simplify(lim_em + isw.subs(ts, to).doit() + obs)
    ta = TA_full(l, ts).subs(ts, to).doit()
    assert tb == 0 and ta == 0, (l, tb, ta)
print("[ok] coincidence t*->to (l=1..4, from the complete expressions): T^B -> 0 = T^A")
print("     (l=1: emitter and observer Doppler cancel; l>=2 trivial since j_l(0)=j_l'(0)=0)")
print("     => T^A_l == T^B_l for all l >= 1, all k, t*, to: first-order gauge invariance verified")

# superhorizon behaviour of the source: k^2 E' = O(k^2) -> T_l = O(k^2 t^2) (no anisotropy from constant zeta)
assert sp.simplify(sp.series(sp.diff(k2E, t), k, 0, 2).removeO()) == 0
print("[ok] k^2 E'(t) = O(k^2): constant-zeta (k->0) mode produces no anisotropy (Weinberg adiabatic mode)")
print("FIRST ORDER VERIFIED")
