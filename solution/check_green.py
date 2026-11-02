"""Retarded Green functions of the unified triangular system (exact).

lapse   nu'' + (K^2/3) nu - 2 nu/t^2 = S :
        G_N(t,t') = (1/c)[(1 + 1/(x x')) sin(x - x') - ((x - x')/(x x')) cos(x - x')],
        c = K/sqrt3, x = c t, x' = c t';  K -> 0:  G_N^0 = (t^3 - t'^3)/(3 t t')
        regular homogeneous solution  u1(x) = sin x/x - cos x = x j1(x) ~ x^2/3
tensor  h'' + 2h'/t + K^2 h = S :  G_h = (t'/t) sin(K(t - t'))/K ;  K = 0:  G_0 = (t'/t)(t - t')
The particular solution int_0^t G S dt' of a source S = O(t^2) is O(t^4) (lapse), so the t^2
coefficient of nu_r is exactly the free regular coefficient, fixed by the Hamiltonian constraint.
"""

import sympy as sp

t, tp, K = sp.symbols("t tp K", positive=True)
c = K / sp.sqrt(3)
x, xp = c * t, c * tp
GN = ((1 + 1 / (x * xp)) * sp.sin(x - xp) - ((x - xp) / (x * xp)) * sp.cos(x - xp)) / c
GN0 = (t**3 - tp**3) / (3 * t * tp)
Gh = tp / t * sp.sin(K * (t - tp)) / K
G0 = tp / t * (t - tp)


def ok(e, lab):
    assert sp.simplify(e) == 0, (lab, sp.simplify(e))
    print("[ok]", lab)


ok(sp.diff(GN, t, 2) + K**2 / 3 * GN - 2 * GN / t**2, "G_N solves the lapse operator in t")
ok(GN.subs(t, tp), "G_N(t,t) = 0")
ok(sp.diff(GN, t).subs(t, tp) - 1, "dG_N/dt (t=t+) = 1")
ok(sp.limit(GN, K, 0) - GN0, "G_N -> (t^3 - t'^3)/(3 t t') as K -> 0")
ok(sp.diff(GN0, t, 2) - 2 * GN0 / t**2, "G_N^0 solves the K=0 lapse operator")
ok(sp.diff(Gh, t, 2) + 2 * sp.diff(Gh, t) / t + K**2 * Gh, "G_h solves the tensor operator")
ok(Gh.subs(t, tp), "G_h(t,t) = 0")
ok(sp.diff(Gh, t).subs(t, tp) - 1, "dG_h/dt (t=t+) = 1")
ok(sp.limit(Gh, K, 0) - G0, "G_h -> G_0 as K -> 0")
u1 = sp.sin(x) / x - sp.cos(x)
ok(sp.diff(u1, t, 2) + K**2 / 3 * u1 - 2 * u1 / t**2, "u1 = x j1(x) is the regular homogeneous lapse mode")
ok(sp.series(9 * u1 / K**2, t, 0, 3).removeO() - t**2, "9 u1/K^2 = t^2 + O(t^4)")
a = sp.Symbol("a")
part = sp.integrate(GN0 * a * tp**2, (tp, 0, t))
ok(part - a * t**4 / 10, "K=0 particular solution of an O(t^2) source is O(t^4) (leading term a t^4/10)")
print("GREEN FUNCTIONS VERIFIED")
