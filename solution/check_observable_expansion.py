"""Generic verification of the observable expansion formulas through third order.

All fields are ARBITRARY (SymPy undefined functions of (t, x, y, z, eps)), constrained only by
the vanishing of the background shear/connection.  Each displayed formula of
the report's section 'The observable through third order' is compared with the direct eps-expansion of its exact
definition.  WLOG nhat = (0, 0, 1) (all formulas are rotation covariant); background ray
Xbar(t) = nhat (t_o - t), ebar = -nhat.

Notation: a quantity Y has  Y = eps Y1 + eps^2 Y2/2 + eps^3 Y3/6.
"""
import sympy as sp

t, to, tb, eps, lam = sp.symbols('t t_o tbar eps lambda', positive=True)
xs = sp.symbols('x y z')
nh = sp.Matrix([0, 0, 1])
order = 3


def taylor_coeffs(expr, n=order):
    """[Y0, Y1, Y2, Y3] with Y = sum eps^r Y_r / r!  (Y_r = d^r Y/d eps^r at 0)."""
    out, e = [], expr
    for r in range(n + 1):
        out.append(sp.expand(e.subs(eps, 0).doit()))
        e = sp.diff(e, eps)
    return out


_monos = {}


def _monomials(deg):
    if deg not in _monos:
        out = []
        for a in range(deg + 1):
            for b in range(deg + 1 - a):
                for c in range(deg + 1 - a - b):
                    out.append((a, b, c))
        _monos[deg] = out
    return _monos[deg]


def field(name):
    """An arbitrary field vanishing on the background, as its general jet about the background
    ray:  sum_{r>=1} eps^r/r! sum_alpha c_{r,alpha}(t) (X - Xbar(t))^alpha / alpha!,
    with independent arbitrary functions c of t and every spatial derivative order that the
    third-order formulas can reach (2 at r=1, 1 at r=2, 0 at r=3)."""
    dev = [xs[0], xs[1], xs[2] - (to - t)]
    out = 0
    for r, deg in ((1, 2), (2, 1), (3, 0)):
        for al in _monomials(deg):
            c = sp.Function(f'{name}_{r}_{al[0]}{al[1]}{al[2]}')(t)
            mon = dev[0]**al[0] * dev[1]**al[1] * dev[2]**al[2] / (
                sp.factorial(al[0]) * sp.factorial(al[1]) * sp.factorial(al[2]))
            out += eps**r / sp.factorial(r) * c * mon
    return out


# ---------------------------------------------------------------- abstract fields
sig = sp.Matrix(3, 3, lambda i, j: field(f's{min(i, j)}{max(i, j)}'))           # sigma_jk (symmetric)
Smix = sp.Matrix(3, 3, lambda i, j: field(f'S{i}{j}'))                           # S^i_j
Gam = [[[field(f'G{a}{min(b, c)}{max(b, c)}') for c in range(3)] for b in range(3)] for a in range(3)]
nu = field('nu')

# ray perturbations: arbitrary functions of t (their equations are checked separately)
Xr = [sp.Matrix(3, 1, lambda i, _: sp.Function(f'X{r}_{i}')(t)) for r in (1, 2)]
Xpert = eps * Xr[0] + eps**2 / 2 * Xr[1]
er = [x_.diff(t) for x_ in Xr]
Xbar = nh * (to - t)
Xfull = Xbar + Xpert
efull = -nh + Xpert.diff(t)
on_ray = {xs[i]: Xfull[i] for i in range(3)}


def at_ray(e):
    return e.subs(on_ray, simultaneous=True)


def bar(e):
    """evaluate an expression at the background ray (eps-independent positions)."""
    return e.subs({xs[i]: Xbar[i] for i in range(3)}, simultaneous=True)


def grad(f):
    return sp.Matrix([sp.diff(f, xx) for xx in xs])


def dX(f, V):
    """(V . d) f, f a scalar field expression in (t, x, y, z, eps)."""
    return sum(V[i] * sp.diff(f, xs[i]) for i in range(3))


def nth(Fexpr, r):
    """r-th eps-Taylor coefficient d^r/deps^r at 0 of an expression (field values off the ray)."""
    return sp.diff(Fexpr, eps, r).subs(eps, 0)


def check(lhs, rhs, label):
    d = sp.expand((lhs - rhs).doit())
    assert d == 0, (label, d)
    print(f'[ok] {label}')


# =============================================== (1) integrand F = sigma_jk e^j e^k
F = at_ray((efull.T * sig * efull)[0, 0])
Fc = taylor_coeffs(F)
s_r = lambda r: sig.applyfunc(lambda e: nth(e, r))       # sigma_r (field, off-ray)
nn = lambda M: (nh.T * M * nh)[0, 0]
bil = lambda M, a, b: (a.T * M * b)[0, 0]
F1 = bar(nn(s_r(1)))
F2 = bar(nn(s_r(2)) + 2 * dX(nn(s_r(1)), Xr[0]) - 4 * bil(s_r(1), nh, er[0]))
F3 = bar(nn(s_r(3)) + 3 * dX(nn(s_r(2)), Xr[0]) - 6 * bil(s_r(2), nh, er[0])
         + 3 * dX(nn(s_r(1)), Xr[1]) + 3 * dX(dX(nn(s_r(1)), Xr[0]), Xr[0])
         - 12 * dX(bil(s_r(1), nh, er[0]), Xr[0]) - 6 * bil(s_r(1), nh, er[1])
         + 6 * bil(s_r(1), er[0], er[0]))
check(Fc[0], 0, 'F0 = 0')
check(Fc[1], F1, 'F1 = sigma1_nn')
check(Fc[2], F2, 'F2 = sigma2_nn + 2 X1.d sigma1_nn - 4 sigma1(n, e1)')
check(Fc[3], F3, 'F3 (eight terms; no X3/e3 because sigma vanishes on the background)')

# =============================================== (2) ray acceleration  A = -2 S e - Gamma(e,e) + e sigma(e,e)
def Avec(Sm, G, sg, X, e):
    out = sp.zeros(3, 1)
    see = (e.T * sg * e)[0, 0]
    for i in range(3):
        out[i] = -2 * sum(Sm[i, j] * e[j] for j in range(3)) \
            - sum(G[i][j][k] * e[j] * e[k] for j in range(3) for k in range(3)) + e[i] * see
    return out


Afull = at_ray(Avec(Smix, Gam, sig, Xfull, efull))
Ac = [taylor_coeffs(a, 2) for a in Afull]
S_r = lambda r: Smix.applyfunc(lambda e: nth(e, r))
G_r = lambda r: [[[nth(Gam[a][b][c], r) for c in range(3)] for b in range(3)] for a in range(3)]
def A_lin(r, e, w):
    """2 S_r w - Gamma_r(e, w) - Gamma_r(w, e) + ... : directional pieces used below."""
    pass
A1 = bar(2 * S_r(1) * nh - sp.Matrix([sum(G_r(1)[i][j][k] * nh[j] * nh[k] for j in range(3) for k in range(3))
                                         for i in range(3)]) - nh * nn(s_r(1)))
def A1field(Sm, G, sg):
    return 2 * Sm * nh - sp.Matrix([sum(G[i][j][k] * nh[j] * nh[k] for j in range(3) for k in range(3))
                                    for i in range(3)]) - nh * nn(sg)
A2 = bar(A1field(S_r(2), G_r(2), s_r(2))
         + 2 * A1field(S_r(1), G_r(1), s_r(1)).applyfunc(lambda f: dX(f, Xr[0]))
         + 2 * (-2 * S_r(1) * er[0]
                + 2 * sp.Matrix([sum(G_r(1)[i][j][k] * nh[j] * er[0][k] for j in range(3) for k in range(3))
                                 for i in range(3)])
                + er[0] * nn(s_r(1)) + 2 * nh * bil(s_r(1), nh, er[0])))
for i in range(3):
    check(Ac[i][0], 0, f'A0_{i} = 0')
    check(Ac[i][1], A1[i], f'A1_{i} = 2 S1 n - Gamma1(n,n) - n sigma1_nn')
    check(Ac[i][2], A2[i], f'A2_{i} (field, transport and direction terms)')

# ray solution with final conditions X_r(t_o) = 0, X_r'(t_o) = v_r
v = sp.Symbol('v')
Aab = sp.Function('Aab')
Xsol = v * (t - to) + sp.Integral((sp.Symbol('tp') - t) * Aab(sp.Symbol('tp')), (sp.Symbol('tp'), t, to))
check(sp.diff(Xsol, t, 2), Aab(t), "X_r'' = A_r  for X_r = v_r (t - t_o) + int_t^{t_o} (t'-t) A_r(t') dt'")
check(Xsol.subs(t, to), 0, 'X_r(t_o) = 0')
check(sp.diff(Xsol, t).subs(t, to), v, "X_r'(t_o) = v_r")

# =============================================== (3) emission root:  ln t_e + phi(t_e) = ln tbar
phi = sp.Function('phi')                       # phi(t, eps) = nu(t, X(t; eps); eps), phi(t,0)=0
ph = lambda tt: eps * phi(tt, eps)
de = sp.symbols('d1 d2 d3')
te = tb + eps * de[0] + eps**2 / 2 * de[1] + eps**3 / 6 * de[2]
root = sp.log(te / tb) + ph(te)
rc = taylor_coeffs(root)
phr = [sp.diff(eps * phi(t, eps), eps, r).subs(eps, 0) for r in range(4)]
d1 = -tb * phr[1].subs(t, tb)
u1 = d1 / tb
dphi = lambda r, k=1: sp.diff(phr[r], t, k).subs(t, tb)
u2 = u1**2 - phr[2].subs(t, tb) - 2 * dphi(1) * d1
d2 = tb * u2
u3 = 3 * u1 * u2 - 2 * u1**3 - phr[3].subs(t, tb) - 3 * dphi(2) * d1 - 3 * dphi(1) * d2 - 3 * dphi(1, 2) * d1**2
d3 = tb * u3
sol = sp.solve(rc[1], de[0])[0]
check(sol, d1, 'd1 = -tbar phi1')
sol2 = sp.solve(rc[2].subs(de[0], d1), de[1])[0]
check(sol2, d2, 'd2 = tbar [u1^2 - phi2 - 2 phi1\' d1]')
sol3 = sp.solve(rc[3].subs({de[0]: d1, de[1]: d2}), de[2])[0]
check(sol3, d3, 'd3 = tbar [3u1u2 - 2u1^3 - phi3 - 3phi2\'d1 - 3phi1\'d2 - 3phi1\'\'d1^2]')
# phi_r in terms of nu and ray displacements
phifull = at_ray(nu)
pc = taylor_coeffs(phifull)
nu_r = lambda r: nth(nu, r)
check(pc[1], bar(nu_r(1)), 'phi1 = nu1')
check(pc[2], bar(nu_r(2) + 2 * dX(nu_r(1), Xr[0])), 'phi2 = nu2 + 2 X1.d nu1')
check(pc[3], bar(nu_r(3) + 3 * dX(nu_r(2), Xr[0]) + 3 * dX(nu_r(1), Xr[1]) + 3 * dX(dX(nu_r(1), Xr[0]), Xr[0])),
      'phi3 = nu3 + 3 X1.d nu2 + 3 X2.d nu1 + 3 (X1.d)^2 nu1')

# =============================================== (4) moving lower limit:  I = int_{t_e}^{t_o} F dt
Fg = sp.Function('Fg')                         # F(t, eps) = eps * Fg(t, eps) along the ray
Ffun = lambda tt: eps * Fg(tt, eps)
Iexact = sp.Integral(Ffun(t), (t, te, to))
Ic = taylor_coeffs(Iexact)
Fr = [sp.diff(eps * Fg(t, eps), eps, r).subs(eps, 0) for r in range(4)]
LOS = lambda r: sp.Integral(Fr[r], (t, tb, to))
Fe = lambda r, k=0: sp.diff(Fr[r], t, k).subs(t, tb)
check(Ic[1], LOS(1), 'I1 = int F1')
check(Ic[2], LOS(2) - 2 * Fe(1) * de[0], 'I2 = int F2 - 2 F1(te) d1')
check(Ic[3], LOS(3) - 3 * Fe(1) * de[1] - 3 * Fe(2) * de[0] - 3 * Fe(1, 1) * de[0]**2,
      "I3 = int F3 - 3 F1 d2 - 3 F2 d1 - 3 F1' d1^2")

# =============================================== (5) observer clock (fixed proper time from T = T_I)
tIb, tob = sp.symbols('tbar_I tbar_o', positive=True)
om = sp.Function('omega')                      # omega(t, eps) = nu(t, X_o; eps)
omf = lambda tt: eps * om(tt, eps)
dI = sp.symbols('dI1 dI2')
do = sp.symbols('do1 do2')
tI = tIb + eps * dI[0] + eps**2 / 2 * dI[1]
tO = tob + eps * do[0] + eps**2 / 2 * do[1]
tp = sp.Symbol('tp', positive=True)
clock = sp.Integral(tp * sp.exp(omf(tp)), (tp, tI, tO)) - (tob**2 - tIb**2) / 2
cc = taylor_coeffs(clock, 2)
omr = [sp.diff(eps * om(t, eps), eps, r).subs(eps, 0) for r in range(3)]
W1 = lambda T: sp.Integral(tp * omr[1].subs(t, tp), (tp, 0, T))
W2 = lambda T: sp.Integral(tp * (omr[2] + omr[1]**2).subs(t, tp), (tp, 0, T))
dI1 = -tIb * omr[1].subs(t, tIb)
uI1 = dI1 / tIb
dI2 = tIb * (uI1**2 - omr[2].subs(t, tIb) - 2 * sp.diff(omr[1], t).subs(t, tIb) * dI1)
do1 = (tIb * dI1 - W1(tob) + W1(tIb)) / tob
do2 = (tIb * dI2 + dI1**2 - W2(tob) - 2 * tob * omr[1].subs(t, tob) * do1
       + W2(tIb) + 2 * tIb * omr[1].subs(t, tIb) * dI1 - do1**2) / tob
# the clock condition, with W written as int_0, is linear in the shifts at each order
c1 = cc[1].subs({dI[0]: dI1, do[0]: do1}).doit()
c2 = cc[2].subs({dI[0]: dI1, do[0]: do1, dI[1]: dI2, do[1]: do2}).doit()
for lab, c_ in (('clock O(eps)', c1), ('clock O(eps^2)', c2)):
    c_ = sp.expand(c_)
    # split int_{tIb}^{tob} into int_0^{tob} - int_0^{tIb}
    c_ = c_.replace(lambda e: isinstance(e, sp.Integral) and e.limits[0][1] == tIb and e.limits[0][2] == tob,
                    lambda e: sp.Integral(e.function, (e.limits[0][0], 0, tob))
                    - sp.Integral(e.function, (e.limits[0][0], 0, tIb)))
    c_ = sp.expand(c_)
    terms = sp.Add.make_args(c_)
    integ = {}
    rest = 0
    for term in terms:
        ints = [a for a in term.atoms(sp.Integral)]
        if ints:
            assert len(ints) == 1
            key = (ints[0].function.subs(ints[0].limits[0][0], tp), ints[0].limits[0][2])
            integ[key] = integ.get(key, 0) + sp.cancel(term / ints[0])
        else:
            rest += term
    # integrals of distinct integrands/limits are independent: their coefficients must cancel
    left = sum(cf * sp.Integral(f, (tp, 0, up)) for (f, up), cf in integ.items() if sp.simplify(cf) != 0)
    tot = sp.simplify(sp.expand(rest + left))
    if tot != 0:
        # integrands that are sums: expand linearly
        tot = sp.simplify(sp.expand(sp.expand(tot).replace(
            lambda e: isinstance(e, sp.Integral),
            lambda e: sum(sp.Integral(a, *e.limits) for a in sp.Add.make_args(sp.expand(e.function))))))
    assert tot == 0, (lab, tot)
    print(f'[ok] {lab}: do_r and dI_r formulas solve the fixed-proper-time condition')

# =============================================== (6) composition with the random observation time
Thc = sp.Function('Thc')                       # coordinate-event map Theta_c(t_o, eps) = eps Thc
dd = sp.symbols('D1 D2')
tO2 = tob + eps * dd[0] + eps**2 / 2 * dd[1]
Th = eps * Thc(tO2, eps)
thc = taylor_coeffs(Th)
Tr = [sp.diff(eps * Thc(t, eps), eps, r).subs(eps, 0) for r in range(4)]
d_ = lambda r, k: sp.diff(Tr[r], t, k).subs(t, tob)
check(thc[2], d_(2, 0) + 2 * dd[0] * d_(1, 1), 'Theta2 = Thc2 + 2 D1 dThc1')
check(thc[3], d_(3, 0) + 3 * dd[1] * d_(1, 1) + 3 * dd[0] * d_(2, 1) + 3 * dd[0]**2 * d_(1, 2),
      'Theta3 = Thc3 + 3 D2 dThc1 + 3 D1 dThc2 + 3 D1^2 d2Thc1')
print('OBSERVABLE EXPANSION FORMULAS VERIFIED FOR ARBITRARY FIELDS')
