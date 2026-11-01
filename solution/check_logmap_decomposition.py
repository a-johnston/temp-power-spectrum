"""Exact check of the log-map decomposition of the fourth-order spectrum.

Claim (report, Thm. "log-map decomposition").  Let L = ln T_obs - <ln T_obs>_Omega (centred), with
L = eps theta + eps^2 L2/2 + eps^3 L3/6, theta Gaussian and isotropic with spectrum C_l, and
Theta = e^L/<e^L>_Omega - 1.  For l >= 1,
    (1/4) C^(22)_l + (1/3) Re C^(13)_l = C^{L,4}_l + N_l[C] + M_l,
    C^{L,4}_l = (1/4) C_l[L2, L2] + (1/3) Re C_l[theta, L3],
    N_l[C]    = (1/2) sum_{l1,l2} (2l1+1)(2l2+1)/(4 pi) (l1 l2 l; 0 0 0)^2 C_l1 C_l2 - C_l^2/(2 pi),
    M_l       = (1/2) Re C_l[L2, theta^2] + Re C_l[theta, theta L2]   (collapsed <theta theta L2>).

Model (fully explicit, exact): theta is a band-limited isotropic Gaussian sky (l = 1, 2 with
independent spectra), chi an independent l = 1 Gaussian sky, and
    L2 = P0(theta chi + theta^2),   L3 = P0(theta chi^2)
(non-Gaussian, non-local in the harmonic sense).  The map Theta is expanded from its exact
definition, every harmonic projection is an exact sphere integral of polynomials, and every
ensemble average is an exact Gaussian moment.  Nothing is numerical.
"""
import itertools
import sympy as sp
from sympy.physics.wigner import wigner_3j

x, y, z = sp.symbols('x y z', real=True)
xp, yp, zp = sp.symbols('xp yp zp', real=True)
eps = sp.Symbol('epsilon')
C1, C2, D1 = sp.symbols('C1 C2 D1', positive=True)
pi = sp.pi


def sphere_monomial(a, b, c):
    if a % 2 or b % 2 or c % 2:
        return sp.S.Zero
    return 2 * sp.gamma(sp.Rational(a + 1, 2)) * sp.gamma(sp.Rational(b + 1, 2)) * \
        sp.gamma(sp.Rational(c + 1, 2)) / sp.gamma(sp.Rational(a + b + c + 3, 2))


def sphere_integral(poly):
    p = sp.Poly(sp.expand(poly), x, y, z)
    return sp.expand(sum(coef * sphere_monomial(*mon) for mon, coef in p.terms()))


Y1 = [sp.sqrt(3 / (4 * pi)) * v for v in (x, y, z)]
Y2 = [sp.sqrt(15 / (4 * pi)) * x * y, sp.sqrt(15 / (4 * pi)) * y * z, sp.sqrt(15 / (4 * pi)) * x * z,
      sp.sqrt(15 / (16 * pi)) * (x**2 - y**2), sp.sqrt(5 / (16 * pi)) * (3 * z**2 - (x**2 + y**2 + z**2))]
basis = Y1 + Y2
# orthonormality on the unit sphere (z^2 kept explicit; integrals use exact monomial moments)
for i, j in itertools.product(range(8), repeat=2):
    assert sp.simplify(sphere_integral(basis[i] * basis[j]) - (1 if i == j else 0)) == 0

a = sp.symbols('a0:8', real=True)          # theta coefficients: var C1 (l=1), C2 (l=2)
b = sp.symbols('b0:3', real=True)          # chi coefficients: var D1
variance = {**{a[i]: C1 for i in range(3)}, **{a[i]: C2 for i in range(3, 8)}, **{v: D1 for v in b}}
gauss_vars = list(a) + list(b)

theta = sum(ai * Yi for ai, Yi in zip(a, basis))
chi = sum(bi * Yi for bi, Yi in zip(b, Y1))


def sky_mean(f):
    return sphere_integral(f) / (4 * pi)


def P0(f):
    return sp.expand(f - sky_mean(f))


L2 = P0(theta * chi + theta**2)
L3 = P0(theta * chi**2)


def gaussian_mean(poly):
    p = sp.Poly(sp.expand(poly), *gauss_vars)
    out = 0
    for mon, coef in p.terms():
        term = coef
        for v, k in zip(gauss_vars, mon):
            if k % 2:
                term = 0
                break
            term *= sp.factorial2(k - 1) * variance[v]**(k // 2)
        out += term
    return sp.expand(out)


def legendre_kernel(l):
    """(2l+1)/(4 pi) P_l(n.n') as a polynomial in n=(x,y,z), n'=(xp,yp,zp)."""
    return sp.expand((2 * l + 1) / (4 * pi) * sp.legendre(l, x * xp + y * yp + z * zp))


def spectrum(F, G, l):
    """C_l[F, G] = (1/(2l+1)) <sum_m F_lm G_lm> = (1/(2l+1)) <(2l+1)/(4pi) int int F(n) G(n') P_l>."""
    K = sp.Poly(legendre_kernel(l), xp, yp, zp)
    total = 0
    for mon, coef in K.terms():
        Gm = sphere_integral(sp.expand(G * x**mon[0] * y**mon[1] * z**mon[2]))   # int G(n') n'^mon
        Fm = sphere_integral(sp.expand(F * coef))
        total += Fm * Gm
    return sp.expand(gaussian_mean(total) / (2 * l + 1))


def series_map():
    """Theta from its exact definition, Taylor coefficients Theta^(2), Theta^(3) (with r! convention)."""
    Lfull = eps * theta + eps**2 * L2 / 2 + eps**3 * L3 / 6
    num = sum(Lfull**n / sp.factorial(n) for n in range(4))
    num = sp.expand(num)
    num = sum(num.coeff(eps, k) * eps**k for k in range(4))
    den = sky_mean(num)
    den = sum(sp.expand(den).coeff(eps, k) * eps**k for k in range(4))
    inv = sp.expand(sum((1 - den)**n for n in range(4)))
    inv = sum(inv.coeff(eps, k) * eps**k for k in range(4))
    Th = sp.expand(num * inv - 1)
    T1 = Th.coeff(eps, 1)
    T2 = 2 * Th.coeff(eps, 2)
    T3 = 6 * Th.coeff(eps, 3)
    return T1, T2, T3


def N_formula(l, Cs):
    tot = 0
    for l1, l2 in itertools.product(Cs, repeat=2):
        w = wigner_3j(l1, l2, l, 0, 0, 0)
        tot += sp.Rational(1, 2) * (2 * l1 + 1) * (2 * l2 + 1) / (4 * pi) * w**2 * Cs[l1] * Cs[l2]
    return sp.expand(tot - Cs.get(l, 0)**2 / (2 * pi))


if __name__ == '__main__':
    T1, T2, T3 = series_map()
    assert sp.expand(T1 - theta) == 0
    Cs = {1: C1, 2: C2}
    for l in (1, 2, 3, 4):
        exact = sp.expand(spectrum(T2, T2, l) / 4 + spectrum(T1, T3, l) / 3)
        CL4 = spectrum(L2, L2, l) / 4 + spectrum(theta, L3, l) / 3
        M = spectrum(L2, theta**2, l) / 2 + spectrum(theta, theta * L2, l)
        N = N_formula(l, Cs)
        assert sp.simplify(exact - (CL4 + N + M)) == 0, l
        # The lognormal term alone: theta-only map, L2 = L3 = 0
        Tn2, Tn3 = theta**2, theta**3 - 3 * theta * sky_mean(theta**2)
        assert sp.simplify(spectrum(Tn2, Tn2, l) / 4 + spectrum(theta, Tn3, l) / 3 - N) == 0, l
        print(f'[ok] l={l}: exact quartic spectrum = C^L4 + N_l[C] + M_l;  N_l = {sp.factor(N)}')
    print('LOG-MAP DECOMPOSITION VERIFIED (exact sphere integrals and Gaussian moments)')
