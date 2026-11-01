"""Superhorizon limit k0 t_o -> 0 of the complete observable, exact in eps, and the resulting
closed forms for C_2^(11), C_2^(22), C_2^(13) at leading order O(k0^4).

Leading-order observable (derived in the report, Sec. 9.2): with tau^2 = (t_o^2 - tbar_*^2)/6,
    Theta(nhat) = -tau^2 e^{-2 zeta} TF[d_i d_j zeta - d_i zeta d_j zeta] nhat^i nhat^j + O(k0^3),
all derivatives at the observer, EXACT in zeta (the O(t^2) gradient-expansion data are exact
in zeta).  Ray bending, time delay, endpoint shifts, the clock and the frame enter at O(k0^3) or
beyond; the normalisation removes nothing at this order (the map is a pure quadrupole).

Statistics (monochromatic, <g^2> = sigma^2): one-point Gaussian covariances of (g, g_i, g_ij):
    <g g> = s2, <g g_ij> = -s2 k^2 d_ij/3, <g_i g_j> = s2 k^2 d_ij/3,
    <g_ij g_kl> = s2 k^4 (d_ij d_kl + d_ik d_jl + d_il d_jk)/15, odd-derivative cross terms vanish.
"""
import itertools
import sympy as sp

eps, s2, k, tau = sp.symbols('epsilon sigma2 k0 tau', positive=True)

# ---------------------------------------------------------------- Gaussian variables at a point
names = ['g'] + [f'g{i}' for i in range(3)] + [f'g{i}{j}' for i in range(3) for j in range(i, 3)]
V = {n: sp.Symbol(n) for n in names}
g = V['g']
gi = [V[f'g{i}'] for i in range(3)]
gij = [[V[f'g{min(i, j)}{max(i, j)}'] for j in range(3)] for i in range(3)]
d = lambda i, j: 1 if i == j else 0


def cov(a, b):
    def idx(n):
        return [int(c) for c in n[1:]]
    ia, ib = idx(a.name), idx(b.name)
    na, nb = len(ia), len(ib)
    if na == 0 and nb == 0:
        return s2
    if (na + nb) % 2 == 1:
        return 0
    if {na, nb} == {0, 2}:
        i, j = ia if na == 2 else ib
        return -s2 * k**2 * d(i, j) / 3
    if na == 1 and nb == 1:
        return s2 * k**2 * d(ia[0], ib[0]) / 3
    i, j = ia
    m, n = ib
    return s2 * k**4 * (d(i, j) * d(m, n) + d(i, m) * d(j, n) + d(i, n) * d(j, m)) / 15


def wick(monomial_factors):
    """<prod of Gaussian variables> by explicit enumeration of pairings (exact)."""
    f = list(monomial_factors)
    if not f:
        return 1
    if len(f) % 2:
        return 0
    a = f[0]
    tot = 0
    for j in range(1, len(f)):
        rest = f[1:j] + f[j + 1:]
        c = cov(a, f[j])
        if c != 0:
            tot += c * wick(rest)
    return tot


def expect(poly):
    poly = sp.expand(poly)
    out = 0
    for term in sp.Add.make_args(poly):
        c, factors = 1, []
        for fct in sp.Mul.make_args(term):
            b, e = fct.as_base_exp()
            if b in V.values():
                factors += [b] * int(e)
            else:
                c *= fct
        out += c * wick(factors)
    return sp.expand(out)


# ---------------------------------------------------------------- the map, exact in eps (to eps^3)
zeta = eps * g
E = sp.series(sp.exp(-2 * zeta), eps, 0, 4).removeO()
Phi = sp.Matrix(3, 3, lambda i, j: eps * gij[i][j] - eps**2 * gi[i] * gi[j])
Psi = (E * Phi).applyfunc(lambda e: sp.expand(sp.series(sp.expand(e), eps, 0, 4).removeO()))
Psi = Psi - sp.eye(3) * Psi.trace() / 3                      # trace-free part
Psi_r = {r: Psi.applyfunc(lambda e: sp.factorial(r) * e.coeff(eps, r)) for r in (1, 2, 3)}
# Theta^(r)(n) = -tau^2 Psi_r[ij] n_i n_j

# ---------------------------------------------------------------- harmonic normalisation, exact
# For a pure quadrupole f_B(n) = B_ij n_i n_j (B trace-free), the addition theorem gives
# int dOmega' P_2(n.n') f_B(n') = (4 pi/5) f_B(n), hence
#   sum_m a_2m[A] a_2m[B]* = (5/4pi) int int P_2(n.n') f_A f_B = int dOmega f_A f_B.
# The last integral is evaluated with exact spherical moments.
nv = sp.symbols('n1 n2 n3', real=True)
nvec = sp.Matrix(nv)
A = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'A{min(i, j)}{max(i, j)}'))
A = A - sp.eye(3) * A.trace() / 3
Bm = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'B{min(i, j)}{max(i, j)}'))
Bm = Bm - sp.eye(3) * Bm.trace() / 3
fA = (nvec.T * A * nvec)[0, 0]
fB = (nvec.T * Bm * nvec)[0, 0]


def sphere_integral(e):
    P = sp.Poly(sp.expand(e), *nv)
    out = 0
    df = lambda m: sp.factorial2(m - 1) if m > 0 else 1
    for (a, b, c), cf in P.terms():
        if a % 2 == 0 and b % 2 == 0 and c % 2 == 0:
            out += cf * 4 * sp.pi * df(a) * df(b) * df(c) / sp.factorial2(a + b + c + 1)
    return sp.expand(out)


S = sphere_integral(fA * fB)
contr = sum(A[i, j] * Bm[i, j] for i in range(3) for j in range(3))
ratio = sp.simplify(S / sp.expand(contr))
assert ratio == 8 * sp.pi / 15, ratio
# the map is a pure quadrupole: its l = 0 part vanishes (trace-free), and no other l occurs
assert sphere_integral(fA) == 0
print('[ok] sum_m a_2m[A] a_2m[B]* = (8 pi/15) A:B for trace-free A, B; only l = 2 is present')

def C2(r, s_):
    """C_2^(rs) = (1/5) sum_m <a^(r) a^(s)*>_c = (1/5)(8 pi/15) tau^4 <Psi_r : Psi_s>  (connected)."""
    val = expect(sum(Psi_r[r][i, j] * Psi_r[s_][i, j] for i in range(3) for j in range(3)))
    disc = sum(expect(Psi_r[r][i, j]) * expect(Psi_r[s_][i, j]) for i in range(3) for j in range(3))
    return sp.factor(sp.Rational(1, 5) * 8 * sp.pi / 15 * tau**4 * (val - disc))

C11 = C2(1, 1)
C22 = C2(2, 2)
C13 = C2(1, 3)
total4 = sp.factor(C22 / 4 + C13 / 3)
print('C_2^(11) =', C11)
print('C_2^(22) =', C22)
print('C_2^(13) =', C13)
print('C_2 eps^4 coefficient  [C22/4 + Re C13/3] =', total4)
print('ratio of fourth-order correction to eps^2 C11:  eps^2 *', sp.factor(total4 / C11))
# direct check: expand <Psi:Psi> to eps^4 as a whole
full = expect(sum(Psi[i, j] * Psi[i, j] for i in range(3) for j in range(3)))
full = sp.expand(sp.Rational(1, 5) * 8 * sp.pi / 15 * tau**4 * full)
assert sp.simplify(full.coeff(eps, 2) - C11) == 0
assert sp.simplify(full.coeff(eps, 4) - total4) == 0
print('[ok] direct eps-expansion of <|a_2m|^2> reproduces eps^2 C11 + eps^4 [C22/4 + C13/3]')
# first-order cross-check against the line-of-sight transfer function: T_2 -> -(2/15) k^2 tau^2
assert sp.simplify(C11 - 4 * sp.pi * s2 * (sp.Rational(2, 15) * k**2 * tau**2)**2) == 0
print('[ok] C_2^(11) = 4 pi sigma^2 T_2^2 with T_2 = -(2/15) k0^2 tau^2 (superhorizon limit of T_2)')
with open('results/superhorizon_C2.txt', 'w') as fh:
    fh.write(f'C2_11 = {C11}\nC2_22 = {C22}\nC2_13 = {C13}\nC2_eps4 = {total4}\n')
