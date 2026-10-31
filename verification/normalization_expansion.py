"""Exact expansion of the normalised map through O(eps^3).

T_obs(n) = T_fluid(o) * exp(-I(n)),  I = eps I1 + eps^2 I2/2 + eps^3 I3/6  (line-of-sight exponent,
including every endpoint/observer/tetrad term).  The direction-independent prefactor T_fluid(o)
cancels exactly in Theta = T_obs / <T_obs>_sky - 1.  Sky averages <.> are linear functionals
(NOT ensemble averages): they are realisation-dependent random numbers.  We model <.> by an
abstract linear operator acting on monomials of I1, I2, I3 and verify  <Theta^(r)> = 0.
"""
import sympy as sp

eps = sp.Symbol('epsilon')
I1, I2, I3 = sp.symbols('I1 I2 I3')
# sky-average symbols for every monomial that can appear through third order
mon = {}
def avg(expr):
    expr = sp.expand(expr)
    out = 0
    for term in sp.Add.make_args(expr):
        c, m = term.as_independent(I1, I2, I3)
        if m == 1:
            out += c
        else:
            if m not in mon:
                mon[m] = sp.Symbol('<' + str(m).replace('**', '^').replace('*', ' ') + '>')
            out += c * mon[m]
    return out

I = eps * I1 + eps**2 * I2 / 2 + eps**3 * I3 / 6
R = sp.series(sp.exp(-I), eps, 0, 4).removeO()
Theta = sp.series(R / avg(R), eps, 0, 4).removeO() - 1
coef = {r: sp.expand(sp.factorial(r) * Theta.coeff(eps, r)) for r in (1, 2, 3)}
for r in (1, 2, 3):
    print(f"Theta^({r}) =", coef[r])
# consistency: sky average of each order vanishes identically (exact normalisation)
for r in (1, 2, 3):
    a = avg(coef[r].subs({v: 0 for v in []}))
    # <c * <m>> = c <m> for random-but-direction-independent sky averages
    expr = 0
    for term in sp.Add.make_args(coef[r]):
        c, m = term.as_independent(I1, I2, I3)
        expr += c * (avg(m) if m != 1 else 1)
    assert sp.simplify(expr) == 0, (r, expr)
    print(f"[ok] sky average of Theta^({r}) vanishes identically")

# ---------------------------------------------------------------------------------------------
# Non-tautological check (audit W3): a discrete "sky" of M directions with exact rational
# weights w_a (sum 1) and exact rational values I1_a, I2_a, I3_a.  We compute the exact
# Taylor coefficients of  exp(-I_a) / sum_b w_b exp(-I_b) - 1  and compare them with the
# printed formulas Theta^(r), in which <.> is replaced by the weighted discrete average.
# A deliberately mutated formula must FAIL (negative control).
import random

def discrete_check(formulas, seed, M=5):
    rnd = random.Random(seed)
    wts = [sp.Rational(rnd.randint(1, 9)) for _ in range(M)]
    tot = sum(wts)
    wts = [w_ / tot for w_ in wts]
    vals = [{I1: sp.Rational(rnd.randint(-9, 9), rnd.randint(1, 5)),
             I2: sp.Rational(rnd.randint(-9, 9), rnd.randint(1, 5)),
             I3: sp.Rational(rnd.randint(-9, 9), rnd.randint(1, 5))} for _ in range(M)]
    Rs = [sp.exp(-(eps * v[I1] + eps**2 * v[I2] / 2 + eps**3 * v[I3] / 6)) for v in vals]
    avgR = sum(w_ * R_ for w_, R_ in zip(wts, Rs))
    def davg(expr):  # discrete average of a monomial expression in I1, I2, I3
        return sum(w_ * expr.subs(v) for w_, v in zip(wts, vals))
    ok = True
    for a in range(M):
        exact = sp.series(Rs[a] / avgR - 1, eps, 0, 4).removeO()
        for r in (1, 2, 3):
            want = sp.factorial(r) * exact.coeff(eps, r)
            f = formulas[r]
            # replace sky-average symbols by discrete averages, then I's by direction values
            sub = {sym: davg(m) for m, sym in mon.items()}
            got = f.xreplace(sub).subs(vals[a])
            if sp.nsimplify(sp.expand(got - want)) != 0:
                ok = False
    return ok

for seed in (11, 12, 13):
    assert discrete_check(coef, seed), seed
print("[ok] printed Theta^(1,2,3) reproduce the exact expansion of T/<T>-1 on random exact discrete skies")
mut = dict(coef)
mut[3] = coef[3] + 3 * I1 * I2 - 3 * I1 * I2 * sp.Rational(1, 2)  # coefficient 3 I1 I2 -> 3/2 I1 I2
assert not discrete_check(mut, 11)
mut2 = dict(coef)
mut2[2] = coef[2] + 2 * mon[I1] * I1  # drop the normalisation interference term -2<I1> I1
assert not discrete_check(mut2, 11)
print("[ok] negative controls: mutated Theta^(3) and Theta^(2) (dropped <I1> I1 term) are rejected")
