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
