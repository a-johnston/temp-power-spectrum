"""Reconstruct the first-order transfer function from its superhorizon series alone (exact).

The kernels are entire in the momenta at fixed times, so the complete small-k0 series determines the
answer at every k0. If, in addition, the function class is known and finite, matching the series is
linear algebra with a unique solution.  This script demonstrates that for T_l:

  input   the small-k0 Taylor series of T_l at fixed times, computed directly from the defining
          line-of-sight integral  T_l = -3 int_{x*}^{x_o} s'(x) w_l(sqrt3 (x_o - x)) dx  (no Prop. Tl);
  ansatz  sum over  phase(x_o) * x_o^j * J_m(omega, y)  with
            y in {x_o, x*, x_o - x*},  omega in {+-(1 - sqrt3), +-(1 + sqrt3)},  m = 1..l+3,
            phases e^{i p x_o} with p in {0, +-1, +-sqrt3},  j in [JMIN, 1],
          and UNKNOWN complex coefficients (a superset of the class of Prop. Tl);
  solve   match the Laurent coefficients in k0 through order N, exactly over Q(sqrt3, i);
  check   (a) the solution predicts the next orders that were not used;
          (b) it equals T_l at finite k0 to 30 digits (mpmath quadrature of the definition, as a
              check of the identity, not as evidence), including deep inside the horizon.

Units: t_o = sqrt3 (so x_o = k0) and t* = sqrt3 rho, i.e. x* = rho x_o, with rho symbolic.

    python reconstruct_transfer.py [l] [N]
"""

import sys
import time
from itertools import product

import mpmath as mp
import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.rings import ring

L = int(sys.argv[1]) if len(sys.argv) > 1 else 2
NORD = int(sys.argv[2]) if len(sys.argv) > 2 else 12
JMIN = -(L + 4)  # Laurent range in x_o: bounded by the pole orders (2 at x = 0, l + 3 at x = x_o)
MD = L + 3  # J order at the x = x_o pole
MO = 2  # J order at the x = 0 pole
SHIFT = -JMIN + MD - 1  # multiply by k^SHIFT to clear the Laurent poles
T0 = time.time()

r3 = sp.sqrt(3)
K = sp.QQ.algebraic_field(r3, sp.I)
R, kk, rr = ring("k rho", K)
k, rho = sp.symbols("k rho")
FIT = SHIFT + NORD  # shifted k-degree used in the fit
DEG = FIT + 8  # shifted k-degree checked; FIT+1..DEG are predictions
c = lambda e: K.from_sympy(sp.nsimplify(e))
IU = c(sp.I)


def tr(P, deg=DEG):
    return R({m: v for m, v in P.items() if m[0] <= deg})


def mul(P, Q):
    return tr(P * Q)


def exp_series(p):
    out, term = R.one, R.one
    for n in range(1, DEG + 2):
        term = tr(term * kk * IU * c(p) * K.from_sympy(sp.Rational(1, n)))
        out += term
    return out


def J_poly(m, om, y):
    """J_m(om, y) = y^(1-m) * P(y),  P(y) = sum_{j>=m} (i om)^j y^j / (j! (j-m+1))."""
    out, ypow = R.zero, R.one
    iom = IU * c(om)
    for j in range(DEG + m + 2):
        if j >= m:
            out += ypow * (iom**j * K.from_sympy(sp.Rational(1, sp.factorial(j) * (j - m + 1))))
        ypow = tr(ypow * y)
    return out


def target():
    """k^SHIFT * T_l from the defining integral, expanded in k (exact, rational)."""
    x, y, u = sp.symbols("x y u")
    jl = sp.expand_func(sp.jn(L, u))
    w = sp.series(jl / 3 + sp.diff(jl, u, 2), u, 0, DEG + 2).removeO()
    sp_ = sp.series(sp.diff(2 * sp.sin(x) / x - sp.cos(x), x), x, 0, DEG + 2).removeO()
    integrand = sp.expand(-3 * k * sp_.subs(x, k * y) * w.subs(u, r3 * k * (1 - y)))
    integrand = sp.Add(*[t for t in sp.Add.make_args(integrand) if sp.degree(t, k) <= DEG - SHIFT])
    T = sp.expand(sp.integrate(integrand, (y, rho, 1)) * k**SHIFT)
    return R({m: K.from_sympy(v) for m, v in sp.Poly(T, k, rho).terms()})


def basis():
    """Function class read off the integrand, not from Prop. Tl.

    s'(x) = sum_a e^{i a x} S_a(x) and w_l(sqrt3 (x_o - x)) = sum_b e^{i b (x_o - x)} V_b(x_o - x), with
    a = +-1, b = +-sqrt3, S_a Laurent with a double pole at 0, V_b Laurent with a pole of order l+3.
    Partial fractions in x give phase e^{i b x_o} times J_m(a - b, .) at x_o and x* (m <= 2), and phase
    e^{i a x_o} times J_m(b - a, x_o - x*) (m <= l + 3).  Laurent coefficients in x_o: all unknown.
    """
    out = []
    for a, b in product((1, -1), (r3, -r3)):
        for an in ("o", "s"):
            for m, j in product(range(1, MO + 1), range(JMIN, 1)):
                out.append((b, a - b, an, m, j))
        for m, j in product(range(1, MD + 1), range(JMIN, 1)):
            out.append((a, b - a, "d", m, j))
    return out


# J_m(om, rho k) carries rho^(1-m) and J_m(om, (1-rho) k) carries (1-rho)^(1-m): clear with COMMON
COMMON = (rr ** (MO - 1)) * ((1 - rr) ** (MD - 1))


_CORE = {}


def column(b):
    p, om, an, m, j = b
    key = (p, om, an, m)
    if key not in _CORE:
        y = {"o": kk, "s": rr * kk, "d": (1 - rr) * kk}[an]
        core = mul(exp_series(p), J_poly(m, om, y))
        # COMMON / yf^(m-1)
        if an == "o":
            fac = COMMON
        elif an == "s":
            fac = rr ** (MO - m) * (1 - rr) ** (MD - 1)
        else:
            fac = rr ** (MO - 1) * (1 - rr) ** (MD - m)
        _CORE[key] = tr(core * fac)
    return tr(_CORE[key] * kk ** (SHIFT + j + 1 - m))


def main():
    tgt = tr(target() * COMMON)
    print(f"[{time.time() - T0:6.1f}s] target series ready", flush=True)
    B = basis()
    cols = [column(b) for b in B]
    print(f"[{time.time() - T0:6.1f}s] {len(B)} basis columns ready", flush=True)
    keys = sorted({mon for P in cols + [tgt] for mon in P if mon[0] <= FIT})
    idx = {mon: i for i, mon in enumerate(keys)}
    M = [[K.zero] * (len(B) + 1) for _ in keys]
    for jc, P in enumerate(cols):
        for mon, v in P.items():
            if mon in idx:
                M[idx[mon]][jc] = v
    for mon, v in tgt.items():
        if mon in idx:
            M[idx[mon]][len(B)] = v
    DM = DomainMatrix(M, (len(keys), len(B) + 1), K)
    print(f"[{time.time() - T0:6.1f}s] system {len(keys)} x {len(B)}; row reduction ...", flush=True)
    Rm, piv = DM.rref()
    assert len(B) not in piv, "inconsistent: the ansatz cannot reproduce the series"
    sol = [K.zero] * len(B)
    for r, pc in enumerate(piv):
        sol[pc] = Rm.rep.to_ddm()[r][len(B)]  # free variables set to zero
    print(
        f"[{time.time() - T0:6.1f}s] solved: rank {len(piv)} of {len(B)}; {sum(1 for v in sol if v != K.zero)} nonzero coefficients",
        flush=True,
    )
    pred = R.zero
    for v, P in zip(sol, cols, strict=True):
        if v != K.zero:
            pred += P * v
    resid = tr(pred - tgt)
    bad = sorted({mon[0] for mon, v in resid.items() if v != K.zero})
    assert not [d for d in bad if d <= FIT], f"fit orders not reproduced: {bad}"
    assert not bad, f"predictions beyond the fit fail at shifted k-degrees {bad}"
    print(f"[ok] fitted k0 orders <= {NORD}; orders {NORD + 1}..{DEG - SHIFT} predicted exactly (not used in the fit)")
    return [K.to_sympy(v) for v in sol], B


def evaluate(sol, B, x_o, x_s):
    """The reconstructed function at finite k0, from its exact coefficients (high precision throughout:
    the terms cancel strongly, and J_m is summed from its entire power series)."""
    x_o, x_s = mp.mpf(x_o), mp.mpf(x_s)
    with mp.workdps(60 + int(2 * float(x_o))):
        num = lambda e: mp.mpc(mp.mpf(str(sp.N(sp.re(e), mp.mp.dps))), mp.mpf(str(sp.N(sp.im(e), mp.mp.dps))))

        def J(m, om, y):
            # J_m = sum_{j>=m} (i om)^j y^(j-m+1) / (j! (j-m+1))
            tot, j, term = mp.mpc(0), m, (1j * om) ** m / mp.factorial(m)
            while True:
                add = term * y ** (j - m + 1) / (j - m + 1)
                tot += add
                if j > m + 10 and abs(add) < mp.mpf(10) ** (-mp.mp.dps):
                    return tot
                j += 1
                term *= 1j * om / j

        tot = mp.mpc(0)
        for cf, (p, om, an, m, j) in zip(sol, B, strict=True):
            if cf == 0:
                continue
            y = {"o": x_o, "s": x_s, "d": x_o - x_s}[an]
            tot += num(cf) * mp.expj(num(p).real * x_o) * x_o**j * J(m, num(om).real, y)
        return +tot


def direct(x_o, x_s):
    mp.mp.dps = 30

    def w(yy):
        if yy == 0:
            return mp.mpf(2) / 15 if L == 2 else mp.mpf(0)
        j = lambda z: mp.sqrt(mp.pi / (2 * z)) * mp.besselj(L + 0.5, z)
        return j(yy) / 3 + mp.diff(j, yy, 2)

    sprime = lambda x: mp.diff(lambda z: 2 * mp.sin(z) / z - mp.cos(z), x)
    return -3 * mp.quad(lambda x: sprime(x) * w(mp.sqrt(3) * (x_o - x)), [x_s, x_o])


if __name__ == "__main__":
    import os
    import pickle

    cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", f"transfer_reconstruction_l{L}_N{NORD}.pkl")
    if os.path.exists(cache):
        with open(cache, "rb") as fh:
            sol = [sp.sympify(v) for v in pickle.load(fh)]
        B = basis()
        print(f"(exact solution loaded from {os.path.basename(cache)}; delete it to re-solve)")
    else:
        sol, B = main()
        with open(cache, "wb") as fh:
            pickle.dump([str(v) for v in sol], fh)
    for x_o, x_s in ((0.7, 0.2), (3.0, 1.0), (12.0, 4.0), (40.0, 10.0)):
        a, b = evaluate(sol, B, x_o, x_s), direct(x_o, x_s)
        print(
            f"  x_o={x_o:5.1f} x*={x_s:4.1f}: reconstructed {mp.nstr(a.real, 12)} (imag {mp.nstr(a.imag, 2)})   direct {mp.nstr(b, 12)}"
        )
