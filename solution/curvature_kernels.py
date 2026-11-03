"""Fast local-curvature map through four gradients, at general (k,q,-q).

No field or ray engine is imported. --compare checks the full normalized map
against the saved independent generator, BEFORE projection or Wick averaging.
The dipole is retained. Default operation is read-only.
"""

import argparse
import pickle
import time
from pathlib import Path

import check_cubic_curvature_hand as c
import sympy as s


def kernels():
    # General momenta, not a monochromatic specialization. The previous hand
    # check uses a unit shell; here no square-root quotient relation is imposed.
    p, b1, b2, kap = s.symbols("p b1 b2 kap", real=True)
    saved = c.SHELL_REDUCTION, c.momenta, c.K
    try:
        c.SHELL_REDUCTION = False
        c.momenta = ((p, 0, 0), (b1, b2, 0), (-b1, -b2, 0))
        c.K = [tuple(sum(c.momenta[a][j] for a in range(3) if mask & (1 << a)) for j in range(3)) for mask in range(8)]
        leading, odd, pieces = c.map_pieces(c.geometry())
        fourth = sum((time * value for time, value in pieces.values()), c.J())
        out = {}
        for mask in (1, 3, 7):
            out[mask] = s.expand(
                sum(
                    kap**degree * (part.coeff(mask) - c.sky(part.coeff(mask)))
                    for degree, part in ((2, leading), (3, odd), (4, fourth))
                )
            )
        return out, (p, b1, b2, kap)
    finally:
        c.SHELL_REDUCTION, c.momenta, c.K = saved


def compare(out, symbols):
    """Clear denominators in a polynomial ring; do not trust a simplified cache."""
    from sympy.polys.domains import QQ_I
    from sympy.polys.rings import ring

    root = Path(__file__).resolve().parent / "results"
    for name, masks in (("2", ((1, (1,)), (3, (1, 2)))), ("3", ((1, (1,)), (7, (1, 2, 3))))):
        # These are trusted local repository artifacts, not downloaded pickles.
        with (root / f"superhorizon_kernels_{name}.pkl").open("rb") as stream:
            reference = pickle.load(stream)
        inverse = ("ip", "u", "v", "its", "itI", "itob")
        names = [g for g in reference["gens"] if g not in inverse]
        R, *gens = ring(names, QQ_I)
        G = dict(zip(names, gens, strict=True))
        inv_indices = [reference["gens"].index(g) for g in inverse]
        reg_indices = [reference["gens"].index(g) for g in names]
        denoms = (
            G["p"],
            (G["p"] + G["b1"]) ** 2 + G["b2"] ** 2,
            (G["p"] - G["b1"]) ** 2 + G["b2"] ** 2,
            G["ts"],
            G["tI"],
            G["tob"],
        )
        p, b1, b2, kap = symbols
        sub = dict(zip((p, b1, b2, kap, *c.n, c.o, c.e, c.ti), s.symbols("p b1 b2 kap n1 n2 n3 tob ts tI"), strict=True))
        for mask, labels in masks:
            terms = reference["kernels"][labels]
            maxima = tuple(max(mon[j] for mon in terms) for j in inv_indices)
            groups = {}
            for mon, coefficient in terms.items():
                powers = tuple(mon[j] for j in inv_indices)
                group = groups.setdefault(powers, {})
                group[tuple(mon[j] for j in reg_indices)] = QQ_I.from_sympy(s.sympify(coefficient))
            total = R.zero
            common = R.one
            for denom, power in zip(denoms, maxima, strict=True):
                common *= denom**power
            for powers, group in groups.items():
                numerator = R.from_dict(group)
                for denom, maximum, power in zip(denoms, maxima, powers, strict=True):
                    numerator *= denom ** (maximum - power)
                total += numerator
            expected = R.from_expr(out[mask].xreplace(sub)) * common
            assert total == expected, labels
            print(f"[ok] RAW kernel {labels}: all sky components, general momenta and times, degrees 2--4", flush=True)


def spectrum(out, symbols):
    """Exact quadrupole coefficients, including general spectral moments."""
    p, b1, b2, kap = symbols
    q, mu, v = s.symbols("q mu v", real=True)
    m = {j: s.Symbol(f"m{j}", real=True) for j in (0, 2, 4, 6)}

    def radial_average(expr, pair):
        expr = s.expand(expr.subs({b1: q * mu, b2: q * v}))
        expr = s.Poly(expr, v).rem(s.Poly(v**2 + mu**2 - 1, v)).as_expr().expand()
        assert not expr.has(v)
        expr = s.integrate(expr, (mu, -1, 1)) / 2
        gens = (p, q) if pair else (p,)
        return s.factor(sum(coef * s.prod(m[j] for j in powers) for powers, coef in s.Poly(expr, *gens).terms()))

    a, b = (out[1].coeff(kap, j) for j in (2, 4))
    u, w = (out[3].coeff(kap, j) for j in (2, 4))
    v2, v4 = (out[7].coeff(kap, j) for j in (2, 4))
    # Every degree-two map is pure l=2. Orthogonality lets us avoid projecting
    # its degree-four partner; degree-three products have only odd multipoles.
    entries = {
        "C11 degree4": (4 * s.pi / 5 * c.sky(a * a), False),
        "C11 degree6": (8 * s.pi / 5 * c.sky(a * b), False),
        "C22/4 degree4": (2 * s.pi / 5 * c.sky(u * u), True),
        "C22/4 degree6": (4 * s.pi / 5 * c.sky(u * w), True),
        "C13/3 degree4": (4 * s.pi / 5 * c.sky(a * v2), True),
        "C13/3 degree6": (4 * s.pi / 5 * c.sky(a * v4 + b * v2), True),
    }
    results = {name: radial_average(expr, pair) for name, (expr, pair) in entries.items()}
    Delta = c.o**2 - c.e**2
    L = c.o**4 / 12 - c.o**2 * c.e**2 / 2 + 2 * c.o * c.e**3 / 3 - c.e**4 / 4
    E = c.o**4 - c.e**4
    betaH = -11 * L / 14 - 3 * E / 40
    betaa = 27 * L / 7 + 7 * E / 20 + (c.o**4 + 3 * c.ti**4) / 4 - c.e**4
    compact = {
        "C11 degree4": 4 * s.pi / 2025 * Delta**2 * m[4],
        "C22/4 degree4": 4 * s.pi / 6075 * Delta**2 * (12 * m[0] * m[4] + 5 * m[2] ** 2),
        "C13/3 degree4": 16 * s.pi / 2025 * Delta**2 * m[0] * m[4],
        "C11 degree6": 16 * s.pi / 6075 * Delta * betaH * m[6],
        "C22/4 degree6": 32 * s.pi / 18225 * Delta * (12 * betaH * m[0] * m[6] + (18 * betaH - betaa) * m[2] * m[4]),
        "C13/3 degree6": 8 * s.pi / 18225 * Delta * (60 * betaH * m[0] * m[6] + (50 * betaH - 9 * betaa) * m[2] * m[4]),
    }
    for name, value in compact.items():
        assert s.expand(results[name] - value) == 0, name

    # Independent leading octupole, retaining then subtracting the measured
    # dipole of the odd map. Do not confuse a raw cubic polynomial with its STF part.
    def octupole(expr):
        return s.expand(expr - 3 * sum(ni * c.sky(ni * expr) for ni in c.n))

    oct1, oct2, oct3 = (octupole(out[mask].coeff(kap, 3)) for mask in (1, 3, 7))
    beta3 = (c.o - c.e) ** 2 * (c.o + 2 * c.e) / 18
    octnorm = 16 * s.pi / 1225 * beta3**2
    octentries = {
        "octupole C11 degree6": (4 * s.pi / 7 * c.sky(oct1 * s.conjugate(oct1)), False, octnorm * m[6]),
        "octupole C22/4 degree6": (
            2 * s.pi / 7 * c.sky(oct2 * s.conjugate(oct2)),
            True,
            octnorm * (9 * m[0] * m[6] + 28 * m[2] * m[4]),
        ),
        "octupole C13/3 degree6": (4 * s.pi / 7 * c.sky(oct1 * s.conjugate(oct3)), True, 9 * octnorm * m[0] * m[6]),
    }
    for name, (expr, pair, expected) in octentries.items():
        results[name] = radial_average(expr, pair)
        assert s.expand(results[name] - expected) == 0, name
    print("[ok] leading octupole: independent 22/13 contributions and general-moment factor", flush=True)
    total = results["C22/4 degree6"] + results["C13/3 degree6"]
    dilation = 18 * m[0] * results["C11 degree6"]
    remainder = 8 * s.pi / 18225 * Delta * (122 * betaH - 13 * betaa) * m[2] * m[4]
    assert s.expand(total - dilation - remainder) == 0
    print("[ok] general-moment betaH/betaa forms; subtracting observer dilation leaves only m2*m4", flush=True)
    for name, value in results.items():
        print(name, "=", value, flush=True)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--spectrum", action="store_true", help="derive general-moment quadrupole through gradient degree six")
    parser.add_argument("--output", type=Path, help="write exact expressions to a NEW text artifact")
    args = parser.parse_args()
    start = time.monotonic()
    out, symbols = kernels()
    if args.compare:
        compare(out, symbols)
    if args.spectrum:
        spectrum(out, symbols)
    if args.output:
        with args.output.open("x") as stream:
            for mask, value in out.items():
                stream.write(f"mask={mask}\n{value}\n\n")
    for mask, value in out.items():
        print(f"kernel mask {mask}: {len(s.Add.make_args(value))} polynomial terms")
    print(f"Completed exact curvature construction in {time.monotonic() - start:.2f} s")


if __name__ == "__main__":
    main()
