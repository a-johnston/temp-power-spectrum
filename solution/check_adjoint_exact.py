"""Exact algebra behind the report's cubic adjoint and its momentum form.

No numerical integration or field-engine imports. Fubini itself is justified
analytically in report/curvature_and_adjoint_derivations.md; these checks establish the generic
operator, diagonal-jump, source-compatibility, and endpoint identities it uses.
"""

import sympy as s


def main():
    t, r, c = s.symbols("t r c", positive=True)
    k2 = 3 * c * c

    def u1(x):
        return s.sin(x) / x - s.cos(x)

    def u2(x):
        return -s.cos(x) / x - s.sin(x)

    def H(f):
        return 6 * f + 18 / k2 * (s.diff(f, t) / t + f / t**2)

    def LN(f):
        return 3 * s.diff(f, t, 2) + k2 * f - 6 * f / t**2

    def zero(f):
        return s.simplify(s.expand_trig(f)) == 0

    G = -(u2(c * t) * u1(c * r) - u1(c * t) * u2(c * r)) / (3 * c)
    assert zero(LN(G))
    assert zero(G.subs(t, r))
    assert zero(s.diff(G, t).subs(t, r) + s.Rational(1, 3))
    # The t derivative of H acting on the retarded convolution generates this
    # diagonal term. Forgetting it would give a wrong Hamiltonian-source adjoint.
    jump = s.simplify((s.Rational(1, 2) * 18 / (k2 * t) * s.diff(G, t)).subs(t, r))
    assert zero(jump + 3 / (k2 * r))
    print("[ok] exact lapse Green equation, zero diagonal, derivative jump -1/3; contact term -3/(k^2 s)")

    h1 = -6 * (2 * s.sin(c * t) / (c * t) - s.cos(c * t) - 1)
    assert zero(H(u1(c * t)) - (6 - h1))
    assert zero(s.diff(H(u1(c * t)), t) / 2 + s.diff(h1, t) / 2)
    for u in (u1(c * t), u2(c * t)):
        assert zero(s.diff(H(u), t) / 2 - 3 * (s.diff(u, t) - u / t))
    c1 = -k2 / 9
    c3 = s.Rational(2, 9) * (-2 * (3 * k2) - k2)  # (k,q,-q), |k|=|q|
    assert s.simplify(c3 / c1) == 14
    print("[ok] homogeneous factor 14, constant 6 drops out, and both partial-transfer identities")

    nu = s.Function("nu")(t)
    cn, ch, cm = (s.Function(x)(t) for x in ("cN", "cHam", "iKdotcM"))
    h = H(nu) + 3 * ch / (2 * k2)
    nudd = (-k2 * nu + 6 * nu / t**2 - cn) / 3
    dh = s.diff(h, t).subs(s.diff(nu, t, 2), nudd)
    momentum_dh = 6 * (s.diff(nu, t) - nu / t) + 3 * cm / k2
    assert zero(dh - momentum_dh - 3 / (2 * k2) * (s.diff(ch, t) - 4 * cn / t - 2 * cm))
    # Thus cHam' - 4 cN/t = 2 i K.cM for consistent lower-order sources.
    # It cancels the contact term before any time integration is done.
    combined = -3 * cn / (k2 * t) + 3 * s.diff(ch, t) / (4 * k2)
    assert zero(combined.subs(s.diff(ch, t), 4 * cn / t + 2 * cm) - 3 * cm / (2 * k2))
    assert not zero((-3 * cn / (k2 * t) + 3 * (2 * cm) / (4 * k2)) - 3 * cm / (2 * k2))
    print("[ok] source compatibility and exact cancellation of the adjoint contact/Hamiltonian terms")

    z = s.Symbol("eps")
    d1, d2, d3, D1, D2, D3 = s.symbols("d1 d2 d3 D1 D2 D3")
    f1, f1p, f1pp, f2, f2p, f3 = s.symbols("f1 f1p f1pp f2 f2p f3")
    shift = z * d1 + z * z * d2 / 2 + z**3 * d3 / 6
    # Integral from the background endpoint to the shifted one, using local
    # arbitrary Taylor jets of F. The background integrand is identically zero.
    endpoint = (
        z * (f1 * shift + f1p * shift**2 / 2 + f1pp * shift**3 / 6)
        + z * z / 2 * (f2 * shift + f2p * shift**2 / 2)
        + z**3 * f3 * shift / 6
    )
    third = s.expand(endpoint).coeff(z, 3) * 6
    assert s.expand(third - (3 * f1 * d2 + 3 * f1p * d1**2 + 3 * f2 * d1)) == 0
    assert not third.has(d3)
    shift = z * D1 + z * z * D2 / 2 + z**3 * D3 / 6
    composed = z * (f1 + f1p * shift + f1pp * shift**2 / 2) + z * z / 2 * (f2 + f2p * shift) + z**3 * f3 / 6
    third = s.expand(composed).coeff(z, 3) * 6
    assert s.expand(third - (f3 + 3 * D2 * f1p + 3 * D1 * f2p + 3 * D1**2 * f1pp)) == 0
    assert not third.has(D3)
    print("[ok] d3 and D3 absent, while d2,d1^2,D2,D1,D1^2 remain at cubic order")
    print("EXACT ADJOINT INGREDIENTS VERIFIED")


if __name__ == "__main__":
    main()
