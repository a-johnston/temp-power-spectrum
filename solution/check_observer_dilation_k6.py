"""Observer-dilation decomposition of the eps^4 k0^6 shell quadrupole coefficient (exact).

Prop. obs-dilation: (1/4)C22 + (1/3)C13 = C^{c,4} + (sigma^2/2) D^2 C11 - 2 D X on the shell, with
  C^{c,4} = centred 22 (|K^c_2|^2) + centred 13 (K_1 K^c_3),  X = (4pi/25) K_1 <K^c_2>.
At leading order this is 29/3 = 8 + 5/3 + 0 + 0.  Here every term is evaluated at the next order
from the independent hand map (check_quadrupole_next_hand), the linear transfer and eq.
stationary-four-gradient, and the sum is checked against the known complete k0^6 coefficient.
Units: k0 = 1, sigma = 1; e, o, ti are the emission, observation and initial times.
Run: .venv/bin/python -B solution/check_observer_dilation_k6.py
"""

import check_quadrupole_next_hand as h
import sympy as s
from check_observer_centred_kappa4 import J3_2, J3_4, K1_2, K1_4, K2_4, Kc4, k

o, e, ti, t, n = h.o, h.e, h.ti, h.t, h.n
pi = s.pi


def linear(P, pn):
    """Even-in-n linear quadrupole kernel at degrees 2 and 4 (centred)."""
    tau2, lam = (o**2 - e**2) / 6, o - t
    two = tau2 * h.nn(P)
    four = -h.time_integral(h.k2 * t**3 * h.nn(P) / 30 + t * lam**2 * pn**2 * h.nn(P) / 6)
    return s.expand(two - h.sky(two)), s.expand(four - h.sky(four))


def main():
    U, V = h.gradient_fields()
    leading, parts = h.hand_map(U, V)
    part = s.expand(sum(parts.values()))
    lp2, lp4 = linear(h.P, h.p.dot(n))
    lq2, lq4 = linear(h.Q, h.q.dot(n))
    # K^c_2 = K_2 + D K_1(p) + D K_1(q); D = 2 on degree 2 and 4 on degree 4.
    lead_c = s.expand(leading + 2 * (lp2 + lq2))
    part_c = s.expand(part + 4 * (lp4 + lq4))
    sub = {k: 1}
    c22_4 = s.factor(pi * h.shell_integral(h.sky(lead_c**2)) / 5).subs(sub)
    c22_6 = s.factor(2 * pi * h.shell_integral(h.sky(lead_c * part_c)) / 5).subs(sub)
    C11_4 = 4 * pi * (o**2 - e**2) ** 2 / 2025
    assert s.simplify(c22_4 - s.Rational(5, 3) * C11_4) == 0, "leading centred 22 must be 5/3 C11"
    J12, J14 = K1_2.subs(sub), K1_4.subs(sub)
    norm = 4 * pi / 25
    centred13 = norm * J12 * Kc4.subs(sub)
    dilation = 36 * norm * J12 * J14  # (1/2) D^2 C11 with C11 = norm K1^2, D^2 = 36 on k0^6
    Kc2_4 = K2_4.subs(sub) + 4 * J14  # <K^c_2(k,q)> at k^4: <K2> + D K1(k), <D K1(q)> = 0
    cross = -12 * norm * J12 * Kc2_4  # -2 D X with D = 6
    total = s.factor(c22_6 + centred13 + dilation + cross)
    # Known complete coefficient: (1/4)C22 (hand map) + (1/3)C13 (eq. stationary-four-gradient).
    tot22 = -8 * pi * (o**2 - e**2) * (105 * ti**4 + 719 * o**4 - 1920 * o**2 * e**2 + 2560 * o * e**3 - 1464 * e**4) / 637875
    tot13 = norm * (J12 * J3_4.subs(sub) + J14 * J3_2.subs(sub))
    assert s.simplify(total - (tot22 + tot13)) == 0, "decomposition must reproduce the known k0^6 coefficient"
    # Why centred 22 is proportional at k0^6: K^c_2^[2] = -deltaK1^[2] (polarisation of the linear
    # kernel), and the overlap of K^c_2^[4] with it equals that of -2 deltaK1^[4] = -(D/2) deltaK1^[4].
    Kv = h.p + h.q
    PK = Kv * Kv.T - h.identity * Kv.dot(Kv) / 3
    tau2, lam = (o**2 - e**2) / 6, o - t
    lK2 = tau2 * h.nn(PK)
    lK4 = -h.time_integral(Kv.dot(Kv) * t**3 * h.nn(PK) / 30 + t * lam**2 * Kv.dot(n) ** 2 * h.nn(PK) / 6)
    lK2, lK4 = s.expand(lK2 - h.sky(lK2)), s.expand(lK4 - h.sky(lK4))
    pol2, pol4 = s.expand(-(lK2 - lp2 - lq2)), s.expand(-(lK4 - lp4 - lq4))
    assert s.expand(lead_c - pol2) == 0, "leading centred kernel is minus the linear polarisation"

    def overlap(A, B):
        return s.expand(2 * pi * h.shell_integral(h.sky(A * B)) / 5).subs(sub)

    assert s.simplify(c22_6 - 2 * overlap(pol2, pol4)) == 0, "centred 22 at k0^6 = twice the polarisation overlap"
    print("[ok] k0^6 centred 22 = 2 <deltaK1^[2] deltaK1^[4]>: fixed by the linear transfer alone")
    # Mechanism: observer-local pieces and the dilation subtraction are orthogonal to the leading
    # structure; the bulk field and the two first-order ray pieces sum to the factor 2.
    ref = overlap(pol2, pol4)
    for name in ("endpoint", "clock", "normalization"):
        assert overlap(pol2, parts[name]) == 0, name
    assert overlap(pol2, s.expand(4 * (lp4 + lq4))) == 0
    bulk = sum(overlap(pol2, parts[x]) for x in ("field", "position", "direction"))
    assert s.simplify(bulk - 2 * ref) == 0
    assert all(s.simplify(overlap(pol2, parts[x]) / ref).has(e) for x in ("field", "position", "direction"))
    print("[ok] mechanism: endpoint/clock/normalisation/dilation are orthogonal; field + position + direction = 2")
    C11_6 = 2 * norm * J12 * J14
    assert s.simplify(c22_6 - s.Rational(28, 3) * C11_6) == 0, "centred 22 at k0^6 is 28/3 C11[k0^6]"
    assert s.simplify(dilation - 18 * C11_6) == 0
    print("[ok] centred 22 = (28/3) C11 and dilation = 18 C11 at k0^6 (5/3 and 8 at k0^4)")
    for name, val in (
        ("centred 22", c22_6),
        ("centred 13", centred13),
        ("dilation (sigma^2/2)D^2C11", dilation),
        ("cross -2DX", cross),
    ):
        print(f"{name:28s} {s.factor(val)}")
    print("total", total)
    print("[ok] k0^6 decomposition reproduces the known complete coefficient; leading centred 22 = 5/3 C11")


if __name__ == "__main__":
    main()
