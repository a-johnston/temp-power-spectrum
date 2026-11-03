"""Exact leading observer-centred cubic shell quadrupole (report Sec. stationary-evaluated).

K^c_3 = K_3 + sum_pairs D K_2 + sum_a D^2 K_1 (exact; Sec. obs-dilation). On the shell, with the
quadrupole weight aligned with k, the q average of the (q,-q) and single-q terms vanishes, so
<K^c_3> = Jbar_3 + 2 D<K_2(k,q)> + D^2 K_1(k), with D = n on a degree-n gradient coefficient.
Inputs: Jbar_3 (eq. stationary-four-gradient), K_1 through k^4 (linear transfer), and the independent
hand map K_2 through four gradients (check_quadrupole_next_hand.hand_map).
Run: .venv/bin/python -B solution/check_observer_centred_kappa4.py
"""

import check_quadrupole_next_hand as h
import sympy as s
from quadrupole_adjoint_expansion import Q, angular_A
from unified_adjoint_source import inner

b, c, o, e, ti, n = h.b, h.c, h.o, h.e, h.ti, h.n
U, V = h.gradient_fields()
leading, parts = h.hand_map(U, V)
k2 = b**2 + c**2
p = s.Matrix([b, 0, c])
w = s.Rational(15, 2) * (n.dot(p)) ** 2 / k2 - s.Rational(5, 2)  # 5 P2(n.phat)


def shellavg(expr):  # <.>_qhat on the shell, as a function of k (=|p|)
    num = s.expand(h.sky(s.expand(w * k2 * expr)))  # multiply by k2 to clear the 1/k2 of w
    k = s.Symbol("k", real=True)
    return s.expand(h.shell_integral(num) / 2 / k**2)


K2_2 = shellavg(leading)
K2_4 = shellavg(sum(parts.values()))
k = s.Symbol("k", real=True)
tau2 = (o**2 - e**2) / 6
A = o**4 / 12 - o**2 * e**2 / 2 + 2 * o * e**3 / 3 - e**4 / 4
D = o**4 - e**4
C = Q(s.Matrix([1, 0, 0]))
kk = s.Matrix([1, 0, 0])
K1_2 = s.Rational(2, 3) * tau2 * k**2
K1_4 = (-A * inner(angular_A(C, kk, kk), Q(kk)) / 6 - D * inner(C, Q(kk)) / 120) * k**4
J3_2 = 4 * (o**2 - e**2) / 9 * k**2
J3_4 = (3984 * e**4 - 6256 * e**3 * o + 4692 * e**2 * o**2 - 1853 * o**4 - 567 * ti**4) / 3402 * k**4
lead = s.factor(J3_2 + 2 * 2 * K2_2 + 4 * K1_2)
assert lead == 0, "leading observer-centred cubic kernel must vanish"
Kc4 = s.factor(J3_4 + 2 * 4 * K2_4 + 16 * K1_4)
expected = k**4 * (4 * (e - o) * (12 * e**3 - 8 * e**2 * o + 7 * e * o**2 + 7 * o**3) + 5 * (o**4 - ti**4)) / 270
assert s.expand(Kc4 - expected) == 0
print("[ok] <K^c_3(k,q,-q)> = 0 at kappa^2 and =", s.factor(expected), "at kappa^4")
