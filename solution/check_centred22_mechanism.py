"""Mechanism of the exact k0^6 centred-22 ratio (exact).

Projected on the leading structure deltaK1^[2], the k0^4 centred quadratic kernel splits by time
structure into a local part (time function D = o^4 - e^4) and a ray part (A = int_e^o t(o-t)^2 dt).
Each equals exactly twice the matching part of the polarised linear kernel deltaK1^[4]:
local four-gradient shear <-> local t^3 transfer correction; field phase + ray position + direction
<-> free-streaming phase. The dilation-subtraction terms are orthogonal to both.
Run: .venv/bin/python -B solution/check_centred22_mechanism.py
"""

import check_quadrupole_next_hand as h
import sympy as s
from check_observer_dilation_k6 import linear

o, e, t, n = h.o, h.e, h.t, h.n
U, V = h.gradient_fields()
P, Q, p, q = h.P, h.Q, h.p, h.q
K = p + q
lam = o - t
sig2 = t * U - 2 * t * (P + Q) / 3
sig24 = 2 * t**3 * V + t**3 * (P * Q + Q * P) / 9 + h.k2 * t**3 * (P + Q) / 15
cen = lambda x: s.expand(x - h.sky(x))
local = cen(-h.time_integral(h.nn(sig24)))
phase = cen(-h.time_integral(-(K.dot(n) ** 2) * lam**2 * h.nn(sig2) / 2))
leading, parts = h.hand_map(U, V)
assert s.expand(local + phase - parts["field"]) == 0
lp2, _ = linear(P, p.dot(n))
lq2, _ = linear(Q, q.dot(n))
Kv = K
PK = Kv * Kv.T - h.identity * Kv.dot(Kv) / 3
tau2 = (o**2 - e**2) / 6
pol2 = s.expand(-(cen(tau2 * h.nn(PK)) - lp2 - lq2))
ov = lambda B: s.factor(2 * s.pi * h.shell_integral(h.sky(pol2 * B)) / 5)
k = s.Symbol("k", real=True)


def lin_split(Pm, kk2, kn):
    loc = cen(-h.time_integral(kk2 * t**3 * h.nn(Pm) / 30))
    ph = cen(-h.time_integral(t * lam**2 * kn**2 * h.nn(Pm) / 6))
    return loc, ph


lKl, lKp = lin_split(PK, s.expand(K.dot(K)), K.dot(n))
lpl, lpp = lin_split(P, h.k2, p.dot(n))
lql, lqp = lin_split(Q, h.k2, q.dot(n))
pol_loc = s.expand(-(lKl - lpl - lql))
pol_ph = s.expand(-(lKp - lpp - lqp))
ray = s.expand(phase + parts["position"] + parts["direction"])
assert s.simplify(ov(local) / ov(pol_loc)) == 2, "local four-gradient shear = 2 x polarised local transfer"
assert s.simplify(ov(ray) / ov(pol_ph)) == 2, "phase + position + direction = 2 x polarised free-streaming phase"
assert ov(s.expand(4 * (lpl + lql))) == 0 and ov(s.expand(4 * (lpp + lqp))) == 0
print("[ok] k0^6 centred-22 relation = two separate factor-2 identities (local field; ray kinematics)")
