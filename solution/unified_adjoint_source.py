"""First field-adjoint source: exact off-shell field and measurement coefficients.

A fixed first direction has (zeta, nu, H)=(z,a,h), momentum p. An independent
STATE probe has (0,m,j), momentum q: its primordial input is held fixed.
No field equations, shell restriction, time expansion or momentum division is
used in these coefficients. Primed field amplitudes are partial time derivatives;
ray responses and adjoints supplied to measurement_source are pulled back.

The first-adjoint weak source is measurement_source MINUS field_variation.
The global sky/clock integrals and regular Green transpose are specified in
report/first_field_adjoint_source.md. In particular, clock terms are NOT to be
integrated over sky a second time. This module does not evaluate a spectrum.
"""

from collections.abc import Callable
from dataclasses import dataclass

import sympy as s


def inner(A, B):
    # sympy.trace(Matrix(...)) calls Matrix.doit(deep=True), inadvertently
    # attempting every definite Green integral stored in the coefficients.
    return sum(A[i, j] * B[i, j] for i in range(A.rows) for j in range(A.cols))


def stf(A):
    return (A + A.T) / 2 - s.eye(3) * sum(A[i, i] for i in range(3)) / 3


def sym_outer(a, b):
    return (a * b.T + b * a.T) / 2


def christoffel(k, G):
    """Linear connection, returned as three matrices Gamma^l_ij."""
    return [s.Matrix(3, 3, lambda i, j, l=l: s.I * (k[i] * G[l, j] + k[j] * G[l, i] - k[l] * G[i, j]) / 2) for l in range(3)]


def ricci_linear(k, G):
    return (k.dot(k) * G + sum(G[i, i] for i in range(3)) * k * k.T - k * (G * k).T - (G * k) * k.T) / 2


def ricci_mixed(p, q, G, j):
    """Mixed lower Ricci in the symmetric orthonormal frame of exp(alpha G+beta j).

    This is an off-shell matrix-exponential variation, including the mixed
    metric coefficient {G,j}/2 and both orthonormal-frame conversion terms.
    """
    P = p + q
    ga, gb = christoffel(p, G), christoffel(q, j)
    metric_ab = (G * j + j * G) / 2
    direct = christoffel(P, metric_ab)
    gab = [direct[l] - sum((G[l, r] * gb[r] + j[l, r] * ga[r] for r in range(3)), s.zeros(3)) for l in range(3)]

    def component(i, k):
        value = s.I * sum(P[l] * gab[l][i, k] - P[k] * gab[l][i, l] for l in range(3))
        for l in range(3):
            for r in range(3):
                value += ga[l][l, r] * gb[r][i, k] + gb[l][l, r] * ga[r][i, k]
                value -= ga[l][k, r] * gb[r][i, l] + gb[l][k, r] * ga[r][i, l]
        return value

    Rab = s.Matrix(3, 3, component)
    Ra, Rb = ricci_linear(p, G), ricci_linear(q, j)
    return Rab - (G * Rb + Rb * G + j * Ra + Ra * j) / 2


def field_variation(p, q, z, a, ap, h, hp, m, mp, j, jp, t):
    """D_a F_y[(m,j)]: mixed F_N and 2 e^(H/2) F_T e^(-H/2)."""
    P = p + q
    G = 2 * z * s.eye(3) + h
    pq = p.dot(q)
    AN = -q.dot(h * P) - z * (2 * q.dot(q) - pq) + 4 * pq * a + 12 * a / t**2
    CN = AN * m - a * p.dot(j * P) + inner(hp, jp) / 2
    Hess = (m * (P * (G * q).T + (G * q) * P.T - pq * G) + a * (P * (j * p).T + (j * p) * P.T - pq * j)) / 2
    gradient = -a * m * (p * q.T + q * p.T)
    CH = 2 * ap * jp + 2 * mp * hp + 2 * stf(ricci_mixed(p, q, G, j) - 2 * Hess + 2 * gradient)
    return CN, CH


def stf_dual(action):
    """Trace-pairing dual of a linear functional on symmetric trace-free matrices."""
    basis = [s.diag(1, 0, -1), s.diag(0, 1, -1)]
    for i, j in ((0, 1), (0, 2), (1, 2)):
        E = s.zeros(3)
        E[i, j] = E[j, i] = 1
        basis.append(E)
    c = [action(E) for E in basis]
    return s.Matrix(
        [
            [(2 * c[0] - c[1]) / 3, c[2] / 2, c[3] / 2],
            [c[2] / 2, (2 * c[1] - c[0]) / 3, c[4] / 2],
            [c[3] / 2, c[4] / 2, -(c[0] + c[1]) / 3],
        ]
    )


@dataclass(frozen=True)
class FieldVariationCoefficients:
    """Coefficients in <lambda_0, D_a F_y[probe]>, before time IBP."""

    nu: object
    nu_prime: object
    H: object
    H_prime: object

    def action(self, m, mp, j, jp):
        return self.nu * m + self.nu_prime * mp + inner(self.H, j) + inner(self.H_prime, jp)


def field_coefficients(p, q, z, a, ap, h, hp, ln, lh, t):
    """Off-shell transpose coefficients; SUBTRACT these in the first source."""
    P = p + q
    G = 2 * z * s.eye(3) + h
    pq = p.dot(q)
    AN = -q.dot(h * P) - z * (2 * q.dot(q) - pq) + 4 * pq * a + 12 * a / t**2
    Vn = ln * AN - 4 * P.dot(lh * G * q) + 2 * pq * inner(lh, G) - 8 * a * p.dot(lh * q)
    Vh = -ln * a * stf(sym_outer(p, P))
    Vh += 2 * stf_dual(lambda j: inner(lh, ricci_mixed(p, q, G, j)))
    Vh -= 2 * a * stf(p * (lh * P).T + (lh * P) * p.T - pq * lh)
    return FieldVariationCoefficients(Vn, 2 * inner(lh, hp), Vh, ln * hp / 2 + 2 * ap * lh)


@dataclass(frozen=True)
class MeasurementSource:
    """Coefficients of -D_a J_y after eliminating measurement/ray STATE probes.

    H,H_prime: integrate on [e,o] and over dOmega.
    nu_emission,H_prime_emission,H_prime_observer: integrate over dOmega only.
    nu_clock_bulk: integrate on [initial,o], NO sky integral.
    nu_clock_initial: one initial-event coefficient, NO sky integral.
    All probe plane-wave phases are already included.
    """

    H: object
    H_prime: object
    nu_emission: object
    H_prime_emission: object
    H_prime_observer: object
    nu_clock_bulk: object
    nu_clock_initial: object

    def ray_bulk(self, j, jp):
        return inner(self.H, j) + inner(self.H_prime, jp)

    def sky_endpoints(self, m_e, jp_e, jp_o):
        return self.nu_emission * m_e + inner(self.H_prime_emission, jp_e) + inner(self.H_prime_observer, jp_o)


def measurement_source(q, n, w, X, v, R, K, lambda_I_a, D_a, a_e, F_a_e, A_a, e, o, initial, t):
    """Complete -D_a J_y coefficients for an independent field STATE probe.

    lambda_I_a = w Theta_a - integral(w Theta_a)/(4 pi).
    D_a is the physical first clock shift C[a]. X,v are coordinate-observer
    responses with X(o)=v(o)=0: do NOT add a clock displacement D_a*n to X.
    Clock composition is already included separately in this source.
    a_e and F_a_e include the first direction's emission phase.
    A_a = integral dOmega w partial_o U_a is ALREADY sky-integrated.
    First direction X,v,R,K include their own pulled-back phases. The probe's
    primordial input is zero; X_probe,v_probe have been eliminated exactly.
    """
    P = s.eye(3) - n * n.T
    N = stf(n * n.T)
    E1 = s.I * q.dot(X) * N / 2 - stf(sym_outer(v, n)) - stf(sym_outer(P * K, n)) / 2
    E0 = -stf(sym_outer(R, n)) / 2 - s.I * q.dot(K) * N / 2 + s.I * q.dot(n) * stf(sym_outer(K, n)) / 2
    phase = s.exp(s.I * q.dot(n) * (o - t))
    phase_e = s.exp(s.I * q.dot(n) * (o - e))
    H = phase * w * E0
    Hp = phase * (w * E1 + lambda_I_a * N / 2 + s.I * w * D_a * q.dot(n) * N / 2)
    return MeasurementSource(
        H, Hp, e * w * F_a_e * phase_e, e * w * a_e * phase_e * N / 2, w * D_a * N / 2, -A_a * t / o, -A_a * initial**2 / o
    )


@dataclass(frozen=True)
class FirstFieldAdjointSource:
    """Complete weak source as specified definite integrals, with no probe ray.

    sky is the FULL dOmega integral (not its normalized average). Field
    coefficients are already Fourier-dual/sky-contracted. Clock coefficients
    are already sky-contracted too, so neither is passed through sky again.
    action can be evaluated on the causal regular Green probes in the note.
    """

    field: FieldVariationCoefficients
    measurement: MeasurementSource
    t: object
    emission: object
    observation: object
    initial: object
    sky: Callable[[object], object]

    def action(self, m, j):
        m = s.sympify(m)
        mp, jp = m.diff(self.t), j.diff(self.t)
        field = s.Integral(self.field.action(m, mp, j, jp), (self.t, 0, self.observation))
        ray = s.Integral(self.sky(self.measurement.ray_bulk(j, jp)), (self.t, self.emission, self.observation))
        ends = self.sky(
            self.measurement.sky_endpoints(
                m.subs(self.t, self.emission), jp.subs(self.t, self.emission), jp.subs(self.t, self.observation)
            )
        )
        clock = s.Integral(self.measurement.nu_clock_bulk * m, (self.t, self.initial, self.observation))
        clock += self.measurement.nu_clock_initial * m.subs(self.t, self.initial)
        return -field + ray + ends + clock


def frame_driver_ab(h, hp, j, jp):
    """Mixed rotation generator in R=e^zeta e^(H/2) E_FW; no linear term."""
    return (h * jp - jp * h + j * hp - hp * j) / 8


def frame_terminal_covector(w, theta_gradient, n):
    """Integrand of K_a = integral skew(w grad_S Theta_a n^T) dOmega.

    The frame first-adjoint is -K_a on [initial,o]. It does not enter the first
    FIELD source: the background derivative of the rotation driver is zero.
    Its pairing with frame_driver_ab must still be retained in the cubic map.
    """
    M = w * theta_gradient * n.T
    return (M - M.T) / 2
