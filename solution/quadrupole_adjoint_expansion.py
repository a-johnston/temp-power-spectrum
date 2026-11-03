"""Closed leading-gradient field adjoints for a normalized quadrupole.

Scale physical/probe momenta by kappa. The first physical field response has
nu=-kappa^2 |p|^2 t^2/9 and H=-kappa^2 t^2 Q_p/3 + O(kappa^4).
For J=<w Theta> and w=(15/2) n.C.n, <w STF(nn)>=C. The functions below
return the kappa^2 coefficient of the FIRST field adjoint, not its full
finite-wavelength value. Times remain arbitrary, with 0<initial<e<o.

No time integration or angular integration remains in these coefficients.
The integrated spatial adjoint contracts the leading (constant-time) physical
quadratic field source, giving part, NOT all, of the kappa^4 cubic observable.
"""

from dataclasses import dataclass

import sympy as s
from unified_adjoint_source import inner, stf, sym_outer


def Q(p):
    return stf(p * p.T)


def angular_A(C, u, v):
    """<w (u.n)(v.n) STF(nn)>, exactly for quadrupolar w."""
    return stf(u.dot(v) * C + 2 * sym_outer(u, C * v) + 2 * sym_outer(v, C * u)) / 7


def angular_B(C, u, v):
    """<w (v.n) STF(sym(u n^T))>."""
    return stf(sym_outer(u, C * v))


def conformal_ricci_dual(p, q, C):
    """STF dual of ricci_mixed(p,q,2I,j), without curvature construction."""
    return stf((-q.dot(q) + p.dot(q) / 2) * C + 2 * sym_outer(q, C * q) - sym_outer(p, C * p) - sym_outer(p, C * q))


def physical_pair_source(u, v):
    """Leading kappa^2 coefficient C_H,uv (not a STATE-probe derivative)."""
    return -4 * (Q(u) + Q(v)) - 4 * stf(sym_outer(u, v))


def initial_pair(u, v):
    return (2 * (u.dot(u) + v.dot(v)) - u.dot(v)) / 9


def leading_cubic_source(a, b, c):
    return 8 * (Q(a) + Q(b) + Q(c)) + 8 * stf(sym_outer(a, b) + sym_outer(a, c) + sym_outer(b, c))


@dataclass(frozen=True)
class QuadrupoleFirstAdjoint:
    """Coefficient of kappa^2, at first-input p and STATE-probe momentum q.

    C is symmetric trace-free. Branch arguments are explicit: early means
    0<s<e; before_initial additionally means s<initial. At discontinuities use
    the one-sided prescription from the first-source note.
    """

    p: s.MatrixBase
    q: s.MatrixBase
    C: s.MatrixBase
    e: object
    o: object
    initial: object

    def __post_init__(self):
        for name in ("p", "q", "C"):
            object.__setattr__(self, name, s.Matrix(getattr(self, name)))
        for name in ("e", "o", "initial"):
            object.__setattr__(self, name, s.sympify(getattr(self, name)))
        assert self.p.shape == self.q.shape == (3, 1)
        assert self.C.shape == (3, 3) and self.C == self.C.T and s.expand(s.trace(self.C)) == 0

    @property
    def tau2(self):
        return (self.o**2 - self.e**2) / 6

    @property
    def U(self):
        p, q, C = self.p, self.q, self.C
        return (
            -p.dot(q) * C / 4
            + 3 * angular_A(C, p, q) / 2
            + angular_A(C, q, q) / 2
            - angular_B(C, p, q)
            + (angular_A(C, p, p) - angular_B(C, p, p)) / 2
        )

    @property
    def V(self):
        return 3 * angular_A(self.C, self.p, self.p) / 2 - angular_B(self.C, self.p, self.p) - self.p.dot(self.p) * self.C / 6

    @property
    def rho(self):
        return conformal_ricci_dual(self.p, self.q, self.C)

    def measurement_tensor(self, t):
        """Sky-integrated first measurement source multiplying j'(t)."""
        return (self.o - t) ** 2 * self.U + self.tau2 * self.V

    def moment(self, r, x, early):
        """Integral_x^o t^r W0(t) dt, explicitly evaluated; r=-1,0,2."""
        assert r in (-1, 0, 2)
        e, o = self.e, self.o
        lo = e if early else x
        value = (o ** (r + 2) - lo ** (r + 2)) / (2 * (r + 2)) - (o ** (r + 3) - lo ** (r + 3)) / (2 * o * (r + 3))
        if early:
            value += (1 / e - 1 / o) * (e ** (r + 3) - x ** (r + 3)) / (2 * (r + 3))
        return value

    def spatial(self, x, early):
        """Lambda_H,a^[2](x), including observer trace and emission jump."""
        p2, e, o, ti = self.p.dot(self.p), self.e, self.o, self.initial
        lo = e if early else x
        F = o**2 / lo - lo - 2 * o * s.log(o / lo)
        G = 1 / lo - 1 / o
        out = x**2 * (F * self.U + self.tau2 * G * self.V)
        out -= 2 * x * (self.moment(0, x, early) - x * self.moment(-1, x, early)) * self.rho
        out += 4 * p2 * x**2 * self.moment(-1, x, early) * self.C / 9
        if early:
            out -= p2 * e * x**2 * self.C / 18
        out += p2 * x**2 * (o**4 + 3 * ti**4) * self.C / (72 * o**3)
        return out

    def lapse(self, x, early, before_initial):
        """Lambda_N,a^[2](x); its nonzero 1/x origin is retained."""
        assert not before_initial or early
        p, q, C, e, o, ti = self.p, self.q, self.C, self.e, self.o, self.initial
        Cp, D = inner(C, Q(p)), (p + q).dot(C * q)
        lo = ti if before_initial else x
        clock = (o**4 - lo**4) / (4 * x) - x**2 * (o - lo)
        if before_initial:
            clock += ti**4 / x - x**2 * ti
        if early:
            clock -= e**4 / x - x**2 * e
        J2, Jminus1 = self.moment(2, x, early), self.moment(-1, x, early)
        return Cp * clock / 27 + 8 * D * (J2 / x - x**2 * Jminus1) / 9 + 4 * Cp * (2 * J2 / x + x**2 * Jminus1) / 27

    def origin(self):
        """Lambda_c,a^[2] = 9 lim_{x->0} x Lambda_N,a^[2](x)."""
        p, q, C, e, o, ti = self.p, self.q, self.C, self.e, self.o, self.initial
        return (3 * o**4 + 5 * ti**4 - 8 * e**4) * inner(C, Q(p)) / 20 + (o**4 - e**4) * (p + q).dot(C * q) / 5

    def integrated_spatial(self):
        """Integral_0^o Lambda_H,a^[2](x) dx, with all logarithms canceled."""
        e, o, ti = self.e, self.o, self.initial
        A = o**4 / 12 - o**2 * e**2 / 2 + 2 * o * e**3 / 3 - e**4 / 4
        out = A * self.U / 3 + self.tau2**2 * self.V - (o**4 - e**4) * self.rho / 120
        return out + self.p.dot(self.p) * (3 * o**4 + 5 * ti**4 - 8 * e**4) * self.C / 360


def first_pair_correction(p, u, v, C, e, o, initial):
    """The kappa^4 first-adjoint term for ONE single|pair partition."""
    adjoint = QuadrupoleFirstAdjoint(p, u + v, C, e, o, initial)
    return inner(adjoint.integrated_spatial(), physical_pair_source(u, v)) - adjoint.origin() * initial_pair(u, v)


def averaged_shell_first_correction(C, k, e, o, initial):
    """Isotropic |q|=|k| average of ONLY the first-adjoint kappa^4 block."""
    polynomial = 13908 * e**4 - 20960 * e**3 * o + 15720 * e**2 * o**2 - 5833 * o**4 - 2835 * initial**4
    return k.dot(k) * inner(C, Q(k)) * polynomial / 11340


def homogeneous_lapse_correction(C, k, e, o, initial):
    """The nonzero (k|q,-q) lapse contribution at kappa^4 after q averaging."""
    return -k.dot(k) * inner(C, Q(k)) * (3 * o**4 + 5 * initial**4 - 8 * e**4) / 36


def first_cubic_correction(a, b, c, C, e, o, initial):
    """Sum three partitions, each once. Not the complete cubic coefficient."""
    return sum(first_pair_correction(p, u, v, C, e, o, initial) for p, u, v in ((a, b, c), (b, a, c), (c, a, b)))
