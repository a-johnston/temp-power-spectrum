"""Entire Cartesian regular field inverse and its complete weak transpose.

The lapse-to-spatial Green convolution and proper-clock Green pairing are
EVALUATED, not left as nested Volterra definitions. All remaining integrals in
WeakGreenAdjoint are the specified source pairings, with physical endpoints and
origin kept. No projectors or inverse powers of momentum occur.
"""

from dataclasses import dataclass
from functools import cached_property

import sympy as s
from quadrupole_angular_source import ContractedMeasurement, QuadrupoleAngularSource
from unified_adjoint_source import FieldVariationCoefficients, MeasurementSource, field_coefficients, inner, stf


class TimeTail(s.Function):
    r"""T_r(a,d)=sum_m (-a)^m d^(2m+r)/(2m+r)!, entire in a and d."""

    nargs = 3

    @classmethod
    def eval(cls, r, a, d):
        if not isinstance(r, s.Integer) or r < 0:
            raise ValueError("TimeTail index must be a nonnegative integer")
        if a.is_zero:
            return d**r / s.factorial(r)
        if d.is_zero:
            return s.S.One if r == 0 else s.S.Zero

    def fdiff(self, argindex=3):
        r, a, d = self.args
        if argindex == 3:
            return TimeTail(r - 1, a, d) if r else -a * TimeTail(1, a, d)
        if argindex == 2:
            return -(d * TimeTail(r + 1, a, d) - r * TimeTail(r + 2, a, d)) / 2
        raise ValueError("The TimeTail index is discrete")

    def _eval_rewrite_as_hyper(self, r, a, d, **hints):
        return d**r * s.hyper((1,), ((r + 1) / 2, (r + 2) / 2), -a * d**2 / 4) / s.factorial(r)

    def _eval_evalf(self, prec):
        return self.rewrite(s.hyper)._eval_evalf(prec)


def spatial_operator(p, B):
    p, B = s.Matrix(p), s.Matrix(B)
    return stf(p.dot(p) * B - p * (B * p).T - (B * p) * p.T)


def scalar_convolutions(a, d):
    """A=Sh*So, B=Sh*(u So), D=Sh*(So-u Cos)/a; all entire."""
    A = (TimeTail(3, a, d) + TimeTail(3, -a, d)) / 2
    B = (d * TimeTail(3, a, d) + TimeTail(4, -a, d) - TimeTail(4, a, d)) / 2
    D = (d * TimeTail(4, a, d) - 2 * TimeTail(5, a, d) + TimeTail(5, -a, d)) / 2
    return A, B, D


def initial_lapse(a, t):
    return 3 * (TimeTail(2, a, t) - TimeTail(3, a, t) / t)


def clock_origin_primitive(a, t):
    """A primitive of t*v_q(t), vanishing at zero."""
    return 3 * (t * TimeTail(3, a, t) - 2 * TimeTail(4, a, t))


def clock_origin_weight(a, initial, observation):
    return clock_origin_primitive(a, observation) - clock_origin_primitive(a, initial) + initial**2 * initial_lapse(a, initial)


@dataclass(frozen=True)
class RegularFieldGreen:
    p: s.MatrixBase
    t: s.Symbol
    source_time: s.Expr

    @property
    def z(self):
        return self.p.dot(self.p)

    @property
    def Q(self):
        return stf(self.p * self.p.T)

    @property
    def d(self):
        return self.t - self.source_time

    @cached_property
    def spatial_coefficients(self):
        z, t, r, d = self.z, self.t, self.source_time, self.d
        return (
            r * d / t,
            -r * (TimeTail(3, z, d) + 3 * TimeTail(3, -z / 3, d)) / (4 * t),
            r * (3 * TimeTail(5, z, d) + TimeTail(5, -z / 3, d)) / (4 * t),
        )

    def spatial(self, B, derivative=0):
        """Derivative acts on the KERNEL, never on the supplied B(t)."""
        B = s.Matrix(B)
        BB = spatial_operator(self.p, B)
        BBB = spatial_operator(self.p, BB)
        coefficients = [c.diff(self.t, derivative) for c in self.spatial_coefficients]
        return coefficients[0] * B + coefficients[1] * BB + coefficients[2] * BBB

    @cached_property
    def lapse(self):
        a, d, t, r = self.z / 3, self.d, self.t, self.source_time
        return (TimeTail(1, a, d) + (d * TimeTail(2, a, d) - TimeTail(3, a, d)) / (t * r)) / 3

    @cached_property
    def coupled_scalar(self):
        """EVALUATED int_s^t G_spatial,scalar(t,u) G_N(u,s) du."""
        A, B, D = scalar_convolutions(self.z / 3, self.d)
        return (self.source_time * A + B + D / self.source_time) / (3 * self.t)

    def lapse_probe(self):
        m = self.lapse
        j = -4 * self.Q * self.coupled_scalar
        return m, m.diff(self.t), j, j.diff(self.t)

    def origin_probe(self):
        m = initial_lapse(self.z / 3, self.t)
        D = scalar_convolutions(self.z / 3, self.t)[2]
        j = -12 * self.Q * D / self.t
        return m, m.diff(self.t), j, j.diff(self.t)

    def clock_primitive(self, upper):
        """Primitive of t G_N(t,s), zero at t=s."""
        a, r = self.z / 3, self.source_time
        d = upper - r
        return (upper * TimeTail(2, a, d) - TimeTail(3, a, d) + (d * TimeTail(3, a, d) - 2 * TimeTail(4, a, d)) / r) / 3


def integral(expr, var, lo, hi):
    if isinstance(expr, s.MatrixBase):
        return expr.applyfunc(lambda x: integral(x, var, lo, hi))
    if expr == 0 or lo == hi:
        return s.S.Zero
    return s.Integral(expr, (var, lo, hi))


def below(expr, point, event):
    """Strict branch: point=event denotes the right-hand trace."""
    if isinstance(expr, s.MatrixBase):
        return expr.applyfunc(lambda x: below(x, point, event))
    if expr == 0:
        return s.S.Zero
    return s.Piecewise((expr, point < event), (0, True))


@dataclass(frozen=True)
class FieldCovectors:
    nu: s.Expr
    H: s.MatrixBase
    origin: s.Expr


@dataclass(frozen=True)
class WeakGreenAdjoint:
    """Transpose of an ALREADY sky-contracted weak source, 0<s<=o.

    There is no sky callback: common clock and field terms cannot accidentally
    receive a second sky integral. Each pairing uses a fresh integration dummy.
    """

    p: s.MatrixBase
    field: FieldVariationCoefficients
    measurement: ContractedMeasurement

    @cached_property
    def clock_amplitude(self):
        M, t, o, i = self.measurement.coefficients, self.measurement.t, self.measurement.o, self.measurement.initial
        A = -o * s.sympify(M.nu_clock_bulk) / t
        if t in A.free_symbols or s.expand(M.nu_clock_initial + A * i**2 / o) != 0:
            raise ValueError("Expected the complete physical clock source, including its initial event")
        return A

    def _at(self, value, t):
        if isinstance(value, s.MatrixBase):
            return value.subs(self.measurement.t, t)
        return s.sympify(value).subs(self.measurement.t, t)

    def spatial(self, source_time):
        t = s.Dummy("green_time", positive=True)
        e, o = self.measurement.e, self.measurement.o
        G = RegularFieldGreen(self.p, t, source_time)
        F, M = self.field, self.measurement.coefficients
        out = -integral(G.spatial(self._at(F.H, t)) + G.spatial(self._at(F.H_prime, t), 1), t, source_time, o)
        out += integral(G.spatial(self._at(M.H, t)) + G.spatial(self._at(M.H_prime, t), 1), t, s.Max(e, source_time), o)
        out += below(G.spatial(s.Matrix(M.H_prime_emission), 1).subs(t, e), source_time, e)
        out += G.spatial(s.Matrix(M.H_prime_observer), 1).subs(t, o)
        return out

    def clock_lapse(self, source_time):
        A = self.clock_amplitude
        if A == 0:
            return s.S.Zero
        i, o = self.measurement.initial, self.measurement.o
        t = s.Dummy("clock_time", positive=True)
        G = RegularFieldGreen(self.p, t, source_time)
        weight = G.clock_primitive(o) - G.clock_primitive(s.Max(i, source_time))
        weight += below(i**2 * G.lapse.subs(t, i), source_time, i)
        return -A * weight / o

    def _scalar_action(self, source_time=None):
        t = s.Dummy("green_time", positive=True)
        e, o, i = self.measurement.e, self.measurement.o, self.measurement.initial
        G = RegularFieldGreen(self.p, t, s.S.One if source_time is None else source_time)
        F, M = self.field, self.measurement.coefficients
        m, mp, j, jp = G.origin_probe() if source_time is None else G.lapse_probe()
        lo = 0 if source_time is None else source_time
        lower_ray = e if source_time is None else s.Max(e, source_time)
        out = -integral(
            self._at(F.nu, t) * m + self._at(F.nu_prime, t) * mp + inner(self._at(F.H, t), j) + inner(self._at(F.H_prime, t), jp),
            t,
            lo,
            o,
        )
        out += integral(inner(self._at(M.H, t), j) + inner(self._at(M.H_prime, t), jp), t, lower_ray, o)
        emission = M.nu_emission * m.subs(t, e) + inner(M.H_prime_emission, jp.subs(t, e))
        out += emission if source_time is None else below(emission, source_time, e)
        out += inner(M.H_prime_observer, jp.subs(t, o))
        if source_time is None:
            out -= self.clock_amplitude * clock_origin_weight(self.p.dot(self.p) / 3, i, o) / o
        else:
            out += self.clock_lapse(source_time)
        return out

    def lapse(self, source_time):
        return self._scalar_action(source_time)

    @cached_property
    def origin(self):
        return self._scalar_action()

    def at(self, source_time):
        return FieldCovectors(self.lapse(source_time), self.spatial(source_time), self.origin)

    def stationary_pair(self, Cnu, CH, c, source_time):
        covectors = self.at(source_time)
        return integral(covectors.nu * Cnu + inner(covectors.H, CH), source_time, 0, self.measurement.o) - covectors.origin * c


def background_adjoint(p, C, t, e, o, initial):
    """Complete background transpose; BOTH H and nu now use ONE time integral."""
    zero = s.zeros(3)
    coefficients = MeasurementSource(zero, s.Matrix(C) / 2, 0, zero, zero, 0, 0)
    measurement = ContractedMeasurement(coefficients, t, e, o, initial)
    field = FieldVariationCoefficients(0, 0, zero, zero)
    return WeakGreenAdjoint(p, field, measurement)


def physical_first_fields(p, t):
    """Regular scalar mode zeta=1, with no inverse momentum projectors."""
    z = p.dot(p)
    a = -z * initial_lapse(z / 3, t) / 9
    h = -2 * stf(p * p.T) * (TimeTail(2, z / 3, t) - 2 * TimeTail(3, z / 3, t) / t)
    return a, h


def physical_first_clock(p, initial, observation):
    z = p.dot(p)
    return z * clock_origin_weight(z / 3, initial, observation) / (9 * observation)


def physical_quadrupole_adjoints(p, q, C, t, u, e, o, initial):
    """Insert exact angular source in closed regular Green transpose.

    Background momentum is p+q; first reconstruction/probe momentum is q,
    including q=0. Returns background and first field-adjoint functionals, NOT
    evaluated generic time pairings or a complete cubic spectrum.
    """
    angular = QuadrupoleAngularSource(p, q, C, e, o, initial)
    background = background_adjoint(p + q, angular.background_tensor(p + q, t), t, e, o, initial)
    a, h = physical_first_fields(p, u)
    measurement = angular.assemble(1, h, a, t, u, clock_shift=physical_first_clock(p, initial, o))
    at, ht = physical_first_fields(p, t)
    field = field_coefficients(p, q, 1, at, at.diff(t), ht, ht.diff(t), background.lapse(t), background.spatial(t), t)
    return background, WeakGreenAdjoint(q, field, measurement)
