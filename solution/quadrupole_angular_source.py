"""Exact finite-wavelength sky closure of the first quadrupole field source.

J=<w Theta>, w=(15/2)n.C.n, so <w STF(nn)>=C. All returned coefficients
are ALREADY sky integrated. Mixed-time kernels use F_0,...,F_7 and retain
normalization products at separate phases. Time/Green integrals remain.
There is no gradient expansion, momentum projection or division by |phase|.
"""

from dataclasses import dataclass, field

import sympy as s
from entire_sphere import phase_average
from unified_adjoint_source import FirstFieldAdjointSource, MeasurementSource, inner, stf, sym_outer


def integrate(expr, u, lo, hi):
    if isinstance(expr, s.MatrixBase):
        return expr.applyfunc(lambda z: integrate(z, u, lo, hi))
    return s.S.Zero if expr == 0 else s.Integral(expr, (u, lo, hi))


@dataclass(frozen=True)
class RayPolynomials:
    X_upper: s.MatrixBase
    v_upper: s.MatrixBase
    R_lower: s.MatrixBase
    K_lower: s.MatrixBase


def ray_polynomials(p, n, z, h, hp, t, u):
    """v=-int_t^o phi v_upper, X=+int_t^o phi X_upper.

    R=+int_e^t phi R_lower, K=+int_e^t phi K_lower. Field derivatives
    are taken BEFORE pullback. n is a unit vector; v_upper is transverse.
    """
    P = s.eye(3) - n * n.T
    G = 2 * z * s.eye(3) + h
    V = P * hp * n / 2 + s.I * (p * n.dot(G * n) - p.dot(n) * G * n) / 2
    r = s.I * p * n.dot(hp * n) / 2
    return RayPolynomials((u - t) * V - G * n / 2, V, r, -P * hp * n - (t - u) * r)


@dataclass(frozen=True)
class AngularBlocks:
    H_lower: s.MatrixBase
    Hprime_upper: s.MatrixBase
    Hprime_lower: s.MatrixBase
    Hprime_normalization: s.MatrixBase


@dataclass(frozen=True)
class ContractedMeasurement:
    """No angular variable remains. Never apply a sky integral a second time."""

    coefficients: MeasurementSource
    t: object
    e: object
    o: object
    initial: object

    def with_field(self, field_coefficients):
        """Combine with the already Fourier-dual background field variation."""
        return FirstFieldAdjointSource(
            field_coefficients, self.coefficients, self.t, self.e, self.o, self.initial, lambda expr: expr
        )


@dataclass(frozen=True)
class QuadrupoleAngularSource:
    p: s.MatrixBase
    q: s.MatrixBase
    C: s.MatrixBase
    e: object
    o: object
    initial: object
    n: s.MatrixBase = field(default_factory=lambda: s.Matrix(s.symbols("n0:3", real=True)))

    def __post_init__(self):
        for name in ("p", "q", "C", "n"):
            object.__setattr__(self, name, s.Matrix(getattr(self, name)))
        for name in ("e", "o", "initial"):
            object.__setattr__(self, name, s.sympify(getattr(self, name)))
        assert self.p.shape == self.q.shape == self.n.shape == (3, 1)
        assert self.C.shape == (3, 3) and self.C == self.C.T and s.expand(s.trace(self.C)) == 0

    @property
    def weight(self):
        return s.Rational(15, 2) * self.n.dot(self.C * self.n)

    @property
    def N(self):
        return stf(self.n * self.n.T)

    def plain(self, poly, K):
        return phase_average(poly, self.n, K, maximum_rank=7)

    def weighted(self, poly, K):
        return phase_average(poly, self.n, K, self.weight, maximum_rank=7)

    def phases(self, t, u):
        Kp = self.p * (self.o - u)
        Kq = self.q * (self.o - t)
        return Kp, Kq, Kp + Kq

    def polynomial_blocks(self, z, h, hp, t, u):
        """Phase-stripped Volterra atoms, before their exact sky contraction."""
        n, q, N = self.n, self.q, self.N
        P = s.eye(3) - n * n.T
        R = ray_polynomials(self.p, n, z, h, hp, t, u)
        upper = s.I * q.dot(R.X_upper) * N / 2 + stf(sym_outer(R.v_upper, n))
        lower = -stf(sym_outer(P * R.K_lower, n)) / 2
        bulk = (
            -stf(sym_outer(R.R_lower, n)) / 2 - s.I * q.dot(R.K_lower) * N / 2 + s.I * q.dot(n) * stf(sym_outer(R.K_lower, n)) / 2
        )
        return bulk, upper, lower

    def normalization(self, hp, t, u):
        """Complete random-mean normalization atom; includes BOTH products."""
        Kp, Kq, K = self.phases(t, u)
        hpnn = self.n.dot(hp * self.n)
        out = -self.weighted(hpnn * self.N, K) / 4
        out += self.plain(hpnn, Kp) * self.weighted(self.N, Kq) / 4
        out += self.weighted(hpnn, Kp) * self.plain(self.N, Kq) / 4
        return out

    def blocks(self, z, h, hp, t, u):
        bulk, upper, lower = self.polynomial_blocks(z, h, hp, t, u)
        K = self.phases(t, u)[2]
        return AngularBlocks(
            self.weighted(bulk, K), self.weighted(upper, K), self.weighted(lower, K), self.normalization(hp, t, u)
        )

    def background_tensor(self, P, t):
        """Exact sky tensor entering the background spatial adjoint."""
        return self.weighted(self.N, P * (self.o - t))

    def frame_density(self, hp, u):
        """First frame covector K_a=[S_a,C], using only the linear quadrupole.

        Integrate this density on [e,o]; the frame FIRST ADJOINT is -K_a.
        This is separate from the field source, not an extra term to add to
        stationary_cubic's already included frame-exponent boundary.
        """
        S = -s.Rational(15, 4) * self.plain(self.n.dot(hp * self.n) * self.N, self.p * (self.o - u))
        return S * self.C - self.C * S

    def frame_covector(self, h, u):
        return integrate(self.frame_density(h.diff(u), u), u, self.e, self.o)

    def assemble(self, z, h, a, t, u):
        """Complete angularly integrated measurement source, with known h(u),a(u).

        The remaining integrals have explicit fixed/Volterra limits. z is the
        time-independent primordial amplitude; the probe has delta zeta=0.
        Combine through .with_field(...), whose sky operation is IDENTITY.
        """
        assert t != u and s.diff(s.sympify(z), u) == 0
        h = s.Matrix(h)
        a = s.sympify(a)
        hp = h.diff(u)
        b = self.blocks(z, h, hp, t, u)
        D = -(self.initial**2 * a.subs(u, self.initial) + integrate(u * a, u, self.initial, self.o)) / self.o
        A = inner(self.C, hp.subs(u, self.o)) / 2
        A += integrate(
            self.weighted(s.I * self.p.dot(self.n) * self.n.dot(hp * self.n) / 2, self.p * (self.o - u)), u, self.e, self.o
        )
        H = integrate(b.H_lower, u, self.e, t)
        Hp = (
            integrate(b.Hprime_upper, u, t, self.o)
            + integrate(b.Hprime_lower, u, self.e, t)
            + integrate(b.Hprime_normalization, u, self.e, self.o)
        )
        Hp += s.I * D * self.weighted(self.q.dot(self.n) * self.N, self.q * (self.o - t)) / 2
        Ke = (self.p + self.q) * (self.o - self.e)
        lapse_e = self.e * self.weighted(self.n.dot(hp.subs(u, self.e) * self.n), Ke) / 2
        Hpe = self.e * a.subs(u, self.e) * self.weighted(self.N, Ke) / 2
        Hpo = D * self.C / 2
        coeffs = MeasurementSource(H, Hp, lapse_e, Hpe, Hpo, -A * t / self.o, -A * self.initial**2 / self.o)
        return ContractedMeasurement(coeffs, t, self.e, self.o, self.initial)
