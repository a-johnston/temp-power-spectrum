"""Finite exact square-free amplitude algebra; no mode evolution or quadrature.

Mask bits label distinct Fourier inputs. Multiplication kills repeated bits,
which preserves the fully polarized coefficient exactly. Spatial differentiation
multiplies each coefficient by i times its subset's total momentum.
"""
import sympy as s


class AmplitudeJet:
    def __init__(self, coeff=None):
        self.coeff = {m: s.sympify(v) for m, v in (coeff or {}).items() if v != 0}

    @staticmethod
    def lift(x):
        return x if isinstance(x, AmplitudeJet) else AmplitudeJet({0: x})

    def __add__(self, other):
        out = dict(self.coeff)
        for m, v in self.lift(other).coeff.items():
            out[m] = out.get(m, 0) + v
        return AmplitudeJet(out)

    __radd__ = __add__

    def __neg__(self):
        return AmplitudeJet({m: -v for m, v in self.coeff.items()})

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) - self

    def __mul__(self, other):
        out = {}
        for a, x in self.coeff.items():
            for b, y in self.lift(other).coeff.items():
                if not a & b:
                    out[a | b] = out.get(a | b, 0) + x * y
        return AmplitudeJet(out)

    __rmul__ = __mul__

    def __truediv__(self, divisor):
        return self * (s.S.One / divisor)

    def dt(self, time):
        return AmplitudeJet({m: s.diff(v, time) for m, v in self.coeff.items()})

    def dx(self, i, momenta):
        return AmplitudeJet({m: s.I * sum(k[i] for a, k in enumerate(momenta)
                                         if m & (1 << a)) * v
                             for m, v in self.coeff.items() if m})

    def coefficient(self, mask):
        return self.coeff.get(mask, s.S.Zero)

    def normalized(self):
        return AmplitudeJet({m: s.expand(v) for m, v in self.coeff.items()})


def matrix(entries):
    return [[AmplitudeJet.lift(entries[i][j]) for j in range(3)] for i in range(3)]


def identity():
    return matrix([[int(i == j) for j in range(3)] for i in range(3)])


def add(*matrices):
    return [[sum(x[i][j] for x in matrices) for j in range(3)] for i in range(3)]


def scale(value, x):
    return [[value * x[i][j] for j in range(3)] for i in range(3)]


def product(x, y):
    return [[sum(x[i][a] * y[a][j] for a in range(3)).normalized()
             for j in range(3)] for i in range(3)]


def trace(x):
    return sum(x[i][i] for i in range(3))


def derivative(x, time):
    return [[v.dt(time) for v in row] for row in x]


def coefficient(x, mask):
    return s.Matrix(3, 3, lambda i, j: x[i][j].coefficient(mask))


def from_modes(linear, quadratic):
    """h=sum z_a H_a; j=2 sum_{a<b} z_a z_b J_ab, not j/2."""
    h = [[AmplitudeJet({1 << a: mode[i, j] for a, mode in enumerate(linear)})
          for j in range(3)] for i in range(3)]
    j = [[AmplitudeJet({(1 << a) | (1 << b): 2 * mode[i, j]
                       for (a, b), mode in quadratic.items()})
          for j in range(3)] for i in range(3)]
    return h, j


def connection_linear(h, momenta):
    return [[[(h[a][j].dx(i, momenta) + h[a][i].dx(j, momenta)
               - h[i][j].dx(a, momenta)) / 2 for j in range(3)]
             for i in range(3)] for a in range(3)]


def connection_product(matrix_left, gamma):
    return [[[sum(matrix_left[a][b] * gamma[b][i][j] for b in range(3)).normalized()
              for j in range(3)] for i in range(3)] for a in range(3)]


def ricci_derivative(gamma, momenta):
    return [[sum(gamma[a][i][j].dx(a, momenta) - gamma[a][i][a].dx(j, momenta)
                 for a in range(3)) for j in range(3)] for i in range(3)]


def ricci_product(gamma, other):
    return [[sum(gamma[a][a][b] * other[b][i][j] - gamma[a][j][b] * other[b][i][a]
                 for a in range(3) for b in range(3)) for j in range(3)] for i in range(3)]
