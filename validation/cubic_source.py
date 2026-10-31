"""Explicit cubic source, using the finite polarization specified in the derivation.

Inputs are three linear metric matrices and three quadratic pair kernels at the
same symbolic time. The actual pair kernels are the exact retarded expressions
in QUADRATIC_INITIAL_DATA.md; tests may instead supply declared boundary jets.
"""
import sympy as s
from amplitude_jets import (
    AmplitudeJet, add, coefficient, connection_linear, connection_product,
    derivative, from_modes, identity, product, ricci_derivative, ricci_product,
    scale, trace,
)
from initial_data import spatial_operator


def cubic_source(time, momenta, linear, quadratic):
    h, j = from_modes(linear, quadratic)
    I = identity()
    hp, jp = derivative(h, time), derivative(j, time)
    hh = product(h, h)
    zeta = AmplitudeJet({1 << a: 1 for a in range(3)})
    c = trace(h) / 6 - zeta
    d = trace(j) / 6 - trace(hh) / 6
    w = trace(product(h, j)) / 2 - trace(product(hh, h)) / 3
    cp, dp, wp = c.dt(time), d.dt(time), w.dt(time)
    g1 = connection_linear(h, momenta)
    gj = connection_linear(j, momenta)
    hg = connection_product(h, g1)
    g2 = [[[gj[a][i][l] - 2 * hg[a][i][l] for l in range(3)]
           for i in range(3)] for a in range(3)]
    hgj = connection_product(h, gj)
    other = connection_product(add(scale(6, hh), scale(-3, j)), g1)
    g3 = [[[-3 * hgj[a][i][l] + other[a][i][l] for l in range(3)]
           for i in range(3)] for a in range(3)]
    r1 = ricci_derivative(g1, momenta)
    r2 = add(ricci_derivative(g2, momenta), scale(2, ricci_product(g1, g1)))
    r3 = add(ricci_derivative(g3, momenta), scale(3, ricci_product(g1, g2)),
             scale(3, ricci_product(g2, g1)))
    spatial = []
    for i in range(3):
        spatial.append([])
        for l in range(3):
            ci, cl = c.dx(i, momenta), c.dx(l, momenta)
            di, dl = d.dx(i, momenta), d.dx(l, momenta)
            first = cl.dx(i, momenta) - r1[i][l]
            second = (dl.dx(i, momenta)
                      - 2 * sum(g1[a][i][l] * c.dx(a, momenta) for a in range(3))
                      + 2 * ci * cl - r2[i][l])
            spatial[i].append(
                -2 * w.dx(l, momenta).dx(i, momenta)
                -6 * sum(g1[a][i][l] * d.dx(a, momenta)
                         + g2[a][i][l] * c.dx(a, momenta) for a in range(3))
                + 6 * (ci * dl + di * cl) - 2 * r3[i][l]
                + 12 * c * second + (12 * d + 24 * c * c) * first
            )
    source = add(
        scale(3, product(hp, jp)), scale(3, product(jp, hp)),
        scale(-6, product(product(hp, h), hp)),
        scale(-6 * cp, jp), scale(-6 * dp, hp),
        scale(4 * (w / time**2 + wp / time) + 24 * c * d / time**2
              - 16 * c * c * c / time**2, I),
        scale(-12 * d / time**2 + 24 * c * c / time**2 - 12 * dp / time, h),
        scale(-12 * (c / time**2 + cp / time), j), spatial,
    )
    return coefficient(source, 7) / 6


def cubic_regular_source(time, momenta, linear, quadratic):
    total = sum((s.Matrix(k) for k in momenta), s.zeros(3, 1))
    return (cubic_source(time, momenta, linear, quadratic)
            - 16 * s.eye(3) / time**2 - spatial_operator(total, 8 * s.eye(3)))
