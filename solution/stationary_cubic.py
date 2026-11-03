"""Cubic stationary assembly with an explicit first-field-only measurement term.

No second physical fields or second ray responses occur. This module supplies
local ray forcing/driver identities, endpoint/clock composition, normalization,
and the field/initial-data contraction. Integrals and sky operations are passed
explicitly; it does not evaluate the complete spectrum.
"""

from dataclasses import dataclass

import sympy as s
from unified_adjoint_source import inner


def comm(A, B):
    return A * B - B * A


def directional(gradient, v):
    return sum((v[i] * gradient[i] for i in range(3)), s.zeros(3))


def linear_direction_driver(n, hp, gradG):
    P = s.eye(3) - n * n.T
    grad_nn = s.Matrix([inner(n * n.T, g) for g in gradG])
    return P * hp * n / 2 + (grad_nn - directional(gradG, n) * n) / 2


def gradient_linear_driver(n, grad_hp, grad2G):
    """grad2G[i][j] = partial_i partial_j G; n is held fixed."""
    return [linear_direction_driver(n, grad_hp[i], grad2G[i]) for i in range(3)]


def second_ray_forcing(n, G, hp, gradG, grad_driver, X, v):
    """X2'=u2+QX, u2'=QV, v2=u2+n |v1|^2 for first-field-only metrics.

    G=2 zeta_1 I+H_1, hp=partial_t H_1; X,v are fixed-coordinate-observer
    first responses, v transverse. All field jets are pulled back afterwards.
    u2 is transverse and has the physical frame boundary u2(o)=-R_FW,2 n.
    """
    P = s.eye(3) - n * n.T
    Gn = G * n
    nG = directional(gradG, n)
    grad_nn = s.Matrix([inner(n * n.T, g) for g in gradG])
    grad_vGn = s.Matrix([v.dot(g * n) for g in gradG])
    pull_driver = sum((X[i] * grad_driver[i] for i in range(3)), s.zeros(3, 1))
    QX = n * v.dot(v) - G * v - G * G * n / 4 + directional(gradG, X) * n
    raw = 2 * pull_driver - hp * v + v * inner(n * n.T, hp) - comm(G, hp) * n / 4
    raw += directional(gradG, v) * n + nG * v - 2 * grad_vGn
    raw += directional(gradG, Gn) * n / 2 + comm(G, nG) * n / 4 - G * grad_nn / 2
    return QX, P * raw


def first_only_driver_cubic(n, G, hp, grad_hp, grad2_hp, X, v):
    """Known local remainder F3-3(r.X2+s.u2), including the radial v2 jet."""
    nn = lambda A: inner(n * n.T, A)
    hessian = sum((X[i] * X[j] * nn(grad2_hp[i][j]) for i in range(3) for j in range(3)), s.Integer(0))
    return (
        nn(comm(G, comm(G, hp))) / 8
        + 3 * hessian / 2
        - 6 * v.dot(directional(grad_hp, X) * n)
        + 3 * v.dot(hp * v)
        - 3 * v.dot(v) * nn(hp)
    )


def reduced_cubic_bulk(n, G, hp, gradG, grad_driver, grad_hp, grad2_hp, X, v, R, K):
    QX, QV = second_ray_forcing(n, G, hp, gradG, grad_driver, X, v)
    local = first_only_driver_cubic(n, G, hp, grad_hp, grad2_hp, X, v)
    return local - 3 * R.dot(QX) - 3 * K.dot(QV)


def frame_exponent_boundary(K_observer, frame_second, n):
    """Add to U3; frame_second = (1/4) integral_[initial,o] [h,h'] dt."""
    return -3 * K_observer.dot(frame_second * n)


def emission_first_only(e, a_e, total_a_prime_e, X_grad_a_e):
    d1 = -e * a_e
    d2 = e * (a_e**2 - 2 * X_grad_a_e - 2 * d1 * total_a_prime_e)
    return d1, d2


def clock_first_only(initial, o, a_initial, ap_initial, a_observer, int_t_a, int_t_a_squared):
    """Physical initial T event and proper clock for nu=epsilon a only."""
    D1 = -(initial**2 * a_initial + int_t_a) / o
    D2 = (2 * initial**3 * a_initial * ap_initial - int_t_a_squared - 2 * o * a_observer * D1 - D1**2) / o
    return D1, D2


@dataclass(frozen=True)
class ExponentJets:
    first: object
    second: object
    third: object


def compose_first_only(U1, U2, U3, F1e, F2e, total_F1prime_e, d1, d2, D1, D2, U1o, U1oo, U2o, F1eo, d1o):
    """Moving emission + physical clock, with all first-only cross products.

    o-derivatives are taken at fixed coordinate observer before composition;
    total_F1prime_e is the background-ray derivative, not partial_t F1.
    """
    second = U2 - 2 * F1e * d1 + 2 * D1 * U1o
    third = U3 - 3 * F1e * d2 - 3 * F2e * d1 - 3 * total_F1prime_e * d1**2
    third += 3 * D1 * (U2o - 2 * (F1eo * d1 + F1e * d1o)) + 3 * D2 * U1o + 3 * D1**2 * U1oo
    return ExponentJets(U1, second, third)


def normalized_map_cubic(jets, average):
    """Exact third derivative of exp(-I)/average(exp(-I))-1.

    No intermediate dipole projection. average is normalized to average(1)=1.
    """
    u = jets.first - average(jets.first)
    v = jets.second - average(jets.second)
    raw = -jets.third + 3 * u * v - u**3 + 3 * average(u**2) * u
    return raw - average(raw)


def cubic_field_integrand(ln0, lh0, N3, H3, ln1, lh1, N2, H2):
    """Diagonal cubic: lambda0.C3+3 lambda1.C2; physical-input sources."""
    return ln0 * N3 + inner(lh0, H3) + 3 * (ln1 * N2 + inner(lh1, H2))


def assemble_polarized_cubic(measurement, ln0, lh0, N3, H3, partitions, lc0, c3, initial_partitions, t, o):
    """Arbitrary-triple stationary functional with exactly three partitions.

    partitions = [(Lambda_N,a,Lambda_H,a,C_N,bc,C_H,bc), ...].
    initial_partitions = [(Lambda_c,a,c_bc), ...], in matching partition order.
    All sources are PHYSICAL-input polarized kernels, not diagonal expressions.
    The field/initial Fourier pairings have already been taken; time remains.
    """
    assert len(partitions) == len(initial_partitions) == 3
    bulk = ln0 * N3 + inner(lh0, H3)
    for ln, lh, N2, H2 in partitions:
        bulk += ln * N2 + inner(lh, H2)
    initial = -lc0 * c3 - sum(lc * c2 for lc, c2 in initial_partitions)
    return measurement + s.Integral(bulk, (t, 0, o)) + initial


def initial_cubic_pairing(lc0, c3, lc1, c2):
    """Residual is coefficient(nu/t^2)-c[g], fixing both minus signs."""
    return -lc0 * c3 - 3 * lc1 * c2


def polarize_diagonal_cubic(expr, labels):
    """Extract the arbitrary-triple kernel from a diagonal third derivative.

    A cubic diagonal T[g,g,g] gives 6 T[a,b,c] as the coefficient of three
    distinct input labels. This differs from extracting labels from the FULL
    Taylor observable, whose 1/3! has already canceled the permutations.
    """
    assert len(labels) == 3
    out = s.expand(expr)
    for label in labels:
        out = out.coeff(label, 1)
    return out / 6


def initial_diagonal_sources(g, gradg, lapg):
    """c2[g,g], c3[g,g,g] from the nonlinear regular Hamiltonian matching."""
    square = sum(x * x for x in gradg)
    return (square - 4 * g * lapg) / 9, (4 * g * g * lapg - 2 * g * square) / 3


def direct_field_residuals(nu, zeta, H):
    """Finite exact coefficient prescription for PHYSICAL-input C2,C3 and c.

    Inputs are engine.Field objects containing FIRST singleton labels only.
    Unlike first_field_adjoint_source, every primordial input label is present.
    Query degree-2/3 label coefficients of the returned fields directly: their
    Taylor factorials are already canceled, so do NOT divide those by 2!/3!.
    The engine context must have maxdeg=3 and the input momenta configured.
    """
    import engine as eng

    assert eng.Ctx.maxdeg == 3
    for f in [nu, zeta, *[H[i][j] for i in range(3) for j in range(3)]]:
        assert all(len(labels) == 1 for labels in f.d)
    assert all(s.diff(x, eng.t) == 0 for x in zeta.d.values())
    eq = eng.equations(nu, zeta, H)
    plus = eng.mexp(eng.mscale(H, s.Rational(1, 2)))
    minus = eng.mexp(eng.mscale(H, -s.Rational(1, 2)))
    EH = eng.mscale(eng.mmul(eng.mmul(plus, eq["T"]), minus), 2)
    lap = sum((zeta.d_(i).d_(i) for i in range(3)), eng.Field())
    grad_square = sum((zeta.d_(i) * zeta.d_(i) for i in range(3)), eng.Field())
    initial = eng.fexp(-2 * zeta) * (2 * lap + grad_square) / 18
    return eq["N"], EH, initial


def direct_field_kinetic_cubic(h, hp, hpp, t):
    """Noncommuting temporal part of D^3 E_H along a first-field-only family.

    D^3 tr(S^2)=0 in E_N; the remaining direct cubic lapse terms are spatial
    geometry and -24 a^3/t^2. No on-shell/monochromatic assumption is used.
    """
    return comm(h, comm(h, hpp)) / 4 + comm(hp, comm(h, hp)) + comm(h, comm(h, hp)) / (2 * t)
