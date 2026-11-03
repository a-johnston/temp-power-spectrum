"""Numerical regular field profiles through third order at arbitrary momenta (VALIDATION ONLY).

Independent of the stationary/constraint-row machinery: every order is solved from the lapse
equation (N) and the five Cartesian trace-free rows (T), with the regular lapse coefficient from
the initial Hamiltonian row.  Sources come from the exact Christoffel engine (engine.py) with
symbolic profile functions, lambdified and fed with numerical profiles.  Profiles are stored as
Chebyshev interpolants on [0, T]; Green integrals use Gauss-Legendre quadrature of the entire
regular kernels (no projectors, no momentum division, exact K = 0 limit).

Conventions (report, Sec. 4): L_K nu + c_N = 0 with L_K = 3 d_t^2 + K^2 - 6/t^2, and
A'' + 2A'/t + B_K A = -4 TF(KK) nu - 2 c_T.  A field "label set" S carries exp(i K_S . X).
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from functools import lru_cache
from math import factorial
from typing import Any, cast

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from engine import Ctx, Field, equations  # noqa: E402
from engine import t as T_SYM  # noqa: E402

S3 = np.sqrt(3.0)


# ----------------------------------------------------------------------------- entire tails
def tails(a, d, jmax=6):
    """T_j(a, d) = sum_m (-a)^m d^(2m+j)/(2m+j)!, j = 0..jmax, entire and numerically stable."""
    d = np.asarray(d, dtype=float)
    out = []
    x2 = a * d * d
    small = np.abs(x2) < 4.0
    for j in range(jmax + 1):
        # series (used when |a| d^2 is moderate, and always for a <= 0)
        ser = np.zeros_like(d)
        term = d**j / factorial(j)
        m = 0
        while True:
            ser = ser + term
            m += 1
            term = term * (-x2) / ((2 * m + j - 1) * (2 * m + j))
            if m > 6 and np.all(np.abs(term) <= 1e-18 * np.maximum(np.abs(ser), 1e-300)):
                break
            if m > 400:
                break
        out.append(ser)
    if a > 0 and not np.all(small):
        c = np.sqrt(a)
        closed = [np.cos(c * d), np.sin(c * d) / c]
        for j in range(2, jmax + 1):
            closed.append((d ** (j - 2) / factorial(j - 2) - closed[j - 2]) / a)
        for j in range(jmax + 1):
            out[j] = np.where(small, out[j], closed[j])
    return out


def first_profiles(kmag, t):
    """nu = a = cos x - sin x / x, H = b (k^k^ - 1/3), b = -6(2 sin x/x - cos x - 1); derivatives to 3."""
    if kmag == 0:
        z = np.zeros_like(np.asarray(t, dtype=float))
        return dict(a=[z] * 4, b=[z] * 4)
    a2 = kmag**2 / 3.0  # x = sqrt(a2) t
    T = tails(a2, t, 6)
    # a = cos x - sin x/x = T0 - T1/t ; with T1 = t - a2 T3 -> a = T0 - 1 + a2 T3/t = -a2 T2 + a2 T3/t
    a = -a2 * (T[2] - T[3] / t)
    # b = -6(2 sin x/x - cos x - 1) = -6[(1 - cos x) - 2(1 - sin x/x)] = -6 a2 [T2 - 2 T3/t]
    b = -6 * a2 * (T[2] - 2 * T[3] / t)
    # derivatives: dT_j/dt = T_{j-1} (j>=1), dT0/dt = -a2 T1
    ap = -a2 * (T[1] - T[2] / t + T[3] / t**2)
    bp = 6 * (ap - a / t)  # first Hamiltonian/momentum row
    app = -(kmag**2) * a / 3 + 2 * a / t**2
    bpp = 6 * (app - ap / t + a / t**2)
    appp = -(kmag**2) * ap / 3 + 2 * ap / t**2 - 4 * a / t**3
    bppp = 6 * (appp - app / t + 2 * ap / t**2 - 2 * a / t**3)
    return dict(a=[a, ap, app, appp], b=[b, bp, bpp, bppp])


# ----------------------------------------------------------------------------- Chebyshev store
class Cheb:
    """Vector-valued Chebyshev interpolant of a function on [0, T] (trailing value shape)."""

    def __init__(self, T, values_fn, n=72, shape=()):
        self.T, self.n = T, n
        k = np.arange(n)
        x = np.cos(np.pi * (k + 0.5) / n)  # Chebyshev points of the first kind (no endpoints)
        tt = 0.5 * T * (x + 1)
        vals = values_fn(tt)  # (n, *shape)
        vals = np.asarray(vals, dtype=complex).reshape(n, -1)
        self.shape = shape
        self.coef = np.stack(
            [
                np.polynomial.chebyshev.chebfit(x, vals[:, j].real, n - 1)
                + 1j * np.polynomial.chebyshev.chebfit(x, vals[:, j].imag, n - 1)
                for j in range(vals.shape[1])
            ],
            axis=-1,
        )

    def __call__(self, t, deriv=0):
        t = np.asarray(t, dtype=float)
        x = 2 * t / self.T - 1
        c = self.coef
        for _ in range(deriv):
            c = np.polynomial.chebyshev.chebder(c, axis=0) * (2 / self.T)
        out = np.polynomial.chebyshev.chebval(x, c, tensor=True)  # (ncomp, *t.shape)
        out = np.moveaxis(out, 0, -1)
        return out.reshape(t.shape + self.shape)


GL_X, GL_W = np.polynomial.legendre.leggauss(64)


# ----------------------------------------------------------------------------- linear operators
def stf_basis():
    B = []
    for i, j in ((0, 0), (1, 1), (0, 1), (0, 2), (1, 2)):
        M = np.zeros((3, 3))
        if i == j:
            M[i, i], M[2, 2] = 1.0, -1.0
        else:
            M[i, j] = M[j, i] = 1.0
        B.append(M)
    return B


def bop_matrix(K):
    """Matrix of B_K on symmetric 3x3 (acting on flattened 9-vectors, TF-projected)."""
    K = np.asarray(K, dtype=float)
    z = K @ K
    out = np.zeros((9, 9))
    for col in range(9):
        X = np.zeros(9)
        X[col] = 1.0
        X = X.reshape(3, 3)
        Y = z * X - np.outer(K, X @ K) - np.outer(X @ K, K)
        Y = Y - np.eye(3) * np.trace(Y) / 3
        out[:, col] = Y.reshape(9)
    return out


def green_spatial(K):
    """Return functions G(t,s) and dG/dt(t,s) as (..., 9, 9) operator arrays (Cartesian, entire)."""
    K = np.asarray(K, dtype=float)
    z = float(K @ K)
    B = bop_matrix(K)
    B2 = B @ B
    eye = np.eye(9)

    def parts(d):
        Tp = tails(z, d, 6)
        Tm = tails(-z / 3, d, 6)
        P = (
            d[..., None, None] * eye
            - ((Tp[3] + 3 * Tm[3]) / 4)[..., None, None] * B
            + ((3 * Tp[5] + Tm[5]) / 4)[..., None, None] * B2
        )
        dP = eye - ((Tp[2] + 3 * Tm[2]) / 4)[..., None, None] * B + ((3 * Tp[4] + Tm[4]) / 4)[..., None, None] * B2
        return P, dP

    def G(t, s):
        d = t - s
        P, _ = parts(d)
        return (s / t)[..., None, None] * P

    def dG(t, s):
        d = t - s
        P, dP = parts(d)
        return (s / t)[..., None, None] * dP - (s / t**2)[..., None, None] * P

    return G, dG


def green_lapse(K):
    """G_N^+(t,s) (positive-forcing inverse of L_K), its t derivative, and v_K(t)."""
    a = float(np.dot(K, K)) / 3.0

    def G(t, s):
        d = t - s
        T = tails(a, d, 4)
        return (T[1] + (d * T[2] - T[3]) / (t * s)) / 3

    def dG(t, s):
        d = t - s
        T = tails(a, d, 4)
        return (T[0] + d * T[1] / (t * s) - (d * T[2] - T[3]) / (t**2 * s)) / 3

    def v(t):
        T = tails(a, t, 4)
        return 3 * (T[2] - T[3] / t)

    def vp(t):
        T = tails(a, t, 4)
        return 3 * (T[1] - T[2] / t + T[3] / t**2)

    return G, dG, v, vp


# ----------------------------------------------------------------------------- symbolic sources
@dataclass
class Profile:
    """A numerical field at a label set: lapse and spatial (3x3) interpolants, or first-order closed form."""

    labels: frozenset
    K: np.ndarray
    kind: str  # "first" or "solved"
    nu: Cheb | None = None
    nup: Cheb | None = None
    A: Cheb | None = None
    Ap: Cheb | None = None
    kmag: float = 0.0
    Qhat: np.ndarray | None = None


def _fn(name):
    """An undefined SymPy function of time (typed as callable)."""
    return cast(Any, sp.Function(name))(T_SYM)


def _sym_name(S):
    return "".join(str(x) for x in sorted(S))


@lru_cache(None)
def symbolic_sources(momenta_key, target, members):
    """Lambdified sources at label set `target` built from fields on `members` (tuple of label sets).

    momenta_key: tuple of (label, (kx,ky,kz)) with sympy-rational components.
    Returns (functions, argument names).  Singletons carry zeta = 1, nu = a_S(t), H = b_S(t) Qhat;
    larger sets carry nu = n_S(t), H = A_S(t) with five independent STF components.
    """
    momenta = {lab: tuple(sp.nsimplify(c) for c in k) for lab, k in momenta_key}
    Ctx.maxdeg = len(target)
    Ctx.momenta = momenta
    Ctx.trunc = None
    t = T_SYM
    nu, zeta = Field(), Field()
    H = [[Field() for _ in range(3)] for _ in range(3)]
    funcs = []
    for S in members:
        S = frozenset(S)
        nm = _sym_name(S)
        if len(S) == 1:
            (lab,) = tuple(S)
            k = sp.Matrix(momenta[lab])
            Qh = k * k.T / k.dot(k) - sp.eye(3) / 3 if k.dot(k) != 0 else sp.zeros(3, 3)
            fa, fb = _fn(f"a{nm}"), _fn(f"b{nm}")
            funcs += [fa, fb]
            nu += Field({S: fa})
            zeta += Field({S: sp.Integer(1)})
            for i in range(3):
                for j in range(3):
                    H[i][j] += Field({S: fb * Qh[i, j]})
        else:
            fn = _fn(f"n{nm}")
            comps = {(i, j): _fn(f"A{nm}_{i}{j}") for i, j in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2))}
            funcs += [fn] + list(comps.values())
            Am = sp.zeros(3, 3)
            for (i, j), f in comps.items():
                Am[i, j] = Am[j, i] = f
            Am[2, 2] = -Am[0, 0] - Am[1, 1]
            nu += Field({S: fn})
            for i in range(3):
                for j in range(3):
                    H[i][j] += Field({S: Am[i, j]})
    E = equations(nu, zeta, H)
    L = frozenset(target)
    exprs = [E["N"].coeff(L), E["Ham"].coeff(L)] + [E["M"][i].coeff(L) for i in range(3)]
    exprs += [E["T"][i][j].coeff(L) for i in range(3) for j in range(3)]
    rep, names = {}, []
    for f in funcs:
        for d in (3, 2, 1):
            sym = sp.Symbol(f"{f.func.__name__}_d{d}")
            rep[sp.Derivative(f, (t, d))] = sym
            names.append(sym)
        sym = sp.Symbol(f"{f.func.__name__}_d0")
        rep[f] = sym
        names.append(sym)
    lam = sp.lambdify([t] + names, [sp.expand(e).xreplace(rep) for e in exprs], "numpy")
    return lam, tuple(str(n) for n in names)


def evaluate_sources(momenta, target, profiles, t):
    """Numerical (c_N, c_Ham, c_M[3], c_T[3x3]) at times t from the given member profiles."""
    key = tuple((lab, tuple(sp.nsimplify(c) for c in k)) for lab, k in sorted(momenta.items()))
    members = tuple(sorted((tuple(sorted(S)) for S in profiles), key=lambda x: (len(x), x)))
    lam, names = symbolic_sources(key, tuple(sorted(target)), members)
    vals = {}
    for S, prof in profiles.items():
        nm = _sym_name(S)
        if prof.kind == "first":
            f = first_profiles(prof.kmag, t)
            for d in range(4):
                vals[f"a{nm}_d{d}"] = f["a"][d]
                vals[f"b{nm}_d{d}"] = f["b"][d]
        else:
            nus = [prof.nu(t), prof.nup(t), prof.nup(t, 1), prof.nup(t, 2)]
            As = [prof.A(t), prof.Ap(t), prof.Ap(t, 1), prof.Ap(t, 2)]
            for d in range(4):
                vals[f"n{nm}_d{d}"] = nus[d]
                for i, j in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2)):
                    vals[f"A{nm}_{i}{j}_d{d}"] = As[d][..., i, j]
    out = lam(t, *[vals[n] for n in names])
    out = [np.broadcast_to(np.asarray(o, dtype=complex), np.shape(t)) for o in out]
    cN, cH = out[0], out[1]
    cM = np.stack(out[2:5], axis=-1)
    cT = np.stack(out[5:14], axis=-1).reshape(np.shape(t) + (3, 3))
    return cN, cH, cM, cT


# ----------------------------------------------------------------------------- the regular solve
def solve_order(momenta, target, profiles, c_reg, T, n=72):
    """Regular solution (nu, A) at label set `target` with lapse coefficient c_reg (initial row)."""
    K = sum((np.asarray(momenta[lab], dtype=float) for lab in target), np.zeros(3))
    GN, dGN, vK, vKp = green_lapse(K)
    src = lambda s: evaluate_sources(momenta, target, profiles, s)  # noqa: E731

    def nu_vals(t):
        cN = lambda s: -src(s)[0]  # noqa: E731
        s = 0.5 * t[:, None] * (GL_X[None, :] + 1)
        w = 0.5 * t[:, None] * GL_W[None, :]
        f = cN(s)
        return c_reg * vK(t) + np.sum(GN(t[:, None], s) * f * w, axis=1), c_reg * vKp(t) + np.sum(
            dGN(t[:, None], s) * f * w, axis=1
        )

    cache = {}

    def nu_pair(t):
        key = id(t)
        if key not in cache:
            cache[key] = nu_vals(t)
        return cache[key]

    nu = Cheb(T, lambda t: nu_pair(t)[0], n)
    nup = Cheb(T, lambda t: nu_pair(t)[1], n)
    G, dG = green_spatial(K)
    QK = np.outer(K, K) - np.eye(3) * (K @ K) / 3

    def rhs(s):
        cT = src(s)[3]
        return (-4 * QK * nu(s)[..., None, None] - 2 * cT).reshape(s.shape + (9,))

    def A_vals(t, deriv):
        s = 0.5 * t[:, None] * (GL_X[None, :] + 1)
        w = 0.5 * t[:, None] * GL_W[None, :]
        Gm = (dG if deriv else G)(t[:, None], s)  # (nt, ng, 9, 9)
        f = rhs(s)  # (nt, ng, 9)
        return np.einsum("tgij,tgj,tg->ti", Gm, f, w).reshape(t.shape + (3, 3))

    A = Cheb(T, lambda t: A_vals(t, 0), n, (3, 3))
    Ap = Cheb(T, lambda t: A_vals(t, 1), n, (3, 3))
    return Profile(frozenset(target), K, "solved", nu=nu, nup=nup, A=A, Ap=Ap)


def initial_row(momenta, S):
    """Regular lapse t^2 coefficient (report eq. grad): c_2 and c_3 kernels."""
    ks = [np.asarray(momenta[lab], dtype=float) for lab in sorted(S)]
    if len(ks) == 2:
        return (2 * (ks[0] @ ks[0] + ks[1] @ ks[1]) - ks[0] @ ks[1]) / 9
    if len(ks) == 3:
        return 2 * (-2 * sum(k @ k for k in ks) + ks[0] @ ks[1] + ks[0] @ ks[2] + ks[1] @ ks[2]) / 9
    raise ValueError


def build_profiles(momenta, T, n=72, third=True):
    """All profiles for labels 1..3: first order closed form, pairs and the triple solved."""
    labels = sorted(momenta)
    prof = {}
    for lab in labels:
        k = np.asarray(momenta[lab], dtype=float)
        km = float(np.sqrt(k @ k))
        Qh = np.outer(k, k) / km**2 - np.eye(3) / 3 if km > 0 else np.zeros((3, 3))
        prof[frozenset([lab])] = Profile(frozenset([lab]), k, "first", kmag=km, Qhat=Qh)
    from itertools import combinations

    for a, b in combinations(labels, 2):
        sub = {frozenset([a]): prof[frozenset([a])], frozenset([b]): prof[frozenset([b])]}
        prof[frozenset([a, b])] = solve_order(momenta, (a, b), sub, initial_row(momenta, (a, b)), T, n)
    if third and len(labels) == 3:
        prof[frozenset(labels)] = solve_order(momenta, tuple(labels), dict(prof), initial_row(momenta, labels), T, n)
    return prof
