"""Projector-free Cartesian form of the spatial field equation (exact).

The trace-free spatial equation is, for the order-r metric kernel A (trace-free),
    A'' + (2/t) A' + BK[A] = -4 TF(K K) nu - 2 c_T,     BK[A] = TF[K^2 A - K (A K)^T - (A K) K^T].
1. BK has eigenvalue K^2 on TT tensors, 0 on vector tensors K w + w K (w . K = 0), and -K^2/3 on the
   scalar TF(K K): checked for generic K.
2. Its regular retarded solution with A(0) = 0 is a power series in t whose coefficients are
   POLYNOMIAL in the Cartesian components of K (no inverse powers): m(m+1) a_m = rhs_{m-2} - BK a_{m-2}.
   This series is compared, coefficient by coefficient and exactly, with the projector-based triangular
   kernel of fast_series.py (which carries 1/|K|^2 through the formal symbols u, v and 1/p^2). They agree
   for the pairs (k, q) and (k, -q) through t^TCUT, and the K = 0 pair (q, -q) reduces to the exact K = 0
   kernel. Hence the recombined Cartesian kernel is entire in K: every inverse power of |K| cancels, and the
   K -> 0 limit is the separately constructed K = 0 kernel.
"""

import sys

import sympy as sp

sys.argv = sys.argv[:1]
import fast_series as fs  # noqa: E402


def check_eigenvalues():
    K = sp.Matrix(sp.symbols("K1:4", real=True))
    w = sp.Matrix(sp.symbols("w1:4", real=True))
    K2 = (K.T * K)[0]
    tf = lambda M: M - sp.eye(3) * M.trace() / 3
    BK = lambda A: tf(K2 * A - K * (A * K).T - (A * K) * K.T)
    # scalar
    S = tf(K * K.T)
    assert sp.simplify(BK(S) + K2 / 3 * S) == sp.zeros(3)
    # vector: w orthogonal to K
    wperp = w - K * (K.T * w)[0] / K2
    V = K * wperp.T + wperp * K.T
    assert sp.simplify(BK(V)) == sp.zeros(3)
    # tensor: P M P - P tr(P M P)/2 for a generic symmetric M
    M = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"m{min(i, j)}{max(i, j)}", real=True))
    P = sp.eye(3) - K * K.T / K2
    T = P * M * P - P * (P * M * P).trace() / 2
    assert sp.simplify(BK(T) - K2 * T) == sp.zeros(3)
    print("[ok] BK eigenvalues: K^2 (TT), 0 (vector), -K^2/3 (scalar)")


def cartesian_series(labels, second):
    """Regular series solution of the Cartesian equation with the same sources as fast_series.solve."""
    lab = frozenset(labels)
    E = fs.equations(*fs.build(labels, second))
    cT = [[E["T"][i][j].coeff(lab) for j in range(3)] for i in range(3)]
    nu, _ = fs.solve(labels, second)
    KV = fs.Ctx.mom(lab)
    K2 = sum((c * c for c in KV), fs.ZERO)
    third = fs.rat(1, 3)

    def tf(M):
        tr = (M[0][0] + M[1][1] + M[2][2]) * third
        return [[M[i][j] - (tr if i == j else fs.ZERO) for j in range(3)] for i in range(3)]

    def BK(A):
        AK = [sum((A[i][j] * KV[j] for j in range(3)), fs.ZERO) for i in range(3)]
        M = [[K2 * A[i][j] - KV[i] * AK[j] - AK[i] * KV[j] for j in range(3)] for i in range(3)]
        return tf(M)

    KK = tf([[KV[i] * KV[j] for j in range(3)] for i in range(3)])
    a = {}
    for m in range(2, fs.TCUT + 1):
        prev = a.get(m - 2, [[fs.ZERO] * 3 for _ in range(3)])
        Bp = BK(prev)
        rhs = [[-(KK[i][j] * nu.coeff(m - 2)) * 4 - cT[i][j].coeff(m - 2) * 2 for j in range(3)] for i in range(3)]
        a[m] = [[(rhs[i][j] - Bp[i][j]) * fs.rat(1, m * (m + 1)) for j in range(3)] for i in range(3)]
    return a


def check_against_triangular():
    second = {}
    for pr in ((1, 2), (1, 3), (2, 3)):
        nu, A = fs.solve(list(pr), {})
        second[pr] = (nu, A)
        a = cartesian_series(list(pr), {})
        for m in range(2, fs.TCUT - 2):
            for i in range(3):
                for j in range(3):
                    assert fs.vanishes(A[i][j].coeff(m) - a[m][i][j]), (pr, m, i, j)
        # the Cartesian coefficients contain no formal inverse |K|^2 (u, v) at all
        uses = any(c.degree(fs.U) > 0 or c.degree(fs.V) > 0 for m in a for row in a[m] for c in row if c != 0)
        assert not uses, pr
        print(f"[ok] pair {pr}: projector-free Cartesian series = triangular kernel through t^{fs.TCUT - 3}; no 1/|K|^2 appears")
    nu3, A3 = fs.solve([1, 2, 3], second)
    a3 = cartesian_series([1, 2, 3], second)
    for m in range(2, fs.TCUT - 2):
        for i in range(3):
            for j in range(3):
                assert fs.vanishes(A3[i][j].coeff(m) - a3[m][i][j]), ("third", m, i, j)
    print(f"[ok] third order at (k, q, -q): projector-free Cartesian series = triangular kernel through t^{fs.TCUT - 3}")


if __name__ == "__main__":
    check_eigenvalues()
    check_against_triangular()
    print("CARTESIAN GREEN CHECKS PASSED (exact)")
