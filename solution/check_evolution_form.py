"""Exact equivalence of the unified evolution/constraint system with the 4D Einstein equations.

Unified form (n = ln N, S = gamma^{-1} sigma the mixed shear, all indices raised with gamma):

  (N)  R_00 = 3 a1^2/N^2        <=>  3 n'' = Delta n + 2|Dn|^2 - tr(S S) - 3 a1^2/N^2
  (M)  R_0i = 0                 <=>  D_j S^j_i - 2 d_i n' + 2 S^j_i d_j n + 2 n' d_i n = 0
  (T)  TF[gamma^{-1} R_(ij)] = 0 <=>  TF[ S' + 2 n' S + R^k_j - 2 D^k D_j n + 2 D^k n D_j n ] = 0

(T) is new relative to both earlier branches: in mixed form the sigma.sigma terms of the
spatial equation cancel exactly, because (gamma^{-1})' sigma = -2 S S.  The physical source
R_ij = a1^2 gamma_ij / N^2 is pure trace, so only the trace-free (TF) part is an evolution
equation; its trace combines with (N) into the Hamiltonian constraint (proved upstream in
verification/check_field_equations_generic.py).

Proof method (as upstream): every function value and derivative is replaced by an independent
symbol (2-jet), for a generic six-function gamma with t-independent det and arbitrary N.
"""
import os
import sys
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'verification'))
sys.argv = sys.argv[:1]  # upstream module parses argv for its mode; force 'jet'
import check_field_equations_generic as up  # noqa: E402

ev = up.make_ev('jet')
X, X3, t = up.X, up.X3, up.t
gam, gi, N, G3, G4 = up.gam, up.gi, up.N, up.G3, up.G4

R4 = up.ricci_at(G4, X, ev)
R3 = up.ricci_at(G3, X3, ev)
n = sp.log(N)
dn = [sp.diff(n, q) for q in X3]
ndot = sp.diff(n, t)
sig = sp.Matrix(3, 3, lambda i, j: sp.diff(gam[i, j], t) / 2)
S = gi * sig                                   # S^k_j
DDn = sp.Matrix(3, 3, lambda i, j: sp.diff(n, X3[i], X3[j]) - sum(G3[k][i][j] * dn[k] for k in range(3)))
dnn = sp.Matrix(3, 3, lambda i, j: dn[i] * dn[j])

gip = gi.applyfunc(ev)
lhs = gip * sp.Matrix(3, 3, lambda i, j: R4[i + 1, j + 1])          # gamma^{-1} R4_ij (mixed)
rhs = (S.diff(t) + 2 * ndot * S).applyfunc(ev) + gip * R3 \
    + (gi * (-2 * DDn + 2 * dnn)).applyfunc(ev)
D = (lhs - rhs).applyfunc(sp.cancel)
tr = sp.cancel(sum(D[i, i] for i in range(3)))
for i in range(3):
    for j in range(3):
        r = sp.cancel(D[i, j] - (tr / 3 if i == j else 0))
        assert r == 0, (i, j, r)
print("[ok] (T): trace-free part of gamma^{-1} R4_ij equals TF[S' + 2n'S + R + -2DDn + 2DnDn] (generic 2-jets)")

# sigma.sigma cancellation identity used in the derivation
Sd = (gi * sig.diff(t)).applyfunc(ev) - (S.diff(t)).applyfunc(ev) - 2 * (S * S).applyfunc(ev)
assert all(sp.cancel(e) == 0 for e in Sd)
print("[ok] gamma^{-1} sigma' = S' + 2 S S  (generic 2-jets)")

# trace of S vanishes exactly because det gamma is t-independent
assert sp.cancel(ev(sum(S[i, i] for i in range(3)))) == 0
print("[ok] tr S = 0 exactly")
print("UNIFIED EVOLUTION FORM PROVED")
