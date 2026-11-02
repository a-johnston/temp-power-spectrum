"""Exact checks for the shell-locality proposition of the report.

1. On the shell, (1/2) Q'(t) D_n e^{ik.X} = -sigma_{1,nn} for every |k| = k0 mode, with
   D_n = (n.d)^2 + k0^2/3; hence theta is a local line integral of g.
2. H_1 = -Q (d_i d_j + delta_ij k0^2/3) g reproduces the Fourier kernel Q (k_i k_j - delta k^2/3).
3. The second-order shear from gamma = e^{2 zeta} e^H is
   sigma^(2) = (1/2) A_2' + (1/2) (H_1^2)' + 2 g H_1'     (Taylor convention, generic matrices),
   so A_2 is the only non-local second-order field in the redshift integrand.
"""

import sympy as sp

t, k0, eps = sp.symbols("t k0 epsilon", positive=True)
X = sp.symbols("X1:4", real=True)
n = sp.symbols("n1:4", real=True)
kh = sp.symbols("khat1:4", real=True)
Q = sp.Function("Q")(t)

# 1-2: a single shell mode, |k| = k0 imposed by khat.khat = 1
k = [k0 * c for c in kh]
unit = {kh[2]: sp.sqrt(1 - kh[0] ** 2 - kh[1] ** 2), n[2]: sp.sqrt(1 - n[0] ** 2 - n[1] ** 2)}
mode = sp.exp(sp.I * sum(ki * xi for ki, xi in zip(k, X, strict=True)))
H1 = sp.Matrix(3, 3, lambda i, j: -Q * (sp.diff(mode, X[i], X[j]) + (k0**2 / 3 if i == j else 0) * mode))
kernel = sp.Matrix(3, 3, lambda i, j: Q * (k[i] * k[j] - (k0**2 / 3 if i == j else 0)) * mode)
assert sp.simplify((H1 - kernel).subs(unit)) == sp.zeros(3)
nv = sp.Matrix(n)
sigma_nn = (nv.T * kernel.diff(t) * nv)[0] / 2
Dn = (
    sum(ni * sp.diff(sum(nj * sp.diff(mode, xj) for nj, xj in zip(n, X, strict=True)), xi) for ni, xi in zip(n, X, strict=True))
    + k0**2 / 3 * mode
)
assert sp.simplify((sp.diff(Q, t) * Dn / 2 + sigma_nn).subs(unit)) == 0
print("[ok] shell first-order fields are local; theta = P0[(1/2) int Q' D_n g]")


# 3: second-order shear with generic symmetric time-dependent matrices
def symmat(name):
    f = {(i, j): sp.Function(f"{name}{i}{j}")(t) for i in range(3) for j in range(i, 3)}
    return sp.Matrix(3, 3, lambda i, j: f[min(i, j), max(i, j)])


A1, A2 = symmat("H"), symmat("A")
g = sp.Symbol("g")  # time independent
H = eps * A1 + eps**2 * A2 / 2
expH = sp.eye(3) + H + H * H / 2
gamma = sp.exp(2 * eps * g) * expH
sigma = (gamma.diff(t) / 2).applyfunc(lambda e: sp.series(e, eps, 0, 3).removeO())
sigma2 = 2 * sigma.applyfunc(lambda e: sp.expand(e).coeff(eps, 2))
target = A2.diff(t) / 2 + (A1 * A1).diff(t) / 2 + 2 * g * A1.diff(t)
assert (sigma2 - target).applyfunc(sp.expand) == sp.zeros(3)
print("[ok] sigma^(2) = A2'/2 + (H1^2)'/2 + 2 g H1' (generic matrices)")
print("SHELL-LOCALITY CHECKS PASSED (exact symbolic algebra)")
