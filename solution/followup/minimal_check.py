"""Stand-alone: the K->0 pair fields solve the plane-symmetric equations (q=1, x=t/sqrt3)."""

# fmt: off
import sympy as sp

t, y, e, z1, z2 = sp.symbols("t y e z1 z2", positive=True)
x = t / sp.sqrt(3)
j0, c = sp.sin(x) / x, sp.cos(x)
C2 = sp.Function("C2")(t)                       # Cin(2x); C2'(t) = (1-cos 2x)/t
n1, h1 = c - j0, 6 * (1 + c - 2 * j0)
n2 = -sp.cos(2 * x) - 6 * c + 6 * j0 + 3 * sp.sin(2 * x) / x - 5 * j0**2
a2 = -78 - 36 * c - 6 * sp.cos(2 * x) + 72 * j0 + 36 * sp.sin(2 * x) / x - 24 * j0**2 + 24 * C2
g = z1 * sp.exp(sp.I * y) + z2 * sp.exp(-sp.I * y)
Z = e * g
f = e * h1 * g + e**2 * a2 * z1 * z2             # H = f diag(2/3,-1/3,-1/3); gamma = e^{2Z} exp H
nu = e * n1 * g + e**2 * n2 * z1 * z2
D = lambda u, v: sp.diff(u, v).subs(sp.Derivative(C2, t), (1 - sp.cos(2 * x)) / t)
fy, ny, zy = D(f, y), D(nu, y), D(Z, y)
ex = sp.exp(-2 * f / 3 - 2 * Z)
Ff = ex * (-fy**2 / 18 + fy * ny / 3 + fy * zy / 6 + D(fy, y) / 6 + 2 * ny**2
           + 4 * ny * zy - 2 * D(ny, y) + zy**2 - D(zy, y))
ft, nt = D(f, t), D(nu, t)
lap = ex * (D(ny, y) - (zy + fy / 3) * ny + 2 * (zy - fy / 6) * ny)
T_ = D(ft, t) + 2 * (1 / t + nt) * ft + 2 * Ff                                # (T)
N_ = (3 * D(nt, t) - lap - 2 * ex * ny**2 + ft**2 / 6
      + 3 / t**2 * (sp.exp(-2 * nu) - 1))                                     # (N)
for name, eq in (("T", T_), ("N", N_)):
    first = sp.simplify(sp.diff(eq, e).subs(e, 0))
    second = sp.expand(sp.diff(eq, e, 2).subs(e, 0) / 2).coeff(z1, 1).coeff(z2, 1)
    print(name, "first:", first, "pair:", sp.simplify(sp.expand(sp.expand_trig(second))))
