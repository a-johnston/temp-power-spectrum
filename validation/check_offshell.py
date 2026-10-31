"""Exact local offshell ADM/Bianchi checks, with no gravitational constraints imposed."""
import unittest
import sympy as s


class OffshellTests(unittest.TestCase):
    def zero(self, value):
        self.assertEqual(s.simplify(s.expand(value)), 0)

    def test_mixed_adm_to_metric_reduction(self):
        # At a spatial orthonormal point, arbitrary symmetric M and covariant
        # V=Hessian(chi)+dchi*dchi-R suffice to check the matrix identity.
        H, chi = s.symbols('H chi', real=True)
        M = s.Matrix(3, 3, lambda i, j: s.Symbol(f'M{min(i,j)}{max(i,j)}'))
        V = s.Matrix(3, 3, lambda i, j: s.Symbol(f'V{min(i,j)}{max(i,j)}'))
        I = s.eye(3)
        cp = s.trace(M) / 3  # Fluid determinant identity, not a gravity constraint.
        B = H * I + M
        # Direct mixed ADM with H'=-H² and 8*pi*G*N²*rho/3=H² exp(-2chi).
        Bp = (H + cp - s.trace(B)) * B + s.exp(2 * chi) * V + H**2 * s.exp(-2 * chi) * I
        Mp = Bp + H**2 * I
        qpp = 4 * M * M + 2 * Mp
        proposed = (4 * M * M - 4 * (H + cp) * M
                    + 2 * (H**2 * (s.exp(-2 * chi) - 1) - 2 * H * cp) * I
                    + 2 * s.exp(2 * chi) * V)
        for value in qpp - proposed:
            self.zero(value)

    def test_offshell_spatial_einstein_residual(self):
        # F_ij=E_ij-(gamma_ij/2) tr_4 E=0 is what ADM evolution enforces,
        # not E_ij=0 when the Hamiltonian constraint has not been imposed.
        En, spatial_trace = s.symbols('En spatial_trace')
        solution = s.solve(spatial_trace - s.Rational(3, 2) * (-En + spatial_trace),
                           spatial_trace)[0]
        self.zero(solution - 3 * En)
        self.zero((-En + solution) / 2 - En)

    def test_coordinate_bianchi_produces_constraint_propagation(self):
        # Spatial Riemann normal coordinates at a point: gamma_ij=delta_ij,
        # spatial dgamma=0, but gamma'_ij=-2 N K_ij is arbitrary.
        N = s.Symbol('N', positive=True)
        Nt = s.Symbol('Nt')
        Ni = s.symbols('N1:4')
        K = s.Matrix(3, 3, lambda i, j: s.Symbol(f'K{min(i,j)}{max(i,j)}'))
        Hc, Ht = s.symbols('Hc Ht')
        Hi = s.symbols('H1:4')
        M = s.symbols('M1:4')
        Mt = s.symbols('Mt1:4')
        divM = s.Symbol('divM')
        G = [[[s.S.Zero for _ in range(4)] for _ in range(4)] for _ in range(4)]
        G[0][0][0] = Nt / N
        for i in range(3):
            G[0][0][i + 1] = G[0][i + 1][0] = Ni[i] / N
            G[i + 1][0][0] = N * Ni[i]
            for j in range(3):
                G[0][i + 1][j + 1] = -K[i, j] / N
                G[i + 1][0][j + 1] = G[i + 1][j + 1][0] = -N * K[i, j]
        # E^mu_nu after ADM evolution, but with arbitrary constraint residuals.
        E = s.zeros(4)
        E[0, 0] = -Hc / 2
        for i in range(3):
            E[0, i + 1] = M[i] / N
            E[i + 1, 0] = -N * M[i]
            E[i + 1, i + 1] = Hc / 2
        partial = [-Ht / 2 - N * divM - sum(Ni[i] * M[i] for i in range(3))]
        partial += [Mt[i] / N - M[i] * Nt / N**2 + Hi[i] / 2 for i in range(3)]
        divergences = []
        for nu in range(4):
            divergences.append(partial[nu] + sum(
                G[mu][mu][a] * E[a, nu] - G[a][mu][nu] * E[mu, a]
                for mu in range(4) for a in range(4)))
        expected_Ht = -2 * N * divM + 2 * N * s.trace(K) * Hc - 4 * sum(
            M[i] * Ni[i] for i in range(3))
        self.zero(divergences[0] + (Ht - expected_Ht) / 2)
        for i in range(3):
            expected_Mt = -N * Hi[i] / 2 + N * s.trace(K) * M[i] - Hc * Ni[i]
            self.zero(divergences[i + 1] - (Mt[i] - expected_Mt) / N)


if __name__ == '__main__':
    unittest.main(verbosity=2)
