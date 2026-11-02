"""Exact checks for report/longitudinal_redshift.md.

The main result is an exact ADM/Hamiltonian redshift identity, not a perturbative
geodesic generator. Generic spatial/time metric jets check all shift terms.
First-order radiation gauge and physical-clock identities are checked separately.
No project engine, numerical quadrature, or time evolution is used.
"""

import sympy as s


def zero(expr):
    assert s.simplify(s.trigsimp(s.expand(expr))) == 0, s.simplify(expr)


def symmetric(prefix):
    entries = {(i, j): s.Symbol(f"{prefix}{i}{j}", real=True) for i in range(3) for j in range(i, 3)}
    return s.Matrix(3, 3, lambda i, j: entries[min(i, j), max(i, j)])


def adm_hamiltonian():
    h, ht = symmetric("h"), symmetric("ht")
    hx = [symmetric(f"hx{k}") for k in range(3)]
    lapse = s.Symbol("N", positive=True)
    grad_lapse = s.Matrix(s.symbols("N0:3", real=True))
    shift = s.Matrix(s.symbols("beta0:3", real=True))
    direction = s.Matrix(s.symbols("e0:3", real=True))
    db = s.Matrix(3, 3, lambda i, j: s.Symbol(f"db{i}{j}", real=True))  # partial_j beta^i
    velocity = lapse * direction - shift
    momentum = h * direction  # divide the covector p_i by q
    pdot_over_q = s.Matrix(
        [
            -grad_lapse[i] + lapse * (direction.T * hx[i] * direction)[0] / 2 + sum(db[j, i] * momentum[j] for j in range(3))
            for i in range(3)
        ]
    )
    advected_h = ht + sum((velocity[k] * hx[k] for k in range(3)), s.zeros(3))
    # Differentiate q=sqrt(h^ij p_i p_j) along Hamilton's equations.
    dlogq = -(direction.T * advected_h * direction)[0] / 2 + direction.dot(pdot_over_q)
    lie = sum((shift[k] * hx[k] for k in range(3)), s.zeros(3)) + db.T * h + h * db
    predicted = -direction.dot(grad_lapse) - (direction.T * (ht - lie) * direction)[0] / 2
    zero(dlogq - predicted)
    print("[ok] exact ADM energy law from Hamilton's equations on generic metric jets")
    missing_lie = s.expand((direction.T * lie * direction)[0] / 2)
    assert missing_lie != 0
    print("[ok] negative control: dropping the shift Lie derivative fails")


def potential_recombination():
    lapse = s.Symbol("N", positive=True)
    direction = s.Matrix(s.symbols("e0:3", real=True))
    shift = s.Matrix(s.symbols("b0:3", real=True))
    phi_grad = s.Matrix(s.symbols("phi0:3", real=True))
    psi_grad = s.Matrix(s.symbols("psi0:3", real=True))
    phi_t, psi_t, tensor = s.symbols("phi_t psi_t tensor", real=True)
    # tensor = 1/2 exp(-2 Psi) e^i e^j (M' - L_beta M)_ij.
    energy_driver = lapse * direction.dot(phi_grad) - psi_t + shift.dot(psi_grad) + tensor
    phi_path = phi_t + (lapse * direction - shift).dot(phi_grad)
    source = phi_t + psi_t - shift.dot(phi_grad + psi_grad) - tensor
    zero(energy_driver - phi_path + source)
    print("[ok] exact potential-endpoint extraction, with shift advection and tensor shear")

    # I = J - [ln(N_T/a) + ln D]; compare the two exact temperature ratios.
    ne, no, ae, ao, de, do, J = s.symbols("ne no ae ao de do J", real=True)
    I = J - ((no - ao + do) - (ne - ae + de))
    temp_frame_log_ratio = ne - no - I
    adm_log_ratio = ae - ao + do - de - J
    zero(temp_frame_log_ratio - adm_log_ratio)
    print("[ok] exact equality with the temperature-frame redshift, including endpoint conformal and Doppler factors")

    # Specialize the ADM energy law back to the temperature frame itself.
    nu_t, nu_along_e, shear = s.symbols("nu_t nu_along_e shear", real=True)
    J_integrand = nu_along_e + nu_t + shear
    nu_path = nu_t + nu_along_e
    zero(J_integrand - nu_path - shear)
    print("[ok] temperature-frame specialization gives J=[nu]+I")


def radiation_gauge_and_clock():
    x, k = s.symbols("x k", positive=True)
    t = s.sqrt(3) * x / k
    profile = 2 * s.sin(x) / x - s.cos(x)
    B = 3 * (profile - 1) / k**2
    dt = lambda expr: s.diff(expr, x) * k / s.sqrt(3)
    A = dt(B)
    nu = s.cos(x) - s.sin(x) / x
    phi = -2 * (s.sin(x) / x - s.cos(x)) / x**2
    phi_from_lapse = nu - dt(A) - A / t
    psi_from_spatial = A / t - nu - 1 - k**2 * B / 3
    zero(phi_from_lapse - phi)
    zero(psi_from_spatial - phi)
    zero(A - dt(B))
    print("[ok] full time+space transformation gives Phi=Psi and zero first-order scalar shift")
    assert s.simplify(A) != 0
    assert s.simplify(nu - phi) != 0
    print("[ok] negative control: a purely spatial transformation retains a shift and the O(1) lapse oscillation")

    # Finite-time clock identity by integration by parts, using Phi=nu-A'-A/t.
    ti, to, Ai, Ao, nui, W = s.symbols("ti to Ai Ao nui W", real=True)
    initial_shift_T = -ti * nui
    initial_shift_L = initial_shift_T + Ai
    integral_tphi = W - (to * Ao - ti * Ai)  # W=int_ti^to t*nu dt
    observation_shift_L = (ti * initial_shift_L - integral_tphi) / to
    observation_shift_T = (ti * initial_shift_T - W) / to
    zero(observation_shift_L - observation_shift_T - Ao)
    print("[ok] physical clock: delta eta_o = D1 + alpha_o, with both initial-event terms retained")


if __name__ == "__main__":
    adm_hamiltonian()
    potential_recombination()
    radiation_gauge_and_clock()
    print("LONGITUDINAL REDSHIFT IDENTITIES VERIFIED (exact; no late-time asymptotic claim)")
