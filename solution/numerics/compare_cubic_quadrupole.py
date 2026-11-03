"""Finite-wavelength validation of the assembled cubic quadrupole (NUMERICAL VALIDATION, not evidence).

Forward (independent): forward_cubic.ForwardMap with fields solved from (N) and the Cartesian (T)
rows, Fermi-Walker frame, physical roots, clock and exact normalisation.
Assembled (report formulas): J_111 (repository's second-ray-free direct measurement)
  - sum_{a|bc} S_meas,a[y_bc]  (no-probe-ray Hessian)  + J_lin[y_abc]  (constraint rows).

Also reported: the dilation channel (p,0,0) against (p d_p)^2 of the linear quadrupole, and
J_lin in constraint rows against the Cartesian (hyperbolic-sector) third-order field.

Run: .venv/bin/python -B solution/numerics/compare_cubic_quadrupole.py [--quick]
"""

from __future__ import annotations

import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import assembled_cubic as ac  # noqa: E402
import cubic_fields as cf  # noqa: E402
import forward_cubic as fc  # noqa: E402

E, O, IOTA = 1.2, 3.0, 0.5
C = np.array([[1.0, 0.3, 0.0], [0.3, -0.5, -0.2], [0.0, -0.2, -0.5]])
TRIPLES = {
    "noncoplanar unequal": {1: (0.9, 0.3, -0.4), 2: (-0.2, 0.7, 0.5), 3: (0.5, -0.5, 0.6)},
    "zero pair momentum": {1: (0.9, 0.3, -0.4), 2: (-0.2, 0.7, 0.5), 3: (0.2, -0.7, -0.5)},
    "near-zero pair": {1: (0.9, 0.3, -0.4), 2: (-0.2, 0.7, 0.5), 3: (0.2, -0.7, -0.45)},
}
START = time.time()


def log(msg):
    print(f"[{time.time() - START:6.0f}s] {msg}", flush=True)


def qhat(k):
    k = np.asarray(k, dtype=float)
    return np.outer(k, k) / (k @ k) - np.eye(3) / 3


def compare(name, ks):
    prof = cf.build_profiles(ks, 4.0)
    full = fc.ForwardMap(ks, prof, E, O, IOTA)
    J_fwd = full.quadrupole(full.map(), C, [1, 2, 3])
    singles = {frozenset([a]) for a in ks}
    first = fc.ForwardMap(ks, prof, E, O, IOTA, labels_on=singles)
    J111_fwd = first.quadrupole(first.map(), C, [1, 2, 3])
    J111 = ac.direct_measurement(ks, E, O, IOTA, C)
    hess = 0.0
    for a in (1, 2, 3):
        b, c = (x for x in (1, 2, 3) if x != a)
        y = prof[frozenset([b, c])]
        p = np.asarray(ks[a], dtype=float)
        q = np.asarray(ks[b], dtype=float) + np.asarray(ks[c], dtype=float)
        hess -= ac.measurement_hessian(p, q, float(np.linalg.norm(p)), qhat(p), (y.nu, y.A, y.Ap), E, O, IOTA, C)
    lin, parts = ac.linear_constraint_rows(ks, prof, (1, 2, 3), E, O, C)
    lin_cart = ac.linear_cartesian(ks, prof, (1, 2, 3), E, O, C)
    assembled = J111 + hess + lin
    rel = lambda x, y: abs(x - y) / max(abs(y), 1e-300)  # noqa: E731
    log(f"{name}: |k_b+k_c| = {np.linalg.norm(np.add(ks[2], ks[3])):.3f}")
    print(f"    forward cubic quadrupole          {J_fwd.real:+.12f}  (imag {J_fwd.imag:+.1e})")
    print(f"    assembled J111 - S_meas + J_lin   {assembled.real:+.12f}  rel diff {rel(assembled, J_fwd):.1e}")
    print(f"    J111: repository formula vs forward first-field-only   rel diff {rel(J111, J111_fwd):.1e}")
    print(f"    J_lin: constraint rows vs Cartesian third-order field   rel diff {rel(lin, lin_cart):.1e}")
    print(f"    pieces: J111 {J111.real:+.6f}, -sum S_meas {hess.real:+.6f}, lapse {parts['lapse'].real:+.6f},")
    print(f"            momentum {parts['momentum'].real:+.6f}, tensor {parts['tensor'].real:+.6f}")
    if name == "noncoplanar unequal":
        norot = fc.ForwardMap(ks, prof, E, O, IOTA, rotate=False)
        J_norot = norot.quadrupole(norot.map(), C, [1, 2, 3])
        J_noclock = full.quadrupole(full.map(clock=False), C, [1, 2, 3])
        r_rot, r_clock = rel(J_norot, J_fwd), rel(J_noclock, J_fwd)
        print(f"    control: symmetric triad instead of Fermi-Walker changes the forward value by {r_rot:.1e}")
        print(f"    control: dropping the physical clock changes it by {r_clock:.1e}")
        assert r_rot > 1e-6 and r_clock > 1e-3, "controls must be resolved"
    return rel(assembled, J_fwd), rel(J111, J111_fwd), rel(lin, lin_cart)


def dilation():
    p = np.array([0.9, 0.3, -0.4])
    ks = {1: tuple(p), 2: (0.0, 0.0, 0.0), 3: (0.0, 0.0, 0.0)}
    prof = cf.build_profiles(ks, 4.0)
    fm = fc.ForwardMap(ks, prof, E, O, IOTA)
    J3 = fm.quadrupole(fm.map(), C, [1, 2, 3])

    def J1(scale):
        m = {1: tuple(scale * p)}
        f = fc.ForwardMap(m, cf.build_profiles(m, 4.0, third=False), E, O, IOTA)
        return f.quadrupole(f.map(), C, [1])

    h = 0.01
    v = {j: J1(np.exp(j * h)) for j in (-2, -1, 0, 1, 2)}
    D2 = (-v[-2] + 16 * v[-1] - 30 * v[0] + 16 * v[1] - v[2]) / (12 * h * h)
    r = abs(J3 - D2) / abs(D2)
    log(f"dilation channel (p,0,0): forward {J3.real:+.10f}, (p d_p)^2 J1 {D2.real:+.10f}, rel diff {r:.1e}")
    return r


if __name__ == "__main__":
    results = [dilation()]
    names = list(TRIPLES)[:1] if "--quick" in sys.argv else list(TRIPLES)
    for nm in names:
        results += list(compare(nm, TRIPLES[nm]))
    worst = max(results)
    print(f"worst relative difference {worst:.1e}")
    assert worst < 1e-6, "finite-wavelength validation failed"
    print("FINITE-WAVELENGTH VALIDATION PASSED (numerical validation, not an analytical proof)")
