"""Exact plane-symmetric evolution: <f>/eps^2 against a(x), and the failure of amplitude expansion at late time.

<f> is the label average, the exact K=0 component to all orders in eps for this seed. Spectral tail monitors steepening.
"""

import numpy as np
from forward_plane import Evol
from k0_closed import a_exact

s3 = np.sqrt(3.0)
XS = (3.0, 10.0, 30.0, 100.0, 300.0)

if __name__ == "__main__":
    print("eps; then for each x: <f>/eps^2 - a(x) (a(x) in header); tail = max high-harmonic |nu_m|/|nu_1| over the run")
    print("a(x):", {x: round(float(a_exact(s3 * x)), 4) for x in XS})
    for eps in (1e-3, 3e-3, 1e-2, 3e-2):
        ev = Evol(eps, 0.05, s3 * max(XS) + 1, N=96, dt=0.01, mmax=40)
        row = []
        for x in XS:
            i = int(round((s3 * x - ev.ts[0]) / ev.dt))
            favg = ev.c[i, 0, 0].real
            if not np.isfinite(favg):
                row.append("nan")
                continue
            row.append(f"x={x:g}: {favg / eps**2 - a_exact(ev.ts[i]):+.4f}")
        finite = np.isfinite(ev.c[:, 2, 1].real)
        mag = np.abs(ev.c[:, 2])
        tail = mag[:, 30:].max(axis=1) / (mag[:, 1] + 1e-30)
        bad = np.argmax((tail > 1e-3) | ~finite) if ((tail > 1e-3) | ~finite).any() else -1
        tb = ev.ts[bad] / s3 if bad >= 0 else np.inf
        print(f"eps={eps:g}: " + "  ".join(row) + f"  | spectral tail>1e-3 first at x={tb:.1f}, eps*x={eps * tb:.2f}", flush=True)
