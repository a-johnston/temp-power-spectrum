"""How much of the temperature white-noise power is generated, using the complete Q rather than two pieces.

Averages Q^2 over a band q in (1 +- 0.3) q0, flat weight, at fixed t_o/t_* and t_I = t_*/2. The cross terms mean that the
fractions need not sum to one; 'gen' is -[W]/3 and 'em' is the emission piece.
"""

import numpy as np
from k0_closed import W_R
from map_k0 import s3, theta_pair


def W(x):
    return W_R(np.array([float(x)]))[0]


if __name__ == "__main__":
    print("x*0 r | rms Q, rms gen, rms em | <gen^2>/<Q^2>, <em^2>/<Q^2>")
    for r in (4, 32):
        for x0 in (1, 3, 5, 10, 20):
            if r * x0 > 300:
                continue
            xs = np.linspace(0.7 * x0, 1.3 * x0, 41)
            Q, G, E = [], [], []
            for x in xs:
                xo = r * x
                lam = s3 * (xo - x)
                res = theta_pair(s3 * x, s3 * xo, s3 * x / 2, nmu=int(1.5 * lam + 80), nt=int(24 * lam) + 2000, chunk=4)
                Q.append(res["Q2"])
                G.append(-(W(xo) - W(x)) / 3)
                E.append(res["Q2pieces"]["em"])
            Q, G, E = map(np.array, (Q, G, E))
            q2 = np.mean(Q**2)
            print(
                f"{x0:3} {r:3} | {np.sqrt(q2):7.3f} {np.sqrt(np.mean(G**2)):7.3f} {np.sqrt(np.mean(E**2)):7.3f}"
                f" | {np.mean(G**2) / q2:6.3f} {np.mean(E**2) / q2:6.3f}",
                flush=True,
            )
