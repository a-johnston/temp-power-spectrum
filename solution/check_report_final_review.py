"""Exact regressions for the final manuscript's signs, weights and fixed-label identities.

Reads the printed labelled equations, rather than duplicating their coefficients in a test.
The independent four-gradient hand map supplies the degree-six spectrum reference.
Also checks the off-shell linear Bianchi identities, sky-weight conversion, clock derivatives,
and the report's input/reference graph. No quadrature or numerical evolution is used.
"""

import re
from pathlib import Path

import sympy as s
from sympy.parsing.sympy_parser import implicit_multiplication_application, parse_expr, standard_transformations

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "report"


def equation(text, label):
    marker = rf"\label{{{label}}}"
    assert text.count(marker) == 1, label
    return text.split(marker, 1)[1].split(r"\end{equation}", 1)[0]


def printed_six(text, label, e, o, ti):
    """Read the fraction and time polynomial of the two printed degree-six formulas."""
    rhs = equation(text, label).split("=", 1)[1].strip()
    match = re.fullmatch(
        r"([+-]?)\\frac\{(\d+)\\pi\\delta\}\{(\d+)\}\s*\\big\((.*?)\\big\)[,.]?\s*",
        rhs,
        re.S,
    )
    assert match, f"Unrecognized manuscript formula: {label}"
    sign, numerator, denominator, polynomial = match.groups()
    polynomial = polynomial.replace(r"\iota", "ti").replace("^", "**")
    poly = parse_expr(
        polynomial,
        local_dict={"e": e, "o": o, "ti": ti},
        transformations=standard_transformations + (implicit_multiplication_application,),
    )
    return (-1 if sign == "-" else 1) * s.Rational(numerator, denominator) * s.pi * (o**2 - e**2) * poly


def check_printed_signs():
    from check_observer_centred_kappa4 import J3_2, J3_4, K1_2, K1_4, K2_4, e, k, o, ti

    text = (REPORT / "temperature_spectrum_report.tex").read_text()
    assert r"$\delta=o^2-e^2=-(e-o)(e+o)$" in text
    norm = 4 * s.pi / 25
    j12, j14 = K1_2.subs(k, 1), K1_4.subs(k, 1)
    expected_cross = -12 * norm * j12 * (K2_4.subs(k, 1) + 4 * j14)
    quadratic = (
        -8 * s.pi * (o**2 - e**2) * (105 * ti**4 + 719 * o**4 - 1920 * o**2 * e**2 + 2560 * o * e**3 - 1464 * e**4) / 637875
    )
    expected_total = quadratic + norm * (j12 * J3_4.subs(k, 1) + j14 * J3_2.subs(k, 1))
    for label, expected in [("eq:obs-cross-six", expected_cross), ("eq:obs-total-six", expected_total)]:
        actual = printed_six(text, label, e, o, ti)
        assert s.expand(actual - expected) == 0, label
        assert s.expand(-actual - expected) != 0, "sign-flip control must fail"
    print("[ok] printed cross and total degree-six formulas match the hand-map reference; sign controls")


def sphere_mean(expr, n):
    result = s.S.Zero
    for powers, coefficient in s.Poly(s.expand(expr), *n).terms():
        if any(p % 2 for p in powers):
            continue
        result += coefficient * s.prod(s.factorial2(p - 1) for p in powers) / s.factorial2(sum(powers) + 1)
    return s.expand(result)


def stf_symbols(prefix):
    a, b, c, d, e = s.symbols(f"{prefix}1:6")
    return s.Matrix([[a, b, c], [b, d, e], [c, e, -a - d]])


def check_weights_and_clock():
    n = s.Matrix(s.symbols("nx ny nz"))
    C, S = stf_symbols("c"), stf_symbols("s")
    N = n * n.T - s.eye(3) / 3
    w = s.Rational(15, 2) * (n.T * C * n)[0]
    W = w / (4 * s.pi)
    theta = (n.T * S * n)[0]
    mean = lambda expr: sphere_mean(expr, n)
    assert (w * N).applyfunc(mean) == C
    assert (4 * s.pi * W * N).applyfunc(mean) == C
    J = mean(w * theta)
    assert s.expand(J - s.trace(C * S)) == 0
    assert (s.Rational(15, 2) * theta * N).applyfunc(mean) == S
    multiplier = W * theta - J / (4 * s.pi)
    assert mean(multiplier) == 0
    assert (4 * s.pi * w * N).applyfunc(mean) != C, "unconverted full-angle weight must fail"

    technical = (REPORT / "cubic_technical_appendix.tex").read_text()
    for correct in [
        r"C_K(t)=\int\dd\Omega\,W\,\TF(nn)",
        r"\lambda_I=(1+\Theta)\Big[W-\frac{J}{4\pi}\Big]",
        r"\Lambda_{H,a}(o^-)=\frac{D_a}2\int WN_n\,\dd\Omega=\frac{D_a}{2}\mathsf C",
    ]:
        assert correct in technical, f"Manuscript weight convention changed: {correct}"
    # Differentiate the physical roots, treating lapse values and derivatives as independent jets.
    bar_i, t_i, t_o, nu_i, nu_o, nu_prime_i = s.symbols("bar_i t_i t_o nu_i nu_o nu_prime_i")
    dt_i = s.exp(-nu_i) / (1 + t_i * nu_prime_i)
    at_fixed_o = (t_i * s.exp(nu_i) * dt_i - bar_i) / (t_o * s.exp(nu_o))
    at_fixed_tau = t_i * s.exp(nu_i) * dt_i / (t_o * s.exp(nu_o))
    # Use bar_i = t_i exp(nu_i), without changing the independent derivative jet.
    reported = bar_i / (t_o * s.exp(nu_o)) * (s.exp(-nu_i) / (1 + t_i * nu_prime_i) - 1)
    assert s.simplify((at_fixed_o - reported).subs(bar_i, t_i * s.exp(nu_i))) == 0
    assert s.simplify(at_fixed_tau - at_fixed_o - bar_i / (t_o * s.exp(nu_o))) == 0
    print("[ok] averaged/full-angle weights, quadrupole contraction, normalization adjoint and clock labels")


def check_off_shell_bianchi():
    t = s.Symbol("t", positive=True)
    K = s.Matrix(s.symbols("kx ky kz"))
    z = K.dot(K)
    nu = s.Function("nu")(t)
    a, b, c, d, e = [s.Function(name)(t) for name in "abcde"]
    A = s.Matrix([[a, b, c], [b, d, e], [c, e, -a - d]])

    def tf(M):
        return M - s.eye(3) * s.trace(M) / 3

    LN = 3 * s.diff(nu, t, 2) + z * nu - 6 * nu / t**2
    LM = s.I * A.diff(t) * K / 2 - 2 * s.I * K * (s.diff(nu, t) - nu / t)
    LT = (A.diff(t, 2) + 2 * A.diff(t) / t) / 2 + tf(z * A - K * (A * K).T - (A * K) * K.T) / 2 + 2 * tf(K * K.T) * nu
    LH = 4 * z * nu + 12 * s.diff(nu, t) / t + 12 * nu / t**2 - K.dot(A * K)
    assert s.simplify(s.diff(LH, t) - 4 * LN / t - 2 * s.I * K.dot(LM)) == 0
    assert s.simplify(LM.diff(t) + 2 * LM / t - s.I * LT * K + s.I * K * (LH + 4 * LN) / 6) == s.zeros(3, 1)
    assert s.simplify(s.diff(LH, t) - 3 * LN / t - 2 * s.I * K.dot(LM)) != 0
    print("[ok] both linear Bianchi identities annihilate arbitrary off-shell fields; coefficient control")


def check_report_graph():
    visited = set()

    def load(path):
        assert path not in visited, f"Repeated/cyclic input: {path.name}"
        visited.add(path)
        text = re.sub(r"(?<!\\)%[^\n]*", "", path.read_text())
        return re.sub(r"\\input\{([^}]+)\}", lambda m: load(REPORT / m[1]), text)

    text = load(REPORT / "temperature_spectrum_report.tex")
    labels = re.findall(r"\\label\{([^}]+)\}", text)
    assert len(labels) == len(set(labels)), "duplicate labels"
    refs = set(re.findall(r"\\(?:eqref|ref|pageref)\{([^}]+)\}", text))
    assert not refs - set(labels), f"Missing labels: {sorted(refs - set(labels))}"
    assert text.index(r"\subsection*{Results at a glance}") < text.index(r"\tableofcontents")
    assert r"\subsection{Next-gradient quadrupole and general spectra}\label{sec:next-gradient}" in text
    print(f"[ok] report graph: {len(visited)} inputs, {len(labels)} unique labels, all references resolved")


if __name__ == "__main__":
    check_printed_signs()
    check_weights_and_clock()
    check_off_shell_bianchi()
    check_report_graph()
    print("ALL FINAL-REVIEW CHECKS PASSED")
