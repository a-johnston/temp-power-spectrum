# Handoff: reviewed report and next analytical targets

## Current delivery state

The final-review corrections and structural revision have been applied. The report is
`report/temperature_spectrum_report.pdf`, rebuilt from the current TeX: **74 pages**.
The abstract/problem are on page 1, the linked results summary and reading paths on page 2,
and the contents follow. “What remains unevaluated” starts on page 46 and ends at the top of
page 47, where the appendices begin. Length was not reduced by deleting exploratory mathematics:
most of the reduction in the main text comes from moving detailed machinery and verification.

### Pre-presentation review (2026-11-03)

No mathematical error was found in the printed results. The 17 exact checks below were rerun
after the edits and pass. Printed formulas were also retyped from the TeX and tested separately
from the repository code: a multilinear Wick calculation from the three-tensor map (all `k0^6`
polynomials, the general-moment formulas, `28/3`, the four-piece split, and the linear relation
between the stationary, centred and cross formulas); `t^2`/`t^4` data for the pair sources and
averaged cubic sources against the exact-in-zeta curvature tensors; the acoustic `F_S` and
Kohri–Terada `F_T`; `T_l` table values and node positions; and the curvature identities on a
generic seed. These are now in `solution/check_printed_formulas.py`, also called by `run_all.sh`.
Eq. (1.8) of Mirbabayi–Zaldarriaga matches the quoted `phi_2`.

Edits: amplitudes renamed `S_2,T_2 -> F_S,F_T` (they collided with the moments and the transfer
function), propagation-adjoint vector `W -> Y`, the Hankel sign `sigma -> varsigma` in the `T_l`
proposition (shear and variance keep the standard `sigma`), and `D -> E = o^4-e^4` in the technical appendix.
Corrected a stale `n = ln N`, the false claim that §5.3.1 follows the observable, the Cotton
attribution, the ledger's undefined “L1–L5” and a wrong numerical percentage range. Added the
explicit `R_l`–`J_abc` link, a short derivation of the three-tensor coefficients, the named
verification rows in the requirements table, and two entries in the results list.

The saved kernels `solution/results/superhorizon_kernels_*.pkl` are now tracked, so
`curvature_kernels.py --compare` works from a fresh clone without the 13-hour regeneration.

Not verified in this review: arbitrary-time cubic kernels beyond the existing scope, the numerical
runs (retained results, apart from `C11` and node values), and rendered page appearance.

This handoff describes the revised state. The earlier session history and handoff are in git
(the preceding tip was `de9ac7d`).

### Corrections applied

- **Degree-six signs:** the cross term and total in the observer-dilation subsection now use
  `delta = o^2 - e^2`, not `(e-o)(e+o)`, with the corresponding positive prefactors. They agree
  with the independently derived limits section. New equation labels are `eq:obs-cross-six`
  and `eq:obs-total-six`.
- **Sky weights:** the canonical convention is `w = (15/2) n.C.n`, `W = w/(4pi)`, and
  `J = <w Theta> = integral W Theta dOmega`. Full-solid-angle adjoints use **W**, not w.
  The source, normalization multiplier, frame covector and boundary traces were changed together.
  Older working notes may use different measures; do not copy their weights implicitly.
- **Initial-event identity:** the derivative with a `-1` term holds at fixed background
  observation time `bar t_o`; at fixed elapsed proper time it has no `-1`. The shell
  Fermi–Walker map agrees with the symmetric-frame clock-composition identity through cubic
  amplitude, not to all orders. Multipole norms remain invariant under frame rotation exactly.
- **Generated infrared power:** the remark now assumes regular short-mode input (smooth compact
  spectral support away from zero) and uses bounded convolution integrals at fixed finite time.
  It does not infer compact real-space support from entireness or arbitrary-order regularity
  from `integral xi^2 < infinity`. The monochromatic shell-overlap exception is retained,
  explicitly exemplified by the lapse rather than every observable.
- **Higher-order projection claims:** the factor-two identities are exact through the displayed
  orders. Numerical fits suggesting failure beyond them are confined to the numerical appendix
  and are not promoted to analytical conclusions.
- **Numerical comparison:** the shell recomputation is compatible with the older table's rounding
  intervals, not equal to every printed digit. The newer `x_o=3` mu quadrature has no separate
  refinement test recorded.

### Report layout and authoritative sources

| Content | Source |
|---|---|
| Experiment, fields, map, spectra, curvature multipoles and limits | `temperature_spectrum_report.tex` |
| Contracted cubic sources and preferred acoustic scalar amplitude | `cubic_source_reduction.tex`, `shell_acoustic.tex` |
| Top-order nine-moment adjoint | `cubic_adjoint_reduction.tex` |
| Combined cubic formula first, constraints, stationary identity, evaluated consequences | `stationary_cubic_reduction.tex` |
| Requirements/status checklist | `requirements_appendix.tex` |
| Alternative entire-function Cartesian moment formulas | `scalar_moment_appendix.tex` |
| Detailed field adjoints, normalization, frame, boundary traces, Green blocks | `cubic_technical_appendix.tex`, including `propagation_adjoint.tex` |
| Prospective longitudinal-gauge route | `research_directions.tex` |
| Exact verification ledger | Appendix in `temperature_spectrum_report.tex` |
| Unified numerical appendix | `numerical_appendix.tex`, including `validation_appendix.tex` |

The next-gradient quadrupole now has its own subsection, `sec:next-gradient`, rather than being
nested under the leading octupole. Existing mathematical labels have been retained where possible.

## Validation of the revision

Added `solution/check_report_final_review.py`, also called by `solution/run_all.sh`. It:

1. Reads the two printed degree-six formulas from TeX and compares their actual signs and
   coefficients with the hand-map reference; sign-flip controls fail.
2. Checks generic STF sphere moments, conversion between w and W, the quadrupole objective and
   the normalization multiplier; an unconverted full-angle weight fails.
3. Differentiates both fixed-label clock conventions exactly.
4. Verifies both background Bianchi identities on arbitrary off-shell lapse and STF spatial
   fields. This resolves the specific off-shell linear-operator concern in the earlier handoff.
5. Checks the report input graph, duplicate/missing labels, the summary-before-contents ordering,
   and the independently navigable next-gradient subsection.

The following **17 exact checks passed** in the revision session:

```
check_report_final_review.py
check_report_identities.py
check_observer_centred_kappa4.py
check_observer_dilation_k6.py
check_centred22_mechanism.py
check_cubic_curvature_hand.py
curvature_kernels.py --spectrum
check_constraint_rows.py
check_propagation_adjoint.py
check_unified_adjoint.py
check_first_field_adjoint_source.py
check_stationary_cubic.py
check_stationary_quadrupole_four_gradients.py --symbolic
check_shell_cubic_sources.py --compare-engine
check_observable_expansion.py
check_gradient_octupole.py
check_logmap_decomposition.py
```

The new script passes Ruff and ty checks. The PDF builds in three `pdflatex` passes with no
LaTeX warnings, undefined references or overfull boxes. Opening pages were rendered and inspected;
internal PDF link destinations were checked. Temporary build/check logs are under
`/tmp/david-report-revision/` and are not repository artifacts.

**Not rerun:** the full long-duration `run_all.sh`, legacy cubic generation, finite-wavelength
numerical evolution or quadrature-refinement studies. Existing numerical values in the appendix
are retained results, not new runs. No finite-wavelength spectrum evaluation was added by this
editorial revision.

## Established results relevant to the next investigation

### Combined cubic response

For any triple with nonzero total momentum,

```
J_abc = J_111,abc - sum_{a|bc} S_meas,a[y_bc] + J_lin[y_abc].
```

The second-order fields are closed form. `J_lin` is an acoustic lapse-source pairing, a local
momentum-source term and a TT wave-source pairing. Third-order field symbols denote explicit
source integrals, not an independent field solve. No first field adjoint or second ray is
needed. Every distributed term has at most two independent times except `J_111`, which has three.
At total zero momentum use the regular Cartesian zero-momentum kernel.

Source compatibility follows from the Bianchi identity with lower orders on shell:

```
c_Ham' - 4 c_N/t = 2 i K.c_M
c_M' + 2 c_M/t = i K.c_T - (i K/6)(c_Ham + 4 c_N).
```

The scalar/vector response uses the momentum constraint; the initial Hamiltonian row is essential.
The shell scalar amplitude is a single acoustic Green integral, elementary terms and one Cin.
The forced homogeneous tensor cancels from the entire isotropically contracted cubic response;
the homogeneous lapse and its initial row do not.

### Observer-centred shell coefficient and degree-six split

The exact identity is

```
K3c = K3 + sum_pairs D K2 + sum_a D^2 K1.
```

Its shell quadrupole average starts at degree four:

```
<K3c(k,q,-q)> = k^4 [4(e-o)(12e^3 - 8e^2 o + 7e o^2 + 7o^3)
                    + 5(o^4-iota^4)]/270 + O(k^6).
```

The degree-two term vanishes identically. The first bracket is propagation; the second is the
observer-clock contribution. The full degree-six spectrum splits into dilation (`18` times the
linear degree-six coefficient), centred 22 (`28/3` times it), centred 13, and the cross term.
The initial-event dependence occurs only in centred 13 and the cross term.

The projected degree-four centred quadratic kernel equals the projection of
`-(D/2) deltaK1`, where `deltaK1 = K1(k+q)-K1(k)-K1(q)`. Its mechanism is exact:

- endpoint, clock, normalization and dilation-subtraction pieces are orthogonal to the
  degree-two polarization template;
- local four-gradient shear and the local linear-transfer correction give one factor-two identity;
- field phase, ray position and direction, in overlap ratio **7:1:14**, give the other for
  free streaming.

These identities are independently checked by `check_centred22_mechanism.py`. A conceptual
geometric/polarization derivation is still missing. Numerical fits suggest that neither the
proportionality nor the projected relation extends indefinitely, but their next analytical
failure has not been computed. Do not assume the squared four-gradient map supplies degree eight:
it also needs the two-/six-gradient interference.

## Recommended next analytical targets

1. **Explain the two factor-two identities structurally.** Start from leading-gradient locality
   (`K2c` is pure gradient-squared), the polarized degree-four linear transfer, and the split into
   local shear versus free streaming. A Leibniz/symmetric-pairing identity is a possible mechanism,
   not an established explanation. Preserve both time structures separately.
2. **Reduce the finite-wavelength observer-centred cubic shell average.** Insert the combined
   formula into the exact dilation identity. Compare phases and boundary terms before expanding
   radial functions. Derivatives of the linear/quadratic response are resolvent derivatives,
   which may combine with `J_lin` and `S_meas`. The exact dilation channel with two zero momenta
   supplies controls, not the full nonzero-radius shell contraction.
3. **Evaluate the shell-averaged direct measurement `J_111`.** The internal angular average of
   two opposite shell phases is `sin(k(w-v))/(k(w-v))`. The shell's temporal commutators and frame
   term vanish, leaving a concrete low-dimensional integral target. Keep emission, clock and
   separate random-monopole products when combining it with the Hessian sector.
4. **Compute degree-eight quadrupole / next octupole coefficients.** Extend the curvature map
   to six/five gradients, respectively, or the stationary gradient expansion. Use these to test
   the projection mechanism analytically rather than inferring coefficients from fits.
5. **Longitudinal-gauge late-time route and noncompact-spectrum bounds.** These remain larger
   projects; see the report's open section and the prospective-route appendix. Resonant angles
   require uniform estimates or matched asymptotics before taking late-time limits.

### Numerical targets, supplementary only

The forward/assembled comparison reported agreement at the 3e-9 level for three sampled triples
(noncoplanar unequal, exact zero pair, near-zero pair) at `|k| o` around 3. The shell ratios
`(C13/3)/(sigma^2 C11)` are about `2.8955`, `-6.156`, `-67.29` at `x_o=0.3,1,3` on the appendix's
fixed-ratio trajectory. The centred parts at the first two points are `-0.195` and `-2.96`.
Use these as numerical comparisons after deriving an expression, not substitutes for derivation.

### External coefficient discrepancy

A comparison of an early copy of this report by another group's agent reported a 2.2% difference
in the observer-local degree-six 13 coefficient. The problem author already knows. The cause is
unknown because the other report has not been compared term by term. Three symbolic routes and
`check_printed_formulas.py` all reproduce the printed coefficient. If it is ever compared, look at
the `iota^4` term, both gradient orderings, the physical initial event, clock composition,
normalization and dilation conventions.

## Reproduction

```bash
uv sync
.venv/bin/python -B solution/check_report_final_review.py
.venv/bin/python -B solution/check_observer_dilation_k6.py
.venv/bin/python -B solution/check_centred22_mechanism.py
.venv/bin/python -B solution/curvature_kernels.py --spectrum
.venv/bin/python -B solution/check_constraint_rows.py
solution/run_all.sh                  # full exact suite; long
cd report
pdflatex -interaction=nonstopmode -halt-on-error temperature_spectrum_report.tex
pdflatex -interaction=nonstopmode -halt-on-error temperature_spectrum_report.tex
pdflatex -interaction=nonstopmode -halt-on-error temperature_spectrum_report.tex
```

Optional numerical runs (typically 10–25 minutes, depending on grids):

```bash
.venv/bin/python -B solution/numerics/compare_cubic_quadrupole.py
.venv/bin/python -B solution/numerics/observer_centred.py
.venv/bin/python -B solution/numerics/shell_c13.py 1.0 12
.venv/bin/python -B solution/numerics/check_kappa4_shell.py 0.5
```
