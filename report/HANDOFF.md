# Handoff: reviewed report and next analytical targets

## Current delivery state

The final-review corrections and structural revision have been applied. The report is
`report/temperature_spectrum_report.pdf`, rebuilt from the current TeX: **73 pages**.
The abstract/problem are on page 1, the linked results summary and reading paths on page 2,
and the contents follow. The main text ends with “What remains unevaluated” on page 46;
the appendices begin on page 47. Length was not reduced by deleting exploratory mathematics:
most of the reduction in the main text comes from moving detailed machinery and verification.

`report/FINAL_REVIEW.md` records the findings against the pre-revision source; its old line/page
numbers are historical. This handoff describes the revised state. The earlier session history
and handoff are available in git (the preceding tip was `de9ac7d`).

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

An earlier external calculation reported a 2.2% difference in the observer-local degree-six
13 coefficient. The reference remains the exact expression in the limits section and in
`check_stationary_quadrupole_four_gradients.py`, including both gradient orderings. Compare the
`iota^4` term, the physical initial event, clock composition, normalization and dilation
conventions term by term. The initial-event remark now makes the held-fixed labels explicit;
it should not be used as an all-orders equality of the Fermi–Walker and symmetric maps.

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
