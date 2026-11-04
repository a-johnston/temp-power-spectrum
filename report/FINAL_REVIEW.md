# Critical final review

> **Revision status:** the findings below describe the pre-revision report. The corrections and restructuring have since been applied; see [HANDOFF.md](HANDOFF.md) for the current layout, validation and open research targets. Line/page references below are historical.

## Verdict

**The report contains a strong analytical answer and several genuinely elegant side results, but I would not send the current files unchanged.** The main problems are a pair of definite sign errors in the prose, inconsistent angular-weight conventions in the stationary section, overextended claims in two remarks, and a stale delivery PDF. None of these, on the checks performed here, overturns the central long-wavelength results.

Length is not the main problem. The report repeatedly takes the reader through intermediate representations, their verification history, and their limitations before reaching the cleanest result. It should retain its exploratory mathematics while moving implementation history and superseded machinery off the main reading path.

This is a review, not a revision: the report sources and tracked PDF were left unchanged. Source line numbers below refer to the files reviewed. Page numbers refer to a fresh build of the current TeX, not the repository PDF.

## 1. Corrections needed before delivery

### 1.1 High priority: two signs are wrong in the next-order observer-dilation decomposition

**Locations:** `temperature_spectrum_report.tex:1300` and `:1339`, §8.3, “Next gradient order”; compare `:2049` in §9.2.

The subsection defines

\[
\delta=o^2-e^2,
\]

but two expressions copied from the exact check retain the sign appropriate to `(e-o)(e+o) = -delta`.

In the subsection's units, `sigma = k0 = 1`, the correct cross term is

\[
-2D_{k_0}X_2=
\boxed{+\frac{8\pi\delta}{91125}
(204e^4-280e^3o+210e^2o^2-89o^4-45\iota^4)}.
\]

The correct total is

\[
\boxed{+\frac{2\pi\delta}{1913625}
(38856e^4-64640e^3o+48480e^2o^2-18601o^4-4095\iota^4)}.
\]

Both currently have a minus sign. Restore `sigma^4 k0^6` when leaving those units. The corrected total is identical to the expression in the limits section, whose parenthesis has the opposite sign.

**Verification:** reran `solution/check_observer_dilation_k6.py`. It passes and prints the expressions with `(e-o)(e+o)`, not `delta`. The independent curvature check also passes. This is a manuscript transcription error, not a failure of those calculations. The passing checks do not currently protect the printed formulas against this error.

**Action:** fix both signs and add an assertion or small manuscript-formula check using the report's actual `delta` convention.

### 1.2 High priority: the stationary section switches between sky means and full solid-angle integrals

**Locations:** `stationary_cubic_reduction.tex:18–21`, `:81–100`, `:237–241`, and the evaluated adjoint traces beginning at `:282`.

The section defines

\[
w=\frac{15}{2}n^T\mathsf Cn,
\qquad \langle w\,\mathrm{TF}(nn)\rangle_\Omega=\mathsf C.
\]

That convention gives the averaged objective

\[
J=\langle w\Theta\rangle_\Omega=\mathsf C:S.
\]

But the background source is subsequently defined as

\[
C_K(t)=\int d\Omega\,w\,\mathrm{TF}(nn)e^{iK\cdot n(o-t)},
\]

and the normalization adjoint uses the full-integral convention. At zero phase this `C_K` is `4 pi C`, not `C`. Later the text again uses `J_a = C:S_a`, and the displayed evaluated boundary trace has the averaged normalization.

For example, the general boundary formula gives

\[
\Lambda_{H,a}(o^-)=\frac{D_a}{2}\int w\,\mathrm{TF}(nn)\,d\Omega
=2\pi D_a\mathsf C,
\]

whereas the evaluated trace is `D_a C / 2`. The mismatch is a factor of `4 pi`.

`solution/numerics/assembled_cubic.py:14` explicitly uses averaged weights. Thus agreement between the numerical implementations does not resolve the inconsistency in the printed notation.

**Action:** introduce two distinct symbols once:

\[
w=\frac{15}{2}n^T\mathsf Cn,\qquad W=\frac{w}{4\pi},
\qquad J=\langle w\Theta\rangle_\Omega=\int W\Theta\,d\Omega.
\]

Use `W` in all full-`dOmega` objectives, source pairings and multipliers, and `w` only inside sky means. Audit `C_K`, `A_a`, the normalization multiplier, the frame covector and the boundary traces together. The alternative is to convert the entire section to averaged pairings, but that requires changing the displayed normalization-adjoint formulas too.

### 1.3 Delivery blocker: the tracked PDF is not the reviewed report

**Location:** `report/temperature_spectrum_report.pdf`.

The repository PDF has **59 pages**, no “Results at a glance,” and no stationary-reduction section. A fresh three-pass build of the current source has **72 pages**, includes both, and has no undefined references or LaTeX errors.

The fresh build has two overfull boxes: about `0.33 pt` and `12.48 pt`; the latter is at source lines `2676–2679`. It is not accurate to describe the current build as having only the negligible box.

**Action:** after correcting the source, regenerate the tracked PDF and check its opening pages, links and page count. Update the stale README statement at the same time. Do not distribute the current PDF as the final result.

### 1.4 The initial-event remark needs a fixed-label convention and a narrower map-level claim

**Location:** `stationary_cubic_reduction.tex:345–363`, Remark “Where the initial event enters, exactly.”

There are two separate issues.

1. The displayed derivative of `t_o` holds at **fixed background observation time** `bar t_o`. The `-1` in the bracket comes from differentiating `(bar t_o^2 - bar t_I^2)/2`. At fixed elapsed proper time `tau_o`, that right-hand side is constant, and the `-1` is absent. Either derivative is meaningful, but the held-fixed label must be explicit.
2. The argument removes Fermi–Walker dependence only through the amplitude order needed for the cubic map. Equal first-order shell profiles imply `[H1,H1'] = 0`; they do not imply that the full nonlinear frame rotation vanishes to all orders. Mixed commutators involving `H1` and `H2` can enter the frame generator at cubic order and the map at fourth order. Therefore the asserted all-orders identity for the specified Fermi–Walker **map** does not follow.

**Action:** state the map identity through `O(epsilon^3)` and specify the fixed labels. An exact statement can instead be made for the symmetric-frame coordinate map. There is also a stronger spectrum-level observation worth preserving: a global sky rotation leaves each realization's summed multipole norm invariant, so frame dependence can be removed from `C_l` without claiming that the maps agree to all orders.

This correction does not invalidate the displayed clock dependence of the fourth-order spectrum.

### 1.5 The causal-power remark overstates what analyticity and its correlation assumption prove

**Location:** `temperature_spectrum_report.tex:2250–2272`, Remark “Causal generation of superhorizon power.”

The distinction between a finite-width short-mode spectrum and an exactly monochromatic shell is useful and should stay. The proof as written is not sufficient for its general claim:

- An entire Fourier kernel is **not equivalent** to a compactly supported real-space kernel. Appropriate exponential-type and growth bounds are needed; a Gaussian Fourier transform is an elementary counterexample to the unqualified equivalence.
- `integral xi^2 < infinity` is not the same as absolute integrability of `xi`. Even when it makes a quadratic correlation integrable, it gives continuity of its Fourier transform, not automatically smoothness or an `O(K^2)` expansion. Higher spatial moments or stronger spectral regularity are needed for that expansion.
- The claim is made for an arbitrary order-`r` field. Higher-order Gaussian products can contain lower-order Wick components. For example, the connected covariance of `g^3` contains `9 sigma^4 xi`, as well as `6 xi^3`. A condition on `xi^2` alone does not control that first term. A genuinely short-mode input must be specified rather than inferred from the displayed assumption.

**Action:** restrict the statement to suitably regular short-mode inputs, for example smooth compactly supported spectra bounded away from zero, and prove boundedness/continuity of the relevant convolution using the recombined finite-time Cartesian kernels. A finite `P_Y(0)` gives

\[
\Delta_Y^2(K)=\frac{K^3}{2\pi^2}[P_Y(0)+o(1)],
\]

which is already the intended conclusion. Claim the stronger `O(K^2)` remainder only under explicit differentiability assumptions. If invoking finite propagation, provide the appropriate argument for the physical constrained system rather than deriving it from entireness alone.

Keep the shell-overlap explanation, but identify the field in the example as the lapse response; do not silently promote it to a theorem about every observable. Calling the paragraph an argument rather than a theorem does not repair the unsupported implications.

### 1.6 Numerical counterexamples are being used as analytical conclusions in the main text

**Location:** `temperature_spectrum_report.tex:1304–1336`, §8.3.

The text says that the centered-22 proportionality “does not continue,” that it “is lost” at `k0^8`, and that the projected relation holds “through O(k0^4) only.” The supporting finite-wavelength evidence here is numerical fitting. The analytical calculation establishes the relation through the displayed order, not its failure at the next order.

The presence of `|K2c^[4]|^2` at degree eight is not itself a proof of failure: the degree-eight spectrum also contains the interference between the degree-two and degree-six kernels. The open-problems section correctly recognizes that the six-gradient map is required.

**Action:** retain the exact factor-two identities and their field/ray explanation in the main text. Move the fitted ratios and mismatches to the numerical appendix, using wording such as “suggests that the identity does not extend.” In the main text say that no higher-order proportionality has been established. Repeatedly labeling numbers “not evidence” does not prevent them from functioning as evidence when they are the reason for a categorical claim.

## 2. Structure: preserve the discoveries, remove the work log

### What is already working

The new abstract and linked results list are substantial improvements. The complete observable is built before the main monochromatic specialization; physical clocks, emission surfaces, random-monopole normalization and tensor terms are treated seriously. The distinction between exact integral representations and closed-form evaluations is generally careful. The report should not be shortened by deleting its interesting mathematical consequences.

In particular, keep prominent:

- the temperature-frame redshift identity and cancellation of the observer temperature;
- the closed linear transfer and its acoustic interpretation;
- the local curvature quadrupole and curvature-gradient octupole;
- the Gaussian resummations, `29/3 = 8 + 5/3`, `46 = 18 + 28`, and general-spectrum bounds;
- the three-tensor/four-gradient map and the single spectral moment left after observer dilation;
- the exact observer-dilation and log-map decompositions;
- the shell two-amplitude reduction, homogeneous-tensor cancellation and induced-tensor comparison;
- the acoustic constraint-row reduction and combined cubic response;
- the exact projected factor-two identities, presented as results rather than as a sequence of guesses and tests.

These are exploratory results in the problem's own mathematical space, not extraneous process.

### The first few pages should sell the results, not certify the workflow

In the fresh PDF, page 1 contains the abstract and the beginning of the contents; page 2 continues the contents and begins §1; page 3 contains “Results at a glance”; much of page 4 is a dense requirement/status table. This is usable, but the introduction spends too much of its limited space on completion bookkeeping.

Recommended minimal change:

1. Put a self-contained results summary before the full contents, ideally on pages 1–2.
2. Define `v = epsilon^2 sigma^2` next to the displayed primary spectrum. It is currently undefined in the summary itself. Explicitly label that formula as the monochromatic specialization.
3. Show the resummed leading quadrupole as `e^(8v)(1 + 5v/3)` and then its fourth-order factor `1 + 29v/3`; do not use the truncated factor as the expression for the resummation.
4. Give the octupole's `46` explicitly. Replace the unexplained list `18, 28/3` with a short statement of what those coefficients multiply.
5. Give one compact scope sentence: exact finite-time integral answer; evaluated linear transfer and gradient coefficients; no closed-form general fourth-order spectrum or fourth-order late-time limit.
6. Move the requirement/status table to the verification appendix. Retain its links, but do not make it the introduction's second summary.

Offer two linked reading paths: **results and interpretation** versus **derivation and verification**. The author should not have to infer that most of the evaluated nonlinear physics begins on page 37.

### A lecture-like main path

A good main line is:

> physical experiment → exact redshift → regular constrained fields → complete map → spectrum → evaluated curvature multipoles and limits → reduced cubic representation → open questions.

This mostly preserves the present organization and voice. The changes needed are local and selective:

- **§5:** retain the triangular solution, regularity, shell amplitudes and the clean acoustic scalar formula. Move the superseded hyperbolic scalar-moment construction and two-variable kernel to a technical appendix. Introduce the simplest representation first, not several pages before explaining how it can be bypassed.
- **§6:** retain the complete map and boundary formulas. Present the propagation-adjoint result and its proof concisely; place its detailed coefficient bookkeeping with the cubic-reduction machinery.
- **§8.3:** retain the exact decomposition and the next-gradient mechanism. Remove the fit narrative from between the theorem and its interpretation.
- **§9:** create a separately numbered subsection called “Next-gradient quadrupole and general spectra.” The complete `k0^6` quadrupole, three-tensor formula and general-spectrum bounds are currently buried under “Leading octupole” beginning at source line 1957. This is a genuine navigation defect, not merely a length preference.
- **§10:** lead with the combined cubic formula and explain its three terms. Then give the constraint compatibility and the short stationary identity as the explanation. Move the full first-adjoint source, Green reconstruction and lengthy evaluated Green blocks to a technical appendix. Preserve the boundary terms there; moving them is not permission to omit them.
- **§11:** reduce to a short list of genuinely unevaluated mathematical tasks. It currently repeats most of the report's representations, verification scopes and proposed routes. The detailed research plan belongs in the handoff/working notes.
- **Numerics:** merge the two numerical appendices into one with distinct “implementation comparisons” and “exploratory observations” subsections. Keep the useful discoveries and enough reproducibility information; avoid retelling the chronology of every attempted grouping.

Do not impose an arbitrary total-page target. A long technical appendix is compatible with a clear, rewarding main text.

### Editorial rule for the whole report

A main-text paragraph should do at least one of these: state a result, explain why it is true, explain what it means, or delimit its applicability. Runtime comparisons, fixture descriptions, sampling grids and lists of checks belong in the ledger or numerical appendix.

Examples to move rather than delete include the 252-unknown reconstruction details (§8.1), the 13-hour versus 3-second implementation history (§9.2), off-shell fixture narratives (§10), and the detailed longitudinal-gauge research plan (§9.5). Preserve the underlying reconstruction idea and exact longitudinal-redshift identity.

## 3. Smaller corrections and consistency issues

- `stationary_cubic_reduction.tex:126`: “a homogeneous probe” must be **a homogeneous scalar/lapse probe, `q=0, j=0`**. A general homogeneous tensor probe still contributes `-a p^T j p + h':j'/2` to `C_N`. The supporting working note states the restriction correctly.
- `temperature_spectrum_report.tex:2696`: the two numerical implementations do not agree “to the printed precision” in every row: `-6.156` versus `-6.13`, and `-67.29` versus `-67.5`. They may be compatible with the uncertainty induced by the older table's rounded inputs. Say that and give an uncertainty/rounding bound, rather than claiming digit agreement. The newer `x_o=3` angular quadrature also lacks its own reported refinement test.
- `temperature_spectrum_report.tex:1600`: “up to five nested times” describes the unreduced representation. Explicitly point to the reduced three-time bound, otherwise the summary and later result appear inconsistent.
- Distinguish “no independent third-order field solve” from “no third-order field.” The combined formula still displays the lapse and TT responses as explicit source integrals. The computational elimination is valuable without describing those response functionals as absent.
- Standardize notation for the physical initial time (`i`, `iota`, `bar t_I`), and avoid reusing `A_r` for both metric kernels and ray acceleration without a local warning.
- Replace “numerical, not evidence” by **“numerical validation; not used in the analytical derivation”** where appropriate. The latter accurately states the intended role and avoids sounding as though a validation is both being cited and disavowed.
- The README's main-result list is no longer the clearest account of the latest cubic reduction; the newer result is buried in “Further cancellation research.” Promote the final representation and keep research chronology in `HANDOFF.md`.

## 4. Verification performed and limits of this review

Read the problem statement, handoff, README, the complete main TeX source, and all five included TeX files. Examined the generated opening pages and compared the current PDF with a fresh build.

Reran these exact checks; all passed:

| Check | Approximate runtime |
|---|---:|
| `check_observer_centred_kappa4.py` | 2 s |
| `check_observer_dilation_k6.py` | 4 s |
| `check_centred22_mechanism.py` | 3 s |
| `check_cubic_curvature_hand.py` | 7 s |
| `check_report_identities.py` | 1 s |
| `check_constraint_rows.py` | 61 s |

Also independently substituted arbitrary time-dependent `nu` and a generic symmetric trace-free `A` into the two **linear off-shell Bianchi identities**. Both vanish identically. This addresses the particular handoff concern that the background Bianchi operator's annihilation of the linear operator had not been checked separately; I did not find a defect in that step of the compatibility proof.

The combined-response regrouping is consistent with the third derivative chain rule and the stationary identity as written. Its numerical comparisons are supplementary, not needed to justify that algebraic regrouping. The angular normalization still needs the explicit repair described above.

Built the report with three `pdflatex` passes. Review build and extracted/rendered pages are under `/tmp/david-final-review-build/`; targeted check logs are under `/tmp/david-final-review-checks/`.

**Not performed:** the full `solution/run_all.sh`, the long legacy generator runs, new finite-wavelength numerical evolution/convergence studies, or an independent literature audit. Passing the selected checks is not certification of every stored kernel or every claim in the report.

## Recommended final pass

1. Correct the two signs and unify the angular normalization.
2. Repair the initial-event and causal-power remarks; separate numerical suggestions from analytical conclusions.
3. Make the results summary self-contained and move the status table out of the opening.
4. Give the next-gradient quadrupole its own navigable subsection and foreground the final cubic representation.
5. Move process detail and superseded machinery to appendices without deleting the useful side results.
6. Regenerate the delivery PDF, check its links and opening pages, and synchronize the README.

The desired endpoint is not a shorter research diary. It is a mathematical report whose main story is easy to follow, with a substantial reference layer underneath it.
