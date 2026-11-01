The new perspectives tidied up the hard integrals and made them cheap to compute numerically, but they didn't turn them into closed forms. The numerical exploration also argues against an elegant closed form at general k₀ for the part I could compute.

## What I built

`solution/explore_fourth_order.py` evaluates the second-order part of the fourth-order spectrum numerically at finite k₀. It reuses the exact second-order sources and the report's Green-function and map formulas, and takes about 10 s per x_o. Thanks to the log-map split, it covers everything except ⅓Re C[θ, L₃], which needs third-order fields. The superhorizon split was 17/3 for that part and 4 for the missing term, so this is the larger part there. It is exploratory only and not evidence under the problem's rules; I haven't put any of these numbers in the report.

Checks that it is computing the right thing:
- **First order.** It reproduces the closed-form T₂ (to every printed digit at moderate x_o, about 0.4% off at x_o = 12), and the internal check of the lognormal piece comes out at exactly 1.
- **Dilation Ward identity at finite k (x_o = 1.5).** It holds with an error that falls 4× when q halves, as expected. Dropping the ray terms or A₂ breaks it badly, so the test has teeth.
- **Clock and endpoint terms.** The Ward identity doesn't test these (both vanish for a constant long mode), so I checked the clock shift against a direct root of the exact clock equation (agreement to 8 digits). The endpoint term isn't independently tested.
- **Convergence.** Results are stable to 4–5 digits in every grid.
- **Superhorizon limit.** ¼C⁽²²⁾/σ²C⁽¹¹⁾ trends to 17/3 as x_o → 0.

## What it found

1. **Inside the horizon, the fourth-order terms dominate.** ¼C⁽²²⁾ stays around 2–8 σ⁴ out to x_o = 12, while C⁽¹¹⁾ decays and oscillates through zeros of T₂ (about 0.0004–1 σ²). Relative to C⁽¹¹⁾, the fourth order is 50 to 10⁴ times larger once x_o ≳ 2. So the "small renormalisation" picture only holds above the horizon, and the expansion needs ε²σ² well below 10⁻²–10⁻⁴. This is probably the most important physics result of the exercise.
2. **The measurement nonlinearity is a minor piece.** The closed-form lognormal term 𝒩 and the mixed term ℳ together are a few percent of ¼C⁽²²⁾ at x_o ≲ 4, rising to roughly 40% at x_o = 8. The dynamical part ¼C[L₂, L₂] dominates.
3. **The local/resonant split doesn't help.** At superhorizon scales almost everything sits in the resonant A₂ part. At late times the local part and the A₂ part are each large and cancel against each other: at x_o = 12 they are 12.5 and 17.4 against a cross term of −27.5, leaving 2.4. Each piece belongs to one coordinate choice and only the sum is physical, so evaluating them separately is not a simpler route.
4. **No simple structure in the angle between the two modes.** The integrand is smooth, with nothing singular at the sound-cone resonance. It peaks where the two modes sum to K ≈ 0.8–1.0 k₀, not at K → 0, so the "two short modes beating into a long mode" idea is ruled out. The second-order mode also doesn't behave like a first-order mode of momentum K: density/T₂(K)² changes by orders of magnitude with x_o.

## So, more palatable?

- **Yes, modestly.** One piece (𝒩) is now closed form. The only term still needing third-order fields is isolated as ⅓Re C[θ, L₃]. Everything else can be computed reliably and quickly.
- **No, for analytic evaluation.** Nothing I tested suggests the dominant ¼C[L₂, L₂] collapses at finite k₀. Together with point 1, this makes me doubt that an elegant closed form at general k₀ is what the author intends. The elegant results are more likely the ones already in the report: the exact log-map and lognormal structure, the entire-function T_ℓ, and the superhorizon factorisation.

## Options from here

- **Third order.** The missing ⅓Re C[θ, L₃] needs a numerical third-order solver. That is a larger build, with the Ward identity again available for validation.
- **Late-time limit.** If a closed form exists anywhere at finite k₀, the late-time limit is the most plausible place, since ¼C⁽²²⁾ appears to approach a slowly varying O(1)σ⁴. Testing that needs better resolution at large x_o.
- **Reporting.** Decide whether any of these numerical findings belong in the report, clearly labelled as exploratory. I left the report unchanged because the problem forbids numerical evidence.

To rerun it: `cd solution && uv run --no-project --with numpy --with sympy --with mpmath python explore_fourth_order.py ward` (or `scan`). It caches the sources in `solution/results/second_order_sources_numeric.pkl`.
