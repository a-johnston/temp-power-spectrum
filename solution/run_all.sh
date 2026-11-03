#!/usr/bin/env bash
# Run every exact check behind the report. Exits non-zero on failure.
# PYTHON defaults to the project venv (../.venv, created by `uv sync`), else python3. Slow steps (minutes to ~1 h):
# closed_second_order.py --fresh (rebuild ~20 min), observable_generator.py SH (~5 min),
# check_logmap_decomposition.py (~2 min), observable_fast.py SH3 (~12 min).
# Optional slow cross-checks (hours, not completed): closed_second_order.py --check all,
# third_order_consistency.py 4, observable_generator.py SH3 2 (superseded by observable_fast.py SH3).
# DEPTH sets the second-order Taylor depth (6 -> O(t^4); 10 -> O(t^8)).
# SLOW=1 also regenerates the cubic superhorizon kernels (superhorizon_next.py kernels3, ~13 h on a
# laptop). Without it, `superhorizon_next.py spectrum` assembles from the saved kernels if present
# (~25 min the first time, seconds once cached). curvature_kernels.py regenerates the FULL
# map through four gradients in seconds, without caches, and computes the general-moment
# spectrum. The 13 h run is now only an independent legacy cross-check.
set -euo pipefail
cd "$(dirname "$0")"
VENV_PY="$(cd .. && pwd)/.venv/bin/python"
if [ -n "${PYTHON:-}" ]; then PY="$PYTHON"; elif [ -x "$VENV_PY" ]; then PY="$VENV_PY"; else PY=python3; fi
# Foundations: exact lemmas, generic-jet field equations, first order in two gauges, normalisation.
for s in foundations/check_exact_lemmas.py foundations/check_fw_conformal.py \
         foundations/check_exact_field_equations.py "foundations/check_field_equations_generic.py jet" \
         foundations/check_first_order.py foundations/check_initial_data_gradient.py \
         foundations/normalization_expansion.py; do
  echo "=== $s"
  # shellcheck disable=SC2086
  $PY -B $s
done
echo "=== foundations/check_algebra.py"
$PY -B -m unittest discover -s foundations -p check_algebra.py
# Field solution, observable, spectrum and limits.
for s in check_evolution_form.py first_order_engine_check.py check_green.py second_order.py \
         check_ward_dilation.py check_ward_dilation3.py check_observable_expansion.py \
         check_report_identities.py check_gradient_modes.py check_dilation_factorisation.py \
         "transfer_entire.py 1 2 3" check_shell_reduction.py check_tensor_newtonian.py check_longitudinal_redshift.py check_orthonormal_driver.py check_orthonormal_map.py check_adjoint_exact.py "check_shell_cubic_sources.py --compare-engine" check_propagation_adjoint.py check_unified_adjoint.py check_first_field_adjoint_source.py check_stationary_cubic.py check_quadrupole_adjoint_expansion.py "check_stationary_quadrupole_four_gradients.py --symbolic" check_entire_sphere.py check_quadrupole_angular_source.py check_phase_integral_contraction.py check_direct_cubic_measurement.py check_dilation_ward.py check_regular_field_green.py check_quadrupole_green_source.py check_constraint_rows.py check_observer_centred_kappa4.py check_observer_dilation_k6.py check_cartesian_green.py check_bessel_translation.py check_shell_locality.py check_logmap_decomposition.py check_gradient_octupole.py third_order.py \
         "second_order_consistency.py ${DEPTH:-6}" "fast_series.py 8" \
         closed_second_order.py "observable_generator.py W2 4" "observable_generator.py W3 4" \
         "observable_generator.py SH 2" "observable_fast.py SH 2" "observable_fast.py W2 4" \
         "observable_fast.py W3 4" "observable_fast.py SH3 2" superhorizon.py \
         "superhorizon_next.py kernels2" "check_quadrupole_next_hand.py --compare" \
         check_cubic_curvature_hand.py check_cubic_product_identity.py "curvature_kernels.py --spectrum" "reconstruct_transfer.py 2 64" \
         "transfer_subhorizon.py 2 3" "transfer_subhorizon.py 3 2" "transfer_subhorizon.py 4 1"; do
  echo "=== $s"
  # shellcheck disable=SC2086
  $PY -B $s
done
if [ "${SLOW:-0}" = 1 ]; then
  echo "=== superhorizon_next.py kernels3 (slow)"
  $PY -B superhorizon_next.py kernels3
fi
if [ -f results/superhorizon_kernels_2.pkl ] && [ -f results/superhorizon_kernels_3.pkl ]; then
  echo "=== curvature_kernels.py --compare (raw independent kernels)"
  $PY -B curvature_kernels.py --compare
else
  echo "SKIP raw generator/curvature comparison: saved kernels absent; curvature spectrum above is cache-free"
fi
echo "=== superhorizon_next.py spectrum"
$PY -B superhorizon_next.py spectrum
echo "ALL CHECKS PASSED"
