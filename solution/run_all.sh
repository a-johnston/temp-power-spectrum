#!/usr/bin/env bash
# Run every exact check behind the report. Exits non-zero on failure.
# PYTHON defaults to the project venv (../.venv, created by `uv sync`), else python3. Slow steps (minutes to ~1 h):
# closed_second_order.py --fresh (rebuild ~20 min), observable_generator.py SH (~5 min),
# check_logmap_decomposition.py (~2 min), observable_fast.py SH3 (~12 min).
# Optional slow cross-checks (hours, not completed): closed_second_order.py --check all,
# third_order_consistency.py 4, observable_generator.py SH3 2 (superseded by observable_fast.py SH3).
# DEPTH sets the second-order Taylor depth (6 -> O(t^4); 10 -> O(t^8)).
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
         "transfer_entire.py 1 2 3" check_shell_reduction.py check_tensor_newtonian.py check_longitudinal_redshift.py check_shell_locality.py check_logmap_decomposition.py check_gradient_octupole.py third_order.py \
         "second_order_consistency.py ${DEPTH:-6}" "fast_series.py 8" \
         closed_second_order.py "observable_generator.py W2 4" "observable_generator.py W3 4" \
         "observable_generator.py SH 2" "observable_fast.py SH 2" "observable_fast.py W2 4" \
         "observable_fast.py W3 4" "observable_fast.py SH3 2" superhorizon.py \
         "superhorizon_next.py kernels2" "check_quadrupole_next_hand.py --compare"; do
  echo "=== $s"
  # shellcheck disable=SC2086
  $PY -B $s
done
echo "ALL CHECKS PASSED"
