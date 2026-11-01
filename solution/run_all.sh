#!/usr/bin/env bash
# Run every exact check of the unified derivation and the final report. Exits non-zero on failure.
# PYTHON defaults to python3 (provided by the devbox image). Slow steps (minutes to ~1 h):
# closed_second_order.py --fresh (rebuild ~20 min), observable_generator.py SH (~5 min).
# Optional slow cross-checks (hours, not completed): closed_second_order.py --check all,
# third_order_consistency.py 4, observable_generator.py SH3 2.
# DEPTH sets the second-order Taylor depth (6 -> O(t^4); 10 -> O(t^8)).
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
for s in check_evolution_form.py first_order_engine_check.py check_green.py second_order.py \
         check_ward_dilation.py check_ward_dilation3.py check_observable_expansion.py \
         check_report_identities.py check_gradient_modes.py third_order.py \
         "second_order_consistency.py ${DEPTH:-6}" "fast_series.py 8" \
         closed_second_order.py "observable_generator.py W2 4" "observable_generator.py W3 4" \
         "observable_generator.py SH 2" superhorizon.py; do
  echo "=== $s"
  # shellcheck disable=SC2086
  $PY -B $s
done
echo "ALL UNIFIED CHECKS PASSED"
