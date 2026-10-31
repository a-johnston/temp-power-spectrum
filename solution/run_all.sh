#!/usr/bin/env bash
# Run every exact check of the unified derivation. Exits non-zero on the first failure.
# PYTHON defaults to python3 (provided by the devbox image).
# The second-order Taylor check at depth 6 takes about 7 minutes; pass DEPTH=10 for O(t^8).
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
for s in check_evolution_form.py first_order_engine_check.py check_green.py second_order.py \
         check_ward_dilation.py check_observable_expansion.py third_order.py \
         "second_order_consistency.py ${DEPTH:-6}"; do
  echo "=== $s"
  # shellcheck disable=SC2086
  $PY -B $s
done
echo "ALL UNIFIED CHECKS PASSED"
