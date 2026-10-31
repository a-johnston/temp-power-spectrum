#!/usr/bin/env bash
# Run every exact verification script. Exits non-zero on the first failure.
# PYTHON defaults to python3 (provided by the devbox image); override with PYTHON=/path/to/python3.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
for s in check_exact_lemmas.py check_fw_conformal.py check_exact_field_equations.py \
         "check_field_equations_generic.py jet" check_first_order.py \
         check_initial_data_gradient.py normalization_expansion.py "transfer_closed_form.py 2 3"; do
  echo "=== $s"
  # shellcheck disable=SC2086
  $PY $s
done
echo "ALL VERIFICATION SCRIPTS PASSED"
