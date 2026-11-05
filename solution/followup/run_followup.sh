#!/usr/bin/env bash
# Checks behind report/white_noise_followup.pdf. Steps run one at a time.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-../../.venv/bin/python}"
$PY -B minimal_check.py      # stand-alone closed-form substitution
$PY -B k0_fields.py          # first-order and pair residuals, printed I(x), W(x), R1
$PY -B k_series.py           # K-channel t-series against closed K=0 forms and gradient data
$PY -B check_k0_odes.py      # K=0 source equations for n2 and a
$PY -B check_temperature_theorem.py  # exact ray identities and the printed remainder against the map
$PY -B check_followup.py     # negative controls, Kohri-Terada limit, stress, helicity, omission controls
$PY -B check_chart_R.py      # off-shell chart invariance of R at generic momenta (~3 min)
$PY -B forward_plane.py      # exact PDE and rays against the second-order kernel
$PY -B forward_log.py        # the same in the logarithmic regime (~5 min)
$PY -B forward_K.py          # the same at K != 0, convergence to the K->0 limit (~20 min)
$PY -B map_k0.py             # table of quadrupole pieces
$PY -B late_time.py          # late-time regression (numerical)
$PY -B late_log_bound.py     # bound on a further logarithm (numerical, ~7 min)
$PY -B power_split.py        # division of the temperature power (~1 min)
$PY -B nonuniform.py         # amplitude nonuniformity (numerical)
$PY -B finite_K.py           # finite-K limit and P22(K)/A0 (numerical, ~6 min)
