#!/usr/bin/env bash
# Checks behind report/white_noise_followup.pdf. Total about 10 min; steps run one at a time, peak memory about 2.7 GB.
set -euo pipefail
cd "$(dirname "$0")"
PY="${PYTHON:-../../.venv/bin/python}"
$PY -B minimal_check.py      # stand-alone closed-form substitution
$PY -B k0_fields.py          # first-order and pair residuals, printed I(x), W(x), R1
$PY -B check_followup.py     # negative controls, Kohri-Terada limit, stress, helicity, omission controls
$PY -B check_chart_R.py      # off-shell chart invariance of R at generic momenta (~3 min)
$PY -B forward_plane.py      # exact PDE and rays against the second-order kernel
$PY -B map_k0.py             # table of quadrupole pieces
$PY -B late_time.py          # late-time regression (numerical)
$PY -B nonuniform.py         # amplitude nonuniformity (numerical)
$PY -B finite_K.py           # finite-K limit and P22(K)/A0 (numerical, ~6 min)
