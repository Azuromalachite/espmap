#!/bin/bash
# Environment for the review-3 reproduction scripts. Override any of these.
#   SRC   a checkout of the commit under review (v0.1.0 = 4bddf7b)
#   WORK  this folder (fixtures in fx/, outputs written to out/)
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export PY=${PY:-python}
export SRC=${SRC:-$(cd "$HERE/../.." && pwd)}
export WORK=${WORK:-$HERE}
export PYSCF_TMPDIR=${PYSCF_TMPDIR:-$HOME/pyscf_tmp_r3}
mkdir -p "$PYSCF_TMPDIR" "$WORK/out"
export PYTHONDONTWRITEBYTECODE=1
