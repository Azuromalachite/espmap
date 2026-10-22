export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4
export PYSCF_TMPDIR=$HOME/pyscf_tmp
mkdir -p $PYSCF_TMPDIR
export S=/mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/fe5c9cfb-6f3c-4b9b-b67a-3ad1fd947a5d/scratchpad
export REPO=$S/espmap-review5
export P3=$S/espmap-p3-read
export W=$S/review5-work
export PY=~/.venvs/espmap/bin/python
