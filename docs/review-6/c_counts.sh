. "$(dirname "$0")/env.sh"
cd $S/r6_v021 && echo "v0.2.1 collect:" && $PY -m pytest -q -p no:cacheprovider --collect-only --runslow 2>&1 | tail -1 && $PY -m pytest -q -p no:cacheprovider --collect-only --runslow -m slow 2>&1 | tail -1 && grep -c "xfail(strict" tests/*.py; grep -n "xfail" tests/*.py | wc -l
