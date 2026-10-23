. "$(dirname "$0")/env.sh"
for d in $MASTER $PR; do echo "== $d"; cd $d && /usr/bin/time -f "%e s %M KB" $PY -m pytest -q -p no:cacheprovider 2>&1 | tail -4; done
