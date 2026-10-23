. "$(dirname "$0")/env.sh"
cd $R6 && /usr/bin/time -f "%e s %M KB" $PY e_select.py $PR "$@" 2>&1
