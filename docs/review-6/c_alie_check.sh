. "$(dirname "$0")/env.sh"
cd $R6 && /usr/bin/time -f "%e s %M KB" $PY c_alie_check.py $MASTER
