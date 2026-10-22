source /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/fe5c9cfb-6f3c-4b9b-b67a-3ad1fd947a5d/scratchpad/review5-work/env.sh
cd $HOME
for w in mapping glycine; do
  /usr/bin/time -f "peak %M kB, %e s" $PY $W/e_forced.py $REPO $W/out $w > $W/e_forced_$w.out 2>&1
  grep -A12 "^RECORD\|peak\|WARNING: after\|WARNING: before\|Traceback\|Error" $W/e_forced_$w.out | head -30
done
