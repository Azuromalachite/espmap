source /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/fe5c9cfb-6f3c-4b9b-b67a-3ad1fd947a5d/scratchpad/review5-work/env.sh
cd $W
/usr/bin/time -v $PY e_fast.py $REPO $P3 $W/fx 2> $W/e_fast.time > $W/e_fast.out; echo exit $?
grep -E "Maximum resident|Elapsed" $W/e_fast.time
