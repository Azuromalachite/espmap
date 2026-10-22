source /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/fe5c9cfb-6f3c-4b9b-b67a-3ad1fd947a5d/scratchpad/review5-work/env.sh
cd $HOME
free -m | head -2
$PY $W/e_records.py $REPO $W/out $W/fx "$@" > $W/e_records.out 2>&1; echo exit $?
