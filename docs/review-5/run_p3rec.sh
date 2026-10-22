source /mnt/c/Users/user/AppData/Local/Temp/claude/C--Users-user-Documents-espmap/fe5c9cfb-6f3c-4b9b-b67a-3ad1fd947a5d/scratchpad/review5-work/env.sh
export PYTHONDONTWRITEBYTECODE=1
rm -rf $HOME/p3copy && cp -r $P3 $HOME/p3copy && rm -rf $HOME/p3copy/.git
cd $HOME/p3copy
$PY $W/e_p3records.py $HOME/p3copy $REPO > $W/e_p3records.out 2>&1; echo exit $?
$PY -m pytest -q -p no:cacheprovider tests/test_fixtures.py 2>&1 | tail -3 >> $W/e_p3records.out
ls ~/.venvs/espmap/lib/python3.12/site-packages | grep -i espmap >> $W/e_p3records.out
cat ~/.venvs/espmap/lib/python3.12/site-packages/*espmap*/direct_url.json 2>/dev/null >> $W/e_p3records.out
cat $W/e_p3records.out
