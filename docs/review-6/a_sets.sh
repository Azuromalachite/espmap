. "$(dirname "$0")/env.sh"
cd $MASTER
echo "== committed sets (pytest)"
$PY -m pytest -q -p no:cacheprovider tests/test_r11_set.py tests/test_stereo_set.py 2>&1 | tail -3
echo "== r11 held-out (scratchpad/r11_probe/extra) + committed, by run_probe.py"
$PY $S/r11_probe/run_probe.py $MASTER $MASTER/tests/data/r11 $S/r11_probe/extra 2>&1 | tail -25
echo "== stereo held-out, no flag"
$PY $S/stereo_probe/run_probe.py $MASTER $S/stereo_probe/heldout 2>&1 | tail -2
echo "== stereo held-out, --require-stereo"
$PY $S/stereo_probe/run_probe.py $MASTER $S/stereo_probe/heldout --require-stereo 2>&1 | grep -v "^ok " 
