. "$(dirname "$0")/env.sh"
P=$S/r11_probe
echo "#### master dc9d9ef"; $PY -I $P/run_probe.py $MASTER $MASTER/tests/data/r11 $P/extra 2>&1 | grep -v -i warn
echo "#### v0.2.0 6f8b8eb esp_map.py, master set"; $PY -I $P/run_probe.py $S/r6_v020 $MASTER/tests/data/r11 $P/extra 2>&1 | grep -v -i warn | grep "^==\|esp_map from"
