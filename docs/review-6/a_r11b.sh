. "$(dirname "$0")/env.sh"
cd $R6 && $PY a_r11b.py $MASTER $R6/fx/r11
echo "== CLI --check on the 2D-header z drawing"
cd $MASTER && $PY esp_map.py --check --file $R6/fx/r11/draw_threonine_2Dheader_zwedge.mol 2>&1 | tail -6
