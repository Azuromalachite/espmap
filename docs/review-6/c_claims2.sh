. "$(dirname "$0")/env.sh"
cd $R6 && $PY c_claims2.py $MASTER
cd $MASTER && echo "master collect (fast + --runslow):" && $PY -m pytest -q -p no:cacheprovider --collect-only 2>&1 | tail -1 && $PY -m pytest -q -p no:cacheprovider --collect-only --runslow 2>&1 | tail -1 && $PY -m pytest -q -p no:cacheprovider --collect-only --runslow -m slow 2>&1 | tail -1
rm -rf ~/r6_v021 && mkdir ~/r6_v021 && cd $MASTER && git archive 41a989d | tar -x -C ~/r6_v021 && cd ~/r6_v021 && echo "v0.2.1 41a989d collect --runslow:" && $PY -m pytest -q -p no:cacheprovider --collect-only --runslow 2>&1 | tail -1 && $PY -m pytest -q -p no:cacheprovider --collect-only --runslow -m slow 2>&1 | tail -1 && grep -c "xfail" tests/*.py | head; wc -l esp_map.py
cd $MASTER && wc -l esp_map.py
