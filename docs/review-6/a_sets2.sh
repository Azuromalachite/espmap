. "$(dirname "$0")/env.sh"
cd $MASTER
for f in tests/test_r11_set.py tests/test_stereo_set.py; do echo $f; $PY -m pytest -q -p no:cacheprovider $f 2>&1 | tail -1; done
grep -c "\"name\"" tests/data/stereo/manifest.json; grep -c "\"file\"" tests/data/r11/manifest.json
