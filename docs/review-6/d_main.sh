# run after d_revert.sh has finished (one process at a time)
. "$(dirname "$0")/env.sh"
COPY=$S/review6-work/revert/tree
{ echo "### master worktree"; $PY $R6/d_version.py $MASTER
  echo "### copy without .git (git archive of dc9d9ef)"; $PY $R6/d_version.py $COPY
  echo "### copy, pyproject version vs module"; grep -n "version" $COPY/pyproject.toml | head -5
} > $R6/out/d_version.txt 2>&1
{ $PY $R6/d_keys.py $MASTER; echo; echo "########## with B3LYP_WITH_VWN5"; $PY $R6/d_keys.py $MASTER vwn5; } > $R6/out/d_keys.txt 2>&1
$PY $R6/d_records.py $MASTER $S/review6-work/records > $R6/out/d_records.txt 2>&1
