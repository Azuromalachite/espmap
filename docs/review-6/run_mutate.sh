. "$(dirname "$0")/env.sh"
mkdir -p ~/r6mut && cd $R6 && $PY e_mutate.py $PR ~/r6mut 2>&1
