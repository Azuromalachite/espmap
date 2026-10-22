#!/bin/bash
source "$(dirname "$0")/env.sh"
cd $WORK && $PY e_conf.py 2>&1 && $PY e_scf.py 2>&1 | grep -v "double-zeta\|Adequate for a picture"
