#!/bin/bash
source "$(dirname "$0")/env.sh"
cd $WORK && $PY e_fix3.py 2>&1
