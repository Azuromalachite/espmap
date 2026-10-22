#!/bin/bash
source "$(dirname "$0")/env.sh"
cd $WORK && $PY e_bench.py
