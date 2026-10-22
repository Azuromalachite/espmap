"""Write bench/env.json: the environment every benchmark row was measured in.

Kept separate from timing.py so it can be captured without editing a
script that is mid-run (timing.py re-reads itself for every row).
"""
import json
import os
import platform
import subprocess
from pathlib import Path

import matplotlib
import numpy
import pyscf
import rdkit
import scipy
import skimage

HERE = Path(__file__).resolve().parent


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True,
                          cwd=HERE.parent).stdout.strip()


cpu = next((l.split(":", 1)[1].strip() for l in sh("lscpu").splitlines()
            if l.startswith("Model name")), platform.processor())
env = {
    "git_commit": sh("git rev-parse HEAD"),
    "git_branch": sh("git branch --show-current"),
    "python": platform.python_version(),
    "platform": platform.platform(),
    "cpu": cpu,
    "logical_cpus": os.cpu_count(),
    "mem_total_GB": round(os.sysconf("SC_PAGE_SIZE")
                          * os.sysconf("SC_PHYS_PAGES") / 1024 ** 3, 1),
    "packages": {m.__name__: m.__version__
                 for m in (numpy, scipy, pyscf, rdkit, skimage, matplotlib)},
    "blas": numpy.show_config(mode="dicts")["Build Dependencies"]["blas"]["name"],
}
(HERE / "env.json").write_text(json.dumps(env, indent=2) + "\n")
print(json.dumps(env, indent=2))
