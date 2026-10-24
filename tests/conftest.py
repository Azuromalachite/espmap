import sys
from pathlib import Path

import pytest

# the module lives at the repo root, beside this tests/ directory
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def pytest_addoption(parser):
    parser.addoption("--runslow", action="store_true", default=False,
                     help="also run the tests that perform real SCF calculations")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: performs a real SCF calculation")
    config.addinivalue_line("markers", "release: hours of SCF (ensembles); run for a release by "
                                       "bench/acceptance.py, excluded from the monthly smoke")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--runslow") or config.getoption("-m") == "slow":
        return
    skip = pytest.mark.skip(reason="needs --runslow (performs real SCF work)")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip)
