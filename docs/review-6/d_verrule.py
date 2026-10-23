"""Does a record made from master dc9d9ef satisfy the project's own version
rule (tests/test_fixtures.py::version_agrees_with_code)? usage: d_verrule.py CODE_DIR DESCRIBE"""
import sys
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[1] + "/tests")
import esp_map as em
from test_fixtures import version_agrees_with_code
print("__version__", em.__version__, "| git describe", sys.argv[2],
      "| version_agrees_with_code:", version_agrees_with_code(em.__version__, sys.argv[2]))
