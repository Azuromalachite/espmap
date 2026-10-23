"""Review 6, R18 under a stale install: which version a record names when the
module that runs is not the installed distribution. usage: d_version.py CODE_DIR"""
import json, sys
from importlib import metadata
sys.path.insert(0, sys.argv[1])
import esp_map as em
d = metadata.distribution("espmap") if True else None
print("module file:", em.__file__)
print("module __version__:", em.__version__)
print("installed espmap dist:", d.version, "direct_url:", d.read_text("direct_url.json"))
print("code_identity():", json.dumps(em.code_identity()))
