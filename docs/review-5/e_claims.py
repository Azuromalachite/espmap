"""Review 5: which per-site numbers does test_claude_md_quotes_the_per_site_budget
actually pin? It checks `value in CLAUDE.md` (whole file, substring). For each
value: how many times it occurs as a substring, as a standalone token, and
whether deleting it from the per-site paragraph still passes the test.
Also: the C18 genus split behind "8 extended".
    python e_claims.py <repo>
"""
import json
import re
import sys
from pathlib import Path

REPO = Path(sys.argv[1])
sys.path.insert(0, str(REPO / "tests"))
sys.path.insert(0, str(REPO))
import test_claims as tc  # noqa: E402

rows = tc._site_rows()
text = tc.TEXT["CLAUDE.md"]
para = text[text.index("**Per-site statistics have their own budget**"):text.index("## Deliberately absent")]
keys = ("partition small anions", "partition C18", "conformer C18", "grid mean",
        "conformer C18 extended", "conformer C18 extended max", "partition atoms",
        "partition atoms C18", "atoms buried C18")
print(f"{'row':28} {'value':7} substr(file) token(file) token(para)  passes if removed from para?")
for k in keys:
    v = rows[k]
    tok = re.compile(r"(?<![\d.])" + re.escape(v) + r"(?![\d])")
    stripped = text.replace(para, tok.sub("X", para))
    print(f"{k:28} {v:7} {text.count(v):^12} {len(tok.findall(text)):^11} {len(tok.findall(para)):^11}"
          f"  {v in stripped}")

c18 = [json.loads(ln) for ln in (REPO / "bench" / "site_partition.jsonl").read_text().splitlines()
       if ln.strip() and json.loads(ln)["part"] == "c18"]
print("\nC18 rows by genus:", sorted((r.get("genus"), r.get("conformer", r.get("key"))) for r in c18))
