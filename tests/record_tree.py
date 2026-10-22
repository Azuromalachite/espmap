"""Compare a record against docs/record-schema.json, in both directions.

Review 4 N5: the spec declared keys the code never wrote, and the code wrote
keys the spec renamed. A one-way check (every schema key present) misses the
second; a key the code writes and the schema does not declare fails here too.

Schema notation: an object lists its keys exactly ("$"-prefixed keys are
notation, "$nullable": true lets the whole object be null); a one-element
list is a list whose items all match that element, and {"$list": item,
"$nullable": true} is the same list or null; a string is a leaf type,
"str" | "num" | "int" | "bool" | "xyz" (three numbers) | "any", with a
trailing "?" for nullable.
"""
import json
from pathlib import Path

SCHEMA_PATH = Path(__file__).resolve().parents[1] / "docs" / "record-schema.json"


def load_schema():
    """The schema with every "#name" (or "#name?", nullable) replaced by
    $defs[name], so a shape used four times is written once."""
    raw = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    defs = raw.get("$defs", {})

    def resolve(node):
        if isinstance(node, str) and node.startswith("#"):
            name = node[1:].rstrip("?")
            out = dict(resolve(defs[name]))
            if node.endswith("?"):
                out["$nullable"] = True
            return out
        if isinstance(node, list):
            return [resolve(x) for x in node]
        if isinstance(node, dict):
            return {k: resolve(v) for k, v in node.items()}
        return node
    return {k: resolve(v) for k, v in raw.items() if k != "$defs"}


_LEAF = {"str": (str,), "int": (int,), "bool": (bool,), "num": (int, float)}


def tree_problems(schema, value, path="record"):
    """Every mismatch between schema and value, as 'path: what' strings."""
    if isinstance(schema, str):
        kind, nullable = schema.rstrip("?"), schema.endswith("?")
        if value is None:
            return [] if nullable or kind == "any" else [f"{path}: null, expected {kind}"]
        if kind == "any":
            return []
        if kind == "xyz":
            ok = (isinstance(value, list) and len(value) == 3
                  and all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in value))
            return [] if ok else [f"{path}: {value!r}, expected three numbers"]
        ok =isinstance(value, _LEAF[kind]) and not (kind != "bool" and isinstance(value, bool))
        return [] if ok else [f"{path}: {type(value).__name__}, expected {kind}"]
    if isinstance(schema, list):
        if not isinstance(value, list):
            return [f"{path}: {type(value).__name__}, expected a list"]
        return [p for i, v in enumerate(value) for p in tree_problems(schema[0], v, f"{path}[{i}]")]
    if value is None:
        return [] if schema.get("$nullable") else [f"{path}: null, expected an object"]
    if "$list" in schema:
        return tree_problems([schema["$list"]], value, path)
    if not isinstance(value, dict):
        return [f"{path}: {type(value).__name__}, expected an object"]
    keys = {k for k in schema if not k.startswith("$")}
    out = [f"{path}.{k}: declared by the schema, absent from the record"
           for k in sorted(keys - set(value))]
    out += [f"{path}.{k}: written by the code, absent from the schema"
            for k in sorted(set(value) - keys)]
    for k in sorted(keys & set(value)):
        out += tree_problems(schema[k], value[k], f"{path}.{k}")
    return out
