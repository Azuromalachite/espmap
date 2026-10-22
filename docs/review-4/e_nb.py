import os, sys, json, subprocess, hashlib, ast, builtins
REPO = os.environ["REPO"]; W = os.environ["W"]
out = W + "/nb_regen.ipynb"
subprocess.run([sys.executable, REPO + "/notebooks/build_notebook.py", out], check=True)
a = open(out, "rb").read(); b = open(REPO + "/notebooks/espmap_colab.ipynb", "rb").read()
print("regenerated == shipped:", a == b, hashlib.sha1(b).hexdigest()[:10])
nb = json.loads(b)
bound = set(dir(builtins)) | {"get_ipython"}
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code": continue
    src = "".join(c["source"])
    t = ast.parse(src)
    loads, stores = set(), set()
    for n in ast.walk(t):
        if isinstance(n, ast.Name):
            (stores if isinstance(n.ctx, (ast.Store, ast.Del)) else loads).add(n.id)
        elif isinstance(n, (ast.FunctionDef, ast.ClassDef)): stores.add(n.name)
        elif isinstance(n, ast.arg): stores.add(n.arg)
        elif isinstance(n, (ast.Import, ast.ImportFrom)):
            for al in n.names: stores.add((al.asname or al.name).split(".")[0])
    free = loads - stores - bound
    unbound = sorted(x for x in free if x not in bound)
    print(f"cell {i}: exec_count={c['execution_count']} outputs={len(c['outputs'])} free-not-yet-bound={unbound}")
    bound |= stores
