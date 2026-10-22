#!/bin/bash
source "$(dirname "$0")/env.sh"
cd $WORK
$PY -c "import sys,numpy,scipy,pyscf,rdkit,skimage,matplotlib,plotly; print(sys.version); print(numpy.__version__,scipy.__version__,pyscf.__version__,rdkit.__version__,skimage.__version__,matplotlib.__version__,plotly.__version__)"
$PY -c "import cmcrameri, colorspacious; print('cmcrameri ok')"
# run tests against a copy so no cache dirs land in the release tree
rm -rf $WORK/tree && cp -r $SRC $WORK/tree
cd $WORK/tree && $PY -m pytest -q -p no:cacheprovider 2>&1 | tail -5
