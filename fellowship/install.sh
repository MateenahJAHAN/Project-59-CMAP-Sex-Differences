#!/bin/sh
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
envdir=${1:?Usage: sh fellowship/install.sh /path/to/venv}
"${PYTHON:-python3.12}" -m venv "$envdir"
"$envdir/bin/python" -m pip install -r "$here/requirements.txt"
"$envdir/bin/python" -m pip install --upgrade --no-deps --target "$envdir/phase6" gseapy==1.3.2
