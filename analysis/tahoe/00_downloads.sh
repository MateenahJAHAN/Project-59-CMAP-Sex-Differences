#!/bin/sh
# Fetch and fingerprint the primary caller's four external annotation files.
set -eu
: "${P59_DATA:?Set P59_DATA to the run's data directory}"
code=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
out=${P59_PARTA_OUT:-"$P59_DATA/tahoe/partA_local"}
py=${P59_PYTHON:-"$P59_DATA/venv/bin/python"}
if [ ! -x "$py" ]; then py=python3; fi
"$py" - "$code" "$P59_DATA" "$out" <<'PY'
from pathlib import Path
import sys
sys.path.insert(0, str(Path(sys.argv[1]).parent / "part_a"))
from run import annotations
annotations(Path(sys.argv[2]).resolve(), Path(sys.argv[3]).resolve())
PY
