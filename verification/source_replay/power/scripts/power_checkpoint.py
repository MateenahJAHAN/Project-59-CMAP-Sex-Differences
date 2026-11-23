"""Run the released power study, saving each complete replicate before returning it to the original driver."""
import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile

sys.path.insert(0, str(Path(__file__).resolve().parent))
import power
import numpy as np
import scipy
import pandas as pd

K = power.K
original_task = power._task

def encode_scalar(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(type(value).__name__)

def checkpointed_task(args):
    rep, labels, eligible, permutations = args
    path = Path(os.environ['P59_POWER_CHECKPOINT']) / f'rep_{rep:04d}.npz'
    fingerprint = hashlib.sha256(
        labels.tobytes() + np.asarray(eligible).tobytes() + str(permutations).encode()
    ).hexdigest()
    if path.is_file():
        with np.load(path, allow_pickle=False) as z:
            if str(z['fingerprint'].item()) != fingerprint:
                raise ValueError(f'Assignment mismatch in {path}')
            rows, plants = json.loads(z['records'].item())
            cube = z['intervals']
            estimates = {(t, s): cube[i, j] for i,t in enumerate(K.TERMS) for j,s in enumerate(K.E2_SCEN)}
            return rep, (rows, plants, estimates), z['elapsed'].item()
    rep, payload, elapsed = original_task(args)
    rows, plants, estimates = payload
    assert len(rows) == len(K.TERMS) * len(K.SCENARIOS)
    cube = np.stack([[estimates[(t, s)] for s in K.E2_SCEN] for t in K.TERMS])
    with NamedTemporaryFile(dir=path.parent, prefix='.rep_', suffix='.tmp', delete=False) as f:
        np.savez_compressed(f, fingerprint=fingerprint, elapsed=elapsed, intervals=cube,
                            records=json.dumps([rows, plants],default=encode_scalar,separators=(',', ':')))
        f.flush()
        os.fsync(f.fileno())
    os.replace(f.name, path)
    return rep, payload, elapsed

def main():
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument('--B', type=int, required=True)
    p.add_argument('--R', type=int, required=True)
    p.add_argument('--tag', default='power')
    args, _ = p.parse_known_args()
    # Use the frozen engine, numerical library versions and input manifests as the run identity.
    repo = Path(K.REPO)
    files = [Path(__file__), Path(power.__file__), Path(K.__file__), repo / 'src/G/SHA256SUMS.txt',
             repo / 'config/G_defaults.json', repo / 'config/G_inputs_pinned.tsv',
             Path(K.GRES) / 'G_sets_lines.tsv', Path(K.GRES) / 'G_sets_counts.tsv']
    manifest = {
        'files': {f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files},
        'numpy':np.__version__, 'pandas':pd.__version__, 'scipy':scipy.__version__,
        'python':sys.version, 'B':args.B, 'R':args.R, 'tag':args.tag
    }
    key = hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest()
    directory = Path(K.RAW) / 'checkpoints' / key
    directory.mkdir(parents=True,exist_ok=True)
    (directory/'run.json').write_text(json.dumps(manifest,indent=2)+'\n')
    os.environ['P59_POWER_CHECKPOINT']=str(directory)
    print(f'{len(list(directory.glob("rep_*.npz")))} of {args.R} completed replicates restored', flush=True)
    power._task = checkpointed_task
    power.main()

if __name__=='__main__':
    main()
