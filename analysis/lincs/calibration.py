"""Save each completed assignment from the frozen Part B calibration routine for restart."""
import hashlib
import fcntl
import json
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import ggates as gates

CACHE = Path(os.environ["P59_CALIBRATION_CACHE"])


def array_digest(array):
    value = np.ascontiguousarray(array)
    return value.dtype.str, value.shape, hashlib.sha256(memoryview(value).cast("B")).hexdigest()


def run_id(panel, tp, codes, elig, family, term, permutations, rule):
    code = {}
    for file in ("engine.py", "gcore.py", "ggates.py"):
        path = Path(gates.__file__).parent / file
        code[file] = hashlib.sha256(path.read_bytes()).hexdigest()
    spec = dict(code=code, tensor=array_digest(panel.R), genes=list(panel.genes),
                compounds=list(panel.comps), lines=list(panel.lines), lineage=array_digest(panel.lin),
                permutation_keys=array_digest(panel.keys), tp=array_digest(tp),
                codes=array_digest(codes), elig=array_digest(elig), family=family,
                term=term, permutations=permutations, rule=rule)
    encoded = json.dumps(spec, sort_keys=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def number(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value))


def compute_chunk(args):
    directory, tp, codes, elig, family, term, permutations, rule, indices = args
    for index in indices:
        finished = directory / f'{index}.json'
        if finished.is_file():
            continue
        with (directory / f'.{index}.lock').open('a') as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            if finished.is_file():
                continue
            values = gates._calib_task((tp,codes,elig,family,term,permutations,rule,[index]))
            assert len(values) == 1 and values[0]['pseudo'] == index
            with tempfile.NamedTemporaryFile('w', dir=directory, delete=False, prefix='.writing-') as out:
                json.dump(values[0], out, default=number)
                temp = out.name
            os.replace(temp, finished)


def run_calibration(panel_args, panel, tp, codes, elig, family, term, B, rule, s_idx, workers, chunk=10):
    indices = list(s_idx)
    key = run_id(panel,tp,codes,elig,family,term,B,rule)
    directory = CACHE / key
    directory.mkdir(parents=True, exist_ok=True)
    pending = [i for i in indices if not (directory / f'{i}.json').is_file()]
    if pending:
        workers = min(workers, int(os.environ.get('P59_CAL_WORKERS', workers)))
        with ProcessPoolExecutor(max_workers=max(workers, 1), initializer=gates._init,
                                 initargs=(panel_args,)) as pool:
            tasks = [(directory,tp,codes,elig,family,term,B,rule,pending[k:k+5])
                     for k in range(0,len(pending),5)]
            list(pool.map(compute_chunk, tasks))
    result = [json.loads((directory / f'{i}.json').read_text()) for i in sorted(indices)]
    assert [r['pseudo'] for r in result] == sorted(indices)
    return result
