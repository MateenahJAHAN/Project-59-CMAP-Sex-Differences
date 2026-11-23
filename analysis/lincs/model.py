"""Run the original Part B model after publication with verified settings, real inputs and full calibration."""
from argparse import ArgumentParser, Namespace
import csv
import fcntl
import json
import os
from pathlib import Path
import sys
import tempfile

from reproduce import HERE, file_hash, link

sys.path.insert(0, str(HERE / 'src/G'))


def verify_replay(info, args, cfg, reduced):
    if info['status'] != 'real' or info['n_rows_changed'] or args.mode != 'real' or reduced or args.stop_after:
        raise ValueError('The public replay requires the complete model and the registered labels.')
    for line in (HERE / 'src/G/SHA256SUMS.txt').read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        if file_hash(HERE / name) != digest:
            raise ValueError(f'Released model checksum differs: {name}')
    return {
        'execution': 'post-publication reproduction of registered analysis',
        'historical_plan_tag': 'partB-prereg-v1',
        'historical_tag_object': 'e8899f5216da25ee39760b1ab32506b803bf0f65',
        'historical_plan_commit': '572b97498d8d555e6f3616c99cb3289351223c77',
        'historical_plan_sha256': '32bdc666e9c85a43368a8f6f70269155b24fb0fde0095326c5d2178588a8bd1f',
        'release_code_manifest_sha256': file_hash(HERE / 'src/G/SHA256SUMS.txt'),
    }


def run_model(work, workers):
    with (work / "registered-worker.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        complete_model(work, workers)


def complete_model(work, workers):
    import capture
    import run_g
    import ggates
    os.environ['P59_CALIBRATION_CACHE'] = str(work / 'calibration_checkpoints')
    from calibration import run_calibration
    ggates.run_calibration = run_calibration
    bundle = work / 'registered'
    receipt = bundle / 'receipt.json'
    signature = file_hash(HERE / 'src/G/SHA256SUMS.txt')

    def file_checks(folder):
        return {str(path.relative_to(folder)): file_hash(path)
                for name in ('G', 'capture') for path in sorted((folder / name).rglob('*')) if path.is_file()}

    def expose():
        for item, name in [('G', 'results/G'), ('capture', 'capture'), ('receipt.json', 'registered_model.json')]:
            link(bundle / item, work / name)

    if receipt.exists():
        previous = json.loads(receipt.read_text())
        if previous['release_manifest'] != signature or previous['files'] != file_checks(bundle):
            raise ValueError('Completed model differs from its receipt. Use a fresh work directory for a changed analysis.')
        for row in previous['inputs']:
            source = (work if row['base'] == 'work' else HERE) / row['path']
            if file_hash(source) != row['sha256']:
                raise ValueError(f"Completed model input changed: {source}")
        expose()
        print(f"Verified {len(previous['inputs'])} model inputs and {len(previous['files'])} result files.", flush=True)
        return
    for target in (bundle, work / 'results/G', work / 'capture', work / 'registered_model.json'):
        if target.exists() or target.is_symlink():
            raise FileExistsError(f"Existing output has no matching receipt: {target}; choose a fresh work directory.")
    with tempfile.TemporaryDirectory(prefix='part-b-model-', dir=work) as temporary:
        stage = Path(temporary) / 'completed'
        stage.mkdir()
        opts = Namespace(config='config/G_defaults.json', mode='real', line_table=str(work / 'model_inputs/taskD/D_line_table.tsv'),
                         out=str(stage / 'G'), log=str(work / 'logs/G.log'), plan=None, plan_tag=None,
                         osf_id=None, stop_after=None, workers=workers, inputs=str(work / 'model_inputs'),
                         tensors_D=str(work / 'data/D/tensors'), tensors_E=str(work / 'data/E/tensors'),
                         persig_E=str(work / 'data/E/tensors'), siginfo=str(work / 'sources/siginfo_beta.txt'),
                         prism=str(work / 'sources/prism_secondary_dose_response_curve_parameters.csv'),
                         B=None, calib_n=None, calib_B=None, replay_n=None, relabel_n=None, random_n=None)
        capture.record(stage / 'capture')
        run_g.main(args=opts, replay_check=verify_replay)
        capture.finish()
        if not (stage / 'G/G_H1.tsv').is_file():
            raise RuntimeError(f'The full model stopped before the primary results. See {work / "logs/G.log"}')
        inputs = []
        with (stage / 'G/G_hash_check.tsv').open() as stream:
            for row in csv.DictReader(stream, delimiter='\t'):
                original = HERE / row['path']
                base, prefix = ('work', work) if original.is_relative_to(work) else ('release', HERE)
                inputs.append({'base':base, 'path':str(original.relative_to(prefix)), 'sha256':row['observed_sha256']})
        data = {'release_manifest':signature, 'inputs':inputs, 'files':file_checks(stage)}
        (stage / receipt.name).write_text(json.dumps(data, indent=2) + '\n')
        stage.rename(bundle)
    expose()


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--work', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    a = parser.parse_args()
    run_model(a.work.resolve(), a.workers)
