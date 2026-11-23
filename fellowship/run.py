"""Run the numbered fellowship analyses in an external output directory."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from sources import REPO, prepare

HERE = Path(__file__).resolve().parent
SCREEN = REPO / 'verification/source_replay/phases2_5'


def outside_checkout(path):
    path = Path(path).expanduser().resolve()
    if path == REPO or REPO in path.parents:
        raise argparse.ArgumentTypeError('Choose a dedicated directory outside the source checkout.')
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=outside_checkout, required=True, help='Input cache, populated on first run')
    parser.add_argument('--out', type=outside_checkout, required=True, help='Separate directory for this run')
    parser.add_argument('--workers', type=int, default=2, help='Parallel compound processes in Phases 2–5')
    parser.add_argument('--phases', type=int, nargs='+', choices=range(1, 9), default=list(range(1, 9)))
    parser.add_argument('--prepare-only', action='store_true', help='Fetch and verify inputs without running the screens')
    parser.add_argument('--phase6-python', default=sys.executable,
                        help='Python with GSEApy 1.3.2 for the Phase 6 replay')
    args = parser.parse_args()
    if args.workers < 1 or args.data == args.out:
        parser.error('Use at least one worker and distinct data and output directories.')
    phases = set(args.phases)
    if phases != {8} or args.prepare_only:
        prepare(args.data)
    if args.prepare_only:
        return
    args.out.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, P59_PUBLIC=str(REPO), P59_GTEX=str(args.data), P59_META=str(args.data),
               P59_MATRIX=str(args.data), P59_PREPARED=str(args.out/'prepared'), P59_OUTPUT=str(args.out/'outputs'),
               P59_GENE_TESTS=str(args.out/'phase2_gene_tests'), P59_COMPARISON=str(args.out/'comparison'), P59_PHASE_RUN=str(args.out), PYTHONHASHSEED='0', MPLBACKEND='Agg',
               RAYON_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    if phases - {8}:
        import gseapy, numpy, pandas, scipy
        identity = {'numpy':numpy.__version__, 'pandas':pandas.__version__, 'scipy':scipy.__version__,
                    'gseapy':gseapy.__version__, 'python':sys.version.split()[0], 'pythonhashseed':0,
                    'inputs':hashlib.sha256((SCREEN/'source_manifest.json').read_bytes()).hexdigest(),
                    'source':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in
                              [HERE/'run.py', HERE/'sources.py', HERE/'other_phases.py', SCREEN/'replay.py', SCREEN/'compare.py', SCREEN/'figures.py', SCREEN/'verify_outputs.py',
                               *sorted(REPO.glob('Phase[2-5]*/*.ipynb'))]}}
        marker = args.out/'run-identity.json'
        if marker.exists() and json.loads(marker.read_text()) != identity:
            parser.error('The software or inputs changed. Use a new output directory.')
        if not marker.exists() and (args.out/'outputs').exists():
            parser.error('Existing results lack a run identity. Use a new output directory.')
        marker.write_text(json.dumps(identity, indent=2)+'\n')
    def run(python, script, *items, extra_env=None):
        subprocess.run([python,'-B',str(script),*map(str,items)], env=extra_env or env, check=True)
    if phases - {8}:
        run(sys.executable, SCREEN/'replay.py', '--prepare')
    if 1 in phases:
        run(sys.executable, HERE/'other_phases.py', 1, args.out)
    main_phases = sorted(phases & {2,3,4,5})
    if main_phases:
        if gseapy.__version__ != '1.2.1':
            parser.error('Phases 2–5 use GSEApy 1.2.1. Install fellowship/requirements.txt.')
        run(sys.executable, SCREEN/'replay.py', '--workers',args.workers,'--phases',*main_phases)
        run(sys.executable, SCREEN/'compare.py', *main_phases)
        run(sys.executable, HERE/'compare_release.py', '--out',args.out,'--phases',*main_phases)
        if {3,4,5}.issubset(phases):
            run(sys.executable, SCREEN/'figures.py')
        if {2,3,4,5}.issubset(phases):
            run(sys.executable, SCREEN/'verify_outputs.py')
    if 6 in phases or (7 in phases and not all((args.out/'phase6'/f).is_file() for f in ['tissue_results.csv.gz','comparison.json'])):
        extra_env = env.copy()
        alternate = Path(sys.prefix)/'phase6'
        if args.phase6_python == sys.executable and alternate.is_dir():
            extra_env['PYTHONPATH'] = str(alternate)
        run(args.phase6_python, HERE/'other_phases.py', 6, args.out, extra_env=extra_env)
        if 6 in phases:
            run(args.phase6_python, REPO/'Phase6_NoSampleFilter/code/02_case_study_plots.py',
                '--data',args.data,'--out',args.out/'phase6/cases',extra_env=extra_env)
    for phase in sorted(phases & {7,8}):
        run(sys.executable, HERE/'other_phases.py', phase, args.out)
    notebook_plots = sorted(phases & {1,2,6})
    if notebook_plots:
        run(sys.executable, HERE/'notebook_figures.py', '--data',args.data,'--out',args.out,'--phases',*notebook_plots)


if __name__ == '__main__':
    main()
