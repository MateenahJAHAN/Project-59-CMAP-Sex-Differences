"""Exercise resumption for one already-saved assignment in a fresh process."""
import argparse
import io
import json
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'scripts'))
import power_checkpoint as checkpoint
K=checkpoint.K

def compare(records, path):
    # Read both serialized tables with the same null and round-trip float handling as the full audit.
    frame=pd.DataFrame(records)
    left=pd.read_csv(io.StringIO(frame.to_csv(sep='\t',index=False)),sep='\t',float_precision='round_trip')
    right=pd.read_csv(path,sep='\t',nrows=len(frame),float_precision='round_trip')
    assert left.shape==right.shape and set(left.columns)==set(right.columns)
    assert list(left.rep)==list(right.rep)
    fields=[c for c in left if c!='time_s']
    dif=left[fields].compare(right[fields])
    return sorted(set(dif.columns.get_level_values(0)))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint_dir',type=Path)
    parser.add_argument('archived_raw',type=Path)
    opts=parser.parse_args()
    config=json.loads((opts.checkpoint_dir/'run.json').read_text())
    assert (opts.checkpoint_dir/'rep_0000.npz').is_file()
    os.environ['P59_POWER_CHECKPOINT']=str(opts.checkpoint_dir)
    D=K.Data(with_level3=True)
    pairs,counts=K.factorial_set(D.P)
    labels,eligible=K.draw_codes(D.P,pairs,counts,config['R'],K.SEED_PSEUDO)
    rep,(calls,plants,arrays),_=checkpoint.checkpointed_task((0,labels[:,0].copy(),eligible,config['B']))
    with np.load(opts.archived_raw/'power_e2.npz',allow_pickle=False) as old:
        cube=np.stack([[arrays[(t,s)] for s in K.E2_SCEN] for t in K.TERMS])
        identical=bool(np.array_equal(cube,old['E2'][rep],equal_nan=True))
    result=dict(replicate=rep,restored_in_new_process=True,scenarios=len(calls),planted_rows=len(plants),
                changed_scenario_fields=compare(calls,opts.archived_raw/'power_runs.tsv.gz'),
                changed_planted_fields=compare(plants,opts.archived_raw/'power_planted.tsv.gz'),
                interval_array_value_count=int(cube.size),intervals_exact=identical)
    (ROOT/'results/checkpoint_restore.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)
    assert identical and not result['changed_scenario_fields'] and not result['changed_planted_fields']
if __name__=='__main__':
    main()
