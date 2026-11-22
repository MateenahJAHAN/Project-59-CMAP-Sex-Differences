"""Re-read each entire Broad CEL with Biopython and check the summaries from the released tables."""
import argparse
import gzip
import json
from pathlib import Path

from Bio.Affy import CelFile
import Bio
import numpy as np
import pandas as pd

from broad_cel import CELLS, OUT, cel


def main(cache):
    checked={}
    for gsm in CELLS.values():
        src=cache/(gsm+".CEL.gz")
        with gzip.open(src,"rb") as handle:
            reference=CelFile.read(handle).intensities.ravel()
        values,_,_=cel(src)
        assert np.array_equal(reference,values),gsm
        checked[gsm]=int(reference.size)
    (OUT/"cel_reader_check.json").write_text(json.dumps(
        {"reader":"Bio.Affy.CelFile, Biopython "+Bio.__version__,"exactly_equal_intensities":checked},indent=2)+"\n")
    cols=["study","cell_key","arm"]
    study=pd.read_csv(OUT/"complete_arm_medians.tsv",sep="\t").set_index(cols).sort_index()
    old=pd.read_csv(OUT/"y_probe_intensities.tsv.gz",sep="\t").query('study=="GSE79185"')
    raw=pd.read_csv(OUT/"idat_probes.tsv",sep="\t").query('study=="GSE86078"')
    values=pd.concat([old,raw]).groupby(cols).log2_centered.agg(["count","median"]).sort_index()
    assert study.index.equals(values.index) and np.array_equal(study["count"],values["count"])
    assert np.allclose(study["median"],values["median"],atol=2e-7,rtol=0)
    vals=pd.read_csv(OUT/"broad_y_probes.tsv.gz",sep="\t").groupby(["cell","arm"]).centered.agg(["count","median"])
    stats=pd.read_csv(OUT/"broad_summary.tsv",sep="\t").query('arm!="autosome" and probes=="all"').set_index(["cell","arm"]).sort_index()
    assert vals.index.equals(stats.index) and (vals["count"]==stats.n_sets).all()
    assert np.allclose(vals["median"],stats["median"],atol=2e-7,rtol=0)
    print(f"{sum(checked.values()):,} Broad feature values matched the independent reader; {len(study)+len(stats)} arm summaries reproduced.")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("cache",type=Path)
    main(p.parse_args().cache)
