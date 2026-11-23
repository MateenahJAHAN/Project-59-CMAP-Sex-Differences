"""Step 10: CONSORT table and sha256 of every result file and script (F2_hashes.tsv).
Usage: s10_consort_hashes.py [hashes-only]"""
import glob, hashlib, json, os, sys
import pandas as pd
R, D = "results/F2", "data/F2"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


if len(sys.argv) < 2 or sys.argv[1] != "hashes-only":
    obs = pd.read_parquet(f"{D}/obs_dmso.parquet")
    sh = pd.read_csv(f"{R}/F2_shard_summary.tsv", sep="\t")
    iv = pd.read_csv(f"{R}/F2_index_vs_obs.tsv", sep="\t")
    cs = json.load(open(f"{D}/consort_states.json"))
    st = pd.read_csv(f"{R}/F2_states.tsv", sep="\t")
    rows = [
        ("obs_metadata", "obs_metadata rows (all cells)", 100648790, "", "s01_obs.py assertion"),
        ("obs_metadata", "rows with drug == DMSO_TF", len(obs), "", "data/F2/obs_dmso.parquet"),
        ("obs_metadata", "DMSO_TF rows with pass_filter full", int((obs.pass_filter == "full").sum()), int(obs[obs.pass_filter == "full"].cell_name.nunique()), "lines in column 4"),
        ("shards", "shards read (footer)", len(sh), "", "F2_shard_summary.tsv"),
        ("shards", "shards holding any vehicle row", int((sh.vehicle_rows > 0).sum()), "", "F2_shard_summary.tsv"),
        ("shards", "row groups holding vehicle rows", int(sh.vehicle_row_groups.sum()), "", "F2_dmso_index.tsv"),
        ("shards", "vehicle rows in shards", int(sh.vehicle_rows.sum()), "", "F2_dmso_index.tsv"),
        ("join", "shard vehicle rows matched to a full-pass obs row (analysis set)", int(iv.shard_and_obs_full.sum()), "", "F2_index_vs_obs.tsv"),
        ("join", "shard vehicle rows matched to a minimal-pass obs row (dropped)", int(iv.shard_and_obs_minimal.sum()), "", "F2_index_vs_obs.tsv"),
        ("join", "shard vehicle rows without an obs row (dropped)", int(iv.shard_not_in_obs.sum()), "", "F2_index_vs_obs.tsv"),
        ("join", "full-pass obs vehicle rows without a shard row", int(iv.obs_full_not_in_shard.sum()), "", "F2_index_vs_obs.tsv"),
        ("reduce", "analysis-set cells reduced", cs["n_analysis"], cs["n_lines"], "lines in column 4"),
        ("reduce", "cells whose UMIs differ from tscp_count", cs["umi_mismatch"], "", "F2_umi_mismatch_by_shard.tsv"),
        ("states", "foreign-marker cells removed before calling", cs["n_foreign"], "", "F2_foreign_cells.tsv"),
        ("states", "cells used for the final calls", cs["n_called"], cs["n_eligible"], "eligible lines (>=200 cells) in column 4"),
        ("states", "lines called", int(st.eligible.sum()), "", "F2_states.tsv"),
    ]
    if os.path.exists(f"{R}/F2_mosaic.tsv"):
        M = pd.read_csv(f"{R}/F2_mosaic.tsv", sep="\t")
        rows += [("mosaic", "lines with >=1,000 vehicle cells", int((st.n_vehicle_cells >= 1000).sum()), "", "F2_states.tsv"),
                 ("mosaic", "line x marker tests (Y)", int((M.marker == "Y").sum()), "", "F2_mosaic.tsv"),
                 ("mosaic", "line x marker tests (XIST)", int((M.marker == "XIST").sum()), "", "F2_mosaic.tsv"),
                 ("mosaic", "tests passing (b), (d) and (e)", int(M.flagged_b_d_e.sum()), "", "F2_mosaic.tsv"),
                 ("mosaic", "tests meeting the full mosaic rule", int(M.mosaic.sum()), "", "F2_mosaic.tsv")]
    pd.DataFrame(rows, columns=["stage", "step", "cells_or_count", "lines", "source"]).to_csv(f"{R}/F2_consort.tsv", sep="\t", index=False)

files = sorted(set(glob.glob(f"{R}/*")) - {f"{R}/F2_hashes.tsv"})
pd.DataFrame([dict(file=f, bytes=os.path.getsize(f), sha256=sha(f)) for f in files]).to_csv(f"{R}/F2_hashes.tsv", sep="\t", index=False)
print("hashed", len(files), "result files")
