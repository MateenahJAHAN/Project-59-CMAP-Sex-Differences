"""Step 1c: index tables and the completeness check of shard vehicle rows against obs_metadata."""
import glob, json, sys, time
import pandas as pd, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import runlog
t0 = time.time()
R = "results/F2"
js = [json.load(open(p)) for p in sorted(glob.glob("data/F2/index_parts/shard_*.json"))]
assert len(js) == 3388, len(js)
rg = pd.DataFrame([r for d in js for r in d["row_groups"]])
sh = pd.DataFrame([dict(shard=d["shard"], path=d["path"], bytes=d["size"], rows=d["num_rows"], row_groups=d["num_row_groups"],
                        footer_bytes=d["footer_bytes"]) for d in js])
g = rg.groupby("shard")
sh = sh.merge(pd.DataFrame(dict(
    plates=g.apply(lambda x: ",".join(sorted(set(x.plate_min) | set(x.plate_max), key=lambda p: int(p[5:])))),
    candidate_row_groups=g.candidate.sum(), vehicle_row_groups=g.vehicle_rows.apply(lambda v: int((v > 0).sum())),
    vehicle_rows=g.vehicle_rows.sum(), no_stats_row_groups=g.reason.apply(lambda v: int((v == "no_stats").sum())))).reset_index(), on="shard", validate="1:1")
assert len(sh) == 3388 and sh.shard.is_unique
sh["n_plates"] = sh.plates.str.count(",") + 1
sh.to_csv(f"{R}/F2_shard_summary.tsv", sep="\t", index=False)
idx = rg[rg.vehicle_rows > 0][["shard", "row_group", "rows", "vehicle_rows", "compressed_bytes", "plate_min", "plate_max"]]
idx.to_csv(f"{R}/F2_dmso_index.tsv", sep="\t", index=False)

veh = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob("data/F2/index_parts/veh_*.parquet"))], ignore_index=True)
assert len(veh) == int(rg.vehicle_rows.sum())
assert veh.BARCODE_SUB_LIB_ID.is_unique, "barcode repeated across shard rows"
obs = pd.read_parquet("data/F2/obs_dmso.parquet")
assert obs.BARCODE_SUB_LIB_ID.is_unique
full = obs[obs.pass_filter == "full"]
m = veh.merge(obs[["BARCODE_SUB_LIB_ID", "pass_filter", "cell_name", "cell_line", "plate", "sample"]], on="BARCODE_SUB_LIB_ID",
              how="outer", suffixes=("_shard", "_obs"), indicator=True, validate="1:1")
assert len(m) == len(set(veh.BARCODE_SUB_LIB_ID) | set(obs.BARCODE_SUB_LIB_ID))
both = m[m._merge == "both"]
cons = (both.plate_shard == both.plate_obs) & (both.sample_shard == both.sample_obs) & (both.cell_line_id == both.cell_line)
print("matched rows with inconsistent plate/sample/line:", int((~cons).sum()))
m["line"] = m.cell_name.where(m.cell_name.notna(), None)
m["plate"] = m.plate_obs.fillna(m.plate_shard)
cl = obs[["cell_line", "cell_name"]].drop_duplicates()
m.loc[m.line.isna(), "line"] = m.loc[m.line.isna(), "cell_line_id"].map(dict(zip(cl.cell_line, cl.cell_name)))
m["cat"] = np.select([(m._merge == "both") & (m.pass_filter == "full"), (m._merge == "both") & (m.pass_filter != "full"),
                      m._merge == "left_only", (m._merge == "right_only") & (m.pass_filter == "full")],
                     ["shard_and_obs_full", "shard_and_obs_minimal", "shard_not_in_obs", "obs_full_not_in_shard"], "obs_minimal_not_in_shard")
tab = m.groupby(["line", "plate", "cat"]).size().unstack(fill_value=0).reset_index()
for c in ["shard_and_obs_full", "shard_and_obs_minimal", "shard_not_in_obs", "obs_full_not_in_shard", "obs_minimal_not_in_shard"]:
    if c not in tab:
        tab[c] = 0
tab["obs_full"] = tab.shard_and_obs_full + tab.obs_full_not_in_shard
tab["shard_vehicle"] = tab.shard_and_obs_full + tab.shard_and_obs_minimal + tab.shard_not_in_obs
tab["inconsistent"] = 0
tab.to_csv(f"{R}/F2_index_vs_obs.tsv", sep="\t", index=False)
print(m.cat.value_counts().to_string())
print("shards:", len(sh), "with vehicle rows:", int((sh.vehicle_rows > 0).sum()), "vehicle row groups:", len(idx),
      "vehicle rows:", int(idx.vehicle_rows.sum()), "bytes of vehicle row groups:", int(idx.compressed_bytes.sum()))
print("rows per shard:", sh.rows.value_counts().head().to_dict(), "row groups:", sh.row_groups.value_counts().to_dict(),
      "multi-plate shards:", int((sh.n_plates > 1).sum()), "no-stats row groups:", int(sh.no_stats_row_groups.sum()))
runlog("s02b_index_check", t0, f"vehicle row groups {len(idx)}")
