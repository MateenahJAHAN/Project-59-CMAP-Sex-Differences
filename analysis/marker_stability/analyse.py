"""Analyse a sealed Tahoe count stream. Requires numpy, pandas, scipy, matplotlib, pyarrow."""
import argparse
import ast
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ALPHA = 0.05


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(root):
    stream = root / "stream"
    hashes = dict(line.strip().split(maxsplit=1)[::-1]
                  for line in (root / "stream_SHA256SUMS.txt").read_text().splitlines())
    state_path = root / "partA_local/partA_local_states.tsv"
    assert sha(state_path) == (root / "partA_local/partA_local_states.sha256").read_text().split()[0]
    samples_path = root / "index/meta/sample_metadata.parquet"
    drivers_path = root / "index/meta/cell_line_metadata.parquet"
    provenance = {str(f.relative_to(root)): sha(f) for f in [state_path, samples_path, drivers_path]}
    data = []
    totals = {}
    for plate in range(1, 15):
        rel = f"plates/plate{plate}/pseudobulk_groups.tsv"
        path = stream / rel
        assert sha(path) == hashes[rel]
        provenance["stream/" + rel] = hashes[rel]
        check = json.loads((path.parent / "MERGE_OK.json").read_text())
        assert check["passed"] and not check["partial_ok"]
        table = pd.read_csv(path, sep="\t")
        assert not table.duplicated(["plate", "sample", "cell_line_id", "stratum"]).any()
        assert set(table.plate) == {f"plate{plate}"}
        cells = table.groupby(["sample", "cell_line_id"]).n_cells.sum()
        assert cells.sum() == check["cells"] == check["rows_expected_from_index"]
        totals[str(plate)] = int(cells.sum())
        data.append(cells)
    counts = pd.concat(data).unstack("cell_line_id", fill_value=0).sort_index()
    states = pd.read_csv(state_path, sep="\t").set_index("cellosaurus_id").sort_index()
    assert counts.columns.tolist() == states.index.tolist() and len(states) == 50 and counts.sum().sum() == sum(totals.values())
    meta = pd.read_parquet(samples_path).set_index("sample").loc[counts.index]
    assert meta.index.is_unique and meta.groupby("plate").size().eq(96).all() and len(meta) == 1344
    dose = meta.drugname_drugconc.map(ast.literal_eval)
    assert all(len(x) == 1 and x[0][2] == "uM" for x in dose)
    assert all(x[0][0] == d for x, d in zip(dose, meta.drug))
    meta["dose"] = [x[0][1] for x in dose]
    meta["plate_id"] = meta.plate.str.extract(r'(\d+)').astype(int)
    meta["control"] = meta.drug.eq("DMSO_TF")
    assert meta[meta.control].groupby("plate").size().eq(2).all()
    states["positive"] = states.xist_high | states.y_high
    drivers = pd.read_parquet(drivers_path)
    hits = drivers[(drivers.Driver_Gene_Symbol == "KRAS") & (drivers.Driver_ProtEffect_or_CdnaEffect == "p.G12C")]
    states["g12c"] = states.index.isin(hits.Cell_ID_Cellosaur)
    return counts, states, meta, totals, provenance


def well_responses(counts, meta):
    logshares = np.log2((counts + 0.5).div(counts.sum(axis=1) + 0.5 * counts.shape[1], axis=0))
    ctr = logshares.loc[meta.control]
    ctrl_means = ctr.groupby(meta.loc[ctr.index, "plate"]).mean()
    centering = ctrl_means.loc[meta.plate].set_axis(meta.index)
    return logshares - centering


def design(meta, responses):
    rows = meta.loc[~meta.control]
    avg = responses.loc[rows.index].groupby([rows.drug, rows.dose]).mean()
    per_dose = {float(d): frame.droplevel("dose") for d, frame in avg.groupby(level="dose")}
    drugs = per_dose[5].index.intersection(per_dose[0.05].index)
    return per_dose[5].loc[drugs] - per_dose[0.05].loc[drugs], per_dose


def groups(states, eligible, col="positive", strat=("organ",), minimum=2):
    panel = states.loc[eligible].copy()
    panel["stratum"] = panel[list(strat)].astype(str).agg(" | ".join, axis=1)
    sizes = panel.groupby(["stratum", col]).size().unstack(fill_value=0).reindex(columns=[False, True], fill_value=0)
    panel = panel[panel.stratum.isin(sizes.index[sizes.min(axis=1) >= minimum])]
    return [(name, g.index[g[col]].tolist(), g.index[~g[col]].tolist())
            for name, g in panel.groupby("stratum")]


def contrast(outcome, group):
    if not group:
        return np.nan
    return sum(outcome[pos].mean(axis=1) - outcome[neg].mean(axis=1)
               for _, pos, neg in group) / len(group) if isinstance(outcome, pd.DataFrame) else float(
        np.mean([outcome[pos].mean() - outcome[neg].mean() for _, pos, neg in group]))


def statistics(y, grouping):
    diff = contrast(y, grouping)
    var = []
    df = []
    space = 1
    arrays = []
    for _, pos, neg in grouping:
        for ids in (pos, neg):
            var.append(y[ids].var(ddof=1) / len(ids) / len(grouping) ** 2)
            df.append(len(ids) - 1)
        a = y[pos + neg].to_numpy()
        arrays.append((a, len(pos)))
        space *= math.comb(len(a), len(pos))
    var = np.asarray(var)
    se = np.sqrt(var.sum())
    freedom = var.sum() ** 2 / (var ** 2 / df).sum()
    radius = stats.t.ppf(1 - ALPHA / 2, freedom) * se
    rng = np.random.default_rng(591121)
    if space <= 99999:
        values = []
        for a, k in arrays:
            sums = np.array([sum(a[list(i)]) for i in itertools.combinations(range(len(a)), k)])
            values.append(sums / k - (a.sum() - sums) / (len(a) - k))
        null = np.array(list(itertools.product(*values))).mean(axis=1)
        p = np.mean(np.abs(null) + 1e-12 >= abs(diff))
    else:
        null = np.zeros(99999)
        for a, k in arrays:
            idx = np.argpartition(rng.random((len(null), len(a))), k-1, axis=1)[:, :k]
            sums = a[idx].sum(axis=1)
            null += (sums/k - (a.sum()-sums)/(len(a)-k))/len(arrays)
        p = (1 + np.count_nonzero(np.abs(null) + 1e-12 >= abs(diff))) / (1 + len(null))
    return dict(estimate=diff, ci_low=diff-radius, ci_high=diff+radius, df=freedom,
                permutation_p=float(p), draws=len(null), total_labelings=space,
                n_lines=sum(len(a)+len(b) for _, a, b in grouping), n_strata=len(grouping))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path.home()/"p59_data/tahoe")
    p.add_argument("--out", type=Path, default=Path(__file__).resolve().parent/"results")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    def table(df, name):
        df.to_csv(args.out/(name+".tsv"), sep="\t", index=False, float_format="%.9g")
    def dump(data, name):
        (args.out/(name+".json")).write_text(json.dumps(data, indent=2, allow_nan=False)+"\n")
    counts, states, meta, plate_totals, provenance = load(args.root)
    min_control = counts.loc[meta.control].min(axis=0)
    eligible = min_control.ge(20)
    g = groups(states, eligible)
    response = well_responses(counts, meta)
    delta, by_dose = design(meta, response)
    y = delta.mean()
    prim = statistics(y, g)
    treated = ~meta.control
    primary_ids = [i for _, a, b in g for i in a+b]
    line_output = states[["cell_name","organ","annotation","y_high","xist_high","confidence","positive","g12c"]].copy()
    line_output["lowest_vehicle_count"] = min_control
    line_output["in_primary"] = line_output.index.isin(primary_ids)
    line_output["mean_high_minus_low_log2_share"] = y
    table(line_output.reset_index(), "lines")
    table(delta.rename_axis(index="compound").reset_index(), "high_minus_low_by_drug_and_line")
    group_rows = []
    for s,a,b in g:
        group_rows.append(dict(organ=s, positive=len(a), neither=len(b), estimate=contrast(y,[(s,a,b)])))
    table(pd.DataFrame(group_rows), "organs")
    splits = []
    for plate, d in meta.loc[meta.control].groupby("plate"):
        a,b = sorted(d.index, key=lambda x:int(x.split("_")[1]))
        splits.append((plate,response.loc[b]-response.loc[a]))
    splits = pd.DataFrame(dict(splits)).T
    control = statistics(splits.mean(),g)
    table(contrast(splits,g).rename("contrast").rename_axis("plate").reset_index(), "vehicle_pairs")
    positive_group = groups(states, eligible, "g12c", minimum=1)
    positive_result = {"estimate":contrast(delta.loc["Adagrasib"],positive_group),
                       "groups":[[k, states.loc[a,"cell_name"].tolist(),states.loc[b,"cell_name"].tolist()] for k,a,b in positive_group]}
    sens = []
    def record(label, gg, vector):
        sens.append(dict(analysis=label, estimate=contrast(vector,gg), n_organs_or_strata=len(gg),
                         n_lines=sum(len(a)+len(b) for _,a,b in gg)))
    record("primary",g,y)
    for threshold in (5,50):
        record(f"minimum {threshold} vehicle cells",groups(states,min_control.ge(threshold)),y)
    record("one line per group sufficient", groups(states,eligible,minimum=1),y)
    record("exclude low confidence",groups(states, eligible & states.confidence.isin(["medium","high"])),y)
    joint = groups(states,eligible & states.annotation.isin(["Female","Male"]),strat=("organ","annotation"))
    record("within organ and recorded sex",joint,y)
    for org,_,_ in g:
        record("leave out organ "+org,[grp for grp in g if grp[0]!=org],y)
    for cell in primary_ids:
        gg = [(org,[i for i in a if i!=cell],[i for i in b if i!=cell]) for org,a,b in g]
        record("leave out "+states.loc[cell,"cell_name"],gg,y)
    drop = meta.plate_id.ne(13)
    delta13,_ = design(meta.loc[drop],response.loc[drop])
    record("exclude plate 13",g,delta13.mean())
    tr = meta.loc[treated]
    within = response.loc[tr.index].groupby([tr.drug,tr.plate,tr.dose]).mean()
    wl = within.xs(.05,level="dose"); wh = within.xs(5.,level="dose")
    same = wh.index.intersection(wl.index)
    paired = (wh.loc[same] - wl.loc[same]).groupby(level="drug").mean()
    record("extreme doses on same plate",g,paired.mean())
    within13 = same[same.get_level_values("plate")!="plate13"]
    paired13 = (wh.loc[within13]-wl.loc[within13]).groupby(level="drug").mean()
    if len(paired13):
        record("extreme doses on same plate excluding plate 13",g,paired13.mean())
    table(pd.DataFrame(sens),"sensitivity")
    dose_out = pd.concat({str(d):contrast(a,g) for d,a in by_dose.items()},axis=1)
    dose_out["high_minus_low"] = contrast(delta,g)
    table(dose_out.reset_index(),"compound_contrasts")
    wells = meta[["plate","drug","dose","control"]].copy()
    wells["captured_cells"] = counts.sum(axis=1)
    wells["state_contrast"] = contrast(response,g)
    table(wells.reset_index(),"well_contrasts")
    reps = []
    for (drug,dose), r in wells.loc[treated].groupby(["drug","dose"]):
        means = r.groupby("plate").state_contrast.mean()
        ordered = sorted(means.index,key=lambda k:int(k[5:]))
        if len(means) > 1:
            subset = (ordered == ["plate6","plate14"])
            reps.append(dict(compound=drug,dose=dose,first_plate=ordered[0],other_plates=",".join(ordered[1:]),
                             only_6_14=subset, first=means[ordered[0]],others=means[ordered[1:]].mean()))
    rep=pd.DataFrame(reps)
    table(rep,"repeated_conditions")
    rep_stats={}
    for name,rows in [("all",rep),("only_6_14",rep[rep.only_6_14])]:
        rep_stats[name]=dict(n=len(rows),pearson=float(rows["first"].corr(rows.others)),
                             median_abs_difference=float((rows["first"]-rows.others).abs().median()))
    stats_out=dict(primary=prim,vehicle=control,positive_control=positive_result,
                   counts=dict(plates=14,cell_lines=len(states),positive_lines=int(states.positive.sum()),
                               cells=int(counts.sum().sum()),treated_cells=int(counts.loc[treated].sum().sum()),
                               control_cells=int(counts.loc[meta.control].sum().sum()),wells=len(meta),
                               treated_wells=int(treated.sum()),compounds=len(delta),treated_zero_entries=int(counts.loc[treated].eq(0).sum().sum()),
                               primary_treated_zero_entries=int(counts.loc[treated,primary_ids].eq(0).sum().sum()),
                               eligible_vehicle_lines=int(eligible.sum()), same_plate_compounds=len(paired),
                               same_plate_compounds_excluding_13=len(paired13),compounds_excluding_13=len(delta13),
                               n_primary_group1=sum(map(lambda row:len(row[1]),g)),n_primary_group0=sum(map(lambda row:len(row[2]),g))),
                   repeated_conditions=rep_stats,plate_cells=plate_totals,
                   in_both_sex_and_organ=[(k,len(a),len(b)) for k,a,b in joint])
    dump(stats_out,"summary")
    dump(provenance,"inputs_sha256")
    print(json.dumps(stats_out,indent=2))


if __name__=="__main__":
    main()
