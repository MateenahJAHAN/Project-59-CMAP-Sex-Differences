"""Step 3: one streaming pass over a sci-Plex UMI triplet matrix (gene_index cell_index count, 1-based).

Cell labels (plan section 2):
  sciPlex2: authors' hash call (top_oligo) with qval <= 0.05, top_to_second_best_ratio >= 10, hash_umis >= 10;
            vehicle = hash dose "0".
  sciPlex3: authors' cell_type and vehicle columns.
  sciPlex4: top hash from hashTable.out with >= 10 UMIs and >= 10x the second-best; vehicle = both doses "0".
Depth = sum over human (ENSG) genes.
Outputs (data/F1/sciplex_proc/): <tag>_vehicle_cells.parquet (vehicle cells only), <tag>_acc.npz
(deep-cell profile sums over all conditions by line x group, no condition label), <tag>_labels_summary.json.
"""
import sys, os, json, time, gzip
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import MARKERS, SEED, DEEP_MIN, NBINS, NGROUP

ROOT, TAG = sys.argv[1], sys.argv[2]
D = os.path.join(ROOT, "data/F1/sciplex")
OUT = os.path.join(ROOT, "data/F1/sciplex_proc"); os.makedirs(OUT, exist_ok=True)
PRE = {"sciPlex2": "GSM4150377_sciPlex2_A549_Transcription_Modulators_",
       "sciPlex3": "GSM4150378_sciPlex3_A549_MCF7_K562_screen_",
       "sciPlex4": "GSM4150379_sciPlex4_A549_MCF7_HDACi_"}[TAG]
GSM = PRE.split("_")[0]
summary = {}


def read_genes():
    p = os.path.join(D, PRE + "gene.annotations.txt.gz")
    g = pd.read_csv(p, sep=r"\s+", header=None, names=["id", "symbol"], engine="python")
    if g.iloc[0, 0] == "id":
        g = g.iloc[1:].reset_index(drop=True)
    g["human"] = g.id.str.startswith("ENSG")
    summary["n_genes_file"] = len(g); summary["n_genes_human"] = int(g.human.sum())
    midx = np.full(len(g), -1, np.int64)
    for i, s in enumerate(MARKERS):
        hit = np.where((g.symbol == s) & g.human)[0]
        assert len(hit) == 1, (s, hit)
        midx[hit[0]] = i
    return g, midx


def read_cells():
    ca = pd.read_csv(os.path.join(D, PRE + "cell.annotations.txt.gz"), sep="\t", header=None, names=["cell", "sample_id"])
    assert ca.cell.is_unique
    pd_ = pd.read_csv(os.path.join(D, f"{GSM}_{TAG}_pData.txt.gz"), sep=" ", low_memory=False)
    pd_ = pd_.rename(columns={"Cell": "cell"})
    assert pd_.cell.is_unique and len(pd_) == len(ca), (len(pd_), len(ca))
    n0 = len(ca)
    c = ca.merge(pd_, on="cell", how="left", validate="1:1")
    assert len(c) == n0 and c["n.umi"].notna().all()
    summary["n_cells_file"] = n0; summary["min_n_umi_pdata"] = int(c["n.umi"].min())
    c["line"] = pd.Series([None] * len(c), dtype=object); c["vehicle"] = c.get("vehicle", False); c["batch"] = "all"; c["bg_batch"] = "all"; c["well"] = ""
    if TAG == "sciPlex3":
        c["line"] = c["cell_type"]; c["vehicle"] = c["vehicle"].astype("boolean").fillna(False).astype(bool)
        c["batch"] = c["culture_plate"]; c["bg_batch"] = "rt" + c["rt_plate"].astype(str); c["well"] = c["top_oligo_W"].astype(str)
        assert set(c.line.dropna()) <= {"A549", "MCF7", "K562"}
        assert not (c.vehicle & c.line.isna()).any()
    elif TAG == "sciPlex2":
        ok = (c.qval <= 0.05) & (c.top_to_second_best_ratio >= 10) & (c.hash_umis >= 10) & c.top_oligo.notna()
        parts = c.top_oligo.str.rsplit("_", n=2)
        c["agent"] = parts.str[0]; c["dose"] = parts.str[1]; c["well"] = parts.str[2]
        assert set(c.agent.dropna()) <= {"BMS", "Dex", "Nutlin", "SAHA"}
        c.loc[ok, "line"] = "A549"; c["vehicle"] = ok & (c.dose == "0")
        summary["hash_pass"] = int(ok.sum())
    elif TAG == "sciPlex4":
        h = pd.read_csv(os.path.join(D, PRE + "hashTable.out.txt.gz"), sep="\t", header=None,
                        names=["sample_id", "cell", "hash", "axis", "umi"])
        h = h.sort_values(["cell", "umi"], ascending=[True, False])
        top = h.groupby("cell").head(1).set_index("cell")
        second = h.groupby("cell").umi.apply(lambda s: s.iloc[1] if len(s) > 1 else 0)
        top["second"] = second
        top["ratio"] = np.where(top.second > 0, top.umi / top.second.where(top.second > 0, 1), np.inf)
        ok = (top.umi >= 10) & (top.ratio >= 10)
        top = top[ok]
        p = top.hash.str.split("_")
        assert (p.str.len() == 7).all()
        lab = pd.DataFrame({"cell": top.index, "d1": p.str[0].values, "c1": p.str[1].values, "d2": p.str[2].values,
                            "c2": p.str[3].values, "ct": p.str[4].values, "plate": p.str[5].values, "w": p.str[6].values})
        assert set(lab.ct) <= {"A549", "MCF7"}
        n1 = len(c)
        c = c.merge(lab, on="cell", how="left", validate="1:1"); assert len(c) == n1
        c["line"] = c["ct"]; c["vehicle"] = (c.d1 == "0") & (c.d2 == "0") & c.ct.notna()
        c["batch"] = c["plate"]; c["well"] = c["plate"].astype(str) + "_" + c["w"].astype(str)
        summary["hash_cells_in_table"] = int(h.cell.nunique()); summary["hash_pass"] = int(ok.sum())
    c["vehicle"] = c["vehicle"].astype(bool)
    summary["n_assigned"] = int(c.line.notna().sum()); summary["n_vehicle"] = int(c.vehicle.sum())
    summary["vehicle_by_line"] = c[c.vehicle].line.value_counts().to_dict()
    return c


def main():
    t0 = time.time()
    g, midx = read_genes()
    c = read_cells()
    ncell, ngene = len(c), len(g)
    human = g.human.values
    lines = sorted(c.line.dropna().unique()); lix = {l: i for i, l in enumerate(lines)}
    li = c.line.map(lix).fillna(-1).astype(int).values
    rng = np.random.default_rng([SEED, int(TAG[-1])])
    by = rng.integers(0, NBINS, ncell); bx = rng.integers(0, NBINS, ncell)
    depth = np.zeros(ncell); ngenes = np.zeros(ncell, np.int64); mc = np.zeros((ncell, len(MARKERS)))
    acc = np.zeros((len(lines) * NGROUP, ngene)); cnt = np.zeros(len(lines) * NGROUP, np.int64)
    carry = None; last_c = -1; nrows = 0; maxc = 0
    path = os.path.join(D, PRE + "UMI.count.matrix.gz")
    reader = pd.read_csv(path, sep="\t", header=None, names=["g", "c", "v"], dtype={"g": np.int64, "c": np.int64, "v": np.float64},
                         chunksize=20_000_000, engine="c")

    def handle(ch):
        gi = ch.g.values - 1; ci = ch.c.values - 1; v = ch.v.values
        assert gi.min() >= 0 and gi.max() < ngene and ci.max() < ncell
        hm = human[gi]
        cells = np.unique(ci)
        depth[:] += np.bincount(ci[hm], weights=v[hm], minlength=ncell)
        ngenes[:] += np.bincount(ci[hm], minlength=ncell)
        mk = midx[gi]; m = mk >= 0
        mc[:] += np.bincount(ci[m] * len(MARKERS) + mk[m], weights=v[m], minlength=ncell * len(MARKERS)).reshape(ncell, len(MARKERS))
        # cells in this block are complete: accumulate deep-cell profiles
        dc = cells[(depth[cells] >= DEEP_MIN) & (li[cells] >= 0)]
        if len(dc):
            ys = mc[dc, 1:9].sum(1); xs = mc[dc, 0]
            gy = np.where(ys == 0, 0, 1 + by[dc]); gx = np.where(xs == 0, 1 + NBINS, 2 + NBINS + bx[dc])
            sy = np.full(ncell, -1); sx = np.full(ncell, -1)
            sy[dc] = li[dc] * NGROUP + gy; sx[dc] = li[dc] * NGROUP + gx
            for s_arr in (sy, sx):
                s = s_arr[ci]; k = (s >= 0) & hm
                flat = np.bincount(s[k] * ngene + gi[k], weights=v[k], minlength=acc.size)
                acc[:] += flat.reshape(acc.shape)
            np.add.at(cnt, sy[dc], 1); np.add.at(cnt, sx[dc], 1)

    for ch in reader:
        nrows += len(ch)
        cvals = ch.c.values
        assert cvals[0] >= last_c and np.all(np.diff(cvals) >= 0), "matrix not sorted by cell; single pass invalid"
        last_c = cvals[-1]
        if carry is not None:
            ch = pd.concat([carry, ch], ignore_index=True)
        lastcell = ch.c.values[-1]
        tail = ch.c.values == lastcell
        carry = ch[tail]; body = ch[~tail]
        if len(body):
            handle(body)
        print(f"{TAG} rows {nrows:,} cell {last_c:,}/{ncell:,} {time.time()-t0:.0f}s", flush=True)
    if carry is not None and len(carry):
        handle(carry)
    summary["matrix_rows"] = nrows
    c["depth"] = depth; c["n_genes"] = ngenes
    for j, s in enumerate(MARKERS):
        c[s] = mc[:, j].astype(np.int32)
    keep = ["cell", "line", "vehicle", "batch", "bg_batch", "well", "n.umi", "depth", "n_genes"] + MARKERS
    if TAG == "sciPlex3":
        keep += ["replicate", "time_point"]
    v = c.loc[c.vehicle, keep].copy()
    v["batch"] = v["batch"].astype(str); v["bg_batch"] = v["bg_batch"].astype(str)
    v.to_parquet(os.path.join(OUT, f"{TAG}_vehicle_cells.parquet"), index=False)
    # all-assigned-cell counts per line (no condition label), for CONSORT
    a = c[c.line.notna()]
    summary["assigned_by_line"] = a.line.value_counts().to_dict()
    summary["deep_assigned_by_line"] = a[a.depth >= DEEP_MIN].line.value_counts().to_dict()
    summary["depth_vs_n_umi_corr"] = float(np.corrcoef(c.depth, c["n.umi"])[0, 1])
    summary["seconds"] = round(time.time() - t0, 1)
    np.savez_compressed(os.path.join(OUT, f"{TAG}_acc.npz"), acc=acc, cnt=cnt, lines=np.array(lines),
                        gene_id=g.id.values.astype(str), symbol=g.symbol.values.astype(str), human=human)
    json.dump(summary, open(os.path.join(OUT, f"{TAG}_labels_summary.json"), "w"), indent=1, default=str)
    print(json.dumps(summary, default=str))


main()
