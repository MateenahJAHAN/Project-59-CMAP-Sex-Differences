"""tp_reducer (tahoe-pipeline): the one reducer used by the DMSO pass and by the full stream.

reduce(table) takes rows of a raw Tahoe-100M file (a whole file, or some row groups read remotely) and returns
  (a) per-cell markers for every kept cell: identity, depth (UMIs and genes, special tokens 0-2 excluded), raw counts
      of XIST, TSIX, the 8 Y genes, TXLNGY, PRKY, the 8 X gametologs, ACTB and GAPDH, y8 (sum of the 8 Y genes) and
      stratum = 2*(y8 >= 3) + (XIST >= 3); with file, row group and row in file;
  (b) a partial pseudobulk: summed counts of every gene per group (cell line, sample, stratum), sparse, with the
      number of cells and UMIs per group;
  (c) full rows (original columns) of the six reference lines.
merge_*() sum partial pseudobulks per plate, join cells to obs_metadata and run the checks:
  every cell's key equals the key of the obs row the index maps it to (so the join is 1:1 and order is as expected);
  pseudobulk counts sum to the per-cell depth for every (line, sample); row counts equal the footer counts;
  no duplicate barcodes; tscp_count equals total_umi; line, sample and plate agree with obs.
No sex state is called and no drug-response statistic is computed here.

Group key: gkey = (sample_number * 64 + line_index) * 4 + stratum   (line_index from index/lines.tsv)
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tp_common as C
import numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc

IDX = C.DATA + "/index"
META_COLS = ["plate", "sample", "drug", "cell_line_id", "BARCODE_SUB_LIB_ID"]


def build_lines_table():
    """50 Cellosaurus ids seen in obs, sorted; index used in gkey. Written once to index/lines.tsv."""
    path = IDX + "/lines.tsv"
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t")
    obs = pq.read_table(IDX + "/obs_compact.parquet", columns=["cell_line"])
    ids = sorted(pc.unique(pc.cast(obs.column("cell_line"), pa.string())).to_pylist())
    assert len(ids) == 50
    clm = pd.read_parquet(C.META["cell_line_metadata"])
    first = clm.drop_duplicates("Cell_ID_Cellosaur").set_index("Cell_ID_Cellosaur")
    t = pd.DataFrame(dict(line_index=np.arange(len(ids)), cell_line_id=ids))
    t["cell_name"] = first.reindex(ids)["cell_name"].values
    t["organ"] = first.reindex(ids)["Organ"].values
    t["depmap_id"] = first.reindex(ids)["Cell_ID_DepMap"].values
    six = C.six_line_ids()
    t["six_line"] = t.cell_line_id.isin(list(six.values()))
    assert t.cell_name.notna().all() and t.six_line.sum() == 6
    t.to_csv(path, sep="\t", index=False)
    return t


def gkey_decode(g):
    g = np.asarray(g, dtype=np.int64)
    return g // 256, (g // 4) % 64, g % 4          # sample number, line index, stratum


class Reducer:
    def __init__(self):
        self.lut = C.marker_lut()
        self.lines = build_lines_table()
        self.line_ids = self.lines.cell_line_id.tolist()
        self.six_ids = set(C.six_line_ids().values())

    def reduce(self, tb, file_idx, rg_of_row, row_in_file, keep_samples=None, want_six=True):
        t0 = time.time()
        n = tb.num_rows
        assert len(rg_of_row) == n and len(row_in_file) == n
        st = dict(file=int(file_idx), rows_in=int(n))
        g = tb.column("genes").combine_chunks()
        e = tb.column("expressions").combine_chunks()
        assert g.null_count == 0 and e.null_count == 0
        lens = pc.list_value_length(g).to_numpy(zero_copy_only=False).astype(np.int64)
        lens_e = pc.list_value_length(e).to_numpy(zero_copy_only=False).astype(np.int64)
        assert np.array_equal(lens, lens_e), "genes and expressions differ in length"
        gv = pc.list_flatten(g).to_numpy(zero_copy_only=False)
        ev = pc.list_flatten(e).to_numpy(zero_copy_only=False)
        assert len(gv) == lens.sum() == len(ev)
        assert gv.min() >= 0 and gv.max() < C.N_TOKENS
        starts = np.r_[0, np.cumsum(lens)[:-1]]
        row = np.repeat(np.arange(n), lens)
        special = gv < C.N_SPECIAL
        n_special = np.bincount(row[special], minlength=n)
        # expected: exactly one special entry per cell, the CLS token at position 0 with value -2
        first_ok = (lens > 0) & (gv[np.minimum(starts, len(gv) - 1)] == C.CLS_TOKEN) & (ev[np.minimum(starts, len(gv) - 1)] == -2.0)
        st["cells_without_cls_first"] = int((~first_ok).sum())
        st["cells_with_other_special_tokens"] = int((n_special != 1).sum())
        keep_e = ~special
        val = ev[keep_e]
        if val.size:
            bad = (val < 0) | (val != np.floor(val))
            if bad.any():
                raise ValueError("file %d: %d gene entries are negative or not integers" % (file_idx, int(bad.sum())))
        # explicit zero entries exist (seen: token 62712 with value 0 at the end of a cell). They add nothing to depth
        # or pseudobulk, and obs gene_count does not count them, so they are dropped and counted.
        zero = val == 0
        st["zero_entries"] = int(zero.sum())
        st["cells_with_zero_entries"] = int(len(np.unique(row[keep_e][zero])))
        if zero.any():
            keep_e[np.nonzero(keep_e)[0][zero]] = False
            val = ev[keep_e]
        rowk = row[keep_e]
        tok = gv[keep_e].astype(np.int64)
        vali = val.astype(np.int64)
        total = np.bincount(rowk, weights=vali, minlength=n).astype(np.int64)
        ngenes = np.bincount(rowk, minlength=n).astype(np.int64)
        assert total.max(initial=0) < 2**31
        # duplicated gene within a cell would double count in markers; check it is absent
        same_cell = rowk[1:] == rowk[:-1]
        st["cells_with_unsorted_or_repeated_tokens"] = int(len(np.unique(rowk[1:][same_cell & (tok[1:] <= tok[:-1])])))
        midx = self.lut[tok]
        sel = midx >= 0
        M = np.zeros((n, len(C.MARKERS)), np.int64)
        np.add.at(M, (rowk[sel], midx[sel]), vali[sel])
        iy = [C.MARKERS.index(x) for x in C.Y8]
        y8 = M[:, iy].sum(1)
        xist = M[:, C.MARKERS.index("XIST")]
        stratum = C.stratum_code(y8, xist)
        # metadata
        meta = {c: tb.column(c).combine_chunks() for c in META_COLS}
        for c in META_COLS:
            assert meta[c].null_count == 0, c
        snum = C.sample_num(meta["sample"]).astype(np.int64)
        pnum = C.plate_num(meta["plate"]).astype(np.int64)
        lid_arr = pc.index_in(meta["cell_line_id"], value_set=pa.array(self.line_ids))
        assert lid_arr.null_count == 0, "unknown cell line id"
        lid = lid_arr.to_numpy().astype(np.int64)
        key = C.encode_keys(meta["BARCODE_SUB_LIB_ID"])
        st["duplicate_barcodes_in_rows"] = int(n - len(np.unique(key)))
        # rows kept
        keep = np.ones(n, bool) if keep_samples is None else np.isin(snum, np.asarray(sorted(keep_samples)))
        st["rows_kept"] = int(keep.sum())
        kidx = np.nonzero(keep)[0]
        # (a) per-cell markers
        cols = {c: meta[c].take(pa.array(kidx)) for c in META_COLS}
        cols["key"] = pa.array(key[kidx])
        cols["total_umi"] = pa.array(total[kidx].astype(np.int32))
        cols["n_genes"] = pa.array(ngenes[kidx].astype(np.int32))
        for k, mname in enumerate(C.MARKERS):
            cols[mname] = pa.array(M[kidx, k].astype(np.int32))
        cols["y8"] = pa.array(y8[kidx].astype(np.int32))
        cols["stratum"] = pa.array(stratum[kidx])
        cols["file"] = pa.array(np.full(len(kidx), file_idx, np.int16))
        cols["rg"] = pa.array(np.asarray(rg_of_row)[kidx].astype(np.int16))
        cols["row"] = pa.array(np.asarray(row_in_file)[kidx].astype(np.int32))
        cells = pa.table(cols)
        # (b) partial pseudobulk over all genes
        gk_cell = (snum * 64 + lid) * 4 + stratum.astype(np.int64)
        ent_keep = keep[rowk]
        comb = gk_cell[rowk[ent_keep]] * 65536 + tok[ent_keep]
        u, inv = np.unique(comb, return_inverse=True)
        cnt = np.bincount(inv, weights=vali[ent_keep]).astype(np.int64)
        pb = pa.table({"gkey": pa.array((u // 65536).astype(np.int32)), "token": pa.array((u % 65536).astype(np.int32)),
                       "count": pa.array(cnt.astype(np.int32))})
        assert cnt.max(initial=0) < 2**31
        gku, ginv = np.unique(gk_cell[kidx], return_inverse=True)
        grp = pd.DataFrame(dict(gkey=gku.astype(np.int32),
                                n_cells=np.bincount(ginv).astype(np.int32),
                                umi=np.bincount(ginv, weights=total[kidx]).astype(np.int64)))
        pl_of_g = pd.Series(pnum[kidx]).groupby(ginv).agg(["min", "max"])
        assert (pl_of_g["min"] == pl_of_g["max"]).all(), "a sample spans two plates"
        grp["plate"] = pl_of_g["min"].values.astype(np.int8)
        # check: pseudobulk sums equal per-cell depth for every group, hence for every (line, sample)
        pb_sum = pd.Series(cnt).groupby(u // 65536).sum()
        chk = grp.set_index("gkey").umi
        st["pb_groups"] = int(len(grp))
        st["pb_nnz"] = int(len(u))
        st["pb_sum_mismatch_groups"] = int((pb_sum.reindex(chk.index).fillna(0).astype(np.int64) != chk).sum())
        assert st["pb_sum_mismatch_groups"] == 0
        assert int(cnt.sum()) == int(total[kidx].sum())
        # (c) six lines, original columns
        six = None
        if want_six:
            m6 = keep & pc.is_in(meta["cell_line_id"], value_set=pa.array(sorted(self.six_ids))).to_numpy(zero_copy_only=False)
            six = tb.filter(pa.array(m6))
            st["six_rows"] = int(m6.sum())
        st["umi_kept"] = int(total[kidx].sum())
        st["reduce_s"] = round(time.time() - t0, 2)
        return dict(cells=cells, pb=pb, groups=pa.Table.from_pandas(grp, preserve_index=False), six=six, stats=st)


def rows_meta_for(md, rgs):
    """row group and row-in-file of each row returned by read_row_groups(rgs), from file metadata."""
    starts = np.r_[0, np.cumsum([md.row_group(g).num_rows for g in range(md.num_row_groups)])]
    rg_of = np.concatenate([np.full(md.row_group(g).num_rows, g, np.int64) for g in rgs])
    row_of = np.concatenate([np.arange(starts[g], starts[g + 1]) for g in rgs])
    return rg_of, row_of, int(sum(md.row_group(g).num_rows for g in rgs))


def reduce_local_file(red, path, file_idx, keep_samples=None, want_six=True):
    pf = pq.ParquetFile(path)
    md = pf.metadata
    rgs = list(range(md.num_row_groups))
    tb = pf.read_row_groups(rgs)
    rg_of, row_of, nexp = rows_meta_for(md, rgs)
    assert tb.num_rows == nexp == md.num_rows
    out = red.reduce(tb, file_idx, rg_of, row_of, keep_samples=keep_samples, want_six=want_six)
    out["stats"]["rows_footer"] = int(md.num_rows)
    out["schema"] = pf.schema_arrow
    return out


def write_partials(out, base, six_path=None):
    """Per-file partial outputs. base = path prefix (…/train-XXXXX). Atomic per file."""
    C.atomic_write_table(out["cells"], base + ".cells.parquet", compression="zstd")
    C.atomic_write_table(out["pb"], base + ".pb.parquet", compression="zstd")
    C.atomic_write_table(out["groups"], base + ".groups.parquet", compression="zstd")
    if six_path is not None and out.get("six") is not None:
        C.atomic_write_table(out["six"], six_path, compression="snappy")


# ------------------------------------------------------------------------------------------------ merging
class ObsIndex:
    """obs columns and the shard-row -> obs-row map, for exact positional joins validated by key."""

    def __init__(self):
        o = pq.read_table(IDX + "/obs_compact.parquet",
                          columns=["key", "plate", "sample", "cell_line", "pass_full", "phase", "tscp_count",
                                   "gene_count", "pcnt_mito"])
        self.key = o.column("key").to_numpy()
        self.plate = o.column("plate").to_numpy()
        self.sample = o.column("sample").to_numpy()
        cl = o.column("cell_line").combine_chunks()
        self.cl_dict = np.array(cl.dictionary.to_pylist(), dtype=object)
        self.cl_code = cl.indices.to_numpy()
        self.pass_full = o.column("pass_full").to_numpy()
        self.phase = o.column("phase").to_numpy()
        self.tscp = o.column("tscp_count").to_numpy()
        self.gene_count = o.column("gene_count").to_numpy()
        self.mito = o.column("pcnt_mito").to_numpy()
        self.srow = pq.read_table(IDX + "/shard_obs_rows.parquet").column("obs_row").to_numpy()
        files = pd.read_parquet(IDX + "/files.parquet")
        self.fstart = np.r_[0, np.cumsum(files.rows.values)[:-1]]


def join_cells(cells_df, ob):
    """cells_df: pandas frame of per-cell markers with file,row,key. Adds obs columns and runs the checks."""
    grow = ob.fstart[cells_df["file"].values.astype(np.int64)] + cells_df["row"].values.astype(np.int64)
    orow = ob.srow[grow]
    chk = {}
    chk["cells"] = int(len(cells_df))
    chk["key_mismatch_vs_obs"] = int((ob.key[orow] != cells_df["key"].values).sum())
    chk["duplicate_barcodes"] = int(len(cells_df) - len(np.unique(cells_df["key"].values)))
    chk["not_pass_full"] = int((~ob.pass_full[orow]).sum())
    chk["plate_mismatch"] = int((ob.plate[orow] != C.plate_num(pa.array(cells_df["plate"].astype(str).values))).sum())
    chk["sample_mismatch"] = int((ob.sample[orow] != C.sample_num(pa.array(cells_df["sample"].astype(str).values))).sum())
    chk["line_mismatch"] = int((ob.cl_dict[ob.cl_code[orow]] != cells_df["cell_line_id"].astype(str).values).sum())
    chk["tscp_ne_total_umi"] = int((ob.tscp[orow] != cells_df["total_umi"].values).sum())
    chk["gene_count_ne_n_genes"] = int((ob.gene_count[orow] != cells_df["n_genes"].values).sum())
    cells_df = cells_df.copy()
    cells_df["obs_row"] = orow.astype(np.int64)
    cells_df["phase"] = pd.Categorical.from_codes(ob.phase[orow].astype(np.int8), C.PHASES)
    cells_df["obs_gene_count"] = ob.gene_count[orow]
    cells_df["obs_tscp_count"] = ob.tscp[orow]
    cells_df["pcnt_mito"] = ob.mito[orow]
    return cells_df, chk


class DenseMerge:
    """Sum sparse partial pseudobulks for a fixed set of groups in a dense int64 (small) or int32 (large) array."""

    def __init__(self, gkeys, dtype=np.int64):
        self.gk = np.asarray(sorted(set(int(x) for x in gkeys)), dtype=np.int64)
        self.dense = np.zeros((len(self.gk), C.N_TOKENS), dtype=dtype)
        self.flat = self.dense.reshape(-1)

    def add(self, pb):
        g = pb.column("gkey").to_numpy().astype(np.int64)
        m = np.isin(g, self.gk)
        g = g[m]
        t = pb.column("token").to_numpy().astype(np.int64)[m]
        c = pb.column("count").to_numpy()[m]
        pos = np.searchsorted(self.gk, g)
        idx = pos * C.N_TOKENS + t
        # inside one partial each (gkey, token) appears once, so a fancy-index add is exact
        assert len(np.unique(idx)) == len(idx)
        self.flat[idx] += c.astype(self.dense.dtype)

    def result(self, block=1024):
        """nonzero entries, sorted by gkey then token (row-major order), extracted in blocks of groups"""
        pg, pt, pcnt = [], [], []
        for a in range(0, len(self.gk), block):
            sub = self.dense[a:a + block]
            r, t = np.nonzero(sub)
            c = sub[r, t]
            assert c.max(initial=0) < 2**31 and c.min(initial=1) > 0
            pg.append(self.gk[a + r].astype(np.int32))
            pt.append(t.astype(np.int32))
            pcnt.append(c.astype(np.int32))
        cat = lambda L: np.concatenate(L) if L else np.zeros(0, np.int32)
        return pa.table({"gkey": pa.array(cat(pg)), "token_id": pa.array(cat(pt)), "count": pa.array(cat(pcnt))})


def sorted_runs(a):
    """values, start index and length of the runs of a sorted 1-D array"""
    a = np.asarray(a)
    if len(a) == 0:
        return a[:0], np.zeros(0, np.int64), np.zeros(0, np.int64)
    starts = np.r_[0, np.flatnonzero(a[1:] != a[:-1]) + 1]
    return a[starts], starts, np.diff(np.r_[starts, len(a)])


def finalize_tables(pb_g, groups_df, cells_df, lines):
    """Turn gkey-coded pseudobulk and groups into the documented output tables, and run the checks that need
    both cells and pseudobulk. Returns dict of tables plus check dict."""
    chk = {}
    g_all = pb_g.column("gkey").to_numpy()
    assert bool(np.all(g_all[1:] >= g_all[:-1])), "pseudobulk entries must be sorted by group"
    gv_, gstart, glen = sorted_runs(g_all)
    s, li, stt = gkey_decode(gv_)
    lid = np.array(lines.cell_line_id.tolist(), dtype=object)
    us, sinv = np.unique(s, return_inverse=True)
    cnt64 = pb_g.column("count").to_numpy().astype(np.int64)
    gsum = np.add.reduceat(cnt64, gstart) if len(gstart) else np.zeros(0, np.int64)
    del cnt64
    pbt = pa.table({
        "sample": pa.DictionaryArray.from_arrays(pa.array(np.repeat(sinv, glen).astype(np.int32)),
                                                 pa.array(["smp_%04d" % x for x in us])),
        "cell_line_id": pa.DictionaryArray.from_arrays(pa.array(np.repeat(li, glen).astype(np.int32)), pa.array(list(lid))),
        "stratum": pa.array(np.repeat(stt, glen).astype(np.int8)),
        "token_id": pb_g.column("token_id"),
        "count": pb_g.column("count"),
    })
    gr = groups_df.groupby("gkey", as_index=False).agg(plate=("plate", "first"), n_cells=("n_cells", "sum"),
                                                       umi=("umi", "sum"))
    gs, gli, gst = gkey_decode(gr.gkey.values)
    gr["sample"] = ["smp_%04d" % x for x in gs]
    gr["cell_line_id"] = lid[gli]
    gr["cell_name"] = lines.cell_name.values[gli]
    gr["stratum"] = gst.astype(np.int8)
    gr["y8_ge3"] = (gst & 2) > 0
    gr["xist_ge3"] = (gst & 1) > 0
    gr["plate"] = "plate" + gr.plate.astype(int).astype(str)
    gr["n_genes_detected"] = pd.Series(glen, index=gv_).reindex(gr.gkey.values).fillna(0).astype(np.int64).values
    gr["pb_umi"] = pd.Series(gsum, index=gv_).reindex(gr.gkey.values).fillna(0).astype(np.int64).values
    chk["pb_groups_without_cells"] = int(len(set(gv_.tolist()) - set(gr.gkey.tolist())))
    chk["groups"] = int(len(gr))
    chk["group_pb_sum_ne_umi"] = int((gr.pb_umi != gr.umi).sum())
    # (line, sample): pseudobulk total over strata and genes equals the sum of per-cell depth
    c_ls = cells_df.groupby(["cell_line_id", "sample"], observed=True).total_umi.sum()
    p_ls = gr.groupby(["cell_line_id", "sample"]).pb_umi.sum()
    j = pd.concat([c_ls.rename("cells"), p_ls.rename("pb")], axis=1).fillna(-1)
    chk["line_sample_pairs"] = int(len(j))
    chk["line_sample_pb_sum_ne_cell_depth"] = int((j.cells != j.pb).sum())
    n_ls = cells_df.groupby(["cell_line_id", "sample"], observed=True).size()
    g_ls = gr.groupby(["cell_line_id", "sample"]).n_cells.sum()
    chk["line_sample_cells_ne_group_cells"] = int((n_ls.reindex(g_ls.index).fillna(-1) != g_ls).sum())
    # phase counts per (line, sample, stratum)
    ph = cells_df.groupby(["plate", "sample", "cell_line_id", "stratum", "phase"], observed=True).size()
    ph = ph.rename("n_cells").reset_index()
    ph["phase"] = ph["phase"].astype(str)
    groups = gr[["plate", "sample", "cell_line_id", "cell_name", "stratum", "y8_ge3", "xist_ge3", "n_cells", "umi",
                 "n_genes_detected"]].sort_values(["sample", "cell_line_id", "stratum"]).reset_index(drop=True)
    return dict(pseudobulk=pbt, groups=groups, phase=ph, checks=chk)


if __name__ == "__main__":
    # python -I tp_reducer.py <local raw file> <out prefix> [file index]
    red = Reducer()
    path, base = sys.argv[1], sys.argv[2]
    fi = int(sys.argv[3]) if len(sys.argv) > 3 else int(os.path.basename(path)[6:11])
    out = reduce_local_file(red, path, fi)
    write_partials(out, base, six_path=base + ".six.parquet")
    print(json.dumps(out["stats"]))
