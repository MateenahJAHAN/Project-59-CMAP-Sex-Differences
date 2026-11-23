"""Hash-checked access to every input of task G.

- Pinned git inputs (tasks D, E, D2): config/G_inputs_pinned.tsv (written by 00_fetch_inputs.py).
- Tensors: task D's D_hashes.tsv (level 5) and task E's E_hashes.tsv (level 3).
- Other files: config/G_external_sha256.tsv.
Every check is recorded for results/G/G_hash_check.tsv. A file that differs from its frozen sha256 always stops
the run. A missing optional file stops a real run; in a test run the analyses that need it are skipped and the
skip is recorded. Response-size files of tasks D and E are refused by name.
"""
import os, re, hashlib
import numpy as np, pandas as pd
import gcore as C
from gconfig import path as cpath

# The source replay found roundoff in seven ploidy values (maximum 8.9e-16).
# Both verified byte streams are pinned. See verification/source_replay/part_b/ploidy_control_comparison.json.
REPLAY_HASHES = {
    "taskD/D_line_table.tsv": "01922a9fc7d3fe85ea4cbc7ce6de1b07a11c323186d103c014ac55bb17d41b7a",
    "taskD/D_line_covariates.tsv": "9ef519446b4b2115a47b93aafc5005bed89ff59b3b7049cb19f27a6db1c47310",
}

FORBIDDEN_FILES = ("response_size", "plate_S")


class InputError(SystemExit):
    pass


def norm_name(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


class Inputs:
    def __init__(self, root, cfg, overrides, strict):
        self.root, self.cfg, self.strict = root, cfg, strict
        P = cfg.raw["paths"]
        g = lambda k: cpath(root, overrides.get(k) or P[k])
        self.dir, self.tD, self.tE, self.persig_dir = g("inputs"), g("tensors_D"), g("tensors_E"), g("persig_E")
        self.siginfo_path, self.prism_path = g("siginfo"), g("prism")
        self.repo_inputs = g("repo_inputs")
        self.pins = pd.read_csv(cpath(root, P["pins"]), sep="\t").set_index("input")
        self.ext = pd.read_csv(cpath(root, "config/G_external_sha256.tsv"), sep="\t").set_index("object")
        self.rows, self._done, self._logged = [], {}, set()
        self.D_hashes = self.table("taskD/D_hashes.tsv").set_index("path").sha256
        self.E_hashes = self.table("taskE/E_hashes.tsv", dtype=str).set_index("object")

    # ---------- hashing ----------
    def _sha(self, p):
        if p not in self._done:
            h = hashlib.sha256()
            with open(p, "rb") as f:
                for b in iter(lambda: f.read(1 << 22), b""):
                    h.update(b)
            self._done[p] = h.hexdigest()
        return self._done[p]

    def _check(self, obj, p, expected, source, alternatives=()):
        assert not any(f in p for f in FORBIDDEN_FILES), f"refusing to open {p}"
        if not os.path.exists(p):
            st, got = "missing", ""
        else:
            got = self._sha(p)
            st = "equal" if got == expected else ("equal (alternative)" if got in alternatives else "differs")
        if (obj, p) not in self._logged:
            self._logged.add((obj, p))
            self.rows.append(dict(object=obj, path=os.path.relpath(p, self.root) if p.startswith(self.root) else p,
                                  expected_sha256=expected, observed_sha256=got, status=st, expected_from=source))
        if st == "differs":
            raise InputError(f"INPUT CHECK: {obj} differs from its frozen sha256 ({source}); stop.")
        if st == "missing" and self.strict:
            raise InputError(f"INPUT CHECK: {obj} is missing ({p}); a real run needs every input.")
        return st

    def hash_table(self):
        return pd.DataFrame(self.rows)

    # ---------- pinned git inputs ----------
    def file(self, rel):
        p = os.path.join(self.dir, rel)
        st = self._check(rel, p, self.pins.loc[rel, "sha256"], "config/G_inputs_pinned.tsv", (REPLAY_HASHES.get(rel),))
        if st == "missing":
            raise InputError(f"pinned input {rel} missing; run src/G/00_fetch_inputs.py")
        return p

    def table(self, rel, **kw):
        return pd.read_csv(self.file(rel), sep="\t", **kw)

    def ext_file(self, obj, p, required=True):
        r = self.ext.loc[obj]
        st = self._check(obj, p, r.expected_sha256, "config/G_external_sha256.tsv: " + r.source_of_expected)
        if st == "missing" and required:
            raise InputError(f"{obj} missing ({p})")
        return p if st != "missing" else None

    def repo_table(self, rel, **kw):
        p = os.path.join(self.repo_inputs, rel)
        self.ext_file(f"inputs/{rel}", p)
        return pd.read_csv(p, sep="\t", **kw)

    # ---------- line information ----------
    def covariates(self):
        return C.read_covariates(self.file("taskD/D_line_covariates.tsv"))

    def landmarks(self):
        lm = self.table("taskD/D_landmarks.tsv")
        assert len(lm) == 978 and lm.gene_symbol.is_unique
        return lm

    def real_line_table_path(self):
        return self.file("taskD/D_line_table.tsv")

    # ---------- tensors ----------
    def d_tensor(self, name):
        p = os.path.join(self.tD, name)
        st = self._check(f"tensor {name}", p, self.D_hashes[f"data/D/tensors/{name}"], "task D results/D/D_hashes.tsv")
        return p if st != "missing" else None

    def e_object(self, obj):
        row = self.E_hashes.loc[obj]
        p = os.path.join(self.tE, os.path.basename(row.local_path))
        st = self._check(f"tensor {obj}", p, row.sha256, "task E results/E/E_hashes.tsv")
        return p if st != "missing" else None

    def panel(self, pn, cov, B_full):
        tp = self.d_tensor(f"{pn}_resp.npy")
        if tp is None:
            raise InputError(f"tensor {pn}_resp.npy missing")
        dm = self.d_tensor(f"{pn}_DMSO_resp.npy") if pn in ("P1", "P2") else None
        seed = self.cfg.run("key_seed_base") + self.cfg.const("panel_numbers")[pn]
        return C.Panel(pn, tp, self.file(f"taskD/tensors/{pn}_lines.tsv"), self.file(f"taskD/tensors/{pn}_compounds.tsv"),
                       self.landmarks().gene_symbol.tolist(), cov, seed, B_full, dm)

    def e_genes(self):
        """Task E's landmark order (gctx column order), which its tensors use. Task D's tensors use the order of
        D_landmarks.tsv instead, so E-built tensors are reordered on all three axes before they meet D's."""
        lm = self.table("taskE/E_landmarks.tsv")
        assert len(lm) == 978 and lm.gene_symbol.is_unique
        return lm.gene_symbol.tolist()

    def level3(self, obj, P, order="D"):
        """Level-3 (or E-built level-5) tensor `obj` in the axis order of panel P and task D's gene order (order 'D'),
        or in task E's own order (order 'E'). Returns (array, lines, compounds, genes) or None if missing (tests)."""
        axes = self.cfg.const("level3_axes").get(obj, obj.split("_", 1)[1])
        p = self.e_object(obj)
        pl, pc = self.e_object(f"{axes}_lines"), self.e_object(f"{axes}_compounds")
        if p is None or pl is None or pc is None:
            return None
        R = np.load(p).astype(np.float64)
        el = [x.rstrip("\n") for x in open(pl)]
        ec = [x.rstrip("\n") for x in open(pc)]
        eg = self.e_genes()
        assert R.shape == (len(el), len(ec), 978)
        if order == "E":
            return R, el, ec, eg
        assert set(el) == set(P.lines) and set(ec) == set(P.comps) and set(eg) == set(P.genes), f"{obj} axes differ from panel {P.pn}"
        li = [el.index(x) for x in P.lines]
        ci = [ec.index(x) for x in P.comps]
        gp = {g: i for i, g in enumerate(eg)}
        gi = [gp[g] for g in P.genes]
        return R[np.ix_(li, ci, gi)], list(P.lines), list(P.comps), list(P.genes)

    # ---------- per-signature level 5 (S11) ----------
    def persig_level5(self):
        p = os.path.join(self.persig_dir, "l5_P1.npy")
        pr = os.path.join(self.persig_dir, "l5_P1_rows.tsv")
        exp = self.E_hashes.loc["l5_P1"].sha256
        alt = self.ext.loc["l5_P1.npy (task D2, 1.11 uM only)"].expected_sha256
        st = self._check("persig l5_P1", p, exp, "task E E_hashes.tsv (alternative: task D2's 1.11 uM-only copy)", (alt,))
        if st == "missing":
            return None
        rows = pd.read_csv(pr, sep="\t")
        Z = np.load(p)
        assert Z.shape == (len(rows), 978) and rows.sig_id.is_unique
        return Z, rows.sig_id.tolist(), st, self.e_genes()

    def siginfo(self, cols):
        p = self.ext_file("siginfo_beta.txt", self.siginfo_path, required=self.strict)
        if p is None:
            return None
        return pd.read_csv(p, sep="\t", usecols=cols, dtype=str, low_memory=False)

    # ---------- regulons and gene sets ----------
    def collectri(self):
        return self.table("taskD/D_regulon_collectri_landmark.tsv")

    def panda(self):
        return self.repo_table("D/regulon_pandaLCL_top30_landmark.tsv")

    def lopesramos(self):
        return self.table("taskD/D_lopesramos_tfs.tsv")

    def chromosomes(self):
        t = self.repo_table("G/landmark_chromosomes.tsv", dtype={"gene_id": str})
        assert len(t) == 978
        return t

    def old_list(self):
        return self.repo_table("G/old_list_overlap.tsv")

    def activity(self):
        return self.table("taskD/D_activity_scores.tsv")

    def panel_definitions(self):
        return self.table("taskD/D_panel_definitions.tsv")

    # ---------- PRISM 19Q4 secondary screen (viability step) ----------
    def prism_auc(self, P, pdef):
        """AUC matrix (lines x compounds of panel P) with task C's rules: rows with passed_str_profiling true; per
        compound one screen (MTS010 when present, else the screen with the most lines); batches of one (compound,
        line) averaged. Compound match: the 13-character BRD prefix of a LINCS pert_id equals the PRISM broad_id
        prefix, else the normalized name (task C's normalization). Line match: PRISM depmap_id = BROAD_ID, else
        ccle_name. Returns (AUC (L, C) with NaN, coverage table) or None if the file is missing (test runs)."""
        p = self.ext_file("prism_secondary_dose_response_curve_parameters.csv", self.prism_path, required=self.strict)
        if p is None:
            return None
        df = pd.read_csv(p, usecols=["broad_id", "depmap_id", "ccle_name", "screen_id", "auc", "name", "passed_str_profiling"],
                         dtype=str, low_memory=False)
        df = df[df.passed_str_profiling.str.upper() == "TRUE"].copy()
        df["auc"] = pd.to_numeric(df.auc, errors="coerce")
        df = df[np.isfinite(df.auc)]
        df["brd13"] = df.broad_id.str[:13]
        df["nname"] = df.name.map(norm_name)
        d = pdef[(pdef.panel == P.pn) & pdef.included].set_index("cmap_name")
        bid = [str(x) if pd.notna(x) else "" for x in P.cov.BROAD_ID]
        ccl = [str(x) if pd.notna(x) else "" for x in P.cov.ccle_name]
        line_idx = {}
        for i in range(len(P.lines)):
            if bid[i]:
                assert ("id", bid[i]) not in line_idx, "duplicate BROAD_ID among panel lines"
                line_idx[("id", bid[i])] = i
            if ccl[i]:
                line_idx.setdefault(("cc", ccl[i]), i)
        A = np.full((len(P.lines), len(P.comps)), np.nan)
        cov = []
        for j, c in enumerate(P.comps):
            ids = {x.strip()[:13] for x in str(d.loc[c, "pert_ids"]).split(";") if x.strip()} if c in d.index else set()
            sub = df[df.brd13.isin(ids)]
            key = "BRD prefix"
            if sub.empty:
                sub = df[df.nname == norm_name(c)]
                key = "normalized name"
            if sub.empty:
                cov.append(dict(panel=P.pn, compound=c, prism_match="none", screen="", n_panel_lines_with_auc=0))
                continue
            nl = sub.groupby("screen_id").depmap_id.nunique()
            scr = "MTS010" if "MTS010" in nl.index else nl.sort_values(ascending=False, kind="stable").index[0]
            s2 = sub[sub.screen_id == scr]
            per = s2.groupby("depmap_id", dropna=False).agg(auc=("auc", "mean"), ccle_name=("ccle_name", "first")).reset_index()
            n = 0
            for _, r in per.iterrows():
                i = line_idx.get(("id", str(r.depmap_id)), line_idx.get(("cc", str(r.ccle_name))))
                if i is not None and np.isnan(A[i, j]):
                    A[i, j] = r.auc; n += 1
            cov.append(dict(panel=P.pn, compound=c, prism_match=key, screen=scr, n_panel_lines_with_auc=n))
        return A, pd.DataFrame(cov)
