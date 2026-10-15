# Where every input file comes from

Nothing in this phase is simulated. Every number traces back to one of the two
public datasets below.

---

## 1. CMap LINCS 2020 — the drug side

The Broad Institute added drugs to human cancer cell lines, waited 6, 12 or
24 hours, and measured how 12,328 genes reacted. Each experiment is one
**signature**. The value for a gene is a z-score: how far that gene moved
compared with untreated cells on the same plate. CMap caps it at ±10.

Portal: <https://clue.io/data/CMap2020#LINCS2020>
Build path: `s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/`

| File | Size | What it holds | In this repo? |
|------|------|---------------|---------------|
| `level5_beta_trt_cp_n720216x12328.gctx` | 33 GB | the z-score matrix: 12,328 genes × 720,216 signatures | **no** |
| `siginfo_beta.txt` | ~444 MB | one row per signature: drug, cell line, dose, time, QC flags | **no** |
| `cellinfo_beta.txt` | ~38 KB | one row per cell line: donor sex, lineage | yes → `/cmap/` |
| `geneinfo_beta.txt` | ~1.1 MB | Entrez ↔ Ensembl ↔ gene symbol | yes → `/cmap/` |
| `compoundinfo_beta.txt` | small | drug mechanism reference | yes → `/cmap/` |

The two large files stay in Google Drive. GitHub caps single files at 100 MB.

The .gctx used in this project was fetched from:

```
https://cmapreg.s3.amazonaws.com/galenAI/level5_beta_trt_cp_n720216x12328.gctx
```

Read it with [cmapPy](https://github.com/cmap/cmapPy) — never load all 33 GB.
Pass the signature IDs you want:

```python
from cmapPy.pandasGEXpress.parse import parse
expr = parse(gctx_file, cid=valid_ids).data_df
```

**Paper:** Subramanian A, Narayan R, Corsello SM, et al. *A Next Generation
Connectivity Map: L1000 Platform and the First 1,000,000 Profiles.* Cell. 2017;
171(6):1437–1452.e17. doi:[10.1016/j.cell.2017.10.049](https://doi.org/10.1016/j.cell.2017.10.049)

---

## 2. GTEx v8 — the sex-biased gene side

GTEx measured gene expression in real human tissue from real donors and listed,
per tissue, the genes expressed differently in men and women. This is our
reference list — the thing we test the drug's gene list against.

Portal: <https://gtexportal.org/home/downloads/adult-gtex>
Archive: `GTEx_Analysis_v8_sbgenes.tar.gz`

| File | What it holds | In this repo? |
|------|---------------|---------------|
| `signif.sbgenes.txt` | the sex-biased genes of every tissue, with `effsize` | **no** (see `/gtex/`) |
| per-tissue `*.sbgenes.txt` | the same, split by tissue | no |

Two columns matter:

- **`effsize`** — the sign decides the label. Positive → female-biased.
  Negative → male-biased. This is what `split_sets()` uses to cut each tissue
  into a FEMALE half and a MALE half.
- **`lfsr`** — local false sign rate. GTEx calls a gene sex-biased at
  `lfsr ≤ 0.05`.

We keep only tissues with at least 15 sex-biased genes, which leaves **44
tissues**. Library sizes vary a lot — 473 genes for
`Cells_EBV.transformed_lymphocytes`, 2,739 for `Brain_Hypothalamus` — and that
variation is one reason a pair can end up as Category D.

**Paper:** Oliva M, Muñoz-Aguirre M, Kim-Hellmuth S, et al. *The impact of sex
on gene expression across human tissues.* Science. 2020 Sep 11;369(6509).
doi:[10.1126/science.aba3066](https://doi.org/10.1126/science.aba3066)

How they got the list: a linear model per tissue with surrogate variables (to
absorb differences in cell-type composition), then MASH — multivariate adaptive
shrinkage, which borrows strength across the 44 tissues — then `lfsr ≤ 0.05`.
They report 13,294 sex-biased genes, 473 to 4,558 per tissue, with a median
autosomal fold change of about 1.04.

**Worth quoting in any discussion:** the paper itself warns that cell lines are
not representative of real tissue. Our Category A is exactly that comparison, so
the caveat is ours too.

---

## 3. Derived reference files

Built once from the two sources above and stored as JSON in
`/Phase5_PerSignature_Categories/reference/`.

| File | What it is | Expected |
|------|-----------|----------|
| `cmap_to_gtex.json` | CMap cell lineage → GTEx tissue(s). Defines "home tissue" | **20** lineages |
| `sb_library_sets.json` | GTEx tissue → its sex-biased Ensembl IDs | **44** tissues |
| `gid_to_ensg.json` | CMap gene_id → Ensembl ID. The translator between the two datasets | 12,277 genes |
| `top50_names.json` | the 50 drugs | 50 |

`sb_library_sets` **must be loaded as Python sets**, not lists:

```python
sb_library_sets = {k: set(v) for k, v in json.load(open(path)).items()}
```

JSON has no set type, so a plain `json.load` gives lists. The intersection
inside `run_gsea` then fails, the failure is swallowed by a `try/except`, and
every pair quietly becomes C or D with no error printed. This is the Phase 5
bug. Check the type before every run.

---

## 4. Tools

| | Version | Used for |
|---|---|---|
| gseapy | 1.3.1 | `gseapy.prerank`, 1,000 permutations, `min_size=15`, seed 42 |
| cmapPy | ≥ 4.0.1 | reading the .gctx |
| pandas / numpy | ≥ 2.0 / ≥ 1.24 | everything else |
| matplotlib | ≥ 3.7 | the figures |
| Python | ≥ 3.12 | the plot code uses PEP 701 f-strings |

**GSEA method:** Subramanian A, Tamayo P, Mootha VK, et al. *Gene set enrichment
analysis: a knowledge-based approach for interpreting genome-wide expression
profiles.* PNAS. 2005;102(43):15545–15550.
doi:[10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102)

`min_size=15` is the gseapy and Broad default. It is a property of GSEA, not a
filter chosen for this project.
