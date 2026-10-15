# Project 59 — Sex Differences in Drug Response

**Fatima Fellowship 2026**
Mentor: Dr. Marouen Ben Guebila, Viswanathan Lab, Dana-Farber Cancer Institute / Harvard
Fellow: Mateenah Jahan

---

## The question

Do cells from men and cells from women respond differently to the same drug?

Cell lines grow in a dish. They make no hormones. So if we see a difference, it
comes from the genes, not from hormones.

---

## The idea in one picture

We have two lists and we ask whether they overlap more than chance allows.

```
LIST 1  —  the drug          LIST 2  —  the tissue
(CMap)                       (GTEx)

every gene, ranked from       the genes that differ between
most activated to most        men and women in that tissue,
suppressed by the drug        split into a female half and
                              a male half

            \                /
             \              /
              GSEA asks: do the female-biased genes
              sit bunched at one end of the drug's
              ranked list, or scattered at random?
```

If they are bunched, the drug is acting on sex-biased biology. If scattered,
it is not.

Direction matters and is easy to get backwards. The **drug's ranked list is the
source** (`rnk=`), the **tissue's gene list is the target** (`gene_sets=`). We
enrich the drug against the tissue, never the reverse.

---

## Three words to keep straight

**Cell line** — cells from one person. THP1 came from one male leukemia patient.
One person means one sex.

**Lineage** — the organ that cell line came from. THP1 came from blood, so its
lineage is blood. Many cell lines share a lineage.

**Signature** — one experiment: one cell line, one drug, one dose, one time.
One cell line gives many signatures because dose and time keep changing.

```
signature  =  drug  ×  cell line  ×  dose  ×  time
```

GTEx tissue is **not** in that formula. Tissue only enters when we compare. So a
drug's signature count does not change from tissue to tissue — bortezomib has
7,300 signatures whether we test it against brain or against fat.

Example: bortezomib in the blood lineage covers 15 cell lines and 255 signatures.

---

## The four categories

| Category | What it means |
|----------|---------------|
| **A** | Significant, and the tissue is the cell line's own tissue (home) |
| **B** | Significant, but a different tissue (cross) |
| **C** | We tested it, nothing came out |
| **D** | We could not test it at all |

D is not "skipped". Every drug was attempted. D means the test could not run —
too few signatures, nothing past `|z| > 2`, no matching GTEx tissue, or fewer
than 15 genes left after matching.

---

## Where the data comes from

| Source | What we take from it |
|--------|---------------------|
| **CMap LINCS 2020** (Broad) | 720,216 drug signatures × 12,328 genes, plus each cell line's donor sex |
| **GTEx v8** (NIH) | the genes that differ between men and women, per tissue, 44 tissues |

Full detail, download paths and papers:
[`Phase6_NoSampleFilter/sources/README.md`](Phase6_NoSampleFilter/sources/README.md)

The 33 GB `.gctx` is not in this repo. GitHub caps files at 100 MB.

---

## The phases

| Folder | What it did |
|--------|-------------|
| `Phase1_Vorinostat` | first single-drug test, one drug end to end |
| `Phase2_Fishers_50Drugs` | Fisher's exact test across 50 drugs |
| `Phase3_PerLineage_GSEA` | moved to GSEA, run per lineage |
| `Phase4_SplitLibraries_DirectionFix` | split each tissue's gene list into female and male halves; fixed the direction |
| `Phase5_PerSignature_Categories` | A/B/C/D per signature, **with** a ≥5 M / ≥5 F sample-size filter |
| `Phase6_NoSampleFilter` | Same pipeline with the sample-size filter removed, on Marouen's instruction |
| **`Phase7_SexBalance_Audit`** | **current.** Not a pipeline step — an audit of Phase 6's hits, asking how many male and female donors are actually behind them |
| `cmap`, `gtex` | metadata and source notes |

---

## Where things stand now — Phase 6

The sample-size filter is gone. **The only filter on the data is `|z| > 2`.**

| | |
|---|---|
| signatures across the 50 drugs | **36,523** |
| cell lines · lineages · drug × lineage pairs | 135 · 20 · 909 |
| GTEx tissues | 44 |
| from female donors / male donors | 22,813 (62.5%) / 13,710 (37.5%) |

Categories, counted per signature (best of its 44 tissues):

| A | B | C | D |
|---|---|---|---|
| 2,656 (7.3%) | 5,163 (14.1%) | 4,684 (12.8%) | 24,020 (65.8%) |

17 of the 50 drugs got any hit at all.

![A/B/C/D per drug](Phase6_NoSampleFilter/figures/04_signature_ABCD_2panel.png)

The clearest single result so far — **bortezomib × blood → EBV-transformed
lymphocytes**, the one Category A pair that survives requiring more than one
cell line of each sex (10 M / 5 F, 255 signatures):

| | NES | FDR | |
|---|---|---|---|
| female-biased genes | −1.606 | 0.0000 | **significant** |
| male-biased genes | +0.999 | 0.4874 | not significant |

Numbers, figures and the full write-up:
[`Phase6_NoSampleFilter/README.md`](Phase6_NoSampleFilter/README.md)

---

## What is **not** settled

Read this before quoting any number above as a finding.

**1. Category A is essentially one drug.** Bortezomib supplies 2,593 of the
2,656 Category A signatures — **97.6%**. Only 6 of 50 drugs have any A; 44 have
none. Bortezomib also holds 20% of all signatures, so the drug with the most
tests has the most hits. That is arithmetic before it is biology.

**2. Most A hits have no female cell lines.** Of the 29 pair-level A hits,
20 have zero female lines and 1 has zero male. Requiring more than one cell line
per sex leaves exactly **one** pair standing.

**2b. Category A looks 91.3% female, and that is pseudoreplication.** Against a
62.5% baseline it looks like a sex effect. But Category A draws on 17 female and
18 male cell lines — balanced — and the female lines contributed 143 signatures
each against 13 for the male lines. `bortezomib × breast` alone is 1,875 of the
2,656 Category A signatures with no male cell lines at all. Of the 11 Category A
drug × lineage pairs, one survives a requirement of two cell lines per sex and an
average not dominated by a single line: **bortezomib × blood**, the headline
result. Full numbers in [`Phase7_SexBalance_Audit`](Phase7_SexBalance_Audit/README.md).

**3. The only sex-balance guarantee left is at drug level.** In Phase 5 the
≥5 M / ≥5 F rule ran in two places — once when picking the 50 drugs (grouped by
**drug**) and once inside the pipeline (grouped by **drug × lineage**). Phase 6
removed the second one, on Marouen's instruction, because it was turning 36,124
rows into D. The first is still in force: `top50_names.json` is unchanged.

That leaves a gap. A drug qualifies on its lineages pooled together, but results
are produced per lineage, and a qualifying drug can still have an all-male
lineage — bortezomib × CNS has 2 M and 0 F. This is not a coding bug; the code
does what it was told. It is a reporting gap, and it is why 20 of the 29
Category A hits sit in lineages with no female cell lines.

**4. Brain subregions are counted as independent.** 19 of the 29 A hits come
from two CNS pairs spread across 13 correlated brain subregions.

**5. No project-wide multiple-testing correction.** FDR is computed inside each
GSEA run. Across all 39,996 tests there is no Benjamini–Hochberg pass yet.

**6. No negative control.** DMSO, the solvent, should show nothing. It is not in
the 50, so we cannot prove the hits are drug-specific.

**7. Two gene-set definitions are in use.** The batch A/B/C/D run uses one set
per tissue (44 sets, all sex-biased genes together). The case-study figures split
each tissue into female and male halves (88 sets). These are **different tests**
and their NES and FDR must not be quoted side by side.

**8. Cell lines are not tissue.** The GTEx paper says so itself. Category A
compares a cell line against the sex-biased genes of the organ it came from, and
the two are not the same thing.

---

## Open questions for Marouen

- Does "how many signatures" mean 36,523 unique signatures or 1,607,012
  signature × tissue combinations?
- Keep or drop `|z| > 2`?
- Apply the sex-balance rule at drug × lineage level?
- A minimum number of signatures per pair?
- Count the 13 CNS subregions as 13 hits or as 1?
- Which gene-set definition goes in the paper — mixed or split?

---

## A bug worth writing down

We re-ran the same code in a fresh session and everything came out Category C.

The cause: `sb_library_sets` was rebuilt from JSON as Python **lists** instead of
**sets**. Inside `run_gsea` there is a set intersection, which fails on a list.
The failure was caught by a `try/except`, so nothing printed. Every pair quietly
became C or D.

Fixed by converting back to sets. Verified: bortezomib × blood gives 29
significant tissues again.

**Check the inputs before the run, not after:**

```python
len(sb_library_sets) == 44
type(list(sb_library_sets.values())[0]) is set     # a list fails silently
len(cmap_to_gtex)    == 20
```

---

## Running it

```bash
pip install -r requirements.txt
```

Python **3.12 or newer** (the plotting code uses PEP 701 f-strings; Colab is on
3.13). Download the 33 GB `.gctx` first — see `sources/`. The batch A/B/C/D run
takes about 24 minutes on Colab.
