# Phase 7 — Sex-balance audit

This is not a new step in the pipeline. It is an **audit of Phase 6's output**,
asking one question nobody had asked yet:

> Of the signatures that came out significant, how many male and how many female
> donors are actually behind them?

Nothing is filtered, nothing is recomputed. It reads the two tables Phase 6
already produced and counts.

**It was not requested.** It came out of looking at the Phase 6 result and
noticing that Category A hits kept landing in lineages with no female cell
lines. The numbers below are the reason.

---

## What it found

### 1. Category A looks strongly female — against a 62.5% baseline

| Category | ♀ signatures | ♂ signatures | total | % female |
|----------|---|---|---|---|
| **A** — home tissue, significant | 2,426 | 230 | 2,656 | **91.3%** |
| B — cross tissue, significant | 2,958 | 2,205 | 5,163 | 57.3% |
| **C** — tested, not significant | 1,710 | 2,974 | 4,684 | **36.5%** |
| D — could not be tested | 15,719 | 8,301 | 24,020 | 65.4% |
| **whole dataset** | 22,813 | 13,710 | 36,523 | **62.5%** |

A and C are the interesting pair, because both were *tested*. A got a hit and is
91% female. C got nothing and is 64% male.

That looks like a finding. It is not.

### 2. Counting donors instead of experiments kills it

| Category | ♀ lines | ♂ lines | signatures per ♀ line | per ♂ line | ratio |
|----------|---|---|---|---|---|
| **A** | **17** | **18** | **143** | **13** | **11.0×** |
| B | 59 | 52 | 50 | 42 | 1.2× |
| C | 60 | 53 | 29 | 56 | 0.5× |
| D | 76 | 57 | 207 | 146 | 1.4× |

Category A draws on 17 female and 18 male cell lines — **balanced**. The female
lines simply contributed eleven times more experiments each, at different doses
and timepoints. No other category shows that pattern.

This is **pseudoreplication**: 143 signatures from one cell line are not 143
donors, they are one donor measured 143 times.

### 3. Most of it is one drug in one lineage

| Drug × lineage | ♀ sigs | ♂ sigs | ♀ lines | ♂ lines |
|---|---|---|---|---|
| bortezomib × breast | 1,875 | 0 | 7 | **0** |
| bortezomib × kidney | 450 | 5 | 1 | 1 |
| bortezomib × blood | 65 | 190 | 5 | 10 |
| YM-155 × skin | 25 | 3 | 3 | 3 |
| *(7 more, all ≤ 9 signatures)* | | | | |

`bortezomib × breast` alone is **1,875 of 2,656** Category A signatures — 70.6%
— with **no male cell lines at all**. Breast cancer cell lines come from female
donors, so that comparison can never carry a sex contrast.

Of all 2,656 Category A signatures, **1,895 (71.3%) come from single-sex
lineages**.

### 4. Effective sample size: 11 pairs, 1 survives

Kish effective sample size, (Σn)² ⁄ Σn², over signatures per cell line:

| Pair | lines | signatures | biggest line | effective | verdict |
|---|---|---|---|---|---|
| bortezomib × breast | 7 | 1,875 | 43.4% | 2.89 | no male lines |
| bortezomib × kidney | 2 | 455 | **98.9%** | **1.02** | behaves like one donor |
| **bortezomib × blood** | **15** | **255** | **23.9%** | **6.60** | **passes** |
| YM-155 × skin | 6 | 28 | 75.0% | 1.73 | one line dominates |

Requiring **≥ 2 cell lines of each sex** *and* **an average not dominated by one
line**, ten of the eleven pairs fail.

The one that survives — `bortezomib × haematopoietic_and_lymphoid_tissue →
Cells_EBV.transformed_lymphocytes` — is the same pair Phase 6 already reported as
the headline result (female NES −1.606, FDR 0.0000; male not significant).

Across all 11 pairs: **48 real cell lines, 29.5 effective.**

---

## Why this matters

Three things follow, and none of them require changing the method.

1. **The 91.3% female figure must not be reported as a sex effect.** It is a
   sampling pattern. Writing it up as biology would not survive review.
2. **The headline result got stronger, not weaker.** It was put through the
   hardest test available and came out the other side alone.
3. **The ten failures are a finding too.** They say CMap's sampling cannot
   support a broad claim about sex differences across many drugs.

### What this does *not* do

It does not fix pseudoreplication — it only measures it. The fix is to average
within each cell line first, then across lines, so every donor counts once. That
needs the full pipeline rerun against the 33 GB matrix, which is a **method
change** and is not made without the advisor's agreement.

---

## Files

```
Phase7_SexBalance_Audit/
├── README.md
├── notebook/
│   └── phase7_sex_balance_audit.ipynb
└── results/
    ├── signatures_by_category_and_sex.csv
    ├── cell_lines_by_category_and_sex.csv
    ├── category_A_by_drug_lineage_sex.csv
    └── pseudoreplication_category_A.csv
```

| File | Holds |
|------|-------|
| `signatures_by_category_and_sex.csv` | the first table, plus the dataset baseline row |
| `cell_lines_by_category_and_sex.csv` | distinct donors per category, and signatures per line |
| `category_A_by_drug_lineage_sex.csv` | all 11 Category A pairs, with a `both_sexes_2plus` flag |
| `pseudoreplication_category_A.csv` | effective sample size per pair |

**Note on the notebook:** it is committed as it was run, exploratory cells and
all. Cells 1 and 4 load the same file, and cells 2 and 5 produce overlapping
tables — the later ones are the cleaned versions. Cell 5 and cell 6 are the two
that produced everything above.

---

## Inputs

Both come from Phase 6 and live in Google Drive, not in this repository:

| File | Rows | What it is |
|------|------|-----------|
| `per_signature_v4_ABCD.csv` | 39,996 | one row per drug × lineage × tissue test |
| `signature_level_ABCD.csv` | 36,523 | one row per signature, with donor sex and its category |

No `.gctx` is needed. The whole audit runs in seconds.
