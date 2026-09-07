# Project 59 - Sex Differences in Drug Response

**Fatima Fellowship 2026**
Mentor: Dr. Marouen Ben Guebila, Viswanathan Lab, Dana-Farber Cancer Institute / Harvard
Fellow: Mateenah Jahan

---

## The question

Do cells from men and cells from women respond differently to the same drug?

Cell lines grow in a dish. They make no hormones. So if we see a difference,
it comes from the genes, not from hormones.

---

## Where the data comes from

### CMap LINCS 2020 - Broad Institute

The Broad Institute took human cancer cell lines, added drugs, waited 6, 12 or
24 hours, then measured how about 12,328 genes reacted. They did this for
thousands of drugs across thousands of cell lines.

Each single experiment is called a **signature**.

| File | What it holds |
|------|---------------|
| level5_beta_trt_cp_n720216x12328.gctx | The big one. 33 GB. 12,328 genes by 720,216 experiments. Every number says how much a gene moved. |
| siginfo_beta.txt | One row per experiment: which drug, which cell line, what dose, how long |
| cellinfo_beta.txt | One row per cell line: donor sex, and which tissue it came from |
| geneinfo_beta.txt | Gene IDs and gene names |

### GTEx v8 - NIH

GTEx studied real human tissues and listed, for each tissue, the genes that
behave differently in men and women. This is our reference list.

| File | What it holds |
|------|---------------|
| signif.sbgenes.txt | Sex-biased genes for each tissue |

We keep only tissues with at least 15 such genes. That leaves 44 tissues.

---

## Three words to keep straight

**Cell line** - cells from one person. THP1 came from one male patient with
leukemia. One person means one sex.

**Lineage** - the organ that cell line came from. THP1 came from blood, so its
lineage is blood. Many cell lines share a lineage.

**Signature** - one experiment. One cell line, one drug, one dose, one time.
A single cell line can give many signatures because the dose and the time keep
changing.

Example: bortezomib in the blood lineage covers 15 cell lines and 255 signatures.

---

## What we do

1. Pick 50 drugs that have at least 5 male and 5 female cell lines
2. For each drug and lineage, average all its signatures into one profile, and
   keep only the genes that moved a lot
3. Run GSEA to check that gene list against all 44 GTEx tissues
4. Give every combination a category

### The four categories

| Category | What it means |
|----------|---------------|
| **A** | Significant, and the tissue is the cell line's own tissue |
| **B** | Significant, but a different tissue |
| **C** | We tested it, nothing came out |
| **D** | We could not test it at all |

If NES is positive, the signal leans female. If negative, it leans male.

---

## What we found

We have 39,996 rows in total.

| Category | Count |
|----------|-------|
| A | 1 |
| B | 79 |
| C | 540 |
| D | 39,376 |

Only 19 of 909 drug-lineage pairs could be tested. That is 2 percent.

Only 5 drugs out of 50 gave any hit: bortezomib, emetine, YM-155, mitoxantrone,
and NVP-BEZ235. All five change a large number of genes.

The one Cat A hit is **bortezomib in blood cell lines**, matching the sex-biased
genes of EBV-transformed lymphocytes. NES -1.78, FDR 0.013.

### Why so many D

| Reason | Rows |
|--------|------|
| That lineage had fewer than 5 male or 5 female cell lines | 36,124 |
| GSEA found too few genes to work with | 3,036 |
| That tissue's gene list was smaller than 15 | 216 |

---

## Things that still need checking

**1. The male signal might not be real.**
76 of the 80 hits lean male. But the data itself leans male. Bortezomib in
blood has 190 male signatures and only 65 female. Two male cell lines, THP1 and
JURKAT, make up 47 percent on their own. We take a plain average, so the result
is male partly because the input is male. We cannot yet tell a real drug effect
apart from this.

**2. There is no negative control.**
DMSO is the liquid the drug is dissolved in. It should show nothing. It was
dropped from the 50, so we cannot prove the hits are drug-specific.

**3. The Cat A hit is weaker than the Cat B hits.**
The home tissue gives NES -1.78. Brain tissues give -2.3 to -2.9. If the effect
were really about the tissue, the home tissue should be the strongest one.

**4. Brain tissues may win just because they are big.**
10 of the 30 Cat B hits are brain. GTEx brain gene lists are the longest ones.
Longer lists are easier to hit by chance.

**5. Some comparisons can never work.**
7 of 20 lineages have no matching GTEx tissue, so they can never get Cat A.
Seven tissue names in our mapping do not exist in GTEx at all - Ovary, Prostate,
Uterus, Cervix, Bladder, Kidney_Medulla. These organs exist in only one sex, so
there is no male-female comparison to make. On the CMap side it is the same
problem: breast cell lines are all female, prostate cell lines are all male.

---

## A bug worth writing down

We re-ran the same code in a fresh session and everything came out Cat C.

The reason: sb_library_sets was rebuilt using Python lists instead of sets.
Inside run_gsea there is a line that does an intersection with a set, which
fails on a list. The error was caught by a try/except, so nothing was printed.
Every pair quietly failed and became C or D.

Fixed by turning the lists back into sets. Checked: bortezomib in blood now
gives 29 significant tissues again, same as before.

The lesson is to check the inputs before running the pipeline, not after.

---

## Folders

- Phase1_Vorinostat - first single-drug test
- Phase2_Fishers_50Drugs - Fisher test across 50 drugs
- Phase3_PerLineage_GSEA - GSEA per lineage
- Phase4_SplitLibraries - split gene libraries by direction
- Phase5_PerSignature_Categories - current work
- cmap - CMap metadata
- gtex - GTEx files

The 33 GB gctx file is not stored here. Download it from the CMap S3 bucket.

---

## Before you run anything

Check two things:

- sb_library_sets should have 44 tissues, and each one should be a **set**, not a list
- cmap_to_gtex should have 20 lineages
