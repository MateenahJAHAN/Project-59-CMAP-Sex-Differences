# Part A raw-data replay

Independent extraction of the EmeraldBay expression counts recovered **100,253 endpoint-vehicle cells** from 52 lines and **30,245 time-zero cells**. Every checked field in all 130,498 per-cell rows matched the archived table. All 207 line-by-sample rows were byte-identical. No line changed its state, Y-arm class or confidence. The largest numerical difference in the line table was 3.47 × 10⁻¹⁸.

All three sci-Plex matrices were read in full, covering **1,171,291,549 sparse matrix entries**. The table of 23,570 vehicle cells matched the archived compressed file byte-for-byte. Five of six line-by-experiment marker patterns matched the freshly derived bulk RNA pattern. A549 expressed both Y markers and XIST in sciPlex3 and only Y markers in bulk.

The DNA check re-parsed 1,102 Cellosaurus records. All 19 STR output tables matched their archived checksums. Joining the newly calculated RNA scores directly to the re-parsed STR calls confirmed 47/52 in EmeraldBay and 3/3 distinct sci-Plex lines. New bulk RNA measurements retain all 773 CCLE and 1,082 Cell Model Passports line calls. Agreement with DNA typing was 714/734 for CCLE and 1,007/1,031 for Cell Model Passports. [Calls by accession](tables/str_fresh_bulk.tsv)

These comparisons involve the AMELY locus in reference stocks. The historical Tahoe RNA classification gives 48/50 in the same check.

The pinned Tahoe-100M observation file passed its SHA-256 check. All 100,648,790 rows matched the local metadata projection across the eleven source fields used by the analysis. There are 2,205,786 full-filter endpoint-vehicle cells. Across 1,910 expression shards, 4,212 eligible row groups contained 2,205,786 qualifying cells. All 7,640 per-shard scientific files matched the archived files byte-for-byte. The [source receipt](receipts/tahoe_raw_shards.json) records the read counts and diagnostic differences. [Run status](receipts/run_status.json) tracks the two downstream Tahoe callers and the EmeraldBay all-condition identity analysis.

## Files and commands

The tables directory contains the complete line-level and EmeraldBay sample-level results. The receipts include the 116 EmeraldBay shard records, all 16 GEO file URLs and hashes, STR record identifiers and checksums, and aggregate and per-cell comparisons.

Run `python3 scripts/verify.py` to verify the package and recompute the reported counts. The original sci-Plex sparse-matrix parser and scientific constants are in `scripts/`. Place the 16 source files listed in `receipts/source_verification.json` under `RUN/data/F1/sciplex/`. From this directory run `python3 scripts/process_sciplex.py RUN sciPlex2` and repeat with `sciPlex3` and `sciPlex4`. This recreates the per-cell vehicle records, sample labels and identity accumulators for each full matrix. The parser requires NumPy, pandas, SciPy and PyArrow.
