# Reproducing Project 59 from source data

The reference is the [released manuscript and supplement](../../manuscript/). Results below distinguish the archived analysis from the new runs.

| Analysis | Result from the source-data replay |
| --- | --- |
| LINCS2020 and GTEx | All 450,255,544 expression values across 36,523 signatures matched the provider's Level 5 matrix exactly. The GTEx v8 tissue file also matched. [LINCS](lincs_inputs/comparison.json) · [GTEx](lincs_inputs/gtex_comparison.json) |
| Phase 2 | All 50 original compound classifications reproduced. For AS-605240 in liver, restricting the correction to tissues with positive overlap gives an adjusted p of 0.0247. Including all 44 eligible tissues gives 0.0544. [All 2,200 tissue comparisons](phases2_5/comparison/phase2_all_tissue_tests.csv) |
| Phase 3 | Fostamatinib gained an own-tissue hit in the literal notebook replay (FDR 0.0410). Changing only the order of tied scores produced FDRs of 0.0336, 0.0757 and 0.0746. The compound’s classification depends on rank order. The full panel is under comparison. |
| Phase 4 | The executed notebook ranked all mapped genes. The original README incorrectly described the filter used in Phase 5. The full comparison is in progress. |
| Phase 5 | All 39,996 categories reproduced, including one own-tissue hit and 79 other-tissue hits among 620 completed tests. GSEApy 1.3.1 reproduced every saved score. Version 1.2.1 changed 14 FDR values by at most 0.001. The sign of the pooled enrichment score cannot identify which donor group responded more. [Figure](phases2_5/figures/phase5_testable_tissues.png) · [Results](phases2_5/sensitivity_gseapy_1.3.1/comparison/phase5.json) · [Version comparison](phases2_5/comparison/phase5_software_versions.csv) |
| Bulk RNA | All 1,103 line classifications matched after rebuilding the two source expression tables. Two CCLE Y scores changed by at most 6.94 × 10⁻¹⁸. [Seven output comparisons](bulk_RNA/comparison.json) · [Source files](bulk_RNA/sources.json) |
| Part A vehicle cells | All 130,498 EmeraldBay vehicle and time-zero records and all 23,570 sci-Plex vehicle records reproduced. All 52 EmeraldBay line calls matched. The independent Tahoe metadata count is 2,205,786 full-filter vehicle cells. [Cell and atlas comparisons](part_a/) |
| Part B response data | All 45 rebuilt Level 3 and Level 5 outcome files matched byte for byte. All 30 model-design checks and three 10,000-permutation ploidy controls reproduced. [Result tables](part_b/) |
| Manuscript figures | Re-rendering the seven figures using the archived plotting inputs gave byte-identical PDFs. The new source-data figure rebuild is in progress. [PDF and pixel comparisons](figures/baseline_comparison.tsv) |

The Part A Tahoe caller comparison, the full Part B analysis, and the fresh 1,000-replicate power simulation are in progress. The manuscript and supplement linked above are the existing release. Proposed edits clarify the historical comparisons in [Note S19](manuscript_review/note_s19_proposed.patch) and correct a typographical error in [Note S13](manuscript_review/note_s13_typo.patch).

[Per-compound rerun files and restart instructions](phases2_5/) · [Part B rebuild and completed controls](part_b/)

[Post hoc power study, recount and restart commands](power/)
