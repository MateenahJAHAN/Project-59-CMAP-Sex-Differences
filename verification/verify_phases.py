"""Cross-check the fellowship's archived CSV exports at their original units."""

import csv
from collections import Counter, defaultdict
from itertools import combinations
from math import isclose
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[1]


def rows(path):
    with (ROOT / path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def check(condition, message):
    if not condition:
        raise ValueError(message)


def early_screens():
    pilot = "Phase1_Vorinostat"
    signatures = rows(f"{pilot}/vorinostat_signatures.csv")
    check(len({r["sig_id"] for r in signatures}) == len(signatures) == 2641, "Phase 1 signature IDs")
    check(Counter(r["donor_sex"] for r in signatures) == {"M": 1180, "F": 1461}, "Phase 1 donor counts")
    check(len({r["cell_iname"] for r in signatures}) == 116, "Phase 1 cell lines")
    genes = sorted(rows(f"{pilot}/vorinostat_all_genes_results.csv"), key=lambda r: float(r["p_value"]))
    check(len({r["gene_id"] for r in genes}) == len(genes) == 12328, "Phase 1 feature IDs")
    running_q = 1.0
    for rank, gene in reversed(list(enumerate(genes, 1))):
        running_q = min(running_q, float(gene["p_value"]) * len(genes) / rank)
        check(isclose(running_q, float(gene["p_adj"]), rel_tol=1e-10), f"Phase 1 adjusted p for {gene['gene_id']}")
    significant = [r for r in genes if float(r["p_adj"]) < 0.05]
    top = {r["gene_id"] for r in significant if abs(float(r["difference"])) > 1}
    check(len(significant) == 6877 and len(top) == 20 and
          top == {r["gene_id"] for r in rows(f"{pilot}/vorinostat_top20_sex_genes.csv")}, "Phase 1 thresholds")

    phase2 = rows("Phase2_Fishers_50Drugs/top50_drug_classification.csv")
    phase3 = rows("Phase3_PerLineage_GSEA/phase3_classification_final.csv")
    drugs = {r["drug"] for r in phase2}
    check(len(drugs) == 50 and {r["drug"] for r in phase3} == drugs, "Phases 2 and 3 drug panel")
    check(Counter(r["category"] for r in phase2) == {"A": 1, "C": 49}, "Phase 2 classifications")
    check(Counter(r["category"] for r in phase3) == {"A": 3, "B": 23, "C": 24}, "Phase 3 classifications")
    phase4 = "Phase4_SplitLibraries_DirectionFix/phase4_classifications_final"
    primary = rows(f"{phase4}.csv")
    other = rows(f"{phase4} (1).csv")
    check(Counter(tuple(r.items()) for r in primary) == Counter(tuple(r.items()) for r in other),
          "Phase 4 exports differ")
    check(Counter(r["category"] for r in primary) == {"A": 237, "B": 8942, "C": 81}, "Phase 4 classifications")


def coverage_and_sampling():
    phase5 = "Phase5_PerSignature_Categories/results"
    comparisons = rows(f"{phase5}/per_signature_categories.csv")
    keys = {(r["drug"], r["lineage"], r["tissue"]) for r in comparisons}
    check(len(keys) == len(comparisons) == 39996, "Phase 5 comparison IDs")
    by_pair = defaultdict(Counter)
    for row in comparisons:
        by_pair[row["drug"], row["lineage"]][row["category"]] += 1
    check(Counter(row["category"] for row in comparisons) == {"A": 1, "B": 79, "C": 540, "D": 39376},
          "Phase 5 coverage")
    pairs = rows(f"{phase5}/per_drug_lineage.csv")
    check(len(by_pair) == len(pairs) == 909 and sum(c["D"] < 44 for c in by_pair.values()) == 19,
          "Phase 5 tested pairs")
    for pair in pairs:
        counts = by_pair[pair["drug"], pair["lineage"]]
        check(all(counts[c] == int(pair[f"n_{c}"]) for c in "ABCD"), "Phase 5 pair totals")

    phase6 = "Phase6_NoSampleFilter/results"
    signatures = rows(f"{phase6}/signature_ABCD_per_drug.csv")
    totals = {r["category"]: int(r["count"]) for r in rows(f"{phase6}/category_totals.csv")
              if r["level"] == "signature_best_of_44"}
    check(all(sum(int(r[c]) for c in "ABCD") == int(r["total"]) for r in signatures), "Phase 6 per-drug totals")
    check({c: sum(int(r[c]) for r in signatures) for c in "ABCD"} == totals and sum(totals.values()) == 36523,
          "Phase 6 signature totals")
    phase7 = "Phase7_SexBalance_Audit/results"
    by_sex = rows(f"{phase7}/signatures_by_category_and_sex.csv")
    for row in by_sex:
        denominator = sum(totals.values()) if row["category"] == "ALL" else totals[row["category"]]
        check(int(row["signatures_female"]) + int(row["signatures_male"]) == int(row["total"]) == denominator,
              "Phases 6 and 7 signature totals")
    pairs_a = rows(f"{phase7}/category_A_by_drug_lineage_sex.csv")
    a = next(r for r in by_sex if r["category"] == "A")
    check(len(pairs_a) == 11 and all(sum(int(r[f"sigs_{sex}"]) for r in pairs_a) == int(a[f"signatures_{sex}"])
                                    for sex in ("female", "male")), "Phase 7 own-tissue sampling")


def adverse_event_reporting():
    drugs = rows("Phase8_SideEffects_FAERS/per_drug_events.csv")
    phase6 = {r["drug"]: r for r in rows("Phase6_NoSampleFilter/results/signature_ABCD_per_drug.csv")}
    check(len({r["drug"] for r in drugs}) == len(drugs) == 18 and
          sum(int(r["reaction_terms_tested"]) for r in drugs) == 11049, "Phase 8 coverage")
    positive, other = [], []
    volumes, rates = [], []
    for drug in drugs:
        p6 = phase6[drug["drug"]]
        group = next((label for label in "AB" if int(p6[label]) > 0), "none")
        check(group == drug["phase6_best_enrichment"], f"Phase 8 group for {drug['drug']}")
        rate = 100 * int(drug["reaction_terms_flagged"]) / int(drug["reaction_terms_tested"])
        rates.append(rate)
        volumes.append(int(drug["female_and_male_reports"]))
        (other if group == "none" else positive).append(rate)
    check(len(positive) == 7 and len(other) == 11 and len(set(rates)) == len(rates), "Phase 8 rank-test groups")
    check(round(median(positive), 1) == 11.3 and round(median(other), 1) == 14.9, "Phase 8 median percentages")
    # With no ties each assignment of seven ranks is equally likely under the null.
    group_size = len(positive)
    u = sum(a > b for a in positive for b in other)
    rank_sums = [sum(group) for group in combinations(range(1, len(rates) + 1), group_size)]
    shifted_u = u + group_size * (group_size + 1) // 2
    p = 2 * sum(s <= shifted_u for s in rank_sums) / len(rank_sums)
    check(len(set(volumes)) == len(volumes), "Phase 8 ties in report volumes")
    volume_ranks = {v: rank for rank, v in enumerate(sorted(volumes), 1)}
    rate_ranks = {v: rank for rank, v in enumerate(sorted(rates), 1)}
    squared_distances = sum((volume_ranks[v] - rate_ranks[r]) ** 2 for v, r in zip(volumes, rates))
    n = len(drugs)
    rho = 1 - 6 * squared_distances / (n * (n ** 2 - 1))
    check(u == 27 and round(p, 4) == 0.3283 and round(rho, 3) == 0.853, "Phase 8 rank tests")


def verify_fellowship_phases():
    early_screens()
    coverage_and_sampling()
    adverse_event_reporting()
    print("Phases 1–8: archived gene tests, classifications, signature counts and FAERS drug-level comparisons verified.")
