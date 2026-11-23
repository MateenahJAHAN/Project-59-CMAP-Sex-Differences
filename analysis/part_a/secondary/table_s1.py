"""Write Supplementary Table S1 from the full-cohort state calls and cross-caller comparisons."""
import csv
import os
import sys
import tempfile
from pathlib import Path

from lincs_state_rule import CELLO_RULE_LABEL, RULE_LABEL_NAME, index, task_a_recall

SOURCES = ("Tahoe-100M", "EmeraldBay", "sciPlex2", "sciPlex3", "sciPlex4", "LINCS2020", "GDSC2", "PRISM")
FIELDS = ("cellosaurus_id source cell_name called_by line_id cmp_model_id organ sex_label state confidence "
          "arm_state n_vehicle_cells n_vehicle_cells_called n_units units_agree grid_stable xist y_score "
          "yp_score yq_score mosaic_marker mosaic_verdict negative_share_lo negative_share_hi").split()
SCORES = ("xist", "y_score", "yp_score", "yq_score")


def read(path):
    with path.open(newline="") as file:
        return list(csv.DictReader(file, delimiter="\t"))


def tahoe(res):
    local = index(read(res / "tahoe/latest/partA_local_states.tsv"), "cellosaurus_id")
    compare = index(read(res / "consensus/k9_line_comparison.tsv"), "cellosaurus_id")
    mosaic = index(read(res / "consensus/mosaic_resolution.tsv"), "cell_name", "marker")
    pairs = {(r["line"], r["marker"]) for r in read(res / "depth_matching/A3_pairs.tsv")}
    mosaic_by_name = {line: mosaic[line, marker] for line, marker in pairs}
    if len(mosaic_by_name) != len(pairs):
        raise ValueError("Multiple mosaic markers for one Table S1 line")
    for row in read(res / "consensus/tahoe_partA_states_consensus.tsv"):
        accession = row["cellosaurus_id"]
        line, detail = local[accession], compare[accession]
        result = dict(cellosaurus_id=accession, cell_name=row["cell_name"], source="Tahoe-100M",
                      called_by="local and F2", organ=line["organ"], sex_label=row["cellosaurus_sex"],
                      state=row["state_final"], confidence=row["confidence_final"], arm_state=row["arm_state_final"],
                      n_vehicle_cells=detail["n_cells_before_filter_f2"], n_vehicle_cells_called=row["n_vehicle_cells_called"],
                      n_units=line["n_plates"], units_agree=detail["plates_agree_f2"], grid_stable=detail["grid_stable_f2"],
                      **{col: row[col] for col in SCORES})
        if row["cell_name"] in mosaic_by_name:
            mr = mosaic_by_name[row["cell_name"]]
            result.update(mosaic_marker=mr["marker"],
                          mosaic_verdict="mosaic not established" if mr["verdict_final"].startswith("not established") else mr["verdict_final"],
                          negative_share_lo=mr["share_lo_both"], negative_share_hi=mr["share_hi_both"])
        yield result


def external_atlases(res, labels):
    local = index(read(res / "tahoe/latest/partA_local_states.tsv"), "cellosaurus_id")
    for dataset, file in (("EmeraldBay", "F1_emerald_states.tsv"), ("sci-Plex", "F1_sciplex_states.tsv")):
        for row in read(res / "F1" / file):
            accession = row["cellosaurus_id"]
            sex = row["annotation"]
            if dataset == "EmeraldBay" and row["annotation_source"] == "pilot_pie_sex":
                sex = labels[accession]["cellosaurus_sex"] if accession in labels else local[accession]["annotation"]
            yield dict(cellosaurus_id=accession, source=dataset if dataset == "EmeraldBay" else row["sample"],
                       cell_name=row["cell_name"], called_by="F1", organ=row["organ"], sex_label=sex,
                       n_units=row["n_plates"], units_agree=row["plates_agree"], xist=row["XIST"],
                       **{field: row[field] for field in ("state", "confidence", "arm_state", "grid_stable", "n_vehicle_cells", *SCORES[1:])})


def bulk(res, labels):
    data = {"a_states": read(res / "A/A_cellline_sex_states.tsv"),
            "a_grid": read(res / "A/A_state_sensitivity_per_line.tsv"), "a_cutoffs": read(res / "A/A_cutoffs.tsv")}
    calls = task_a_recall(data)
    lincs = {}
    for a in data["a_states"]:
        _, state, confidence, _ = calls[a["cell_iname"]]["cellosaurus"]
        if state == "no_data":
            continue
        ac = a["cellosaurus_ac"]
        record = dict(cellosaurus_id=ac, source="LINCS2020", cell_name=a["cellosaurus_name"],
                      called_by="task A", cmp_model_id=a["sanger_model_id"],
                      sex_label=RULE_LABEL_NAME[CELLO_RULE_LABEL[a["cellosaurus_sex"]]], state=state, confidence=confidence)
        if ac in lincs and lincs[ac] != record:
            raise ValueError(f"Conflicting LINCS records for {ac}")
        lincs[ac] = record
    yield from lincs.values()
    for ac, c in labels.items():
        for screen, flag, column in (("GDSC2", "in_gdsc2", "model_id"), ("PRISM", "in_prism", "BROAD_ID")):
            if c[flag] == "True":
                yield dict(cellosaurus_id=ac, cell_name=c["model_name"], source=screen, called_by="task C",
                           line_id=c[column], cmp_model_id=c["model_id"], sex_label=c["annotated_sex"],
                           state=c["expressed_state"])


def main(work):
    res = work / "results"
    labels = index(read(res / "C/C_cellline_sex_labels.tsv"), "RRID")
    rows = [*tahoe(res), *external_atlases(res, labels), *bulk(res, labels)]
    rows.sort(key=lambda r: (r["cellosaurus_id"], SOURCES.index(r["source"])))
    unique = index(rows, "cellosaurus_id", "source")
    accessions = {r["cellosaurus_id"] for r in rows}
    if len(unique) != 1694 or len(accessions) != 1123:
        raise ValueError(f"Unexpected membership: {len(unique)} accession/source pairs, {len(accessions)} lines")
    dest = res / "supplement/supplementary_table_S1.tsv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", dir=dest.parent, newline="", delete=False) as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(out.name, dest)
    print(f"{dest}: {len(rows)} accession/source rows, {len(accessions)} distinct cell lines, {len(FIELDS)} columns")


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
