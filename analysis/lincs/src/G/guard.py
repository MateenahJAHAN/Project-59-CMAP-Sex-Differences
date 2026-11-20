"""Blinding guard for task G.

The entry point calls check() before it loads any response tensor. check() compares the state and label
columns of the line table it was given with those of the frozen real line table (task D's D_line_table.tsv,
sha256 5da2cc4b...), row by row after aligning on cell_iname.

- If the given table's label columns equal the real ones, or fewer than `min_share` of the label-carrying rows
  differ, the labels count as REAL. Then the run is refused unless ALL of these hold:
    1. --mode real;
    2. a frozen plan file is given (--plan) and the environment variable P59_UNBLIND_PLAN_SHA256 equals its
       sha256;
    3. every code and config file listed in src/G/SHA256SUMS.txt has that sha256 (7m), the list itself is
       committed at HEAD, and under 7m A the plan text contains the sha256 of SHA256SUMS.txt;
    4. under 8b A or B, --plan-tag names a git tag whose copy of analysis_plan/G_plan.md has the same sha256 as
       --plan (8b B also needs --osf-id);
    5. no reduced run setting (--B, --calib-n, --replay-n, --relabel-n) is given.
- Otherwise (labels scrambled, for tests) the run is allowed only with --mode test.

This module never sets P59_UNBLIND_PLAN_SHA256.
"""
import os, subprocess, hashlib
import numpy as np, pandas as pd

ENV = "P59_UNBLIND_PLAN_SHA256"
# task D's state and label columns (src/D/09_test_merge.py, STATE_COLS)
STATE_COLS = ["taskA_state", "taskA_confidence", "taskA_y_score", "taskA_y_score_source", "taskA_Yhi", "taskA_Xhi",
              "taskA_XIST_ccle", "taskA_XIST_sanger", "taskA_XIST_cmpbroad", "taskC_state", "taskC_match_key", "source_conflict",
              "cmap_donor_sex", "taskA_cmap_donor_sex", "elig_factorial", "fact_Y", "fact_X", "elig_strict", "strict_XY",
              "elig_annotated", "annot_M"]


class GuardRefusal(SystemExit):
    def __init__(self, msg):
        super().__init__(f"BLINDING GUARD: refused. {msg}")
        self.msg = msg


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _norm(v):
    if isinstance(v, float) and np.isnan(v):
        return "<NA>"
    if v is None or (isinstance(v, str) and v == ""):
        return "<NA>"
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (int, float, np.integer, np.floating)):
        return repr(float(v))
    return str(v)


def classify(given, real, min_share):
    """Returns dict(status, n_label_rows, n_rows_changed, share_changed, columns_compared)."""
    assert given.cell_iname.is_unique and real.cell_iname.is_unique
    assert set(given.cell_iname) == set(real.cell_iname), "line tables hold different lines"
    missing = [c for c in STATE_COLS if c not in given.columns]
    assert not missing, f"given line table lacks label columns {missing}"
    g = given.set_index("cell_iname").loc[real.cell_iname, STATE_COLS]
    r = real.set_index("cell_iname")[STATE_COLS]
    gt = g.apply(lambda row: tuple(_norm(v) for v in row), axis=1).values
    rt = r.apply(lambda row: tuple(_norm(v) for v in row), axis=1).values
    carries = np.array([any(v != "<NA>" for v in t) for t in rt])
    changed = np.array([a != b for a, b in zip(gt, rt)])
    n_lab = int(carries.sum())
    n_ch = int((changed & carries).sum())
    share = n_ch / n_lab if n_lab else 0.0
    if n_ch == 0:
        status = "real"
    elif share < min_share:
        status = "near_real"
    else:
        status = "scrambled"
    return dict(status=status, n_label_rows=n_lab, n_rows_changed=n_ch, share_changed=share,
                columns_compared=len(STATE_COLS))


def _git(root, *args):
    return subprocess.run(["git", "-C", root, *args], capture_output=True, text=True)


def code_sums_ok(root, sums_path):
    """All files in SHA256SUMS.txt exist with that sha256; returns (ok, problems)."""
    probs = []
    if not os.path.exists(sums_path):
        return False, ["SHA256SUMS.txt missing"]
    for line in open(sums_path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        h, rel = line.split(None, 1)
        p = os.path.join(root, rel.strip())
        if not os.path.exists(p):
            probs.append(f"missing {rel}")
        elif sha256_file(p) != h:
            probs.append(f"changed {rel}")
    rel_sums = os.path.relpath(sums_path, root)
    r = _git(root, "show", f"HEAD:{rel_sums}")
    if r.returncode != 0:
        probs.append("SHA256SUMS.txt is not committed at HEAD")
    elif r.stdout != open(sums_path).read():
        probs.append("SHA256SUMS.txt differs from its committed version")
    return not probs, probs


def check(status_info, mode, cfg, root, plan_path=None, plan_tag=None, osf_id=None, reduced_settings=(), env=None):
    """Raise GuardRefusal unless the run is allowed. Returns a dict for the run manifest."""
    env = os.environ if env is None else env
    status = status_info["status"]
    out = dict(label_status=status, mode=mode)
    if status in ("real", "near_real"):
        if mode != "real":
            raise GuardRefusal(f"the line table carries real labels ({status_info['n_rows_changed']} of "
                               f"{status_info['n_label_rows']} label rows differ from the frozen table); use --mode real.")
        if not plan_path or not os.path.exists(plan_path):
            raise GuardRefusal("real labels need --plan <frozen plan file>.")
        h = sha256_file(plan_path)
        if env.get(ENV) != h:
            raise GuardRefusal(f"environment variable {ENV} does not equal the sha256 of {plan_path}.")
        if reduced_settings:
            raise GuardRefusal(f"real labels need full settings; reduced settings given: {list(reduced_settings)}.")
        sums = os.path.join(root, cfg.raw["paths"]["code_sums"])
        ok, probs = code_sums_ok(root, sums)
        if not ok:
            raise GuardRefusal(f"code differs from src/G/SHA256SUMS.txt: {probs[:5]}")
        opt7m = cfg.decide("7m")
        if opt7m == "A" and sha256_file(sums) not in open(plan_path).read():
            raise GuardRefusal("7m A: the frozen plan must contain the sha256 of src/G/SHA256SUMS.txt.")
        opt8b = cfg.decide("8b")
        if not plan_tag:
            raise GuardRefusal("8b: --plan-tag (the dated tag of the frozen plan) is required.")
        r = subprocess.run(["git", "-C", root, "show", f"{plan_tag}:analysis_plan/G_plan.md"], capture_output=True)
        if r.returncode != 0 or hashlib.sha256(r.stdout).hexdigest() != h:
            raise GuardRefusal(f"8b: tag {plan_tag} does not hold analysis_plan/G_plan.md equal to --plan.")
        if opt8b == "B" and not osf_id:
            raise GuardRefusal("8b B: --osf-id is required.")
        out.update(plan_sha256=h, plan_tag=plan_tag, osf_id=osf_id or "", rule_8a=cfg.decide("8a"), rule_7m=opt7m)
    else:
        if mode != "test":
            raise GuardRefusal("the line table's labels differ from the frozen real table; such a run is a test (--mode test).")
        out.update(plan_sha256="", plan_tag="", osf_id="", rule_8a=cfg.decide("8a"), rule_7m=cfg.decide("7m"),
                   rule_8b=cfg.decide("8b"))
    return out
