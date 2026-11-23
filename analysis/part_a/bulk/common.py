"""Shared paths and helpers for task C."""
import datetime as dt
import os
import time
import zlib

ROOT = os.environ["P59_WORK"]
DATA = os.path.join(ROOT, "sources")
RES = os.path.join(ROOT, "results", "C")
FIG = os.path.join(ROOT, "figures", "C")
LOG = os.path.join(ROOT, "logs", "C_runlog.txt")

Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y", "NLGN4Y"]

SEX_SPECIFIC_TISSUES = ["Breast", "Prostate", "Ovary", "Endometrium", "Uterus", "Cervix", "Testis",
                        "Vulva", "Placenta"]
SEX_SPECIFIC_CANCER_TYPES = ["Breast Carcinoma", "Prostate Carcinoma", "Ovarian Carcinoma",
                             "Endometrial Carcinoma", "Cervical Carcinoma"]
DROP_CANCER_TYPES = ["Unknown", "Non-Cancerous"]

BASE_SEED = 20261114


def seed_for(*parts):
    return BASE_SEED + zlib.crc32("|".join(str(p) for p in parts).encode())


class Step:
    """Context manager that logs the runtime of a step to logs/C_runlog.txt."""

    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.t0 = time.time()
        return self

    def __exit__(self, *exc):
        el = time.time() - self.t0
        os.makedirs(os.path.dirname(LOG), exist_ok=True)
        stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        status = "ok" if exc[0] is None else f"FAILED: {exc[1]!r}"
        with open(LOG, "a") as fh:
            fh.write(f"{stamp}\t{self.name}\t{el:.1f}s\t{status}\n")
        print(f"[step] {self.name}: {el:.1f}s {status}", flush=True)
        return False


def safe_merge(left, right, on=None, left_on=None, right_on=None, how="left", validate=None, name=""):
    """Merge with assertions: no NaN keys, key uniqueness per `validate`, and row counts logged."""
    lk = left_on if left_on is not None else on
    rk = right_on if right_on is not None else on
    lk = [lk] if isinstance(lk, str) else list(lk)
    rk = [rk] if isinstance(rk, str) else list(rk)
    assert not left[lk].isna().any().any(), f"{name}: NaN in left keys"
    assert not right[rk].isna().any().any(), f"{name}: NaN in right keys"
    n0 = len(left)
    out = left.merge(right, on=on, left_on=left_on, right_on=right_on, how=how, validate=validate)
    if how == "left" and validate in ("m:1", "1:1"):
        assert len(out) == n0, f"{name}: row count changed {n0} -> {len(out)}"
    print(f"[join] {name}: left {n0} rows, right {len(right)} rows, out {len(out)} rows", flush=True)
    return out
