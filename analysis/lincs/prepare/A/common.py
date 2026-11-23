"""Shared paths, constants and helpers for task A."""
import datetime as dt
import os
import re
import time
from contextlib import contextmanager

ROOT = os.environ["P59_LINCS_WORK"]
RAW = os.path.join(ROOT, "data", "A", "raw")
DERIVED = os.path.join(ROOT, "data", "A", "derived")  # intermediate, not committed
RES = os.path.join(ROOT, "results", "A")
FIG = os.path.join(ROOT, "figures", "A")
LOG = os.path.join(ROOT, "logs", "A_runlog.txt")
P59 = os.environ["P59_PUBLIC_REPOSITORY"]
for d in (DERIVED, RES, FIG):
    os.makedirs(d, exist_ok=True)

SEED = 20261114
XIST = "XIST"
Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y", "NLGN4Y"]
MARKERS = [XIST] + Y_GENES
STATES = ["XX-like", "XY", "Y-lost", "XIST-lost", "discordant", "ambiguous", "no_data"]
CONF = ["high", "medium", "low", "none"]
LEVEL3_URL = ("https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level3/"
              "level3_beta_ctl_n188708x12328.gctx")


def norm(s):
    """Upper case and drop everything that is not A-Z or 0-9."""
    if s is None or (isinstance(s, float) and s != s):
        return ""
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


def log(step, msg):
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(LOG, "a") as fh:
        fh.write(f"{stamp}\t{step}\t{msg}\n")
    print(f"{stamp} [{step}] {msg}", flush=True)


@contextmanager
def timed(step):
    t0 = time.time()
    log(step, "start")
    yield
    log(step, f"done in {time.time() - t0:.1f} s")
