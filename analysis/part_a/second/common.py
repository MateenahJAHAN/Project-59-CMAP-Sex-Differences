import time, datetime, os
SEED = 20261115
Y8 = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y", "NLGN4Y"]
YP = ["RPS4Y1", "ZFY"]
YQ = ["DDX3Y", "USP9Y", "UTY", "NLGN4Y", "KDM5D", "EIF1AY"]
X8 = ["RPS4X", "DDX3X", "KDM5C", "KDM6A", "EIF1AX", "ZFX", "USP9X", "NLGN4X"]
MARKERS = ["XIST", "TSIX"] + Y8 + ["TXLNGY", "PRKY"] + X8 + ["ACTB", "GAPDH"]
NAMED_SPARSE_LINES = ["A549", "SW48", "NCI-H460", "LoVo"]
DEEP_SPARSE = 2500

def runlog(step, t0, msg=""):
    os.makedirs("logs", exist_ok=True)
    with open("logs/F2_runlog.txt", "a") as f:
        f.write(f"{datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}\t{step}\t{time.time()-t0:.1f}s\t{msg}\n")
