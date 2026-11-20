"""Configuration for task G (config/G_defaults.json).

Every DECIDE choice of the Part B pre-registration draft (7a to 7n, 8a, 8b) is read through Config.decide(),
which validates the option and records that the code read it. The run manifest lists the options read, and the
end-to-end test asserts that every DECIDE ID was read.
"""
import json, os, hashlib

DECIDE_IDS = ["7a", "7b", "7c", "7d", "7e", "7f", "7g", "7h", "7i", "7j", "7k", "7l", "7m", "7n", "8a", "8b"]


def path(root, p):
    return p if os.path.isabs(p) else os.path.join(root, p)


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Config:
    def __init__(self, raw, file_path, track=True):
        self.raw = raw
        self.file_path = file_path
        self.sha256 = sha256_file(file_path)
        self.track = track
        self.read_ids = []
        self.validate()

    def validate(self):
        dec, opt = self.raw["decide"], self.raw["decide_options"]
        missing = [i for i in DECIDE_IDS if i not in dec]
        extra = [i for i in dec if i not in DECIDE_IDS]
        assert not missing and not extra, f"decide IDs missing {missing}, unknown {extra}"
        for i in DECIDE_IDS:
            v = dec[i]
            if i == "7j":
                assert isinstance(v, dict) and set(v) == {"unit", "interval"}, "7j needs unit and interval"
                for part in ("unit", "interval"):
                    assert v[part] in opt["7j"][part], f"7j {part}: option {v[part]} not in {list(opt['7j'][part])}"
            else:
                assert v in opt[i], f"{i}: option {v!r} not in {list(opt[i])}"

    def decide(self, i):
        assert i in DECIDE_IDS, i
        if self.track and i not in self.read_ids:
            self.read_ids.append(i)
        return self.raw["decide"][i]

    def const(self, name):
        return self.raw["plan_constants"][name]

    def run(self, name):
        return self.raw["run"][name]

    def decide_table(self):
        rows = []
        for i in DECIDE_IDS:
            v = self.raw["decide"][i]
            if i == "7j":
                desc = f"unit {v['unit']}: {self.raw['decide_options']['7j']['unit'][v['unit']]}; interval {v['interval']}: " \
                       f"{self.raw['decide_options']['7j']['interval'][v['interval']]}"
                v = f"unit={v['unit']};interval={v['interval']}"
            else:
                desc = self.raw["decide_options"][i][v]
            rows.append(dict(decide_id=i, option=v, meaning=desc, read_by_code=i in self.read_ids))
        return rows


def load(file_path, track=True):
    with open(file_path) as f:
        raw = json.load(f)
    return Config(raw, file_path, track=track)
