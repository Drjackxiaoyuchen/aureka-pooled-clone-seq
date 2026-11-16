"""Primer sets, index sets, plate geometry, and template checks."""
from __future__ import annotations
import os, itertools
from dataclasses import dataclass, field
from .seqio import rc, mismatches

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

_HERE = os.path.dirname(__file__)
_CANDIDATES = [
    os.path.join(_HERE, "data", "primer_sets.yaml"),
    os.path.join(_HERE, "..", "..", "primers", "primer_sets.yaml"),
]


def default_yaml() -> str:
    for c in _CANDIDATES:
        if os.path.exists(c):
            return os.path.abspath(c)
    raise FileNotFoundError("primer_sets.yaml not found; pass --primer-yaml")


@dataclass
class PrimerSet:
    name: str
    handle_f: str
    handle_r: str
    arm_f: str
    arm_r: str
    blind_spot_bp: int = 0
    blind_spot_note: str = ""
    description: str = ""
    status: str = ""
    validated_on: str = ""

    def primer_f(self, index_seq: str) -> str:
        return self.handle_f + index_seq + self.arm_f

    def primer_r(self, index_seq: str) -> str:
        return self.handle_r + index_seq + self.arm_r


@dataclass
class Config:
    sets: dict
    fwd_index: dict
    rev_index: dict
    rows: list
    columns: list
    row_labels: list
    column_labels: list

    def wells(self):
        """Yield (well_label, fwd_index_name, rev_index_name) for every plate position."""
        for ri, fname in enumerate(self.rows):
            for ci, rname in enumerate(self.columns):
                yield f"{self.row_labels[ri]}{self.column_labels[ci]}", fname, rname

    def well_of(self, fname, rname):
        try:
            return f"{self.row_labels[self.rows.index(fname)]}{self.column_labels[self.columns.index(rname)]}"
        except ValueError:
            return None


def load_config(path: str | None = None) -> Config:
    path = path or default_yaml()
    if yaml is None:
        raise RuntimeError("pyyaml is required: pip install pyyaml")
    with open(path) as fh:
        d = yaml.safe_load(fh)
    sets = {k: PrimerSet(name=k, **v) for k, v in d["sets"].items()}
    p = d["plate"]
    for direction in ("forward", "reverse"):
        ix = d["indexes"][direction]
        if len(ix) < 2 or len(set(ix.values())) != len(ix):
            raise ValueError(f"{direction} indexes must be distinct; at least two are required")
        if any(len(s) != 6 or set(s) - set("ACGT") for s in ix.values()):
            raise ValueError("indexes must be six uppercase A/C/G/T bases")
    for key, labels, direction in (("rows", "row_labels", "forward"), ("columns", "column_labels", "reverse")):
        if (len(p[key]) != len(p[labels]) or not p[key]
                or len(set(p[key])) != len(p[key]) or len(set(map(str, p[labels]))) != len(p[labels])
                or any(x not in d["indexes"][direction] for x in p[key])):
            raise ValueError(f"invalid plate {key} or {labels}")
    if any(len(str(x)) != 1 or not str(x).isalpha() for x in p["row_labels"]):
        raise ValueError("row labels must be single letters")
    if any(not str(x).isdigit() or int(x) < 1 for x in p['column_labels']):
        raise ValueError('column labels must be positive integers')
    for ps in sets.values():
        if any(not s or set(s) - set("ACGT") for s in (ps.arm_f, ps.arm_r, ps.handle_f, ps.handle_r)):
            raise ValueError("primer arms and handles must be uppercase A/C/G/T DNA")
    cfg = Config(sets=sets, fwd_index=d["indexes"]["forward"], rev_index=d["indexes"]["reverse"],
                  rows=p["rows"], columns=p["columns"],
                  row_labels=p["row_labels"], column_labels=[str(c) for c in p["column_labels"]])
    if len({w for w, _, _ in cfg.wells()}) != len(cfg.rows) * len(cfg.columns):
        raise ValueError("plate well labels must be unique")
    return cfg


# ---------------------------------------------------------------- index QC

def index_distances(indexes: dict):
    """Min Hamming and min Levenshtein within an index set, plus the closest pairs."""
    def lev(a, b):
        d = list(range(len(b) + 1))
        for i, ca in enumerate(a, 1):
            prev, d[0] = d[0], i
            for j, cb in enumerate(b, 1):
                cur = d[j]
                d[j] = min(d[j] + 1, d[j - 1] + 1, prev + (ca != cb))
                prev = cur
        return d[-1]
    names = list(indexes)
    pairs = []
    for a, b in itertools.combinations(names, 2):
        pairs.append((mismatches(indexes[a], indexes[b], 99), lev(indexes[a], indexes[b]), a, b))
    mh = min(p[0] for p in pairs)
    ml = min(p[1] for p in pairs)
    return {"min_hamming": mh, "min_levenshtein": ml,
            "closest": [(a, b, h, l) for h, l, a, b in pairs if l == ml]}


def max_safe_error_budget(indexes: dict) -> int:
    """Errors per index that still leave assignment unambiguous: floor((d_min-1)/2)."""
    return (index_distances(indexes)["min_levenshtein"] - 1) // 2


# ------------------------------------------------------------ template check

def find_sites(seq: str, probe: str, max_mm: int = 0):
    """All circular positions where probe matches seq (either strand), within max_mm."""
    n = len(seq)
    if n < len(probe):
        return []
    dbl = seq + seq[: len(probe)]
    hits = []
    for strand, p in (("+", probe), ("-", rc(probe))):
        for i in range(n):
            if mismatches(dbl[i:i + len(probe)], p, max_mm) <= max_mm:
                hits.append((i + 1, strand))
    return hits


def check_template(seq: str, pset: PrimerSet, max_mm: int = 0):
    """Does this plasmid give exactly one clean product with this primer set?"""
    f = find_sites(seq, pset.arm_f, max_mm)
    r = find_sites(seq, pset.arm_r, max_mm)
    out = {"set": pset.name, "length": len(seq), "n_fwd_sites": len(f), "n_rev_sites": len(r),
           "fwd_sites": f, "rev_sites": r, "amplicon_bp": None, "blind_spot_bp": None, "ok": False,
           "reason": ""}
    if len(f) != 1 or len(r) != 1:
        out["reason"] = (f"expected one site per arm, found {len(f)} forward and {len(r)} reverse"
                         " - this template will give more than one product or none")
        return out
    n = len(seq)
    (fp, fs), (rp, rs) = f[0], r[0]
    if fs == rs:
        out["reason"] = "both arms anneal to the same strand - they are not facing outward here"
        return out
    # forward-reading primer start, reverse-reading primer 5' end
    if fs == "+":
        start = fp - 1
        end = (rp - 1 + len(pset.arm_r)) % n
    else:
        start = rp - 1
        end = (fp - 1 + len(pset.arm_f)) % n
    amp = (end - start) % n or n
    full_amp = amp + len(pset.handle_f) + len(pset.handle_r) + 12
    out.update(amplicon_bp=full_amp, template_span_bp=amp, blind_spot_bp=n - amp, ok=True,
               reason=f"one predicted product of {full_amp} bp including tails; "
                      f"{amp} bp template span, {n - amp} bp blind spot; wet-lab specificity untested")
    return out
