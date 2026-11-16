"""Primer-arm anchored demultiplexing with optional handle corroboration.

Both six-base indexes are required. This heuristic can tolerate missing handles
and some arm errors; it cannot guarantee recovery of every error pattern.

Read layout, forward orientation:

    [handle_f][index_f][arm_f] ... plasmid ... [rc(arm_r)][rc(index_r)][rc(handle_r)]
                       ^ anchor                 ^ anchor
"""
from __future__ import annotations
import os
from collections import Counter
from dataclasses import dataclass, asdict
from .seqio import rc, mismatches, read_fastq, write_fastq
from .primersets import max_safe_error_budget

SEED = 10        # exact-match seed used to locate the arm
SEED_STEP = 4    # overlapping seeds, so one indel cannot destroy every seed


def _seed_positions(arm: str):
    """Overlapping seeds covering the arm, with their offsets.

    Overlapping matters: a single deleted base destroys any seed spanning it,
    so non-overlapping seeds can leave an arm unfindable."""
    stops = list(range(0, max(1, len(arm) - SEED + 1), SEED_STEP))
    last = len(arm) - SEED
    if last > 0 and stops[-1] != last:
        stops.append(last)
    return [(i, arm[i:i + SEED]) for i in stops] or [(0, arm[:SEED])]


def _edit(a: str, b: str, limit: int) -> int:
    """Levenshtein with early exit. Both strings are short (<=30)."""
    if abs(len(a) - len(b)) > limit:
        return limit + 1
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(b)
        best = cur[0]
        for j, cb in enumerate(b, 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (ca != cb))
            best = min(best, cur[j])
        if best > limit:
            return limit + 1
        prev = cur
    return prev[-1]


def find_arm(read: str, arm: str, max_edit: int, slack: int = 4):
    """Read position aligned to arm[0], or -1.

    Nanopore reads carry indels, so a Hamming-only check throws away arms that
    are perfectly recognisable: one deleted base shifts every later position.
    We seed on an exact 10-mer, then score a small window with edit distance
    and keep the offset that fits best."""
    L = len(arm)
    cands = set()
    for off, seed in _seed_positions(arm):
        start = 0
        while True:
            j = read.find(seed, start)
            if j < 0:
                break
            p = j - off
            if 0 <= p <= len(read) - L + slack:
                cands.add(p)
            start = j + 1
    best, best_d = -1, max_edit + 1
    for p in sorted(cands):
        for d in range(-slack, slack + 1):
            q = p + d
            if q < 0 or q + L - slack > len(read):
                continue
            dist = _edit(read[q:q + L], arm, max_edit)
            if dist < best_d or (dist == best_d and best >= 0 and q < best):
                best, best_d = q, dist
        if best_d == 0:
            break
    return best if best_d <= max_edit else -1


def call_index(obs: str, index_map: dict, max_mm: int = 1):
    """Match an observed 6-mer to an index set. Returns (name, distance) or (None, reason)."""
    if len(obs) != 6 or "N" in obs:
        return None, "short"
    hits = []
    for name, seq in index_map.items():
        d = mismatches(obs, seq, max_mm)
        if d <= max_mm:
            hits.append((d, name))
    if not hits:
        return None, "nomatch"
    hits.sort()
    if len(hits) > 1 and hits[0][0] == hits[1][0]:
        return None, "ambiguous"
    return hits[0][1], hits[0][0]


@dataclass
class ReadCall:
    name: str
    length: int
    orientation: str = "."
    fwd_index: str | None = None
    rev_index: str | None = None
    well: str | None = None
    status: str = "unassigned"
    handle_f_seen: bool = False
    handle_r_seen: bool = False


def _index_near(s: str, end: int, index_map: dict, index_mm: int, shift: int = 1):
    """Read the index immediately 5' of `end`, tolerating a small indel shift."""
    hits = []
    for d in range(-shift, shift + 1):
        e = end + d
        if e - 6 < 0 or e > len(s):
            continue
        name, dist = call_index(s[e - 6:e], index_map, index_mm)
        if name:
            hits.append((name, dist))
    if not hits:
        return None
    best_d = min(d for _, d in hits)
    names = {n for n, d in hits if d == best_d}
    return next(iter(names)) if len(names) == 1 else None


def _index_near_rev(s: str, start: int, index_map: dict, index_mm: int, shift: int = 1):
    """Read the (reverse-complemented) index immediately 3' of `start`."""
    hits = []
    for d in range(-shift, shift + 1):
        b = start + d
        if b < 0 or b + 6 > len(s):
            continue
        name, dist = call_index(rc(s[b:b + 6]), index_map, index_mm)
        if name:
            hits.append((name, dist))
    if not hits:
        return None
    best_d = min(d for _, d in hits)
    names = {n for n, d in hits if d == best_d}
    return next(iter(names)) if len(names) == 1 else None


def call_read(seq: str, pset, cfg, arm_mm: int = 4, index_mm: int = 1, index_shift: int = 1):
    """Try both orientations; return the best ReadCall."""
    best = None
    complete = []
    for orient, s in (("+", seq), ("-", rc(seq))):
        # The barcodes must belong to the amplicon ends. Internal primer-like
        # sequences in concatamers or unrelated molecules cannot identify a well.
        fp = find_arm(s[:250], pset.arm_f, arm_mm)
        tail_start = max(0, len(s) - 250)
        rp = find_arm(s[tail_start:], rc(pset.arm_r), arm_mm)
        rp = rp + tail_start if rp >= 0 else -1
        fi = ri = None
        hf = hr = False
        if fp >= 6:
            fi = _index_near(s, fp, cfg.fwd_index, index_mm, index_shift)
            hs = fp - 6 - len(pset.handle_f)
            if hs >= 0:
                hf = _edit(s[hs:fp - 6], pset.handle_f, 4) <= 4
        if rp >= 0:
            e = rp + len(pset.arm_r)
            ri = _index_near_rev(s, e, cfg.rev_index, index_mm, index_shift)
            he = e + 6 + len(pset.handle_r)
            if ri and he <= len(s):
                hr = _edit(rc(s[e + 6:he]), pset.handle_r, 4) <= 4
        score = (fi is not None) + (ri is not None)
        cand = ReadCall(name="", length=len(seq), orientation=orient, fwd_index=fi, rev_index=ri,
                        handle_f_seen=hf, handle_r_seen=hr)
        if score == 2 and fp + len(pset.arm_f) < rp:
            complete.append(cand)
        if best is None or score > (best.fwd_index is not None) + (best.rev_index is not None):
            best = cand
    pairs = {(c.fwd_index, c.rev_index) for c in complete}
    if len(pairs) > 1:
        return ReadCall(name="", length=len(seq), status="ambiguous")
    if complete:
        best = complete[0]
    elif best.fwd_index and best.rev_index:
        return ReadCall(name="", length=len(seq), status="ambiguous")
    if best.fwd_index and best.rev_index:
        best.well = cfg.well_of(best.fwd_index, best.rev_index)
        best.status = "assigned" if best.well else "impossible_pair"
    elif best.fwd_index or best.rev_index:
        best.status = "one_index_only"
    else:
        best.status = "unassigned"
    return best


def demux(fastq_path, pset, cfg, outdir, min_len=1000, arm_mm=4, index_mm=1,
          index_shift=1, write_reads=True, used_wells=None):
    """Demultiplex a pooled FASTQ. Returns a dict of per-well and run-level stats."""
    safe_budget = min(max_safe_error_budget(cfg.fwd_index), max_safe_error_budget(cfg.rev_index))
    if min_len < 1 or arm_mm < 0 or not 0 <= index_mm <= safe_budget or not 0 <= index_shift <= 2:
        raise ValueError(f"invalid demux settings; index mismatch budget must be 0..{safe_budget}")
    if write_reads and os.path.isdir(outdir) and os.listdir(outdir):
        raise ValueError("demux output already contains files; choose a fresh output directory")
    os.makedirs(outdir, exist_ok=True)
    buckets, calls = {}, []
    n_total = n_short = 0
    status_counts = Counter()
    for name, seq, qual in read_fastq(fastq_path):
        n_total += 1
        if len(seq) < min_len:
            n_short += 1
            continue
        c = call_read(seq, pset, cfg, arm_mm, index_mm, index_shift)
        c.name = name
        calls.append(c)
        status_counts[c.status] += 1
        key = c.well if c.status == "assigned" else (
            "one_index_only" if c.status == "one_index_only" else "unassigned")
        s, q = (seq, qual) if c.orientation == "+" else (rc(seq), qual[::-1])
        buckets.setdefault(key, []).append((name, s, q))

    if write_reads:
        for key, recs in buckets.items():
            write_fastq(os.path.join(outdir, f"{key}.fastq"), recs)

    all_wells = [w for w, _, _ in cfg.wells()]
    used = set(used_wells) if used_wells is not None else None
    per_well = []
    for w in all_wells:
        recs = buckets.get(w, [])
        lens = sorted(len(s) for _, s, _ in recs)
        per_well.append({
            "well": w,
            "n_reads": len(recs),
            "median_len": lens[len(lens) // 2] if lens else 0,
            "max_len": lens[-1] if lens else 0,
            "in_use": (w in used) if used is not None else True,
        })
    assigned = sum(p["n_reads"] for p in per_well)
    sentinel = [p for p in per_well if not p["in_use"] and p["n_reads"] > 0]
    run = {
        "fastq": os.path.basename(fastq_path),
        "primer_set": pset.name,
        "reads_total": n_total,
        "reads_below_min_len": n_short,
        "reads_considered": n_total - n_short,
        "reads_assigned": assigned,
        "assignment_rate": round(assigned / max(1, n_total - n_short), 4),
        "one_index_only": status_counts["one_index_only"],
        "unassigned": status_counts["unassigned"],
        "impossible_pair": status_counts["impossible_pair"],
        "ambiguous": status_counts["ambiguous"],
        "min_len": min_len,
        "arm_edits_allowed": arm_mm,
        "index_mismatches_allowed": index_mm,
        "index_shift_allowed": index_shift,
        "sentinel_wells_with_reads": [(p["well"], p["n_reads"]) for p in sentinel],
        "sentinel_reads": sum(p["n_reads"] for p in sentinel),
    }
    return {"run": run, "per_well": per_well, "calls": [asdict(c) for c in calls]}
