"""Primer-arm anchored demultiplexing with optional handle corroboration.

Both six-base indexes are required. This heuristic can tolerate missing handles
and some arm errors; it cannot guarantee recovery of every error pattern.

Read layout, forward orientation:

    [handle_f][index_f][arm_f] ... plasmid ... [rc(arm_r)][rc(index_r)][rc(handle_r)]
                       ^ anchor                 ^ anchor
"""
from __future__ import annotations
import os
import csv
from contextlib import ExitStack
import edlib
from collections import Counter
from dataclasses import dataclass, asdict
from .seqio import rc, mismatches, read_fastq
from .primersets import max_safe_error_budget

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


def _arm_locations(read, arm, max_edit):
    result = edlib.align(arm, read, mode='HW', task='locations', k=max_edit)
    if result['editDistance'] < 0:
        return []
    return [(start, end + 1, result['editDistance']) for start, end in result['locations']]


def find_arm(read: str, arm: str, max_edit: int, slack: int = 4):
    """Compatibility helper: first optimal start, with full-query edit matching."""
    hits = _arm_locations(read, arm, max_edit)
    return min((start for start, _, _ in hits), default=-1)


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
    fwd_arm_start: int | None = None
    fwd_arm_end: int | None = None
    rev_arm_start: int | None = None
    rev_arm_end: int | None = None
    fwd_index_mismatches: int | None = None
    rev_index_mismatches: int | None = None
    assignment_edits: int | None = None


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


def _end_candidates(seq, arm, indexes, reverse, arm_mm, index_mm, shift):
    offset = max(0, len(seq) - 250) if reverse else 0
    target = seq[offset:] if reverse else seq[:250]
    hits = []
    for start, end, edits in _arm_locations(target, arm, arm_mm):
        start, end = start + offset, end + offset
        for delta in range(-shift, shift + 1):
            left = end + delta if reverse else start + delta - 6
            if left < 0 or left + 6 > len(seq):
                continue
            observed = rc(seq[left:left + 6]) if reverse else seq[left:left + 6]
            name, distance = call_index(observed, indexes, index_mm)
            if name:
                hits.append(dict(name=name, start=start, end=end, index_left=left,
                                 arm_edits=edits, index_edits=distance, score=edits + distance))
    return hits


def call_read(seq: str, pset, cfg, arm_mm: int = 4, index_mm: int = 0, index_shift: int = 1):
    """Compare valid end/index pairs across both orientations; reject equal-best wells."""
    complete, partial = [], []
    for orientation, s in (('+', seq), ('-', rc(seq))):
        forward = _end_candidates(s, pset.arm_f, cfg.fwd_index, False, arm_mm, index_mm, index_shift)
        reverse = _end_candidates(s, rc(pset.arm_r), cfg.rev_index, True, arm_mm, index_mm, index_shift)
        for f in forward:
            for r in reverse:
                if f['end'] >= r['start']:
                    continue
                left = f['index_left'] - len(pset.handle_f)
                right = r['index_left'] + 6
                complete.append(ReadCall(
                    name='', length=len(seq), orientation=orientation,
                    fwd_index=f['name'], rev_index=r['name'],
                    fwd_arm_start=f['start'], fwd_arm_end=f['end'],
                    rev_arm_start=r['start'], rev_arm_end=r['end'],
                    fwd_index_mismatches=f['index_edits'], rev_index_mismatches=r['index_edits'],
                    assignment_edits=f['score'] + r['score'],
                    handle_f_seen=left >= 0 and _edit(s[left:f['index_left']], pset.handle_f, 4) <= 4,
                    handle_r_seen=right + len(pset.handle_r) <= len(s)
                        and _edit(rc(s[right:right + len(pset.handle_r)]), pset.handle_r, 4) <= 4))
        if forward or reverse:
            partial.append(ReadCall(name='', length=len(seq), orientation=orientation,
                                    fwd_index=forward[0]['name'] if forward else None,
                                    rev_index=reverse[0]['name'] if reverse else None,
                                    status='one_index_only' if bool(forward) != bool(reverse) else 'ambiguous'))
    if not complete:
        return partial[0] if partial else ReadCall(name='', length=len(seq))
    score = min(c.assignment_edits for c in complete)
    best = [c for c in complete if c.assignment_edits == score]
    if len({(c.fwd_index, c.rev_index) for c in best}) != 1:
        return ReadCall(name='', length=len(seq), status='ambiguous')
    call = best[0]
    call.well = cfg.well_of(call.fwd_index, call.rev_index)
    call.status = 'assigned' if call.well else 'impossible_pair'
    return call


def demux(fastq_path, pset, cfg, outdir, min_len=1000, arm_mm=4, index_mm=0,
          index_shift=1, write_reads=True, used_wells=None, collect_calls=False):
    """Demultiplex a pooled FASTQ. Returns a dict of per-well and run-level stats."""
    safe_budget = min(max_safe_error_budget(cfg.fwd_index), max_safe_error_budget(cfg.rev_index))
    if min_len < 1 or arm_mm < 0 or not 0 <= index_mm <= safe_budget or not 0 <= index_shift <= 2:
        raise ValueError(f"invalid demux settings; index mismatch budget must be 0..{safe_budget}")
    if write_reads and os.path.isdir(outdir) and os.listdir(outdir):
        raise ValueError("demux output already contains files; choose a fresh output directory")
    os.makedirs(outdir, exist_ok=True)
    lengths, calls = {}, []
    n_total = n_short = 0
    status_counts = Counter()
    seen = set()
    with ExitStack() as stack:
        handles = {}
        writer = None
        if write_reads:
            audit = stack.enter_context(open(os.path.join(outdir, 'read_assignments.tsv'), 'w', newline=''))
            writer = csv.DictWriter(audit, fieldnames=list(ReadCall.__dataclass_fields__), delimiter='\t')
            writer.writeheader()
        for name, seq, qual in read_fastq(fastq_path):
            if name in seen:
                raise ValueError(f'duplicate FASTQ read ID: {name}; duplicate molecules cannot inflate depth')
            seen.add(name)
            n_total += 1
            if len(seq) < min_len:
                n_short += 1
                c = ReadCall(name=name, length=len(seq), status='below_min_length')
            else:
                c = call_read(seq, pset, cfg, arm_mm, index_mm, index_shift)
                c.name = name
                status_counts[c.status] += 1
            if collect_calls:
                calls.append(asdict(c))
            if writer:
                writer.writerow(asdict(c))
            key = c.well if c.status == 'assigned' else c.status
            if c.status == 'assigned':
                lengths.setdefault(key, []).append(len(seq))
            if write_reads:
                if key not in handles:
                    handles[key] = stack.enter_context(open(os.path.join(outdir, f'{key}.fastq'), 'w'))
                s, q = (rc(seq), qual[::-1]) if c.orientation == '-' else (seq, qual)
                handles[key].write(f'@{name}\n{s}\n+\n{q}\n')

    all_wells = [w for w, _, _ in cfg.wells()]
    used = set(used_wells) if used_wells is not None else None
    per_well = []
    for w in all_wells:
        lens = sorted(lengths.get(w, []))
        per_well.append({
            "well": w,
            "n_reads": len(lens),
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
        "arm_matching": "edlib full-query alignment; variable aligned arm span",
        "audit_file": "demux/read_assignments.tsv",
        "sentinel_wells_with_reads": [(p["well"], p["n_reads"]) for p in sentinel],
        "sentinel_reads": sum(p["n_reads"] for p in sentinel),
    }
    return {"run": run, "per_well": per_well, "calls": calls}
