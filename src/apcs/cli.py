"""apcs command line.

    apcs plan      -n 30                      plate assignment + bench checklist
    apcs order     --set ampr                 oligo order sheets
    apcs check     ref.gb --set ampr          will this plasmid amplify cleanly?
    apcs indexqc                              index set distances and error budget
    apcs gel       gel.jpg --expected-bp 5700 band QC + pooling volumes
    apcs run       reads.fastq --expect ref.gb   demux, genotype, pick plate
"""
from __future__ import annotations
import argparse, csv, json, os, sys, glob

from . import __version__
from .primersets import load_config, check_template, index_distances, max_safe_error_budget
from .seqio import load_reference, read_fastq, rc


def _cfg(a):
    return load_config(a.primer_yaml)


def _set(cfg, name):
    if name not in cfg.sets:
        sys.exit(f"unknown primer set {name!r}; have: {', '.join(cfg.sets)}")
    return cfg.sets[name]


def _read_samples(path, cfg):
    """samples.csv: well,sample[,expected]. Returns ordered list of dicts."""
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames or "well" not in [k.strip().lower() for k in reader.fieldnames]:
            raise ValueError("samples.csv requires a well column")
        headers = [k.strip().lower() for k in reader.fieldnames]
        if len(set(headers)) != len(headers):
            raise ValueError("samples.csv contains duplicate column names")
        for r in reader:
            if None in r:
                raise ValueError("samples.csv row has too many columns")
            r = {k.strip().lower(): (v or "").strip() for k, v in r.items()}
            if not r.get("well"):
                if any(r.values()):
                    raise ValueError("samples.csv contains a nonempty row without a well")
                continue
            rows.append(r)
    known = {w for w, _, _ in cfg.wells()}
    bad = [r["well"] for r in rows if r["well"] not in known]
    if bad:
        sys.exit(f"samples.csv names wells that are not on this plate: {', '.join(bad)}")
    if not rows or len({r['well'] for r in rows}) != len(rows):
        raise ValueError("samples.csv must contain nonempty, unique well assignments")
    pairs = {w: (f, r) for w, f, r in cfg.wells()}
    for r in rows:
        for key, actual in zip(("fwd_index", "rev_index"), pairs[r["well"]]):
            if r.get(key) and r[key] != actual:
                raise ValueError(f"{r['well']}: {key} disagrees with plate configuration")
    return rows


# --------------------------------------------------------------------- plan

def _block(cfg, n):
    """Smallest rectangle of forward x reverse indexes that holds n clones.

    A rectangle means fewer primer tubes on the bench than a row-major fill:
    16 clones as 4x4 needs 8 tubes, as 3 rows of 6 it needs 9."""
    nr, nc = len(cfg.rows), len(cfg.columns)
    best = None
    for r in range(1, nr + 1):
        for c in range(1, nc + 1):
            if r * c >= n:
                key = (r + c, r * c)
                if best is None or key < best[0]:
                    best = (key, r, c)
    if best is None:
        return nr, nc
    return best[1], best[2]


def cmd_plan(a):
    cfg = _cfg(a)
    ps = _set(cfg, a.set)
    wells = list(cfg.wells())
    if a.n < 1 or a.n > len(wells):
        sys.exit(f"{a.n} clones but only {len(wells)} index combinations exist in this plate")
    if a.wells:
        want = [w.strip() for w in a.wells.split(",") if w.strip()]
        known = {w: (f, r) for w, f, r in wells}
        bad = [w for w in want if w not in known]
        if bad:
            sys.exit(f"not plate positions: {', '.join(bad)}")
        if len(want) != a.n or len(set(want)) != len(want):
            raise ValueError("--wells must list exactly -n distinct positions")
        chosen = [(w, *known[w]) for w in want]
    elif a.layout == "block":
        nr, nc = _block(cfg, a.n)
        grid = [(w, f, r) for w, f, r in wells
                if cfg.row_labels.index(w[0]) < nr and cfg.column_labels.index(w[1:]) < nc]
        chosen = grid[: a.n]
    else:
        chosen = wells[: a.n]
    os.makedirs(a.out, exist_ok=True)
    if os.path.exists(os.path.join(a.out, "samples.csv")):
        raise ValueError("samples.csv already exists; choose a new plan directory")
    with open(os.path.join(a.out, "samples.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["well", "sample", "expected", "fwd_index", "rev_index",
                    "fwd_primer_tube", "rev_primer_tube"])
        for well, f, r in chosen:
            w.writerow([well, "", "", f, r, f"{ps.name}-F-{f}", f"{ps.name}-R-{r}"])
    chosen_wells = {w for w, _, _ in chosen}
    sentinels = [w for w, _, _ in wells if w not in chosen_wells]
    lines = [
        f"# Pooling plan - {a.n} clones, primer set {ps.name}",
        "",
        f"Fill in the `sample` column of `samples.csv` as you pick colonies. One colony per well.",
        "",
        f"- forward primer tubes needed: {', '.join(sorted({f for _, f, _ in chosen}))}",
        f"- reverse primer tubes needed: {', '.join(sorted({r for _, _, r in chosen}))}",
        f"- unused index combinations ({len(sentinels)}): these are your contamination sentinels.",
        f"  Monitor reads in these positions as a contamination/cross-assignment signal.",
        f"  Unused positions: {', '.join(sentinels) or '(none; no sentinel monitoring)'}",
        "",
        "## Bench order",
        "",
        "1. Colony PCR, one colony per well, with that well's primer pair. 25 uL is plenty.",
        "2. Run every well on a gel. `apcs gel` sizes the bands and tells you which wells to drop",
        "   and how many uL of each to pool.",
        "3. Pool only the wells that passed. Clean up the pool once.",
        "4. Submit the pool as ONE sample. Note the expected amplicon size on the order.",
        "5. `apcs run` on the FASTQ you get back.",
        "",
        "Do not pool a well with no band, the wrong band, or two bands. It cannot give a clean",
        "answer and it takes read depth away from the wells that can.",
    ]
    open(os.path.join(a.out, "README.md"), "w").write("\n".join(lines) + "\n")
    print(f"wrote {a.out}/samples.csv and {a.out}/README.md for {a.n} clones")
    print(f"sentinel (unused) combinations: {len(sentinels)}")


# -------------------------------------------------------------------- order

def cmd_order(a):
    from . import order as order_mod
    cfg = _cfg(a)
    ps = _set(cfg, a.set)
    rows = order_mod.rows_for(ps, cfg, scale=a.scale, purification=a.purification)
    os.makedirs(a.out, exist_ok=True)
    csv_path = os.path.join(a.out, f"order_{ps.name}_primers.csv")
    order_mod.write_csv(csv_path, rows)
    print(f"wrote {csv_path}  ({len(rows)} oligos)")
    try:
        xp = os.path.join(a.out, f"order_{ps.name}_primers.xlsx")
        order_mod.write_vendor_xlsx(xp, rows, title=f"{ps.name} primer set")
        print(f"wrote {xp}")
    except RuntimeError as e:
        print(f"(no xlsx: {e})")
    for r in rows:
        print(f"  {r['name']:<16} {r['sequence']}  {r['length']} nt")


# -------------------------------------------------------------------- check

def cmd_check(a):
    cfg = _cfg(a)
    sets = [cfg.sets[a.set]] if a.set else list(cfg.sets.values())
    paths = []
    for p in a.reference:
        paths.extend(sorted(glob.glob(p)) or [p])
    failed = False
    for p in paths:
        name, seq = load_reference(p)
        print(f"\n{os.path.basename(p)}  ({len(seq)} bp)")
        for ps in sets:
            r = check_template(seq, ps, a.max_mismatch)
            failed |= not r["ok"]
            flag = "OK " if r["ok"] else "NO "
            print(f"  {flag}{ps.name:<6} {r['reason']}")
    if failed:
        raise SystemExit(1)


# ------------------------------------------------------------------ indexqc

def cmd_indexqc(a):
    cfg = _cfg(a)
    for label, ix in (("forward", cfg.fwd_index), ("reverse", cfg.rev_index)):
        d = index_distances(ix)
        print(f"{label}: n={len(ix)}  min Hamming={d['min_hamming']}  "
              f"min Levenshtein={d['min_levenshtein']}  "
              f"safe errors per index={max_safe_error_budget(ix)}")
        for a_, b_, h, l in d["closest"][:4]:
            print(f"    closest: {a_} {ix[a_]} vs {b_} {ix[b_]}  H{h}/L{l}")
    print("\nAn index set with min distance 3 resolves one error and can mis-call at two.")
    print("Keep the demultiplexer at --index-mismatch 1.")


# ---------------------------------------------------------------------- gel

def cmd_gel(a):
    from . import gel as gel_mod
    labels = a.wells.split(",") if a.wells else None
    res = gel_mod.analyse(a.image, expected_bp=a.expected_bp, n_lanes=a.lanes,
                          ladder_lane=a.ladder_lane,
                          ladder_sizes=[float(x) for x in a.ladder_sizes.split(",")]
                          if a.ladder_sizes else None,
                          tolerance=a.tolerance, sample_labels=labels,
                          secondary_frac=a.secondary_frac, prom_frac=a.band_prominence,
                          crop=not a.no_crop)
    os.makedirs(a.out, exist_ok=True)
    json.dump(res, open(os.path.join(a.out, "gel.json"), "w"), indent=1)
    fields = ["lane", "label", "role", "main_bp", "intensity", "saturated", "verdict", "note",
              "pool_uL", "pool_capped"]
    with open(os.path.join(a.out, "gel_qc.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in res["lanes"]:
            w.writerow(r)
    print(f"{res['n_lanes']} lanes, sizing {'calibrated' if res['calibrated'] else 'NOT calibrated'}")
    print(f"{'lane':>4} {'label':>8} {'bp':>7} {'verdict':<16} {'pool uL':>8}  note")
    for r in res["lanes"]:
        print(f"{r['lane']:>4} {str(r.get('label','')):>8} {str(r.get('main_bp') or '-'):>7} "
              f"{r['verdict']:<16} {str(r.get('pool_uL','-')):>8}  {r.get('note','')}")
    drop = [r["label"] for r in res["lanes"] if r["role"] == "sample" and r["verdict"] != "OK"]
    print(f"\ndo not pool: {', '.join(drop) if drop else '(nothing - every lane passed)'}")
    print(f"wrote {a.out}/gel_qc.csv")


# ----------------------------------------------------------------------- run

def cmd_run(a):
    from .demux import demux
    from .genotype import RefPanel, genotype_well
    from .verdict import call_well, rank_and_label, DEFAULTS
    from .report import write_report
    from .seqio import read_fastq as _rf

    cfg = _cfg(a)
    ps = _set(cfg, a.set)
    if a.min_depth < 1 or a.min_len < 1 or a.picks < 0 or a.backups < 0:
        raise ValueError("depth/length must be positive and pick/backup counts nonnegative")
    samples = _read_samples(a.samples, cfg) if a.samples else None
    used = [s["well"] for s in samples] if samples else None
    expected_by_well = {s["well"]: s.get("expected", "") for s in samples} if samples else {}
    label_by_well = {s["well"]: s.get("sample", "") for s in samples} if samples else {}

    expected_paths = [path for p in a.expect for path in (sorted(glob.glob(p)) or [p])]
    other_paths = [path for p in (a.also or []) for path in (sorted(glob.glob(p)) or [p])]
    refs = expected_paths + other_paths
    if not refs:
        sys.exit("give at least one --expect reference (GenBank or FASTA)")
    panel = RefPanel(refs)
    expected_names = panel.names[:len(expected_paths)]
    default_expect = ', '.join(expected_names)
    for well, name in expected_by_well.items():
        if name and name not in panel.names:
            raise ValueError(f"{well}: unknown expected reference {name!r}; use a filename stem")

    thresholds = dict(DEFAULTS)
    thresholds.update(min_depth=a.min_depth, n_picks=a.picks, n_backups=a.backups)
    outdir = a.out
    if os.path.exists(os.path.join(outdir, "report")):
        raise ValueError("report already exists; choose a fresh output directory")
    demux_dir = os.path.join(outdir, "demux")
    res = demux(a.fastq, ps, cfg, demux_dir, min_len=a.min_len, arm_mm=a.arm_mismatch,
                index_mm=a.index_mismatch, index_shift=a.index_shift, used_wells=used)
    run = res["run"]
    print(f"demux: {run['reads_assigned']:,}/{run['reads_considered']:,} assigned "
          f"({run['assignment_rate']*100:.1f}%), {run['sentinel_reads']} into unused combinations")
    import hashlib
    from importlib.metadata import version
    def digest(path):
        h = hashlib.sha256()
        with open(path, 'rb') as fh:
            for block in iter(lambda: fh.read(1024 * 1024), b''):
                h.update(block)
        return h.hexdigest()
    from .primersets import default_yaml
    paths = [a.fastq, a.primer_yaml or default_yaml()] + refs + ([a.samples] if a.samples else [])
    run.update(schema_version=2, software_version=__version__, thresholds=thresholds,
               dependencies={name: version(name) for name in ('mappy', 'edlib', 'numpy')},
               input_sha256=[{'file': os.path.basename(p), 'sha256': digest(p)} for p in paths],
               reference_names=panel.names,
               template_checks={name: check_template(panel.seqs[name], ps) for name in panel.names})

    wells = []
    for p in res["per_well"]:
        w = dict(p)
        w["assigned_reads"] = p["n_reads"]
        w["label"] = label_by_well.get(w["well"], "")
        fq = os.path.join(demux_dir, f"{w['well']}.fastq")
        if w["n_reads"] and os.path.exists(fq):
            g = genotype_well(_rf(fq), panel, min_read_len=a.min_len)
        else:
            g = genotype_well([], panel)
        w.update(g)
        exp = expected_by_well.get(w["well"]) or expected_names
        if not w.get("in_use", True):
            w["verdict"], w["note"] = "NOT USED", "sentinel - no sample loaded"
        else:
            w["verdict"], w["note"] = call_well(g, exp if exp else None, thresholds)
        wells.append(w)
    rank_and_label([w for w in wells if w.get("in_use", True)], thresholds)
    for w in wells:
        w.setdefault("rank", None)

    meta = {
        "title": a.title or f"{os.path.basename(a.fastq)} - pooled clone screen",
        "subtitle": f"primer set {ps.name}; {len(used) if used else len(wells)} wells; "
                    f"expected {default_expect}",
        "blind_spot": ', '.join(f"{name}: {check['blind_spot_bp'] if check['ok'] else 'unknown'}"
                               for name, check in run['template_checks'].items()),
        "headline_note": "This is a screen - confirm the pick with a dedicated whole-plasmid run.",
        "footer": f"Aureka APCS {__version__} · {ps.name} arms · min_len {a.min_len}, "
                  f"index mismatches {a.index_mismatch}",
    }
    picks = write_report(os.path.join(outdir, "report"), cfg, wells, run, meta)
    print(f"\npick: {', '.join(w['well'] for w in picks) or 'nothing cleared the gate'}")
    print(f"report: {outdir}/report/00_PICK_THIS.md  and  plate.html")


# -------------------------------------------------------------------- parser

def build_parser():
    p = argparse.ArgumentParser(prog="apcs", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--primer-yaml", default=None, help="path to primer_sets.yaml")
    p.add_argument("--version", action="version", version=f"apcs {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("plan", help="assign clones to wells and print the bench checklist")
    q.add_argument("-n", type=int, required=True, help="how many clones")
    q.add_argument("--set", default="cole1")
    q.add_argument("--layout", choices=["block", "row"], default="block",
                   help="block (default) uses the fewest primer tubes")
    q.add_argument("--wells", default=None, help="explicit comma separated wells, e.g. a1,a2,b1")
    q.add_argument("-o", "--out", default="plan")
    q.set_defaults(func=cmd_plan)

    q = sub.add_parser("order", help="oligo order sheet for a primer set")
    q.add_argument("--set", default="ampr")
    q.add_argument("--scale", default="10 nmol")
    q.add_argument("--purification", default="PAGE")
    q.add_argument("-o", "--out", default="orders")
    q.set_defaults(func=cmd_order)

    q = sub.add_parser("check", help="does this plasmid give one clean product?")
    q.add_argument("reference", nargs="+")
    q.add_argument("--set", default=None, help="default: every set")
    q.add_argument("--max-mismatch", type=int, default=0)
    q.set_defaults(func=cmd_check)

    q = sub.add_parser("indexqc", help="index set distances and the safe error budget")
    q.set_defaults(func=cmd_indexqc)

    q = sub.add_parser("gel", help="colony-PCR gel QC and pooling volumes")
    q.add_argument("image")
    q.add_argument("--expected-bp", type=float, default=None)
    q.add_argument("--lanes", type=int, default=None)
    q.add_argument("--ladder-lane", type=int, default=None)
    q.add_argument("--ladder-sizes", default=None, help="comma separated, e.g. 10000,5000,3000,1000")
    q.add_argument("--tolerance", type=float, default=0.20)
    q.add_argument("--wells", default=None, help="comma separated labels for the sample lanes")
    q.add_argument("--secondary-frac", type=float, default=0.50,
                   help="a second band this fraction of the main one means do not pool")
    q.add_argument("--band-prominence", type=float, default=0.15,
                   help="peak prominence, as a fraction of the lane maximum")
    q.add_argument("--no-crop", action="store_true", help="skip auto-cropping to the gel")
    q.add_argument("-o", "--out", default="gel_qc")
    q.set_defaults(func=cmd_gel)

    q = sub.add_parser("run", help="demultiplex, genotype and write the pick plate")
    q.add_argument("fastq")
    q.add_argument("--set", default="cole1")
    q.add_argument("--samples", default=None, help="samples.csv from `apcs plan`")
    q.add_argument("--expect", nargs="+", required=True, help="acceptable reference(s)")
    q.add_argument("--also", nargs="*", help="decoys worth telling apart (parental, intermediate)")
    q.add_argument("--min-len", type=int, default=1000)
    q.add_argument("--min-depth", type=int, default=8)
    q.add_argument("--arm-edit", type=int, default=4, dest="arm_mismatch",
                   help="edit distance allowed when matching the primer arm")
    q.add_argument("--index-mismatch", type=int, default=1)
    q.add_argument("--index-shift", type=int, default=1,
                   help="indel slack when locating the index next to the arm")
    q.add_argument("--picks", type=int, default=2)
    q.add_argument("--backups", type=int, default=2)
    q.add_argument("--title", default=None)
    q.add_argument("-o", "--out", default="result")
    q.set_defaults(func=cmd_run)
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    try:
        a.func(a)
    except (ValueError, OSError, RuntimeError) as exc:
        raise SystemExit(f"apcs: {exc}") from exc


if __name__ == "__main__":
    main()
