"""Turn per-well numbers into one word the person at the bench can act on."""
from __future__ import annotations

def _short(name, n=34):
    if not name:
        return "-"
    return name if len(name) <= n else name[:n - 1] + "\u2026"


ORDER = ["PICK", "PICK (low depth)", "BACKUP", "WRONG CONSTRUCT", "MIXED - DISCARD", "NO CALL"]

DEFAULTS = dict(
    min_depth=8,            # reads assigned to the well, after the length filter
    comfortable_depth=15,   # above this, depth is not the limiting factor
    min_purity=0.90,        # fraction of reads agreeing on one reference
    mixed_purity=0.80,      # below this, with real support for a second reference
    mixed_min_second=3,
    min_identity=0.90,      # mean read identity to the chosen reference
    min_breadth=85.0,       # % of the reference the well's reads span
    min_mapped_fraction=0.80,
    min_informative_fraction=0.50,
    n_picks=2,              # how many to call PICK
    n_backups=2,
)


def call_well(g, expected_ref, cfg=None):
    c = dict(DEFAULTS); c.update(cfg or {})
    n = g["n_reads"]
    if n == 0:
        return "NO CALL", "no reads survived the length filter"
    if n < c["min_depth"]:
        return "NO CALL", f"only {n} reads (need {c['min_depth']})"
    if g.get("best_ref") is None:
        return ("NO CALL",
                f"{n} reads map, but none span a region that tells the references apart"
                " - the well may be fine, the reads are just too short or too few here")
    if (g["purity"] < c["mixed_purity"]
            and g.get("second_n", 0) >= c["mixed_min_second"]):
        votes = g.get("read_votes") or {}
        both = len([v for v in votes.values() if v >= c["mixed_min_second"]]) > 1
        return ("MIXED - DISCARD",
                f"{g['purity']:.0%} of {g.get('n_site_calls', n)} informative molecule votes say "
                f"{_short(g['best_ref'])}, the rest say {_short(g['second_ref'])}"
                + (" and whole reads split the same way" if both else "")
                + " - mixed reference evidence; re-streak and confirm")
    informative = g.get("n_informative_reads", n)
    if informative < c["min_depth"]:
        return "NO CALL", f"only {informative} informative molecules (need {c['min_depth']})"
    total = n + g.get("unmapped", 0)
    if n / max(1, total) < c["min_mapped_fraction"]:
        return "NO CALL", "too many assigned molecules fail to map; check unrepresented constructs"
    if informative / max(1, n) < c["min_informative_fraction"]:
        return "NO CALL", "fewer than half of mapped molecules discriminate the reference panel"
    if g.get("has_alternatives") is False:
        return "BACKUP", "reference match only; no alternative construct was tested"
    expected = {expected_ref} if isinstance(expected_ref, str) else set(expected_ref or [])
    if expected and g["best_ref"] not in expected:
        return ("WRONG CONSTRUCT",
                f"{g['purity']:.0%} of informative molecule votes say {_short(g['best_ref'])}; "
                f"expected {', '.join(sorted(expected))}")
    if g["purity"] < c["min_purity"]:
        return "BACKUP", f"purity {g['purity']:.0%} - re-streak before trusting it"
    if g["mean_identity"] < c["min_identity"]:
        return "BACKUP", f"mean read identity {g['mean_identity']:.1%} is low for this depth"
    if g["breadth_pct"] < c["min_breadth"]:
        return "BACKUP", f"reads span only {g['breadth_pct']:.0f}% of the map"
    if informative < c["comfortable_depth"]:
        return "PICK (low depth)", f"reference screen passed with {informative} informative molecules"
    return "PICK", "reference screen passed; confirm the isolated plasmid"


def rank_and_label(wells, cfg=None):
    """wells: list of dicts with 'well', 'verdict', genotype fields. Adds 'rank' and
    promotes the top clean wells to PICK / BACKUP so the bench has a short list."""
    c = dict(DEFAULTS); c.update(cfg or {})
    clean = [w for w in wells if w["verdict"].startswith("PICK")]
    clean.sort(key=lambda w: (-w.get("n_informative_reads", w["n_reads"]),
                              -w["mean_identity"], -w["breadth_pct"], w["well"]))
    for i, w in enumerate(clean):
        w["rank"] = i + 1
        if i >= c["n_picks"] + c["n_backups"]:
            w["verdict"] = "SPARE"
        elif i >= c["n_picks"]:
            w["verdict"] = "BACKUP"
    for w in wells:
        w.setdefault("rank", None)
    return wells
