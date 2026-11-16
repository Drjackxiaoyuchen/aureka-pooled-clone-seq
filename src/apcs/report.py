"""The output a person at the bench actually uses: a plate they can pick from."""
from __future__ import annotations
import os, json, csv, html, datetime
import base64
from importlib.resources import files

GLYPH = {
    "PICK": "**", "PICK (low depth)": "*-", "BACKUP": "+ ", "SPARE": "o ",
    "WRONG CONSTRUCT": "X ", "MIXED - DISCARD": "! ", "NO CALL": ". ", "NOT USED": "  ",
}
COLOR = {
    "PICK": "#1a7f37", "PICK (low depth)": "#4a9c5d", "BACKUP": "#9a6700", "SPARE": "#656d76",
    "WRONG CONSTRUCT": "#cf222e", "MIXED - DISCARD": "#a40e26", "NO CALL": "#8c959f",
    "NOT USED": "#d0d7de",
}


def plate_text(cfg, by_well, title=""):
    rows, cols = cfg.row_labels[:len(cfg.rows)], cfg.column_labels[:len(cfg.columns)]
    out = []
    if title:
        out.append(title)
    out.append("      " + "".join(f"{c:^11}" for c in cols))
    out.append("      " + "".join(f"{cfg.columns[i]:^11}" for i in range(len(cols))))
    out.append("     +" + "----------+" * len(cols))
    for ri, r in enumerate(rows):
        cells = []
        for ci, c in enumerate(cols):
            w = by_well.get(f"{r}{c}")
            if not w or not w.get("in_use", True):
                cells.append(f"{'':^10}|")
            else:
                g = GLYPH.get(w["verdict"], "? ")
                cells.append(f"{g}{w['well']:<3}{w['n_reads']:>4} |")
        out.append(f" {r} {cfg.rows[ri][-3:]}|" + "".join(cells))
        out.append("     +" + "----------+" * len(cols))
    out.append("")
    out.append("  ** PICK      *- PICK, low depth   +  BACKUP     o  spare (clean, not needed)")
    out.append("  X  wrong construct   !  mixed colony, discard   .  no call (too few reads)")
    out.append("  number after the well is mapped reads; see CSV for assigned and informative counts")
    return "\n".join(line.rstrip() for line in out)


def plate_html(cfg, by_well, meta):
    rows, cols = cfg.row_labels[:len(cfg.rows)], cfg.column_labels[:len(cfg.columns)]
    th = "".join(f"<th>{html.escape(c)}<br><span class=ix>{html.escape(cfg.columns[i])}</span></th>" for i, c in enumerate(cols))
    body = []
    for ri, r in enumerate(rows):
        tds = []
        for c in cols:
            w = by_well.get(f"{r}{c}")
            if not w or not w.get("in_use", True):
                count = w.get('assigned_reads', 0) if w else 0
                tds.append(f'<td class="off">{count} sentinel reads</td>')
                continue
            v = w["verdict"]
            tds.append(
                f'<td tabindex="0" role="button" data-well="{html.escape(w["well"], quote=True)}" '
                f'aria-label="{html.escape(w["well"] + ": " + v, quote=True)}" '
                f'style="border-left:6px solid {COLOR.get(v,"#999")}">'
                f'<div class=well>{html.escape(w["well"])}</div>'
                f'<div class=verdict style="color:{COLOR.get(v,"#999")}">{html.escape(v)}</div>'
                f'<div class=num>{w.get("n_informative_reads",0)} informative &middot; {w.get("purity",0)*100:.0f}% support</div>'
                f'<div class=ref>{html.escape(str(w.get("best_ref") or "-"))[:28]}</div></td>')
        body.append(f"<tr><th>{html.escape(r)}<br><span class=ix>{html.escape(cfg.rows[ri])}</span></th>" + "".join(tds) + "</tr>")
    picks = [w for w in by_well.values() if w.get("verdict", "").startswith("PICK")]
    picks.sort(key=lambda w: w.get("rank") or 99)
    headline = ", ".join(w["well"] for w in picks) or "nothing clears the gate"
    data = json.dumps(by_well).replace('<', '\\u003c')
    logo = base64.b64encode(files('apcs').joinpath('data/aureka-logo.png').read_bytes()).decode('ascii')
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(meta.get('title','Pick plate'))}</title>
<style>
 body{{font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
      margin:28px auto;padding:0 20px;color:#33266c;max-width:1100px;background:#f7f6fb}}
 .brand{{height:42px;max-width:220px;object-fit:contain;object-position:left;margin-bottom:24px}}
 h1{{font-size:20px;margin:0 0 2px}} .sub{{color:#656d76;margin-bottom:18px}}
 .headline{{background:#f0eafa;border:1px solid #e6e1ef;border-radius:8px;padding:12px 14px;
            margin:0 0 20px;font-size:16px}}
 .headline b{{color:#4200ac}}
 .scroll{{overflow-x:auto}} .legend{{margin:12px 0;color:#655b79;font-size:12px}}
 #detail{{white-space:pre-wrap;border-left:4px solid #6fc0a9;background:white;padding:16px;margin:18px 0}}
 table{{border-collapse:separate;border-spacing:4px}}
 th{{font-weight:600;color:#656d76;font-size:12px;text-align:center;padding:2px 6px}}
 .ix{{font-weight:400;font-size:10px;color:#8c959f}}
 td{{background:#fff;border:1px solid #d0d7de;border-radius:6px;padding:6px 8px;
     width:118px;height:52px;vertical-align:top}}
 td.off{{background:repeating-linear-gradient(45deg,#fafbfc,#fafbfc 6px,#f0f2f4 6px,#f0f2f4 12px);
         border-style:dashed;color:#77718b;font-size:10px}}
 td[data-well]{{cursor:pointer;min-width:120px}} td:focus{{outline:3px solid #4200ac}}
 .well{{font-weight:700;font-size:13px}} .verdict{{font-size:11px;font-weight:600}}
 .num,.ref{{font-size:10px;color:#656d76}} .ref{{font-family:ui-monospace,Menlo,monospace}}
 .meta{{margin-top:22px;font-size:12px;color:#656d76}}
 .meta code{{background:#f6f8fa;padding:1px 4px;border-radius:3px}}
 @media print{{body{{margin:8px}} .meta{{page-break-inside:avoid}}}}
</style></head><body>
<img class="brand" alt="Aureka" src="data:image/png;base64,{logo}">
<h1>{html.escape(meta.get('title','Pick plate'))}</h1>
<div class=sub>{html.escape(meta.get('subtitle',''))}</div>
<div class=headline>Pick <b>{html.escape(headline)}</b>. {html.escape(meta.get('headline_note',''))}</div>
<div class="scroll"><table aria-label="Clone screening plate"><tr><th></th>{th}</tr>{''.join(body)}</table></div>
<div class="legend">PICK / BACKUP / SPARE: reference screen passed or needs review · red: wrong or mixed evidence · gray: no call. Click a loaded well for details.</div>
<div id="detail" role="status" aria-live="polite">Select a loaded well to inspect its evidence.</div>
<p><a id="genbank-download" hidden>Annotated reference + read evidence (.gb)</a> <a id="amplicon-download" hidden>Predicted indexed amplicon (.gb)</a></p>
<p class="legend"><a href="../demux/read_assignments.tsv">Per-read assignment audit (.tsv)</a> · GenBank exports retain reference sequence; they are not read-derived consensuses.</p>
<div class=meta>{html.escape(meta.get('footer',''))}</div>
<script type="application/json" id="well-data">{data}</script>
<script>
const wells = JSON.parse(document.getElementById('well-data').textContent);
function showWell(el) {{
 const w = wells[el.dataset.well];
 document.getElementById('detail').textContent =
  w.well + ' · ' + w.verdict + '\\n' + (w.label || '') + '\\n' +
  'Reference: ' + (w.best_ref || 'unresolved') + '\\n' +
  'Assigned / mapped / informative molecules: ' + w.assigned_reads + ' / ' + w.n_reads + ' / ' + w.n_informative_reads + '\\n' +
  'Dominant vote support: ' + (w.purity * 100).toFixed(1) + '% · Alignment breadth: ' + w.breadth_pct + '%\\n' + w.note;
 for(const [id,key] of [['genbank-download','genbank_file'],['amplicon-download','amplicon_file']]) {{
  const a=document.getElementById(id); a.hidden=!w[key]; if(w[key]) a.href=w[key];
 }}
}}
document.querySelectorAll('[data-well]').forEach(el => {{
 el.addEventListener('click', () => showWell(el));
 el.addEventListener('keydown', e => {{if(e.key === 'Enter' || e.key === ' ') {{e.preventDefault(); showWell(el);}}}});
}});
</script>
</body></html>"""


def write_report(outdir, cfg, wells, run, meta):
    os.makedirs(outdir, exist_ok=True)
    by_well = {w["well"]: w for w in wells}
    open(os.path.join(outdir, "plate_layout.txt"), "w").write(
        plate_text(cfg, by_well, meta.get("title", "")) + "\n")
    open(os.path.join(outdir, "plate.html"), "w").write(plate_html(cfg, by_well, meta))

    fields = ["well", "label", "in_use", "verdict", "rank", "assigned_reads", "n_reads", "median_len",
              "n_informative_reads", "best_ref", "purity", "n_site_calls", "mean_identity",
              "breadth_pct", "mean_depth", "second_ref", "second_n", "unmapped", "n_ambiguous_reads", "genbank_file", "amplicon_file", "note"]
    with open(os.path.join(outdir, "per_well.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in wells:
            w.writerow(row)
    json.dump(run, open(os.path.join(outdir, "run_stats.json"), "w"), indent=1)

    picks = sorted([w for w in wells if w["verdict"] in ("PICK", "PICK (low depth)")],
                   key=lambda w: w.get("rank") or 99)
    backups = [w for w in wells if w["verdict"] == "BACKUP"]
    bad = [w for w in wells if w["verdict"] in ("WRONG CONSTRUCT", "MIXED - DISCARD")]
    nocall = [w for w in wells if w["verdict"] == "NO CALL" and w.get("in_use", True)]

    L = []
    L.append(f"# {meta.get('title','Pooled clone screen')}\n")
    L.append(f"*{meta.get('subtitle','')}*\n")
    L.append("## Pick this\n")
    if picks:
        L.append(f"**{', '.join(w['well'] for w in picks)}** — candidates for miniprep and independent confirmation. "
                 f"Backups for review: {', '.join(w['well'] for w in backups) or 'none'}.\n")
    else:
        L.append("**Nothing clears the gate.** See the table — then re-pick, or re-pool, "
                 "or fix the step the numbers point at.\n")
    L.append("```\n" + plate_text(cfg, by_well) + "\n```\n")
    L.append("## Every well\n")
    L.append("| Well | Verdict | Construct | Reads | Informative | Identity | Breadth | Note |")
    L.append("|---|---|---|---|---|---|---|---|")
    order = {"PICK": 0, "PICK (low depth)": 1, "BACKUP": 2, "SPARE": 3, "WRONG CONSTRUCT": 4,
             "MIXED - DISCARD": 5, "NO CALL": 6}
    for w in sorted([x for x in wells if x.get("in_use", True)],
                    key=lambda x: (order.get(x["verdict"], 9), -(x["n_reads"]))):
        L.append(f"| {w['well']} | {w['verdict']} | {w.get('best_ref') or '-'} | {w['n_reads']} | "
                 f"{w.get('n_informative_reads', 0)} | "
                 f"{(w.get('mean_identity') or 0)*100:.1f}% | {w.get('breadth_pct',0):.0f}% | "
                 f"{w.get('note','')} |")
    L.append("\n## Run\n")
    L.append(f"- primer set **{run['primer_set']}**, blind spot {meta.get('blind_spot','?')} bp")
    L.append(f"- {run['reads_total']:,} reads in, {run['reads_considered']:,} over "
             f"{run['min_len']} nt, **{run['reads_assigned']:,} assigned "
             f"({run['assignment_rate']*100:.0f}%)**")
    L.append(f"- {run['one_index_only']:,} reads carried only one index, "
             f"{run['unassigned']:,} carried neither")
    sent = run.get("sentinel_reads", 0)
    if run.get("sentinel_wells_with_reads") is not None:
        L.append(f"- **{sent} reads landed in unused index combinations** — a contamination/assignment signal; "
                 f"this does not measure all errors between loaded wells"
                 + (f" ({', '.join(f'{w}:{n}' for w, n in run['sentinel_wells_with_reads'][:6])})"
                    if sent else ""))
    if nocall:
        L.append(f"- {len(nocall)} wells could not be called: "
                 f"{', '.join(w['well'] for w in nocall)} "
                 "(reads present but not spanning a region that separates the references, "
                 "or below min-depth)")
    if bad:
        L.append(f"- do not use: {', '.join(w['well'] + ' (' + w['verdict'] + ')' for w in bad)}")
    L.append("\n## What this does not tell you\n")
    L.append(f"Template blind spots (bp): {meta.get('blind_spot','unknown')}. "
             "Primer-arm bases are supplied by the oligos and do not independently verify the original template. "
             "Dominant vote support is not a validated mixture fraction. No consensus or novel-variant analysis is performed. "
             "Confirm candidate isolates by an independent whole-plasmid sequence check.\n")
    L.append(f"\n<sub>generated {datetime.date.today().isoformat()} by apcs</sub>\n")
    open(os.path.join(outdir, "00_PICK_THIS.md"), "w").write("\n".join(L))
    return picks
