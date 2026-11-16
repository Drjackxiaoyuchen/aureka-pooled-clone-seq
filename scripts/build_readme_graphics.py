"""Build README figures from the committed synthetic report and GenBank output."""
from pathlib import Path
import csv
import html
import json
from Bio import SeqIO

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'assets/graphics'
PURPLE, MINT, INK = '#4200AC', '#6FC0A9', '#33266C'


def frame(title, subtitle, body, height):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="{height}" '
            f'viewBox="0 0 1200 {height}"><rect width="1200" height="{height}" rx="18" fill="#F7F6FB"/>'
            '<g font-family="Arial,sans-serif" fill="#33266C">'
            f'<text x="35" y="43" font-size="27" font-weight="bold">{html.escape(title)}</text>'
            f'<text x="35" y="77" font-size="17">{html.escape(subtitle)}</text>{body}</g></svg>')


def workflow():
    steps = [
        ('Check the map', ['Exact circular template', 'One site per arm', 'Keep colony identity']),
        ('Index each well', ['One row F + column R', 'Individual colony PCR', 'Tails label the product']),
        ('Review and pool', ['Expected single band', 'Quantify representation', 'Pool passing products']),
        ('Sequence the pool', ['Preserve both ends', 'Obtain raw FASTQ', 'Retain the sample map']),
        ('Separate reads', ['Match the primer set', 'Require both indexes', 'Audit every read']),
        ('Inspect results', ['Per-well FASTQ', 'Reference-based maps', 'Confirm an isolate']),
    ]
    parts = ['<defs><marker id="arrow" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto"><path d="M0 0 L7 3.5 L0 7 Z" fill="#6FC0A9"/></marker></defs>']
    for i, (title, lines) in enumerate(steps):
        x = 35 + i * 190
        parts.append(f'<rect x="{x}" y="110" width="175" height="180" rx="12" fill="white" stroke="#E6E1EF"/>'
                     f'<circle cx="{x+25}" cy="139" r="15" fill="{PURPLE}"/>'
                     f'<text x="{x+25}" y="145" text-anchor="middle" fill="white" font-size="17">{i+1}</text>')
        title_lines = ['Sequence', 'the pool'] if title == 'Sequence the pool' else [title]
        for j, line in enumerate(title_lines):
            y = 174 + j * 22 if len(title_lines) > 1 else 182
            parts.append(f'<text x="{x+12}" y="{y}" font-size="18" font-weight="bold">{line}</text>')
        for j, line in enumerate(lines):
            parts.append(f'<text x="{x+12}" y="{218+j*25}" font-size="13">{html.escape(line)}</text>')
        if i < 5:
            parts.append(f'<path d="M{x+176} 200 H{x+187}" stroke="{MINT}" stroke-width="3" marker-end="url(#arrow)"/>')
    parts.append('<text x="35" y="329" font-size="17">Preferred plan: A1–A4, B1–B4, … H1–H4. Read separation needs no reference; construct screening does.</text>')
    return frame('From individual colonies to traceable per-well results',
                 'Pool products after individual PCR and QC. Each well keeps its own forward/reverse index pair.', ''.join(parts), 360)


def gb_location(location):
    parts = [f'{int(p.start)+1}..{int(p.end)}' for p in location.parts]
    text = parts[0] if len(parts) == 1 else 'join(' + ','.join(parts) + ')'
    return f'complement({text})' if location.strand == -1 else text


def output_preview():
    report = ROOT / 'demo/synthetic/result/report'
    with (report / 'per_well.csv').open() as fh:
        wells = {r['well']: r for r in csv.DictReader(fh)}
    stats = json.loads((report / 'run_stats.json').read_text())
    colors = {'PICK': '#247C68', 'PICK (low depth)': '#247C68', 'BACKUP': '#946800', 'SPARE': '#655B79',
              'WRONG CONSTRUCT': '#BE363A', 'MIXED - DISCARD': '#A02C49', 'NO CALL': '#77718B'}
    parts = []
    for ci in range(4):
        parts.append(f'<text x="{180+ci*275}" y="115" text-anchor="middle" font-size="19" font-weight="bold">Column {ci+1}</text>')
    for ri, row in enumerate('abcd'):
        y = 133 + ri * 92
        parts.append(f'<text x="24" y="{y+44}" font-size="20" font-weight="bold">{row.upper()}</text>')
        for ci in range(4):
            w = wells[f'{row}{ci+1}']; x = 50 + ci * 275; color = colors[w['verdict']]
            parts.append(f'<rect x="{x}" y="{y}" width="260" height="80" rx="9" fill="white" stroke="#E6E1EF"/>'
                         f'<rect x="{x}" y="{y}" width="5" height="80" rx="2" fill="{color}"/>'
                         f'<text x="{x+15}" y="{y+25}" font-size="18" font-weight="bold">{w["well"].upper()}</text>'
                         f'<text x="{x+65}" y="{y+25}" font-size="15" fill="{color}">{html.escape(w["verdict"])}</text>'
                         f'<text x="{x+15}" y="{y+51}" font-size="14">{w["assigned_reads"]} assigned · {w["n_informative_reads"]} informative</text>'
                         f'<text x="{x+15}" y="{y+71}" font-size="13">{html.escape(w["best_ref"] or "unresolved")}</text>')
    record = SeqIO.read(ROOT / 'demo/synthetic/result/annotations/wells/a1.reference_evidence.gb', 'genbank')
    features = [f for f in record.features if f.type == 'primer_bind' or f.qualifiers.get('label') == ['OMITTED_FROM_AMPLICON']]
    parts.append('<rect x="35" y="523" width="1130" height="214" rx="12" fill="white" stroke="#E6E1EF"/>'
                 '<text x="55" y="558" font-size="20" font-weight="bold">A1 example files</text>')
    for y, text in [(593, 'a1.fastq → observed, assigned reads'),
                    (625, 'a1.reference_evidence.gb → reference + annotations'),
                    (657, 'a1.predicted_amplicon.gb → expected indexed product')]:
        parts.append(f'<text x="55" y="{y}" font-size="16">{html.escape(text)}</text>')
    parts.append('<text x="665" y="558" font-size="18" font-weight="bold">Features on the 5,200-bp synthetic reference</text>')
    for i, f in enumerate(features):
        text = f'{f.type}  {gb_location(f.location)}'
        parts.append(f'<text x="665" y="{590+i*31}" font-size="15" font-family="monospace">{html.escape(text)}</text>')
    parts.append('<text x="55" y="709" font-size="16" fill="#946800">GenBank files retain reference/predicted bases. They are not read-derived consensuses.</text>')
    return frame('What the generated outputs look like',
                 f'Seed-7 simulation: {stats["reads_total"]} input reads, {stats["reads_assigned"]} assigned. This is not a wet-lab benchmark.',
                 ''.join(parts), 765)


if __name__ == '__main__':
    OUT.mkdir(parents=True, exist_ok=True)
    for name, content in [('workflow_overview', workflow()), ('demultiplex_output_preview', output_preview())]:
        path = OUT / f'{name}.svg'
        path.write_text(content)
        print(path)
