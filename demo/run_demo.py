"""Run and verify both synthetic examples from a fresh checkout."""
from pathlib import Path
import argparse
import csv
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def run(*args):
    subprocess.run([sys.executable, *map(str, args)], check=True)


def verify(report):
    with (report / 'per_well.csv').open() as fh:
        rows = {r['well']: r for r in csv.DictReader(fh)}
    good = ['a1','a2','a3','a4','a5','a6','b1','b2','b3']
    for well in good:
        assert rows[well]['verdict'] in {'PICK', 'PICK (low depth)', 'BACKUP', 'SPARE'}, rows[well]
        assert rows[well]['best_ref'] == 'ref_wanted', rows[well]
    for well in ['b4','b5','b6','c4']:
        assert rows[well]['verdict'] == 'WRONG CONSTRUCT', rows[well]
    assert rows['c1']['verdict'] == 'MIXED - DISCARD', rows['c1']
    for well in ['c2','c3']:
        assert rows[well]['verdict'] == 'NO CALL', rows[well]
    assert all(r['verdict'] == 'NOT USED' for r in rows.values() if r['in_use'] == 'False')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', default='demo-output', type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists():
        parser.error('choose a fresh output directory')
    seq = out / 'sequencing'
    gel = out / 'gel'
    run(ROOT / 'demo/synthetic/make_demo.py', seq)
    run('-m', 'apcs.cli', 'run', seq / 'pool.fastq', '--samples', seq / 'samples.csv',
        '--expect', seq / 'ref_wanted.fa', '--also', seq / 'ref_parent.fa',
        '--title', 'Aureka APCS · synthetic 16-well demonstration', '-o', seq / 'result')
    verify(seq / 'result/report')
    run(ROOT / 'demo/gel/make_demo.py', gel)
    run('-m', 'apcs.cli', 'gel', gel / 'synthetic_gel.png', '--expected-bp', '5000',
        '--ladder-lane', '2', '--ladder-sizes', '10000,5000,2000,500',
        '--wells', 'a1,a2,a3,a4,a5,a6,a7', '--no-crop', '-o', gel / 'result')
    result = json.loads((gel / 'result/gel.json').read_text())
    assert [r['verdict'] for r in result['lanes']] == [
        'OK', 'OK', '-', 'WRONG SIZE', 'MULTIPLE BANDS', 'NO BAND',
        'REVIEW CLIPPING', 'REVIEW VOLUME']
    print(f'Both synthetic examples verified. Open {seq / "result/report/plate.html"}')


if __name__ == '__main__':
    main()
