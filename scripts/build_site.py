"""Bundle the static landing page and generated synthetic reports."""
import argparse
from pathlib import Path
import shutil

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--demo', required=True, type=Path)
parser.add_argument('--out', required=True, type=Path)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=False)
shutil.copy(root / 'site/index.html', args.out / 'index.html')
shutil.copytree(root / 'assets', args.out / 'assets')
shutil.copytree(args.demo / 'sequencing/result/report', args.out / 'demo')
for name in ('samples.csv', 'truth.csv', 'simulation.json', 'ref_wanted.fa', 'ref_parent.fa'):
    shutil.copy(args.demo / 'sequencing' / name, args.out / 'demo' / name)
shutil.copytree(args.demo / 'gel/result', args.out / 'gel')
shutil.copy(args.demo / 'gel/synthetic_gel.png', args.out / 'gel/synthetic_gel.png')
(args.out / '.nojekyll').touch()
print(args.out.resolve())
