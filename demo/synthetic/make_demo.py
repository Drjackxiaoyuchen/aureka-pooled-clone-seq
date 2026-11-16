#!/usr/bin/env python3
"""Build a synthetic pooled run so the whole flow can be tried without wet-lab data.

Produces, in the directory given as the first argument (default: this one):

  ref_wanted.fa     the construct you were trying to build
  ref_parent.fa     the template it came from - differs in one 180 bp block
  pool.fastq        a pooled, indexed, nanopore-flavoured amplicon library
  samples.csv       the plate assignment

The simulated plate deliberately contains every outcome the report has to
handle: clean wells, a well holding the parent, a genuinely mixed colony, a
well with almost no reads, and an empty well.
"""
from __future__ import annotations
import os, random, sys, csv, json

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))
from apcs.primersets import load_config          # noqa: E402
from apcs.seqio import rc                        # noqa: E402

SEED = 7
LEN = 5200
BLOCK_START, BLOCK_LEN = 1800, 180


def random_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def nanopore(seq, rng, sub=0.012, ins=0.006, dele=0.009):
    out = []
    for b in seq:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + ins:
            out.append(b); out.append(rng.choice("ACGT")); continue
        if r < dele + ins + sub:
            out.append(rng.choice([x for x in "ACGT" if x != b])); continue
        out.append(b)
    return "".join(out)


def main(outdir=None):
    outdir = outdir or os.path.dirname(os.path.abspath(__file__))
    os.makedirs(outdir, exist_ok=True)
    rng = random.Random(SEED)
    cfg = load_config()
    ps = cfg.sets["cole1"]

    gap = random_seq(ps.blind_spot_bp, rng)
    body = random_seq(LEN - len(ps.arm_r) - len(ps.arm_f) - len(gap), rng)
    backbone = rc(ps.arm_r) + gap + ps.arm_f + body
    wanted = backbone
    parent = (backbone[:BLOCK_START]
              + random_seq(BLOCK_LEN, rng)
              + backbone[BLOCK_START + BLOCK_LEN:])
    open(os.path.join(outdir, "ref_wanted.fa"), "w").write(">ref_wanted\n" + wanted + "\n")
    open(os.path.join(outdir, "ref_parent.fa"), "w").write(">ref_parent\n" + parent + "\n")

    wells = [(w,f,r) for w,f,r in cfg.wells() if cfg.column_labels.index(w[1:]) < 4][:16]
    plan = {}
    for i, (w, f, r) in enumerate(wells):
        if i < 9:
            plan[w] = ("wanted", 50)
        elif i < 12:
            plan[w] = ("parent", 45)
        elif i == 12:
            plan[w] = ("mixed", 60)
        elif i == 13:
            plan[w] = ("wanted", 4)          # too few reads to call
        elif i == 14:
            plan[w] = ("wanted", 0)          # no band on the gel, nothing pooled
        else:
            plan[w] = ("parent", 20)

    recs = []
    for w, f, r in wells:
        kind, n = plan[w]
        for k in range(n):
            if kind == "mixed":
                tmpl = wanted if k % 2 == 0 else parent
            else:
                tmpl = wanted if kind == "wanted" else parent
            # real premium-PCR libraries are mostly full-length amplicons with a
            # tail of truncated molecules
            core_start = len(ps.arm_r) + len(gap) + len(ps.arm_f)
            amp = ps.primer_f(cfg.fwd_index[f]) + tmpl[core_start:] + rc(ps.primer_r(cfg.rev_index[r]))
            if rng.random() >= 0.8:
                cut = int(len(amp) * rng.uniform(0.3, 0.8))
                amp = amp[:cut] if rng.random() < 0.5 else amp[cut:]
            if rng.random() < 0.25:                      # eroded 5' handle, as in real data
                amp = amp[rng.randrange(1, 18):]
            seq = nanopore(amp, rng)
            if rng.random() < 0.5:
                seq = rc(seq)
            recs.append((f"{w}_read{k}", seq))
    # junk reads that belong to nothing
    for k in range(40):
        recs.append((f"junk{k}", nanopore(random_seq(rng.randrange(400, 3000), rng), rng)))
    rng.shuffle(recs)

    with open(os.path.join(outdir, "pool.fastq"), "w") as fh:
        for name, seq in recs:
            fh.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")

    with open(os.path.join(outdir, "samples.csv"), "w") as fh:
        fh.write("well,sample,expected,fwd_index,rev_index\n")
        for w, f, r in wells:
            fh.write(f"{w},colony_{w},ref_wanted,{f},{r}\n")
    with open(os.path.join(outdir, "truth.csv"), "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(['well', 'source', 'input_molecules'])
        for well, (kind, count) in plan.items():
            writer.writerow([well, kind, count])
    with open(os.path.join(outdir, "simulation.json"), "w") as fh:
        json.dump({'synthetic': True, 'seed': SEED, 'reference_bp': LEN,
                   'blind_spot_bp': len(gap), 'substitution_rate': 0.012,
                   'insertion_rate': 0.006, 'deletion_rate': 0.009,
                   'note': 'Independent toy errors; not an empirical nanopore error model.'}, fh, indent=2)
    print(f"{len(recs)} reads -> {outdir}/pool.fastq")
    print(f"plate: {sum(1 for v in plan.values() if v[0]=='wanted')} wanted, "
          f"{sum(1 for v in plan.values() if v[0]=='parent')} parent, 1 mixed")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else None)
