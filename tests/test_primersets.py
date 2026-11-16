from apcs.primersets import load_config, check_template
from apcs.seqio import rc

cfg = load_config()


def circular_with(pset, gap=200, rest=4000):
    """A plasmid carrying one outward-facing site for each arm."""
    filler = "ACGTTGCAATCGGCTAAGCC"
    body = (filler * ((rest // len(filler)) + 1))[:rest]
    gapseq = (filler * ((gap // len(filler)) + 1))[:gap]
    return rc(pset.arm_r) + gapseq + pset.arm_f + body


def test_outward_pair_geometry():
    for name, ps in cfg.sets.items():
        seq = circular_with(ps)
        r = check_template(seq, ps)
        assert r["ok"], (name, r["reason"])
        assert r["blind_spot_bp"] == 200
        assert r['amplicon_bp'] == len(seq) - 200 + 45


def test_duplicate_site_is_refused():
    ps = cfg.sets["cole1"]
    seq = circular_with(ps) + ps.arm_f
    r = check_template(seq, ps)
    assert not r["ok"] and "one site per arm" in r["reason"]


def test_primer_lengths():
    ps = cfg.sets["ampr"]
    p = ps.primer_f(cfg.fwd_index["M527"])
    assert p.startswith(ps.handle_f) and p.endswith(ps.arm_f) and len(p) == 47
