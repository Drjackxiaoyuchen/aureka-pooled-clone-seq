import pytest
from apcs.primersets import load_config, check_template, index_distances, max_safe_error_budget
from apcs.demux import call_read, find_arm, call_index
from apcs.seqio import rc

cfg = load_config()
ps = cfg.sets["cole1"]
PLASMID = ("ACGTTGCAATCGGCTAAGCCTTAGGCATAC" * 120)


def amplicon(f="M544", r="M547", core=PLASMID):
    return ps.primer_f(cfg.fwd_index[f]) + core + rc(ps.primer_r(cfg.rev_index[r]))


def test_well_geometry():
    assert len(list(cfg.wells())) == 48
    assert cfg.well_of("M544", "M547") == "c3"


def test_round_trip_both_orientations():
    for seq in (amplicon(), rc(amplicon())):
        c = call_read(seq, ps, cfg)
        assert c.status == "assigned" and c.well == "c3"


def test_survives_eroded_handle():
    """The failure mode the shipped 5'-anchored pipeline cannot handle."""
    seq = amplicon()[len(ps.handle_f):]          # handle gone, index intact
    c = call_read(seq, ps, cfg)
    assert c.well == "c3" and c.handle_f_seen is False


def test_survives_one_substitution_in_the_index():
    ix = list(cfg.fwd_index["M544"])
    ix[2] = "A" if ix[2] != "A" else "T"
    seq = ps.handle_f + "".join(ix) + ps.arm_f + PLASMID + rc(ps.primer_r(cfg.rev_index["M547"]))
    assert call_read(seq, ps, cfg).well == "c3"


def test_survives_a_deletion_inside_the_arm():
    arm = ps.arm_f[:10] + ps.arm_f[11:]          # one base deleted
    seq = ps.handle_f + cfg.fwd_index["M544"] + arm + PLASMID + \
        rc(ps.primer_r(cfg.rev_index["M547"]))
    assert call_read(seq, ps, cfg).well == "c3"


def test_two_substitutions_in_the_index_are_not_guessed():
    ix = list(cfg.fwd_index["M544"])
    ix[0], ix[1] = "A", "A"
    seq = ps.handle_f + "".join(ix) + ps.arm_f + PLASMID + \
        rc(ps.primer_r(cfg.rev_index["M547"]))
    c = call_read(seq, ps, cfg)
    assert c.well != "c3"


def test_one_index_only():
    seq = ps.primer_f(cfg.fwd_index["M527"]) + PLASMID
    assert call_read(seq, ps, cfg).status == "one_index_only"


def test_index_set_distance_and_budget():
    for ix in (cfg.fwd_index, cfg.rev_index):
        assert index_distances(ix)["min_levenshtein"] == 3
        assert max_safe_error_budget(ix) == 1


def test_find_arm_and_call_index():
    assert find_arm("N" * 20 + ps.arm_f, ps.arm_f, 4) == 20
    assert find_arm("ACGT" * 30, ps.arm_f, 4) == -1
    assert call_index(cfg.fwd_index["M854"], cfg.fwd_index)[0] == "M854"


def test_shifted_index_ties_are_not_arbitrarily_assigned():
    from apcs.demux import _index_near
    assert _index_near('AACGTAC', 6, {'a': 'AACGTA', 'b': 'ACGTAC'}, 0, 1) is None


def test_empty_used_well_list_means_all_assignments_are_sentinels(tmp_path):
    from apcs.demux import demux
    from apcs.seqio import write_fastq
    seq = amplicon()
    path = tmp_path / 'pool.fastq'
    write_fastq(path, [('r', seq, 'I' * len(seq))])
    result = demux(path, ps, cfg, tmp_path / 'out', used_wells=[])
    assert result['run']['sentinel_reads'] == 1


def test_unsafe_index_error_budget_and_existing_output_rejected(tmp_path):
    from apcs.demux import demux
    with pytest.raises(ValueError):
        demux('missing', ps, cfg, tmp_path / 'out', index_mm=2)
    (tmp_path / 'old.fastq').write_text('old')
    with pytest.raises(ValueError):
        demux('missing', ps, cfg, tmp_path)
