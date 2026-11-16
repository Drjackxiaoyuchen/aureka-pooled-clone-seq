from apcs.verdict import call_well, rank_and_label, DEFAULTS


def g(**kw):
    base = dict(n_reads=20, best_ref="want", purity=1.0, mean_identity=0.97, breadth_pct=99.0,
                second_ref=None, second_n=0, n_site_calls=20, n_informative_reads=20)
    base.update(kw)
    return base


def test_clean_well_is_a_pick():
    assert call_well(g(), "want")[0] == "PICK"


def test_low_depth_is_flagged_not_hidden():
    assert call_well(g(n_reads=10, n_informative_reads=10), "want")[0] == "PICK (low depth)"
    assert call_well(g(n_reads=3), "want")[0] == "NO CALL"


def test_wrong_construct():
    v, why = call_well(g(best_ref="parent"), "want")
    assert v == "WRONG CONSTRUCT" and "parent" in why


def test_mixed_needs_real_support_for_the_second_reference():
    assert call_well(g(purity=0.6, second_ref="parent", second_n=8), "want")[0] == "MIXED - DISCARD"
    # a single stray vote is not a mixed colony
    assert call_well(g(purity=0.6, second_ref="parent", second_n=1), "want")[0] != "MIXED - DISCARD"


def test_no_diagnostic_coverage_is_no_call_not_a_wrong_answer():
    assert call_well(g(best_ref=None), "want")[0] == "NO CALL"


def test_ranking_promotes_the_deepest_clean_wells():
    wells = [dict(well=f"a{i}", verdict="PICK", n_reads=10 * i, mean_identity=0.98,
                  breadth_pct=99.0) for i in range(1, 6)]
    rank_and_label(wells)
    by = {w["well"]: w["verdict"] for w in wells}
    assert by["a5"] == "PICK" and by["a4"] == "PICK"
    assert by["a3"] == "BACKUP" and by["a1"] == "SPARE"


def test_dissimilar_reference_mixture_is_checked_too():
    assert call_well(g(purity=0.5, second_n=10, n_site_calls=0), "want")[0] == "MIXED - DISCARD"


def test_depth_means_informative_molecules():
    assert call_well(g(n_reads=100, n_informative_reads=2), "want")[0] == "NO CALL"


def test_single_reference_does_not_establish_construct_purity():
    assert call_well(g(has_alternatives=False), "want")[0] == "BACKUP"


def test_multiple_expected_constructs_are_accepted():
    assert call_well(g(best_ref="alternate"), {"want", "alternate"})[0] == "PICK"


def test_unrepresented_or_ambiguous_majority_blocks_a_pick():
    assert call_well(g(unmapped=100), "want")[0] == "NO CALL"
    assert call_well(g(n_reads=100, n_informative_reads=20), "want")[0] == "NO CALL"
