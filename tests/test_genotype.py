import random
import pytest
from apcs.genotype import RefPanel, genotype_well, _breadth
from apcs.seqio import rc
from apcs.verdict import call_well


def dna(n, seed):
    r = random.Random(seed)
    return ''.join(r.choice('ACGT') for _ in range(n))


def panel(tmp_path, seqs):
    paths = []
    for name, seq in seqs.items():
        p = tmp_path / f'{name}.fa'
        p.write_text(f'>{name}\n{seq}\n')
        paths.append(p)
    return RefPanel(paths)


def reads(seqs):
    return [(f'r{i}', s, 'I' * len(s)) for i, s in enumerate(seqs)]


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('variant', ['block', 'insert', 'delete'])
def test_related_references_and_indels(tmp_path, reverse, variant):
    want = dna(5000, 1)
    if variant == 'block':
        alt = want[:1800] + dna(120, 2) + want[1920:]
    elif variant == 'insert':
        alt = want[:1800] + dna(25, 2) + want[1800:]
    else:
        alt = want[:1800] + want[1825:]
    stored = rc(alt) if reverse else alt
    p = panel(tmp_path, {'want': want, 'parent': stored})
    g = genotype_well(reads([alt, rc(alt)] * 10), p)
    assert g['best_ref'] == 'parent' and g['purity'] == 1
    assert g['n_informative_reads'] == 20
    assert call_well(g, 'want')[0] == 'WRONG CONSTRUCT'


def test_identical_or_rotated_references_have_no_discriminating_evidence(tmp_path):
    s = dna(5000, 5)
    p = panel(tmp_path, {'want': s, 'same': s[1200:] + s[:1200]})
    g = genotype_well(reads([s] * 20), p)
    assert g['best_ref'] is None
    assert call_well(g, 'want')[0] == 'NO CALL'


def test_partial_reads_missing_diagnostics_do_not_vote(tmp_path):
    s = dna(5000, 5)
    alt = s[:2500] + dna(120, 6) + s[2620:]
    p = panel(tmp_path, {'want': s, 'parent': alt})
    g = genotype_well(reads([s[:1400]] * 20), p)
    assert g['n_informative_reads'] == 0
    assert call_well(g, 'want')[0] == 'NO CALL'


def test_independent_molecules_limit_mixture_support(tmp_path):
    s = dna(5000, 5)
    alt = s[:1200] + dna(50, 8) + s[1250:3000] + dna(50, 9) + s[3050:]
    p = panel(tmp_path, {'want': s, 'parent': alt})
    g = genotype_well(reads([s] * 10 + [alt] * 10), p)
    assert g['n_site_calls'] == 20 and g['second_n'] == 10
    assert call_well(g, 'want')[0] == 'MIXED - DISCARD'


def test_unrelated_mixture_and_winning_reference_breadth(tmp_path):
    a, b = dna(5000, 5), dna(3500, 6)
    p = panel(tmp_path, {'want': a, 'unrelated': b})
    g = genotype_well(reads([a] * 10 + [b] * 10), p)
    assert g['purity'] == 0.5 and g['breadth_pct'] <= 100
    assert call_well(g, 'want')[0] == 'MIXED - DISCARD'


def test_related_plus_unrelated_mixture_is_not_hidden(tmp_path):
    a, b = dna(5000, 5), dna(3500, 6)
    parent = a[:1800] + dna(100, 8) + a[1900:]
    p = panel(tmp_path, {'want': a, 'parent': parent, 'unrelated': b})
    g = genotype_well(reads([a] * 10 + [b] * 10), p)
    assert g['purity'] == 0.5
    assert call_well(g, 'want')[0] == 'MIXED - DISCARD'


def test_shared_alleles_in_three_reference_panel_remain_ambiguous(tmp_path):
    a = dna(5000, 5)
    b = a[:1000] + dna(50, 8) + a[1050:]
    c = a[:3000] + dna(50, 9) + a[3050:]
    p = panel(tmp_path, {'want': a, 'b': b, 'c': c})
    g = genotype_well(reads([a[:1800]] * 20), p)
    assert g['best_ref'] is None


def test_circular_diagnostic_window_crosses_origin(tmp_path):
    a = dna(5000, 5)
    b = dna(25, 8) + a[25:4975] + dna(25, 9)
    p = panel(tmp_path, {'want': a, 'parent': b})
    seq = b[4000:] + b[:4000]
    assert genotype_well(reads([seq] * 20), p)['best_ref'] == 'parent'


def test_coverage_does_not_confuse_zero_length_and_full_circle():
    assert _breadth([(20, 20)], 100) == 0
    assert _breadth([(90, 110), (0, 30)], 100) == 40
    assert _breadth([(20, 120)], 100) == 100
