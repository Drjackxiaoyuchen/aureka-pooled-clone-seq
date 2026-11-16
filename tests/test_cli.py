import csv
import pytest
from apcs.cli import main, _read_samples
from apcs.primersets import load_config


def test_block_plan_sentinels_exclude_actual_chosen_positions(tmp_path):
    main(['plan', '-n', '16', '-o', str(tmp_path)])
    rows = list(csv.DictReader((tmp_path / 'samples.csv').open()))
    assert len(rows) == 16 and len({r['fwd_primer_tube'] for r in rows}) == 4
    sentinel_line = next(s for s in (tmp_path / 'README.md').read_text().splitlines()
                         if 'Unused positions:' in s)
    unused = set(sentinel_line.split(': ', 1)[1].split(', '))
    assert len(unused) == 32 and not unused & {r['well'] for r in rows}


@pytest.mark.parametrize('args', [['-n', '0'], ['-n', '-1'], ['-n', '49'],
                                  ['-n', '2', '--wells', 'a1,a1'],
                                  ['-n', '2', '--wells', 'a1']])
def test_invalid_plan_rejected(tmp_path, args):
    with pytest.raises(SystemExit):
        main(['plan', *args, '-o', str(tmp_path)])


@pytest.mark.parametrize('body', ['well,sample\n', 'well\na1\na1\n',
                                  'well,fwd_index\na1,M545\n',
                                  'well,sample\n,unassigned_colony\n',
                                  'well,well\na1,a2\n'])
def test_empty_duplicate_or_inconsistent_samples_rejected(tmp_path, body):
    p = tmp_path / 'samples.csv'
    p.write_text(body)
    with pytest.raises(ValueError):
        _read_samples(p, load_config())


def test_missing_expected_reference_is_a_cli_error():
    with pytest.raises(SystemExit):
        main(['run', 'not_there.fastq'])
