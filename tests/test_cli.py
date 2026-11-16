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


def test_preferred_four_column_plan_and_uppercase_inputs(tmp_path):
    main(['plan','-n','30','-o',str(tmp_path)])
    rows=list(csv.DictReader((tmp_path/'samples.csv').open()))
    assert [r['well'] for r in rows[:8]]==['a1','a2','a3','a4','b1','b2','b3','b4']
    assert [r['well'] for r in rows[-2:]]==['h1','h2']
    p=tmp_path/'upper.csv';p.write_text('well,sample\nA1,test\n')
    assert _read_samples(p,load_config())[0]['well']=='a1'


def test_standalone_demux_needs_no_reference_and_records_counts(tmp_path):
    from apcs.seqio import rc,write_fastq
    import json
    cfg=load_config();ps=cfg.sets['cole1']
    seq=ps.primer_f(cfg.fwd_index['M544'])+'ACGT'*400+rc(ps.primer_r(cfg.rev_index['M547']))
    path=tmp_path/'pool.fastq';write_fastq(path,[('r1',seq,'I'*len(seq))])
    samples=tmp_path/'samples.csv';samples.write_text('well,sample\nC3,colony\n')
    out=tmp_path/'out'
    main(['demux',str(path),'--samples',str(samples),'-o',str(out)])
    assert (out/'demux/c3.fastq').exists()
    stats=json.loads((out/'demux_stats.json').read_text())
    assert stats['reads_assigned']==1 and stats['sentinel_reads']==0 and stats['input_sha256']
