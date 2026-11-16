from pathlib import Path
from copy import deepcopy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,SimpleLocation
import pytest
from apcs.annotate import annotated_reference,predicted_amplicon,add_read_evidence,depth_runs
from apcs.primersets import load_config
from apcs.seqio import rc

cfg=load_config();ps=cfg.sets['cole1']


def test_public_puc19_primer_coordinates_and_oligo_orientation():
    record=SeqIO.read(Path(__file__).parents[1]/'demo/pUC19/pUC19_M77789.2.gb','genbank')
    annotated,g=annotated_reference(record,ps)
    assert len(annotated.features)==len(record.features)+3
    assert g['fwd_start']==2162 and g['fwd_strand']==-1
    assert g['rev_start']==2458 and g['rev_strand']==1
    assert g['amplicon_bp']==2457 and g['blind_spot_bp']==274
    f,r=annotated.features[-3:-1]
    assert str(f.extract(record.seq))==ps.arm_f
    assert str(r.extract(record.seq))==ps.arm_r
    amp=predicted_amplicon(record,ps,cfg,'c4')
    assert str(amp.seq).startswith(ps.primer_f(cfg.fwd_index['M544']))
    assert str(amp.seq).endswith(rc(ps.primer_r(cfg.rev_index['M881'])))
    for feat in amp.features:
        if feat.type=='primer_bind': assert str(feat.extract(amp.seq))==feat.qualifiers['sequence'][0]


def test_genbank_export_roundtrip_retains_input_features(tmp_path):
    record=SeqIO.read(Path(__file__).parents[1]/'demo/pUC19/pUC19_M77789.2.gb','genbank')
    annotated,_=annotated_reference(record,ps)
    out=tmp_path/'ref.gb';SeqIO.write(annotated,out,'genbank')
    reread=SeqIO.read(out,'genbank')
    assert reread.seq==record.seq
    for a,b in zip(record.features,reread.features):
        assert a.location==b.location and a.qualifiers==b.qualifiers and a.type==b.type


def test_predicted_amplicon_feature_projection_on_both_template_orientations():
    seq=rc(ps.arm_r)+'A'*274+ps.arm_f+'ACGT'*700
    record=SeqRecord(Seq(seq),id='test',annotations={'molecule_type':'DNA','topology':'circular'})
    record.features=[SeqFeature(SimpleLocation(900,1000,strand=1),type='CDS',qualifiers={'label':['demo_gene']})]
    reverse=record.reverse_complement(id=True,name=True,description=True,annotations=True,features=True)
    for source in [record,reverse]:
        amp=predicted_amplicon(source,ps,cfg,'a1')
        projected=next(f for f in amp.features if f.type=='CDS')
        assert projected.extract(amp.seq)==record.features[0].extract(record.seq)


def test_alignment_span_depth_wraps_without_reference_fill():
    assert depth_runs([(90,110),(0,20)],100)==[(0,10,2),(10,20,1),(20,90,0),(90,100,1)]


def test_unresolved_reference_has_unavailable_coverage_not_false_zero():
    r=SeqRecord(Seq('ACGT'*100),id='ref',annotations={'comment':'reference'})
    well=dict(well='a1',verdict='NO CALL',assigned_reads=20,n_reads=20,
              n_informative_reads=0,purity=0,best_ref=None,alignment_intervals=[])
    record,runs=add_read_evidence(r,well,8)
    assert not runs and 'unavailable' in record.annotations['comment']
    assert all(f.qualifiers.get('label')!=['NO_ALIGNED_SPAN'] for f in record.features)
