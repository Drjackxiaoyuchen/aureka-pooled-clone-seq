"""Reference-based GenBank exports. No sequence is inferred from reads here."""
from __future__ import annotations
from copy import deepcopy
from pathlib import Path
import json
import html
import math
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from .seqio import load_reference, _open, rc
from .primersets import check_template, find_sites


def read_record(path):
    name, seq = load_reference(path)
    suffix = str(path).lower().removesuffix('.gz')
    if suffix.endswith(('.gb', '.gbk', '.genbank', '.ape')):
        with _open(path) as fh:
            record = SeqIO.read(fh, 'genbank')
    else:
        record = SeqRecord(Seq(seq), id=name, name=name[:16], description='Input reference; no functional annotations supplied')
    if record.annotations.get('topology') == 'linear':
        raise ValueError('annotation requires a circular plasmid template, not an explicitly linear record')
    record.annotations.update(molecule_type='DNA', topology='circular')
    return record


def circular_location(start, length, n, strand=None):
    start %= n
    if length <= 0 or length > n:
        raise ValueError('invalid circular feature length')
    if start + length <= n:
        return SimpleLocation(start, start + length, strand=strand)
    parts = [SimpleLocation(start, n, strand=strand), SimpleLocation(0, start + length - n, strand=strand)]
    if strand == -1:
        parts.reverse()
    return CompoundLocation(parts)


def feature(location, label, note, kind='misc_feature', **qualifiers):
    return SeqFeature(location, type=kind, qualifiers={'label':[label], 'note':[note], **qualifiers})


def geometry(record, pset):
    seq = str(record.seq).upper()
    check = check_template(seq, pset)
    if not check['ok']:
        raise ValueError(check['reason'])
    fp, fs = check['fwd_sites'][0]
    rp, rs = check['rev_sites'][0]
    n = len(seq)
    fstart = fp - 1 if fs == '+' else (n - (fp - 1 + len(pset.arm_f))) % n
    span = check['template_span_bp']
    gap = check['blind_spot_bp']
    gap_start = (fstart + span) % n
    if fs == '-':
        gap_start = (n - gap_start - gap) % n
    return dict(check, fwd_start=fp-1, fwd_strand=1 if fs=='+' else -1,
                rev_start=rp-1, rev_strand=1 if rs=='+' else -1,
                fwd_arm_length=len(pset.arm_f), rev_arm_length=len(pset.arm_r),
                oriented_fstart=fstart, gap_start=gap_start)


def annotated_reference(record, pset, diagnostic_alleles=()):
    record = deepcopy(record)
    g = geometry(record, pset)
    n = len(record)
    record.description += '; Aureka annotated REFERENCE, not a read-derived consensus'
    record.annotations['comment'] = (str(record.annotations.get('comment','')) +
        '\nREFERENCE SEQUENCE. Added primer and gap annotations are predictions from this input map. '
        'Primer-arm bases originate from oligos and do not independently verify the source template.')
    for direction in ('fwd','rev'):
        arm = pset.arm_f if direction=='fwd' else pset.arm_r
        record.features.append(feature(circular_location(g[direction+'_start'],len(arm),n,g[direction+'_strand']),
            f'{pset.name}_{direction}_annealing_arm',
            'Exact template match. The indexed 5-prime tail is not present in the source template.',
            'primer_bind', sequence=[arm], ApEinfo_fwdcolor=['#4200AC' if direction=='fwd' else '#6FC0A9'],
            ApEinfo_revcolor=['#4200AC' if direction=='fwd' else '#6FC0A9']))
    if g['blind_spot_bp']:
        record.features.append(feature(circular_location(g['gap_start'],g['blind_spot_bp'],n),
            'OMITTED_FROM_AMPLICON', 'Predicted gap between primer arms; no amplicon-read evidence can verify this region.',
            ApEinfo_fwdcolor=['#CF6B4D']))
    for i, allele in enumerate(diagnostic_alleles,1):
        hits = find_sites(str(record.seq).upper(), allele)
        if len(hits)==1:
            pos,strand=hits[0]
            record.features.append(feature(circular_location(pos-1,len(allele),n,1 if strand=='+' else -1),
                f'diagnostic_window_{i}', 'Reference-panel discriminating window; not an observed consensus.',
                ApEinfo_fwdcolor=['#E2B13C']))
    return record,g


def predicted_amplicon(record, pset, cfg, well):
    pairs={w:(f,r) for w,f,r in cfg.wells()}
    well=well.lower()
    if well not in pairs:
        raise ValueError(f'unknown well: {well}')
    g=geometry(record,pset)
    f,r=pairs[well]
    fi,ri=cfg.fwd_index[f],cfg.rev_index[r]
    oriented=record if g['fwd_strand']==1 else record.reverse_complement(
        id=record.id,name=record.name,description=record.description,annotations=True,features=True)
    n,start,span=len(record),g['oriented_fstart'],g['template_span_bp']
    core=(str(oriented.seq)*2)[start:start+span]
    prefix=pset.handle_f+fi
    suffix=rc(ri)+rc(pset.handle_r)
    seq=prefix+core+suffix
    out=SeqRecord(Seq(seq),id=f'{well}_predicted_amp',name=f'{well}_amplicon',
        description=f'PREDICTED indexed amplicon from {record.id}; not read-derived consensus',
        annotations={'molecule_type':'DNA','topology':'linear',
                     'comment':f'Predicted from the supplied map and {pset.name} {f}/{r} oligos. No consensus inference.'})
    # Project inherited features through circular rotation and any reverse complement.
    for original in oriented.features:
        if original.type=='source':
            continue
        parts=[]
        for part in original.location.parts:
            for delta in (-n,0,n):
                lo=max(int(part.start)+delta,start)
                hi=min(int(part.end)+delta,start+span)
                if hi>lo:
                    parts.append(SimpleLocation(lo-start+len(prefix),hi-start+len(prefix),strand=part.strand))
        if not parts:
            continue
        copied=deepcopy(original)
        copied.location=parts[0] if len(parts)==1 else CompoundLocation(parts)
        copied.qualifiers.setdefault('note',[]).append('Projected from supplied reference onto predicted amplicon.')
        if sum(len(p) for p in parts)!=len(original.location):
            copied.qualifiers['note'].append('Partial: original feature intersects omitted sequence.')
            copied.qualifiers.pop('translation',None)
        out.features.append(copied)
    L=len(seq);p=len(prefix);e=p+span
    regions=[(0,len(pset.handle_f),1,'forward_handle'),
             (len(pset.handle_f),p,1,'forward_index_'+f),
             (p,p+len(pset.arm_f),1,'forward_arm_primer_derived'),
             (e-len(pset.arm_r),e,-1,'reverse_arm_primer_derived'),
             (e,e+6,-1,'reverse_index_'+r),
             (e+6,L,-1,'reverse_handle')]
    for lo,hi,strand,label in regions:
        out.features.append(feature(SimpleLocation(lo,hi,strand=strand),label,
            'Primer-derived sequence; predicted oligo architecture. Coordinates refer to this linear amplicon.'))
    out.features += [feature(SimpleLocation(0,len(pset.primer_f(fi)),strand=1),f'{pset.name}-F-{f}',
        'Full forward oligo, including non-template 5-prime tail.','primer_bind',sequence=[pset.primer_f(fi)]),
        feature(SimpleLocation(L-len(pset.primer_r(ri)),L,strand=-1),f'{pset.name}-R-{r}',
        'Full reverse oligo in synthesis orientation.','primer_bind',sequence=[pset.primer_r(ri)])]
    assert len(seq)==g['amplicon_bp'] and seq.startswith(pset.primer_f(fi)) and seq.endswith(rc(pset.primer_r(ri)))
    return out


def depth_runs(intervals,n):
    diff=[0]*(n+1)
    for start,end in intervals:
        span=min(n,end-start)
        if span<=0: continue
        start%=n
        for lo,hi in ((start,min(n,start+span)),(0,max(0,start+span-n))):
            if hi>lo: diff[lo]+=1;diff[hi]-=1
    depths=[];d=0
    for p in range(n): d+=diff[p];depths.append(d)
    runs=[];lo=0
    for p in range(1,n+1):
        if p==n or depths[p]!=depths[lo]: runs.append((lo,p,depths[lo]));lo=p
    return runs


def add_read_evidence(record,well,threshold):
    record=deepcopy(record)
    note=(f"Well {well['well']}; verdict {well['verdict']}; assigned {well['assigned_reads']}; "
          f"mapped {well['n_reads']}; informative {well['n_informative_reads']}; "
          f"dominant vote support {well['purity']:.3f}. REFERENCE-BASED, not read-derived consensus.")
    record.annotations['comment'] += '\n'+note
    record.features.append(feature(SimpleLocation(0,len(record)), 'screen_evidence_'+well['well'],note))
    if well['n_reads'] and well.get('best_ref') is None:
        record.annotations['comment'] += '\nNo reference was resolved: alignment-span coverage track unavailable.'
        return record,[]
    runs=depth_runs(well.get('alignment_intervals',[]),len(record))
    # Coverage classes keep the feature map compact; the exact span-depth track is separate.
    states=[]
    for lo,hi,d in runs:
        state='NO_ALIGNED_SPAN' if d==0 else ('LOW_SPAN_DEPTH' if d<threshold else 'ALIGNED_SPAN_DEPTH')
        if states and states[-1][3]==state:
            a,_,levels,s=states[-1];states[-1]=(a,hi,levels+[d],s)
        else: states.append((lo,hi,[d],state))
    for lo,hi,levels,label in states:
        record.features.append(feature(SimpleLocation(lo,hi),label,
            f'Alignment-span depth {min(levels)}–{max(levels)}. Includes alignment deletions; not a base-call pileup or sequence validation.'))
    return record,runs


def binding_svg(record,g):
    n=len(record);cx,cy,r=270,245,155
    def point(pos,radius=r):
        theta=2*math.pi*pos/n-math.pi/2
        return cx+radius*math.cos(theta),cy+radius*math.sin(theta)
    def arc(start,length,color,width):
        a,b=point(start),point(start+length)
        return f'<path d="M {a[0]:.2f} {a[1]:.2f} A {r} {r} 0 {int(length>n/2)} 1 {b[0]:.2f} {b[1]:.2f}" stroke="{color}" stroke-width="{width}" fill="none"/>'
    gap=g['blind_spot_bp']
    body=arc(g['gap_start']+gap,n-gap,'#6FC0A9',13)+arc(g['gap_start'],gap,'#CF6B4D',13) if gap else ''
    labels=[]
    for direction,color,y in [('fwd','#4200AC',150),('rev','#247C68',245)]:
        pos=g[direction+'_start'];length=g[direction+'_arm_length']
        x1,y1=point(pos+length/2,175)
        labels.append(f'<path d="M{x1:.1f},{y1:.1f} L490,{y}" stroke="{color}" fill="none"/><text x="515" y="{y-6}" fill="{color}" font-size="22" font-weight="bold">{direction.upper()} arm: {pos+1}–{pos+length} ({"+" if g[direction+"_strand"]==1 else "−"})</text>')
        end=pos+length/2+g[direction+'_strand']*45
        a,b=point(pos+length/2,175),point(end,175)
        labels.append(f'<path d="M{a[0]:.1f},{a[1]:.1f} L{b[0]:.1f},{b[1]:.1f}" stroke="{color}" stroke-width="4" marker-end="url(#arrow)"/>')
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="510" viewBox="0 0 1100 510">
<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0 0 L8 4 L0 8 Z" fill="#33266C"/></marker></defs>
<rect width="1100" height="510" rx="16" fill="#F7F6FB"/><g font-family="Arial,sans-serif" fill="#33266C">
<text x="35" y="40" font-size="25" font-weight="bold">Exact-template primer binding</text>
<circle cx="{cx}" cy="{cy}" r="{r}" fill="none" stroke="#E6E1EF" stroke-width="13"/>{body}{''.join(labels)}
<text x="{cx}" y="235" text-anchor="middle" font-size="22">{html.escape(record.id)}</text><text x="{cx}" y="269" text-anchor="middle" font-size="24" font-weight="bold">{n:,} bp</text>
<text x="515" y="335" font-size="20">Omitted gap: {gap} bp</text><text x="515" y="370" font-size="20">Predicted indexed product: {g['amplicon_bp']:,} bp</text>
<text x="35" y="455" font-size="17">Mint: template sequence included in product. Orange: omitted gap. Arrows: extension direction.</text>
<text x="35" y="485" font-size="16">Positions are 1-based on this exact map. Binding-arm sequence is also primer-derived; the map is not a consensus.</text></g></svg>'''


def export_one(path,pset,cfg,outdir,well='a1',diagnostic_alleles=()):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    original=read_record(path)
    record,g=annotated_reference(original,pset,diagnostic_alleles)
    SeqIO.write(record,outdir/'reference_annotated.gb','genbank')
    SeqIO.write(predicted_amplicon(original,pset,cfg,well),outdir/f'{well}_predicted_amplicon.gb','genbank')
    (outdir/'primer_binding.svg').write_text(binding_svg(record,g))
    (outdir/'primer_binding.json').write_text(json.dumps(g,indent=2))
    return record,g


def export_run(paths,panel,pset,cfg,wells,outdir,threshold=8):
    outdir=Path(outdir)
    refs=outdir/'references';well_dir=outdir/'wells'
    refs.mkdir(parents=True,exist_ok=True);well_dir.mkdir(parents=True,exist_ok=True)
    records={};originals={}
    for path,name in zip(paths,panel.names):
        original=read_record(path);originals[name]=original
        alleles=[c['seqs'][name] for c in panel.clusters if name in c['seqs']]
        try:
            record,g=annotated_reference(original,pset,alleles)
            (refs/f'{name}.binding.svg').write_text(binding_svg(record,g))
        except ValueError as exc:
            record=deepcopy(original)
            record.annotations['comment']='REFERENCE, not consensus. Primer overlays unavailable: '+str(exc)
        records[name]=record
        SeqIO.write(record,refs/f'{name}.annotated.gb','genbank')
    for well in wells:
        if not well.get('in_use',True): continue
        name=well.get('best_ref') or panel.primary
        record,runs=add_read_evidence(records[name],well,threshold)
        filename=f"{well['well']}.reference_evidence.gb"
        SeqIO.write(record,well_dir/filename,'genbank')
        well['genbank_file']='../annotations/wells/'+filename
        with (well_dir/f"{well['well']}.alignment_span_depth.bedgraph").open('w') as fh:
            fh.write('# Reference coordinates: 0-based half-open. Alignment-span depth, not base-call pileup.\n')
            if not runs: fh.write('# Track unavailable: no resolved reference.\n')
            for lo,hi,depth in runs: fh.write(f'{record.id}\t{lo}\t{hi}\t{depth}\n')
        try:
            amp=predicted_amplicon(originals[name],pset,cfg,well['well'])
            amp_name=f"{well['well']}.predicted_amplicon.gb"
            SeqIO.write(amp,well_dir/amp_name,'genbank')
            well['amplicon_file']='../annotations/wells/'+amp_name
        except ValueError:
            well['amplicon_file']=''
