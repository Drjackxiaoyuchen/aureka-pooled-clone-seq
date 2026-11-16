"""Generate editable SVG teaching figures; all cartoons are explanatory schematics."""
from pathlib import Path

OUT=Path(__file__).resolve().parents[1]/'assets/graphics'
OUT.mkdir(parents=True,exist_ok=True)
P='#4200AC';M='#6FC0A9';I='#33266C';G='#F7F6FB';O='#BF583B'


def frame(title,body,height=620):
    return f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="{height}" viewBox="0 0 1200 {height}"><rect width="1200" height="{height}" rx="18" fill="{G}"/><g font-family="Arial,sans-serif" fill="{I}"><text x="36" y="48" font-size="27" font-weight="bold">{title}</text>{body}</g></svg>'


body=[]
for y,name,handle,index,arm,alen in [(115,'Forward: 5′ → 3′','ACACTCTTTCCCTACAC','GAAGTC','CGGTAAGACACGACTTATCGCC',22),
                                    (260,'Reverse: 5′ → 3′','GACTGGAGTTCAGACG','ATACGT','GCGTCGATTTTTGTGATGCTCGTC',24)]:
    body.append(f'<text x="36" y="{y}" font-size="21" font-weight="bold">{name}</text>')
    for x,w,color,label,seq in [(36,340,P,f'Handle · {len(handle)} nt',handle),(392,205,M,'Well index · 6 nt',index),(613,548,I,f'Template-binding arm · {alen} nt',arm)]:
        textcolor=I if color==M else '#FFFFFF'
        body.append(f'<rect x="{x}" y="{y+18}" width="{w}" height="78" rx="9" fill="{color}"/><text x="{x+16}" y="{y+48}" fill="{textcolor}" font-size="18">{label}</text><text x="{x+16}" y="{y+77}" fill="{textcolor}" font-size="18" font-family="monospace">{seq}</text>')
body.append('''<text x="36" y="412" font-size="20">First cycles: only the 3′ binding arm anneals to the original template.</text>
<text x="36" y="447" font-size="20">Later cycles: the 5′ handle and index are copied into the double-stranded PCR product.</text>
<text x="36" y="497" font-size="20" font-weight="bold">Forward read layout</text>
<text x="36" y="533" font-size="18" font-family="monospace">F handle | F index | F arm | retained template | rc(R arm) | rc(R index) | rc(R handle)</text>
<text x="36" y="581" font-size="19">Example A1: forward GAAGTC + reverse ATACGT. Both ends are needed to identify a well.</text>''')
(OUT/'primer_architecture.svg').write_text(frame('How an indexed ColE1 primer works',''.join(body)))

f=['M527','M543','M544','M545','M854','M855','M856','M857'];r=['M528','M546','M547','M881']
body=['<text x="36" y="91" font-size="19">One forward oligo per row + one reverse oligo per column = 32 unique combinations from 12 oligos.</text>']
for j,name in enumerate(r):
    x=200+j*240
    body.append(f'<text x="{x+100}" y="140" text-anchor="middle" font-size="21" font-weight="bold">{j+1} · R-{name}</text>')
for i,name in enumerate(f):
    y=163+i*70
    body.append(f'<text x="36" y="{y+32}" font-size="21" font-weight="bold">{chr(65+i)} · F-{name}</text>')
    for j,rname in enumerate(r):
        x=200+j*240
        body.append(f'<rect x="{x}" y="{y}" width="220" height="57" rx="8" fill="white" stroke="#DDD4EB"/><text x="{x+110}" y="{y+24}" text-anchor="middle" font-size="20" font-weight="bold">{chr(65+i)}{j+1}</text><text x="{x+110}" y="{y+46}" text-anchor="middle" font-size="16">{name} + {rname}</text>')
body.append('<text x="36" y="756" font-size="18">Fill A1–A4, then B1–B4, then C1–C4. For 30 colonies, the last two positions are H1 and H2.</text>')
(OUT/'plate_32well.svg').write_text(frame('Your preferred well-to-primer map',''.join(body),790))

body=[]
panels=[('1. Exact template',['Unique binding sites','Opposite extension directions','Predict full indexed size']),
        ('2. Product check',['One band near predicted size','No unexpected second band','Negative control remains clear']),
        ('3. Quantification',['Unsaturated image','Verify all lane / ladder labels','Prefer independent DNA quantification']),
        ('4. Pooling review',['Equal-length products: normalise mass','Different lengths: normalise molarity','Retain original colony-to-well link'])]
for i,(title,lines) in enumerate(panels):
    x=36+i*287
    body.append(f'<rect x="{x}" y="94" width="266" height="210" rx="12" fill="white" stroke="#E6E1EF"/><rect x="{x}" y="94" width="266" height="7" fill="{M}"/><text x="{x+14}" y="134" font-size="20" font-weight="bold">{title}</text>')
    for j,line in enumerate(lines): body.append(f'<text x="{x+14}" y="{173+j*35}" font-size="14">{line}</text>')
body.append(f'<rect x="36" y="338" width="470" height="205" rx="12" fill="#251C44"/>')
for i,(label,bands,color) in enumerate([('Good',[405],'#FFFFFF'),('Blank',[],'#FFFFFF'),('Wrong',[465],'#FFFFFF'),('Multiple',[405,465],'#FFFFFF'),('Clipped',[405],'#FFFFFF')]):
    x=75+i*90
    body.append(f'<text x="{x}" y="575" text-anchor="middle" font-size="17">{label}</text>')
    for y in bands: body.append(f'<rect x="{x-17}" y="{y}" width="34" height="{13 if i==4 else 5}" fill="{color}"/>')
body.append('''<text x="545" y="373" font-size="21" font-weight="bold">Before pooling</text>
<text x="545" y="412" font-size="18">Exclude blank sample lanes, wrong-size and multiple products.</text>
<text x="545" y="447" font-size="18">Re-image clipped lanes before estimating relative abundance.</text>
<text x="545" y="482" font-size="18">Image QC alone cannot prove sequence identity or concentration.</text>
<text x="36" y="615" font-size="17">Gel cartoons are illustrative. Product size comes from the exact map, including the two indexed tails.</text>''')
(OUT/'pcr_qc_guide.svg').write_text(frame('PCR QC before an indexed pool',''.join(body),650))
print(OUT)
