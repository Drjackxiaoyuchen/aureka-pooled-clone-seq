# Aureka APCS · synthetic 16-well demonstration

*primer set cole1; 16 wells; expected ref_wanted*

## Pick this

**a4, c1** — candidates for miniprep and independent confirmation. Backups for review: a3, b3.

```
           1          2          3          4
         M528       M546       M547       M881
     +----------+----------+----------+----------+
 a 527|o a1   30 |o a2   22 |+ a3   31 |**a4   35 |
     +----------+----------+----------+----------+
 b 543|o b1   27 |o b2   31 |+ b3   32 |o b4   29 |
     +----------+----------+----------+----------+
 c 544|**c1   33 |X c2   24 |X c3   27 |X c4   25 |
     +----------+----------+----------+----------+
 d 545|! d1   39 |. d2    2 |. d3    0 |X d4   12 |
     +----------+----------+----------+----------+
 e 854|          |          |          |          |
     +----------+----------+----------+----------+
 f 855|          |          |          |          |
     +----------+----------+----------+----------+
 g 856|          |          |          |          |
     +----------+----------+----------+----------+
 h 857|          |          |          |          |
     +----------+----------+----------+----------+

  ** PICK      *- PICK, low depth   +  BACKUP     o  spare (clean, not needed)
  X  wrong construct   !  mixed colony, discard   .  no call (too few reads)
  number after the well is mapped reads; see CSV for assigned and informative counts
```

## Every well

| Well | Verdict | Construct | Reads | Informative | Identity | Breadth | Note |
|---|---|---|---|---|---|---|---|
| a4 | PICK | ref_wanted | 35 | 35 | 97.3% | 95% | reference screen passed; confirm the isolated plasmid |
| c1 | PICK | ref_wanted | 33 | 33 | 97.4% | 95% | reference screen passed; confirm the isolated plasmid |
| b3 | BACKUP | ref_wanted | 32 | 32 | 97.4% | 95% | reference screen passed; confirm the isolated plasmid |
| a3 | BACKUP | ref_wanted | 31 | 31 | 97.4% | 95% | reference screen passed; confirm the isolated plasmid |
| b2 | SPARE | ref_wanted | 31 | 31 | 97.3% | 95% | reference screen passed; confirm the isolated plasmid |
| a1 | SPARE | ref_wanted | 30 | 30 | 97.3% | 95% | reference screen passed; confirm the isolated plasmid |
| b4 | SPARE | ref_wanted | 29 | 29 | 97.3% | 95% | reference screen passed; confirm the isolated plasmid |
| b1 | SPARE | ref_wanted | 27 | 27 | 97.4% | 95% | reference screen passed; confirm the isolated plasmid |
| a2 | SPARE | ref_wanted | 22 | 22 | 97.4% | 95% | reference screen passed; confirm the isolated plasmid |
| c3 | WRONG CONSTRUCT | ref_parent | 27 | 27 | 97.4% | 95% | 100% of informative molecule votes say ref_parent; expected ref_wanted |
| c4 | WRONG CONSTRUCT | ref_parent | 25 | 25 | 97.3% | 95% | 100% of informative molecule votes say ref_parent; expected ref_wanted |
| c2 | WRONG CONSTRUCT | ref_parent | 24 | 24 | 97.4% | 95% | 100% of informative molecule votes say ref_parent; expected ref_wanted |
| d4 | WRONG CONSTRUCT | ref_parent | 12 | 12 | 97.3% | 95% | 100% of informative molecule votes say ref_parent; expected ref_wanted |
| d1 | MIXED - DISCARD | ref_parent | 39 | 39 | 96.4% | 95% | 56% of 39 informative molecule votes say ref_parent, the rest say ref_wanted and whole reads split the same way - mixed reference evidence; re-streak and confirm |
| d2 | NO CALL | ref_wanted | 2 | 2 | 97.2% | 95% | only 2 reads (need 8) |
| d3 | NO CALL | - | 0 | 0 | 0.0% | 0% | no reads survived the length filter |

## Run

- primer set **cole1**, blind spot ref_wanted: 274, ref_parent: 274 bp
- 709 reads in, 697 over 1000 nt, **399 assigned (57%)**
- 244 reads carried only one index, 54 carried neither
- **0 reads landed in unused index combinations** — a contamination/assignment signal; this does not measure all errors between loaded wells
- 2 wells could not be called: d2, d3 (reads present but not spanning a region that separates the references, or below min-depth)
- do not use: c2 (WRONG CONSTRUCT), c3 (WRONG CONSTRUCT), c4 (WRONG CONSTRUCT), d1 (MIXED - DISCARD), d4 (WRONG CONSTRUCT)

## What this does not tell you

Template blind spots (bp): ref_wanted: 274, ref_parent: 274. Primer-arm bases are supplied by the oligos and do not independently verify the original template. Dominant vote support is not a validated mixture fraction. No consensus or novel-variant analysis is performed. Confirm candidate isolates by an independent whole-plasmid sequence check.


<sub>generated 2026-11-16 by apcs</sub>
