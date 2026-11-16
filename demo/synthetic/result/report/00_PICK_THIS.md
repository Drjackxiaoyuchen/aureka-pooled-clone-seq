# Aureka APCS · synthetic 16-well demonstration

*primer set cole1; 16 wells; expected ref_wanted*

## Pick this

**a4, b3** — candidates for miniprep and independent confirmation. Backups for review: a1, a6.

```
           1          2          3          4          5          6
         M528       M546       M547       M881       M858       M859
     +----------+----------+----------+----------+----------+----------+
 a 527|+ a1   35 |o a2   26 |o a3   34 |**a4   38 |o a5   31 |+ a6   35 |
     +----------+----------+----------+----------+----------+----------+
 b 543|o b1   32 |o b2   32 |**b3   35 |X b4   29 |X b5   33 |X b6   31 |
     +----------+----------+----------+----------+----------+----------+
 c 544|! c1   41 |. c2    2 |. c3    0 |X c4   10 |          |          |
     +----------+----------+----------+----------+----------+----------+
 d 545|          |          |          |          |          |          |
     +----------+----------+----------+----------+----------+----------+
 e 854|          |          |          |          |          |          |
     +----------+----------+----------+----------+----------+----------+
 f 855|          |          |          |          |          |          |
     +----------+----------+----------+----------+----------+----------+
 g 856|          |          |          |          |          |          |
     +----------+----------+----------+----------+----------+----------+
 h 857|          |          |          |          |          |          |
     +----------+----------+----------+----------+----------+----------+

  ** PICK      *- PICK, low depth   +  BACKUP     o  spare (clean, not needed)
  X  wrong construct   !  mixed colony, discard   .  no call (too few reads)
  number after the well is mapped reads; see CSV for assigned and informative counts
```

## Every well

| Well | Verdict | Construct | Reads | Informative | Identity | Breadth | Note |
|---|---|---|---|---|---|---|---|
| a4 | PICK | ref_wanted | 38 | 38 | 97.4% | 95% | reference screen passed; confirm the isolated plasmid |
| b3 | PICK | ref_wanted | 35 | 35 | 97.4% | 95% | reference screen passed; confirm the isolated plasmid |
| a1 | BACKUP | ref_wanted | 35 | 35 | 97.3% | 95% | reference screen passed; confirm the isolated plasmid |
| a6 | BACKUP | ref_wanted | 35 | 35 | 97.3% | 95% | reference screen passed; confirm the isolated plasmid |
| a3 | SPARE | ref_wanted | 34 | 34 | 97.3% | 95% | reference screen passed; confirm the isolated plasmid |
| b1 | SPARE | ref_wanted | 32 | 32 | 97.4% | 95% | reference screen passed; confirm the isolated plasmid |
| b2 | SPARE | ref_wanted | 32 | 32 | 97.3% | 95% | reference screen passed; confirm the isolated plasmid |
| a5 | SPARE | ref_wanted | 31 | 31 | 97.4% | 95% | reference screen passed; confirm the isolated plasmid |
| a2 | SPARE | ref_wanted | 26 | 26 | 97.4% | 95% | reference screen passed; confirm the isolated plasmid |
| b5 | WRONG CONSTRUCT | ref_parent | 33 | 33 | 97.4% | 95% | 100% of informative molecule votes say ref_parent; expected ref_wanted |
| b6 | WRONG CONSTRUCT | ref_parent | 31 | 31 | 97.3% | 95% | 100% of informative molecule votes say ref_parent; expected ref_wanted |
| b4 | WRONG CONSTRUCT | ref_parent | 29 | 29 | 97.3% | 95% | 97% of informative molecule votes say ref_parent; expected ref_wanted |
| c4 | WRONG CONSTRUCT | ref_parent | 10 | 10 | 97.1% | 95% | 90% of informative molecule votes say ref_parent; expected ref_wanted |
| c1 | MIXED - DISCARD | ref_parent | 41 | 41 | 96.4% | 95% | 51% of 41 informative molecule votes say ref_parent, the rest say ref_wanted and whole reads split the same way - mixed reference evidence; re-streak and confirm |
| c2 | NO CALL | ref_wanted | 2 | 2 | 97.2% | 95% | only 2 reads (need 8) |
| c3 | NO CALL | - | 0 | 0 | 0.0% | 0% | no reads survived the length filter |

## Run

- primer set **cole1**, blind spot ref_wanted: 274, ref_parent: 274 bp
- 709 reads in, 697 over 1000 nt, **444 assigned (64%)**
- 212 reads carried only one index, 41 carried neither
- **0 reads landed in unused index combinations** — a contamination/assignment signal; this does not measure all errors between loaded wells
- 2 wells could not be called: c2, c3 (reads present but not spanning a region that separates the references, or below min-depth)
- do not use: b4 (WRONG CONSTRUCT), b5 (WRONG CONSTRUCT), b6 (WRONG CONSTRUCT), c1 (MIXED - DISCARD), c4 (WRONG CONSTRUCT)

## What this does not tell you

Template blind spots (bp): ref_wanted: 274, ref_parent: 274. Primer-arm bases are supplied by the oligos and do not independently verify the original template. Dominant vote support is not a validated mixture fraction. No consensus or novel-variant analysis is performed. Confirm candidate isolates by an independent whole-plasmid sequence check.


<sub>generated 2026-11-16 by apcs</sub>
