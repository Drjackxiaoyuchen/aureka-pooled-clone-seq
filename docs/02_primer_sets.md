# Primer designs and index geometry

Canonical definitions: [`primer_sets.yaml`](../src/apcs/data/primer_sets.yaml). All sequences and oligo orders are written 5′→3′.

| Set | Forward binding arm | Reverse binding arm | Status |
|---|---|---|---|
| ColE1 | `CGGTAAGACACGACTTATCGCC` | `GCGTCGATTTTTGTGATGCTCGTC` | Public upstream design; verify exact template and assay |
| AmpR | `GGTGCCTCACTGATTAAGCATTGG` | `TCGTTCATCCATAGTTGCCTGACTC` | Proposed Aureka design; wet-lab performance unvalidated |

Forward handle: `ACACTCTTTCCCTACAC` (17 nt). Reverse handle: `GACTGGAGTTCAGACG` (16 nt). Each oligo is handle + six-base index + binding arm. ColE1 forward/reverse oligos are 45/46 nt; AmpR oligos are 47/47 nt.

There are eight forward indexes and six reverse indexes. Forward index → row a–h; reverse index → column 1–6. The resulting 48-position index grid is not a claim that every 48-well pool will yield sufficient coverage. Adding combinations requires redesigned configuration, index-distance checks and experimental validation.

The supplied indexes have minimum Hamming and Levenshtein distance 3 within each direction. This permits unique correction of one substitution for a fixed, correctly located six-base observation. The runtime uses Hamming index matching and ±1-base positioning slack, not full index-indel correction. Multiple errors, index boundary errors and synthesis/PCR artifacts can still misassign reads. The all-A/T forward index M545 has 0% GC and should be specifically evaluated in pilots.

`apcs check` scans both strands of a circular map and requires exactly one site per arm on opposite strands. It reports template span, omitted gap and full predicted product length including tails. It does not evaluate melting temperature, primer dimers, host off-targets, PCR efficiency or whether the original bases underneath a primer are correct. It exits nonzero on an incompatible map.

The configuration's nominal 274-bp (ColE1) and 21-bp (AmpR) gaps are design expectations. Only the per-map calculation is applicable to a supplied sequence. No private-map validation counts or matched-temperature claims are included in this release.

`apcs order --set cole1 -o orders` exports the 12 oligos for the preferred four-column layout as CSV and XLSX. Add `--columns 6` for the full 14-oligo catalog. Synthesis scale and purification are suggestions; confirm with the supplier. The [upstream primer consensus](https://github.com/daspacio9/pooledpremiumpcrpipeline2/blob/main/working_directory/ref/cole1-primer-consensus.fasta) and [index definitions](https://github.com/daspacio9/pooledpremiumpcrpipeline2/blob/main/working_directory/ref/cole1-barcodes-48.fasta) were checked on 2026-11-16. Retain their attribution.
