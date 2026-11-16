# Public pUC19 binding example

Source: NCBI GenBank accession [M77789.2](https://www.ncbi.nlm.nih.gov/nuccore/M77789.2), downloaded 2026-11-16. The unmodified source record is retained. Its original broad feature annotations are preserved; functional genes are not invented from a FASTA or an unannotated record.

`apcs annotate demo/pUC19/pUC19_M77789.2.gb --set cole1 --well a1 -o a-fresh-directory` adds exact primer-binding annotations, an omitted-region annotation, a binding SVG and a predicted indexed amplicon GenBank.

On this exact 2686-bp map, the forward arm matches complement(2163..2184), the reverse arm matches 2459..2482, and the omitted interval is 2185..2458 (274 bp). The product is 2686 − 274 + 45 = 2457 bp including indexed tails. The 46 annealing-arm bases are supplied by oligos; only the other 2366 template positions can provide independent amplicon-read evidence. This is a reference illustration, not experimental sequencing validation.
