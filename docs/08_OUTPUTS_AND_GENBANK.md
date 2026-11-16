# What the demultiplexing outputs look like

## Demultiplex only: no reference maps required

```bash
apcs plan -n 16 -o plan
apcs demux pool.fastq.gz --set cole1 --samples plan/samples.csv -o separated
```

```text
separated/
  demux_counts.csv                 all 48 bins, counts and loaded/unloaded status
  demux_stats.json                 counters, settings, input hashes, versions
  demux/
    a1.fastq, a2.fastq, …            reads assigned to each well, orientation normalised
    one_index_only.fastq            insufficient information to name a well
    ambiguous.fastq                 competing equal-best assignments
    unassigned.fastq                no usable index pair
    below_min_length.fastq          reads excluded by the configured length threshold
    read_assignments.tsv            one audit row for EVERY input molecule
```

A bin file is created when it receives reads. The original input remains authoritative. Read IDs and quality strings are preserved; reverse-oriented assigned reads and their qualities are reversed together. Duplicate read IDs cause a failure instead of artificially increasing independent depth.

The audit table includes read ID, length, status, well, orientation, forward/reverse index identities, mismatch counts, aligned arm start/end and combined edit score. Arm positions are **0-based half-open in the orientation-normalised sequence**. Rejected unorientable reads retain their original sequence. The supplied decoder currently requires six-base indexes and the documented arm-adjacent, two-ended architecture; a custom YAML must describe exactly what was used physically. Arbitrary barcode formats are not interchangeable.

## Reference screening and annotated outputs

```bash
apcs run pool.fastq.gz --set cole1 --samples plan/samples.csv \
  --expect wanted.gb --also parent.gb -o screen
```

```text
screen/
  demux/                          FASTQ bins and per-read assignment audit
  annotations/
    references/<name>.annotated.gb   original circular reference + primer/gap/window overlays
    references/<name>.binding.svg    exact-map binding graphic, when sites are compatible
    wells/a1.reference_evidence.gb   reference + well counts/verdict/span-depth overlays
    wells/a1.predicted_amplicon.gb   expected indexed linear product for A1's primer pair
    wells/a1.alignment_span_depth.bedgraph
  report/
    plate.html                    select a well, inspect evidence, download GenBank files
    00_PICK_THIS.md
    per_well.csv
    run_stats.json
```

### GenBank annotation meanings

| Feature | Meaning |
|---|---|
| Original gene/CDS/misc features | Retained from the supplied GenBank map; a FASTA has no such original annotations |
| `primer_bind` | Exact annealing arm on a circular template, or full indexed oligo on a predicted linear product |
| `OMITTED_FROM_AMPLICON` | Template sequence not present in this predicted PCR product |
| `diagnostic_window_*` | Reference-panel region eligible to distinguish supplied constructs, excluding the omitted gap and primer-overwritten arms |
| `screen_evidence_a1` | Assignment, mapping, informative counts and verdict for the well |
| `NO_ALIGNED_SPAN` / `LOW_SPAN_DEPTH` / `ALIGNED_SPAN_DEPTH` | Alignment-span coverage classes, not a pileup of called bases |
| Handle/index/primer-derived features | Expected oligo-derived sequence on the predicted indexed product |

Reference-evidence files **retain reference bases**. They are not read-derived consensuses and must not be treated as proof of a complete sequence. A mixed or no-call well is labelled accordingly. If no reference is resolved, the displayed primary reference is a navigation aid and a coverage track is unavailable; unknown coverage is not filled with reference or falsely reported as zero.

The BEDGraph has 0-based half-open reference coordinates. It describes aligned reference spans, including deletion spans, not a per-base consensus or independent evidence for the original sequence under the primers. GenBank feature coordinates use the format's 1-based convention; complementary strands and origin-crossing joined locations are retained.

To annotate a map without reads:

```bash
apcs annotate template.gb --set cole1 --well a1 -o annotated-map
```

The [public pUC19 example](../demo/pUC19/README.md) includes downloadable files. The [live demo](https://drjackxiaoyuchen.github.io/aureka-pooled-clone-seq/) exposes generated FASTQ bins, the audit table, and per-well GenBank outputs.
