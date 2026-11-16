# Analysis and evidence

## Inputs

One FASTQ/FASTQ.gz with four lines per record, a sample sheet and one or more single-record reference files. FASTA and GenBank references must have unambiguous A/C/G/T sequence; ambiguous bases and multiple records are rejected rather than deleted or silently ignored. All references are treated as circular plasmids. Duplicate filename stems are rejected.

## Demultiplexing

The program aligns the complete primer arm with Edlib within 250 bases of the corresponding read end and uses its actual aligned span, including indels. A complete call needs both indexes in a valid forward/reverse layout. Hamming index matching is exact by default, with ±1-base positional slack. One-substitution recovery is opt-in (`--index-mismatch 1`). Equal-best competing wells across arm/index positions and orientations are rejected. The optional handles are recorded as corroboration; they are not required.

Every considered read enters assigned, one-index-only, unassigned, impossible-pair or ambiguous counters. Reads shorter than `--min-len` (1000 by default) are counted separately. Unused combinations monitor contamination and assignment artifacts, but miss errors between loaded wells; their count is not an unbiased contamination estimate. All 48 combinations loaded means no unused sentinels.

## Reference discrimination

Minimap2 maps reads to doubled circular references. Related maps are aligned with rotated query origins, preventing origin-spanning diagnostic blocks from being silently clipped. Substitutions and indels define windows with 60-bp flanks; nearby sites within 200 bp are grouped. Edlib scores every window allele explicitly in both orientations. A nearly complete window must align within a 15% error budget, and the aggregate winner needs at least a two-edit margin. Shared alleles, ties and unobserved diagnostic regions contribute no genotype vote.

One informative molecule contributes one reference vote, even if it covers multiple windows. Unrelated references use unique best matched-base alignment support with a two-base margin. One group of related backbones plus unrelated single alternatives is supported. Multiple related groups or partial backbone homology may be rejected and should be analysed separately. Large rearrangements, repeats, single-base discrimination and mixtures outside the supplied reference panel require additional methods.

The primary reference must have unique exact primer sites. The omitted gap and primer-overwritten arms are masked from diagnostic windows and their flanks. Differences only in these unobservable regions cannot identify the original template.

Purity in CSV means **dominant-reference fraction among informative molecule votes**, not a calibrated cellular, plasmid or allele mixture fraction. Coverage breadth and mean depth are calculated against the winning reference. A single-reference run establishes alignment compatibility only and cannot receive PICK.

## Reports

- `assigned_reads`: reads with both indexes assigned to a well.
- `n_reads`: assigned reads mapping to at least one reference.
- `n_informative_reads`: mapped reads casting a unique reference vote.
- `n_ambiguous_reads`: mapped reads with no unique vote.
- `unmapped`: assigned reads that do not map.
- `breadth_pct`: union of alignment spans on the winning circular reference.
- `mean_depth`: aligned reference span summed across reads / reference length.
- `n_site_calls`: compatibility field; from version 0.2 onward it counts informative molecules, not windows.

`report/run_stats.json` records software/dependency versions, decision thresholds, exact input SHA-256 hashes, index settings and per-template primer checks. No absolute input paths are embedded automatically. CSV, Markdown, interactive HTML and a text plate are written beside it. Retain raw inputs and this metadata together.

## Gel helper

Image analysis uses heuristic lane detection, band peaks and a log-linear ladder calibration. All supplied ladder sizes must have detected band partners; an inconsistent count is rejected. Without calibration and expected size, no positive pooling volume is produced. Possible clipping is based on pixels near the image range limit; it cannot detect all physical saturation, compressed-image artifacts or 12-bit clipping in a 16-bit file. Pooling volumes assume equal-length products and equal loaded volumes, and require visual verification and preferably independent quantification. `--lanes` limits detected lanes and checks the count; it cannot reconstruct a missing empty lane.
