# Public-release review · 2026-11-16

The supplied prototype was reviewed in a separate workspace; the original download was retained. Version 0.2.0 is a corrected experimental software release, not certification of the wet-lab workflow.

## Corrected issues

- Reject malformed or truncated FASTQ, ambiguous reference bases, multiple reference records and duplicate reference labels. GenBank parsing stops at the record terminator instead of deleting non-ACGT symbols and shifting coordinates.
- Reject duplicate/empty/inconsistent sample assignments, invalid plans, unsafe index mismatch budgets and reuse of existing demultiplex/report outputs.
- Compute unused sentinel positions from the actual chosen layout, rather than a row-major slice. Preserve an explicitly empty used-well set.
- Reject tied shifted index assignments and conflicting orientations, and require primer anchors at valid opposite molecule ends.
- Explicitly compare diagnostic-window alleles instead of treating minimap2's absence of a reported secondary alignment as evidence of a winning allele.
- Include small insertion/deletion diagnostics, reverse-reference coordinates, circular origins and rotated reference maps. Indistinguishable or uncovered reference alternatives receive no genotype vote.
- Count independent molecule votes once per read, including mixtures of unrelated references and mixtures spanning the related and unrelated groups. Gate calls on informative depth, mapping fraction and informative fraction. Single-reference evidence cannot receive PICK.
- Calculate breadth on the winning reference with unwrapped circular spans; do not combine unrelated coordinate systems or confuse zero-length spans with a full circle.
- Reject incomplete ladder calibration, fix sample labels when the ladder is not first, preserve image intensity range, suppress background-only bands, and withhold pooling volume for missing size checks, clipping and volume caps. Retain vertical migration space during auto-cropping.
- Include low-depth picks in HTML headlines, respect requested pick counts, escape HTML inputs, add keyboard-accessible evidence inspection and retain assigned/mapped/informative counts separately.
- Record input hashes, thresholds, software/dependency versions and per-template primer checks. Product lengths include both indexed tails; primer-overwritten bases are explicitly distinguished from verified template sequence.
- Replace unsupported comparative performance, universal-template, precise-mixture and private-map validation claims. Retain upstream primer attribution and MIT notice. Publish only generated sequences and simulated images.

## Validation

- **70 regression tests passed locally on Python 3.12.** Tests cover read parsing, index assignment, template geometry, planning, reference discrimination, mixtures, circular/reverse coordinates, gel failure cases, report escaping and verdict logic.
- The deterministic sequencing example generated 709 reads; 697 met the length filter and 444 were assigned. These are simulation-specific observations, not measured sequencing performance. All 16 loaded-well outcomes passed the independent truth checks: nine covered wanted wells, four parental wells, one mixed well, one low-depth well and one empty well.
- The simulated eight-lane gel produced the expected two passes, ladder, wrong-size, multiple-band, no-band, clipping-review and volume-review results. Its detector/calibration accuracy on real images is unestablished.
- Both 14-oligo CSV exports were regenerated from the canonical definitions, with sequence/length consistency checked. ColE1 handles, arms and indexes were matched against the public upstream files on this date.
- Packaging builds a source distribution and a wheel. CI repeats tests and both demos on Linux with Python 3.10, 3.12 and 3.13; inspect the Actions run for the actual status of your commit.

## Remaining limits

No independent wet-lab accuracy study, assay-specific sensitivity/specificity estimate, or head-to-head upstream benchmark is supplied. AmpR primers remain proposed and experimentally unvalidated. Exact plasmid compatibility, host off-targets, primer chemistry, PCR conditions, sequencing-provider acceptance and read yield require local verification.

Read classification is heuristic. Multiple index errors, index indels, repeats, recombination and unknown constructs can defeat the model. One related-backbone group is supported alongside unrelated single alternatives; more complex panels may be rejected. A two-edit score margin means a single isolated SNP need not be resolvable. Dominant read-vote support is not a calibrated biological mixture fraction. Gel lanes, bands and ladder pairing always require visual review, and clipping detection cannot identify every form of physical saturation.

PICK selects a candidate for independent sequence confirmation. It does not demonstrate base-perfect sequence, absence of low-frequency variants, sequence correctness in the omitted gap or underneath primer arms, or functional activity.
