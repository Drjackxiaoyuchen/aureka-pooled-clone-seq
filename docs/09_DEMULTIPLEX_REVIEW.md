# Demultiplexing review for version 0.3

The review checked source code, adversarial cases and a labelled synthetic pool. This is a software review and simulation exercise, not an experimental performance estimate.

## Changes

1. **Actual aligned arm spans.** Version 0.2 verified an edit distance but still located the reverse index using the original arm length. Two-base insertions/deletions could move the index outside its ±1 search slack. Version 0.3 uses Edlib's actual aligned arm end. It also removes the exact-seed requirement, which could miss an arm with several spaced substitutions.
2. **Assignment ambiguity.** All equal-best valid arm/index pairs across both orientations are compared. Competing wells receive an ambiguous status instead of an arbitrary assignment.
3. **Conservative default.** Exact index matching (`--index-mismatch 0`) is now the default. One-substitution recovery remains an explicit option. A six-base index distance of three protects only the stated one-error, correctly located observation model; multiple errors and boundary artifacts can still cross-assign.
4. **Per-read audit and streaming.** Every input read, including length-rejected reads, receives an audit row. Separate FASTQ buckets retain the rejection reason. Sequence data are written incrementally rather than retained for the whole pool. Read-ID deduplication and read-length summaries still consume memory proportional to read count. Duplicate IDs are rejected as potentially repeated molecules.
5. **Assay-aware diagnostics.** The omitted gap and primer-overwritten annealing-arm bases are excluded from construct-discrimination windows, including their flanks. Differences visible only in those regions cannot produce a confident call. The primary reference must have a unique exact compatible layout.
6. **Expected map versus observed evidence.** GenBank exports preserve supplied annotations and clearly separate circular reference sequence, predicted indexed amplicon and read-span evidence. Unresolved coverage is unavailable, not reference-filled.
7. **Preferred layout.** Planning and the 16-well demo now use A1–A4, B1–B4, C1–C4, D1–D4. Rows A–H support 32 positions with 12 oligos. The full 48-bin decode mapping is preserved, so no existing index pair changes identity.

## Independent read-label check

The same 709-read, seed-7 synthetic pool with the preferred four-column map was evaluated against immutable v0.2.0 and v0.3 code. Its input IDs encode the generating well, so assignment correctness was checked independently of the reported genotype. 697 reads passed the length threshold.

| Decoder | Index mismatch allowance | Assigned | Correct source well | Wrong source well |
|---|---:|---:|---:|---:|
| v0.2.0 | 1 | 443 | 439 | 4 |
| v0.3 | 1, explicit recovery mode | 486 | 483 | 3 |
| v0.3 | 0, new default | 399 | 399 | 0 |

Zero wrong assignments in this one toy dataset is not a measured zero-error guarantee. The stricter default trades read recovery for a more conservative operating point. Compare settings with known controls, unused combinations and sufficient per-well evidence; a higher assignment count alone is not evidence of better accuracy.

Other limits remain: six-base indexes only, primer anchors searched within 250 bases of read ends, inability to identify all PCR chimeras/cross-contamination, bias from amplification and sequencing, unrepresented constructs, and heuristic genotype thresholds. Raw reference inputs, sample identities and read data must be retained.
