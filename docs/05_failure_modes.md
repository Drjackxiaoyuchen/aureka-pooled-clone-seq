# Failure modes

| Observation | Investigation |
|---|---|
| Few reads assigned; many one-index-only | Check physical molecule truncation, sequencing prep and end preservation. One index identifies a row or column, not a well. |
| Neither index found | Check primer set, raw read orientation/layout, read length and arm errors. Seed matching is heuristic. |
| Most reads fail to map | Check exact reference maps and unrepresented constructs; do not interpret a handful of matching reads as a pure sample. |
| Many mapped reads but few informative | Diagnostic regions may be absent from reads or inside the omitted gap. Identical/ambiguous references receive no call. |
| Only one supplied reference | Add parent and relevant intermediates; the tool otherwise reports reference-match backups. |
| Uneven per-well counts | Check amplification, concentration, molecule length and pooling records. Gel intensity is only a relative mass proxy. |
| Reads in unused combinations | Investigate indexing, carryover and sample sheet. Counts do not identify the cause or quantify all loaded-well errors. |
| A row or column looks wrong | Reconcile primer-tube identities and the actual pipetting map with the sample sheet. |
| Gel refuses calibration | Confirm the exact ladder, every detected band, crop and lane identity; do not pair missing ladder bands by guesswork. |
| Gel reports clipping or volume cap | Re-image, quantify, re-amplify or revise the pool manually. No pooling volume is approved automatically. |
| Lane count wrong | Inspect/crop the image. The detector may miss an empty lane; `--lanes` cannot reconstruct its position. |
| Existing output error | Choose a new output directory; stale demultiplexed files must not be mixed with a new run. |

Keep sample sheets, original colony links, reference versions, raw FASTQ/images and run metadata. Reports do not replace the underlying evidence.
