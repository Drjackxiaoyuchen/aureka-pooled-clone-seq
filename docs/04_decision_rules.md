# Decision rules

Thresholds in `src/apcs/verdict.py` are screening heuristics. They have not been calibrated as sensitivity, specificity or mixture detection limits.

1. No mapped reads, fewer than 8 mapped reads, or no unique reference evidence → NO CALL.
2. Dominant informative-molecule fraction <0.80 with at least 3 independent votes for a second reference → MIXED - DISCARD. This label means do not use the well directly; re-streak and confirm. It does not establish the cause of mixed evidence.
3. Fewer than 8 informative molecules, mapped fraction <0.80, or informative fraction <0.50 → NO CALL.
4. Only one reference supplied → BACKUP, because alternatives were not tested.
5. Dominant reference outside the acceptable expected set → WRONG CONSTRUCT.
6. Expected reference with support <0.90, mean aligned identity <0.90, or breadth <85% → BACKUP for review.
7. Otherwise → PICK (low depth) for 8–14 informative molecules or PICK for ≥15.

Passing wells are ranked by informative molecules, identity, breadth and well name. The first `--picks` (default 2) remain picks; the next `--backups` (default 2) become backups; the rest become spares. A BACKUP can therefore mean a lower-ranked passing well or a flagged candidate; always read its note.

`--min-depth`, `--picks` and `--backups` are CLI options. Other heuristics are defined centrally in DEFAULTS; changing them requires documenting and evaluating the new operating point. MIXED is evaluated before the informative-depth gate so strong mixed evidence is not hidden by low depth.

Under independent, equally detectable 50/50 sampling, the chance that eight informative molecules all come from one component is `2 × (1/2)^8 = 0.78125%`. This idealised calculation concerns only that extreme event. It is not the false-negative rate of the configured mixture rule, which also requires three second-reference votes and is affected by assignment, PCR and mapping biases.

An apparently pure reference match does not rule out unrepresented constructs, low-frequency mixtures or new mutations. Independent verification of an isolated plasmid remains necessary.
