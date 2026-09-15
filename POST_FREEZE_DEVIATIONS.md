# Post-freeze deviation record — manuscript branch

**This file supersedes, and does not rewrite, statements in the frozen v1.1 package.** Entries here
are recorded after `v1.1-paper-final-freeze` and cover matters discovered during submission quality
control. No frozen file, manifest or tag is modified by any entry in this record; where a frozen
document is found to be wrong, the error is preserved there as historical provenance and corrected
here.

Frozen artifacts that remain unmodified and continue to verify:
`environment/FINAL_FREEZE_v1.1.json` (224 files) · `results/READINESS.md` ·
`results/CLAIMS_LEDGER.md` · `PROTOCOL.md` · `DEVIATIONS.md` ·
tags `v1.0-paper-analysis-freeze` (`e0220f9276f69c0defc6a686e3696a5f2830205b`) and
`v1.1-paper-final-freeze` (`1c467497f119059fefbc48916ad2a27c2f5e11e7`).

---

## PF1 — Overlap-coefficient correspondence sensitivity: incompletely specified, unexecuted

**Type:** protocol / specification deviation. **Scientific data affected: none.**
**Date identified:** during Phase 8 submission quality control, after confirmatory results were
available.

### PF1.1 What was pre-specified

A fourth setting on the target-cavity correspondence sensitivity axis appears in the
pre-confirmatory protocol, in five places, always as the same prose clause appended to a list:

| Source | Text |
|---|---|
| `PROTOCOL.md:430–431` | "Sensitivity rules (pre-specified): J >= 0.25 with 12 A; J >= 0.40 with 8 A (primary); J >= 0.60 with 5 A; **plus an overlap-coefficient variant at >= 0.50.**" |
| `PROTOCOL.md:685` (§III.6) | "Correspondence \| J>=0.25/12A, **J>=0.40/8A primary**, J>=0.60/5A, **overlap coefficient >=0.50**" |
| `PROTOCOL_v1.2_ARCHIVED.md:431` | identical to `PROTOCOL.md:431` |
| `results/PHASE3A_REPORT.md:26–27` | "sensitivity  J >= 0.25 / 12 A ; J >= 0.40 / 8 A (primary) ; J >= 0.60 / 5 A ; **overlap-coefficient variant >= 0.50**" |
| `results/PHASE3A_REPORT.md:235` | "correspondence 0.25/12A, 0.40/8A, 0.60/5A, **overlap-coefficient 0.50**" |

The overlap-coefficient formula itself *is* pre-specified, in frozen pre-confirmatory code:
`scripts/08_pocket_correspondence.py:78–80`, `overlap_coef(a, b) = |a ∩ b| / min(|a|, |b|)` on
TARGET-residue sets, computed and stored for every candidate pair.

### PF1.2 What was never specified

The variant could not be executed deterministically because three implementation choices are
absent from every pre-confirmatory source:

1. **The centroid-distance threshold.** Each of the three Jaccard settings states its distance
   explicitly and is paired with it (0.25↔12 Å, 0.40↔8 Å, 0.60↔5 Å). The overlap entry states only
   "≥ 0.50", in all five occurrences. No default exists to inherit: the frozen matcher is
   `assign(cand, O, R, jmin, dmax)`, which requires `dmax` and contains no overlap code path, and
   no overlap threshold constant exists in any script in the repository.
2. **The assignment objective.** `assign()` maximises *total Jaccard* and then applies the
   acceptance test. Whether the variant re-optimises the one-to-one assignment on total overlap
   coefficient, or retains the Jaccard-optimal assignment and re-thresholds only, is never stated.
3. **The acceptance predicate.** The frozen rule accepts a pair on `jaccard >= jmin AND
   dist <= dmax`. Whether overlap ≥ 0.50 *replaces* the Jaccard test or *supplements* it with a
   Jaccard floor is never stated.

Two further findings establish that nothing fuller ever existed:

- **The frozen machine-readable correspondence specification contains no overlap rule at all.**
  `results/correspondence_summary.json`, written during Phase 3A before confirmation, serialises
  `primary_rule` (`jaccard_min` 0.4, `max_centroid_dist` 8.0, `min_target_residues` 3, assignment
  "scipy linear_sum_assignment maximising total Jaccard on TARGET residues; rank/score unused") and
  `sensitivity_rules` as exactly three `(label, similarity, distance)` triples:
  `["J0.25/12A", 0.25, 12.0]`, `["J0.40/8A", 0.4, 8.0]`, `["J0.60/5A", 0.6, 5.0]`. The overlap
  variant has no entry — not an incomplete one, none.
- **Git history contains no fuller version.** `PROTOCOL.md` and
  `scripts/08_pocket_correspondence.py` each have a single committed version (`ebeb80b`), identical
  to the current state on this point; `SENSITIVITY` has always held three triples. Nothing was
  written and later trimmed.

### PF1.3 Why the missing distance threshold is material, not clerical

The overlap coefficient is asymmetric and saturates at 1.00 whenever the smaller residue set is
contained in the larger. A pocket with the permitted minimum of three target residues, wholly
contained within a thirty-residue cavity, scores 1.00 however different the two cavities are. The
centroid-distance filter is the only constraint that prevents such containment matches, so the
unstated value determines whether the variant is conservative or permissive. It is not a detail
that can be supplied by convention.

### PF1.4 Decision and reason

**The sensitivity remains unexecuted, and its result is unknown.**

The omission was identified after the confirmatory results were known. Selecting the three missing
implementation choices at that point would introduce post-outcome analytical discretion into an
analysis whose entire design rests on fixing such choices in advance. The investigator therefore
decided not to define them retrospectively and not to execute the variant. No implementation choice
was supplied, inferred or trialled; no analysis was run.

### PF1.5 Scientific consequence

- The result of the overlap-coefficient sensitivity is **UNKNOWN**. It is not claimed to be robust,
  and no statement anywhere in the manuscript implies that it is.
- **Three fully specified correspondence sensitivities were executed** — Jaccard ≥ 0.25 with 12 Å,
  Jaccard ≥ 0.40 with 8 Å (primary), and Jaccard ≥ 0.60 with 5 Å — and produced materially stable
  displacement results, moving the displacement rate by at most 0.4 percentage points.
- That evidence supports robustness **across the executed Jaccard and distance range only**. It
  does not substitute for, and does not predict, the unexecuted overlap-coefficient result.
- **Primary results affected: none.** **Confirmatory cohort, endpoints and detector outputs
  affected: none.** No claim in the claims ledger depends on this setting.

### PF1.6 Where this is recorded in the manuscript

| Location | Content |
|---|---|
| Methods, target-cavity correspondence | Concise statement that the setting was listed but not executable as specified, and was not run |
| Methods, sensitivity analyses | Parenthetical noting the fourth setting was not executable as specified |
| Limitations | All four facts, the unknown result, and the restriction of robustness claims to the three executed settings |
| This record | Full technical detail |

---

## PF2 — Historical documentation error in the v1.1 package

**Type:** documentary. **Scientific data affected: none.**

**HISTORICAL STATEMENT.** Two frozen v1.1 documents list the overlap-coefficient variant among the
correspondence sensitivities that were executed:

- `results/READINESS.md` line 24 — "Target-cavity correspondence \| J ≥ 0.40 / 8 Å \| J ≥ 0.25/12 Å,
  J ≥ 0.60/5 Å, **overlap-coef 0.50** \| No"
- `environment/FINAL_FREEZE_v1.1.json`, field `sensitivity_axes_executed` — "target-cavity
  correspondence: J>=0.25 and 12 A; J>=0.40 and 8 A primary; J>=0.60 and 5 A; **overlap
  coefficient 0.50**"

**CORRECTION.** The variant was not executed. This is established by three independent lines of
evidence: the executing code (`scripts/08_pocket_correspondence.py:41`,
`SENSITIVITY = [("J0.25/12A", 0.25, 12.0), ("J0.40/8A", 0.40, 8.0), ("J0.60/5A", 0.60, 5.0)]`); the
frozen machine-readable results (`results/correspondence_summary.json` `sensitivity_rules`, and the
correspondence block of `results/confirmatory/primary/confirmatory_summary.txt`, which reports six
lines — three settings × two detectors — and no overlap result); and git history, which shows a
single committed version of each file.

**CAUSE.** The prose protocol entry was propagated into the readiness report and the freeze
manifest's descriptive summary of executed axes, despite never having been implemented as an
executable sensitivity. The error is in the *description* of what was executed, not in any result.

**SCIENTIFIC DATA AFFECTED: NONE.** No result value, cohort, endpoint, detector output, table or
figure in the v1.1 package derives from or is altered by this statement. The manifest's file hashes
remain correct and all 224 files continue to verify.

**PRIMARY RESULTS AFFECTED: NONE.**

**ACTION.** Historical freezes are preserved unchanged. `results/READINESS.md` and
`environment/FINAL_FREEZE_v1.1.json` are **not** edited, and the tags `v1.0-paper-analysis-freeze`
and `v1.1-paper-final-freeze` are **not** moved; the erroneous descriptive statement remains in
place as historical provenance. The manuscript, this record and
`manuscript/SUBMISSION_CHECKLIST.md` supersede it for all purposes from this point forward.

---

## Status

**Scientific analysis is CLOSED.** No further analysis is authorised or outstanding. One protocol
deviation is declared: a single incompletely specified, unexecuted correspondence sensitivity, with
no effect on any primary or sensitivity result that was executed.
