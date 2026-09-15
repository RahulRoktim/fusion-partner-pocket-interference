# Protocol Deviations and Amendments

Append-only log. Every entry records: date, amendment ID, whether it was made **before or after**
the analyst saw any outcome data, the change, and the reason.

Protocol v1.0 is preserved verbatim in `PROTOCOL_v1.0_ARCHIVED.md`. `PROTOCOL.md` carries the
current version (v1.1).

---

## Amendment set A — protocol v1.0 → v1.1

**Date:** 2026-09-15
**Outcome status:** **PRE-OUTCOME.** Recorded before any pocket detector was installed, and before
any detector was executed on any structure. No pocket prediction, ranking, or classification result
existed in this project at the time of writing. The only empirical data seen at this point were
RCSB SIFTS annotations, entry-level metadata (resolution, method), and cohort counts — none of
which are outcome variables.
**Requested by:** principal investigator, in writing, as outcome-blind amendments.

---

### A1 — Expression-tag (TAG) residues deleted from BOTH conditions

**Was (v1.0, §4.5):** TAG residues were retained in the ORIGINAL condition and deleted alongside
FUSION and LINKER residues in the FUSION-REMOVED condition.

**Now (v1.1):** TAG residues are deleted from **both** conditions before pocket detection.

The paired contrast is therefore exactly:

| Condition | Residues present |
|---|---|
| ORIGINAL | target + fusion + linker |
| FUSION-REMOVED | target only |

**Reason.** In v1.0 the paired difference confounded two distinct interventions — removing the
fusion domain (the exposure of interest) and removing terminal expression-tag remnants (a nuisance
construct feature that is not the question). Any pocket involving tag residues would have appeared
as a spurious "rescue." Deleting TAG from both arms makes the fusion domain the single manipulated
variable.

**Implementation.** Deleted TAG residue IDs are recorded per structure in the manifest field
`tag_residues_deleted`. TAG residues are, per §3.4, modelled entity positions in no SIFTS-aligned
region located at either terminus.

---

### A2 — Reference-site coverage rule; resolves the §7.2 / §7.4 conflict

**Was (v1.0):** §7.2 named paired restoration of the *biological target reference site* as the
full-study confirmatory primary endpoint, while §7.4 restricted reference-site analyses to the
subset with a defined site and forbade promoting ligand recovery unless coverage exceeded 60%.
These two statements could contradict each other.

**Now (v1.1), frozen before any detector result is viewed:**

- **If ≥ 60%** of the eventual confirmatory cohort has an independently defined biological target
  reference site → paired restoration of the biological target site's rank **may** serve as the
  full-study confirmatory primary endpoint.
- **If < 60%** → the cohort-level fusion-capture / ranking endpoint **remains** primary, and
  biological target-site recovery is a **restricted-subset secondary** endpoint.

**Pilot handling.** The pilot **reports** reference-site coverage but does **not** alter pilot
selection on the basis of ligand availability. Selection remains the seeded, outcome-blind draw of
§9.

**Reference-site definition is outcome-independent.** Reference sites are derived solely from
deposited ligand coordinates plus a pre-declared artifact-ligand exclusion list. Pocket-prediction
output is never consulted when defining a reference site.

---

### A3 — Deletion-boundary artifact control

**New in v1.1.** Excising an internally inserted fusion exposes artificial surface that did not
exist in the deposited structure and could itself generate cavities. An apparent "rescue" of target
ranking could therefore be an artifact of the cut rather than genuine restoration.

**Implementation:**

1. For every fusion-removed structure, record the residues and coordinates **bordering each
   deletion** (`deletion_boundaries`): the last target residue before and the first target residue
   after each excised span, with their CA coordinates.
2. For every predicted pocket in the REMOVED condition, compute the Euclidean distance from the
   pocket centroid to the nearest deletion boundary (`dist_to_deletion_boundary`).
3. **Primary paired analysis: all predictions, no boundary filter.**
4. **Pre-specified sensitivity analysis at 8 Å:** when interpreting apparent target-site rescue,
   disregard *newly appearing* pockets in the REMOVED condition whose centroid lies ≤ 8 Å from a
   deletion boundary.

The 8 Å threshold is fixed now, before results, and will not be re-tuned. Rationale: 8 Å is
approximately the span over which a single excision can nucleate a new alpha-sphere cluster while
remaining small enough not to swallow genuine adjacent target sites.

---

### A4 — Go/no-go language: the 10% figure is a decision threshold, not a statistical null

**Was (v1.0, §10):** framed ≥ 6/24 via a binomial calculation against a 10% rate, phrased in a way
that could read as a hypothesis test.

**Now (v1.1):** the numerical gates **G1/G2/G3 are retained unchanged**. However:

- 10% is **not** an empirically established background rate. No study has measured the rank-1
  fusion-capture rate, which is the entire point of this work.
- The ≥ 6/24 figure is described as a **pre-specified minimum practically interesting pilot signal
  and decision threshold**.
- Any binomial arithmetic against 10% is reported as **descriptive context only** and must never be
  presented as confirmatory hypothesis evidence, a p-value for the scientific hypothesis, or
  evidence against a real null.

The pilot is a feasibility and signal-detection exercise. It generates no confirmatory inference.

---

### A5 — Benchmark / training-set construct audit adopted as a second major aim

**New in v1.1.** Adopted as a **pre-specified second major study aim**, deferred to the phase after
the pocket-prediction pilot passes its gate. Not performed during this pilot.

**Binding language constraints:**

1. Do **not** describe any benchmark dataset as "defective," "broken," or "invalid."
2. Distinguish **TRAINING**, **VALIDATION**, **TEST/BENCHMARK**, and **REFERENCE** datasets; these
   are different objects with different implications and must never be conflated.
3. Do **not** claim any predictor was trained on fusion constructs unless that predictor's actual
   documented training set establishes it. Presence of a structure in a dataset with a similar name
   is not evidence about a given model's training data.
4. Distinguish **presence** of chimeric structures in a dataset from **demonstrated bias** in
   reported performance. Presence is an observation; bias is a claim requiring its own evidence.

---

## Amendment set B — implementation decisions taken during Phase 2A–2C

**Date:** 2026-09-15
**Outcome status:** **PRE-OUTCOME** unless explicitly marked otherwise below.

### B1 — fpocket execution route (pre-outcome)

Docker Desktop's daemon was not running and `sudo` in WSL requires an interactive password, so
neither `apt install` nor a container build was available without user action. fpocket is therefore
**built from pinned source in WSL2 Ubuntu 24.04 with a user-prefix install** (no root). This is a
route decision, not a substitution: fpocket itself is unchanged, and no other detector replaces it.
Exact source tag, commit, build flags and binary hash are recorded in `environment/`.

### B2 — mmCIF parsing library (pre-outcome)

Gemmi is used for coordinate parsing and writing (mmCIF → per-chain PDB). Chosen because it
preserves `label_seq_id` ↔ `auth_seq_id` correspondence, which is required to map SIFTS entity
numbering onto deposited coordinates without hand-built heuristics. Version pinned in
`environment/`.

### B3 — fpocket input format (pre-outcome)

fpocket 4.x consumes PDB format. Both detectors are therefore fed the **same prepared single-chain
PDB file** per structure per condition, so that no format difference can contribute to a
between-detector discrepancy. P2Rank accepts PDB natively.

<!-- Further entries appended below as they arise. Any entry made after outcome data were visible
     MUST be marked POST-OUTCOME in bold and reported both with and without the affected data. -->

---

## Amendment set C — logical impossibility found at implementation

**Date:** 2026-09-15
**Outcome status:** **POST-OUTCOME in timing.** Discovered while implementing the paired
analysis, after detector runs had completed and pockets had been classified. The reasoning is
**structural, not empirical** — it follows from how the REMOVED condition is built, and would have
been true regardless of any result. Both the degenerate and the corrected analyses are reported
so that nothing is concealed.

### C1 — The paired fusion-association endpoint is degenerate; G3 re-read

**The problem.** Amendment A1 defines the REMOVED condition as *target residues only*. A pocket is
classified `FUSION_DOMINATED`, `INTERFACE` or `LINKER` on the basis of FUSION and LINKER residues —
**none of which exist in the REMOVED condition.** Therefore no REMOVED-condition pocket can ever be
fusion-associated. Confirmed empirically: 0 of 86 P2Rank and 0 of 373 fpocket REMOVED pockets are
fusion-associated.

Consequences for the protocol as written:

- Protocol v1.1 §7.2, coverage < 60% branch — "P(rank-1 pocket is fusion-associated), removed vs
  original" — has a structurally fixed cell: the "not-fusion → fusion" discordant count **c ≡ 0**.
  McNemar on it is vacuous; it can only ever report "removal removes fusion pockets," which is a
  restatement of the manipulation, not a finding.
- Gate **G3** — "paired removal flips the rank-1 pocket to the target site in ≥ 5/24 with ≤ 1
  reverse flip" — has a trivially satisfied second clause under that reading.

Amendment A4 permits retaining G1/G2/G3 "unless implementation reveals a logical impossibility."
This is that case, and it is handled by **re-reading G3 rather than weakening it**.

**Resolution.** The paired outcome variable is **"the rank-1 pocket is `TARGET_DOMINATED`"** — the
literal reading of "flips the rank-1 pocket to the target site." Under this definition:

- b (improved) = rank-1 not TARGET_DOMINATED in ORIGINAL → TARGET_DOMINATED in REMOVED
- c (worsened) = rank-1 TARGET_DOMINATED in ORIGINAL → not TARGET_DOMINATED in REMOVED

Both cells can be non-zero, so the test is meaningful. The numeric thresholds of G3 (≥ 5 and ≤ 1)
are **unchanged**. The degenerate fusion-association version is reported alongside, labelled as
structurally determined.

### C2 — Empty prediction sets are not rescues

**Observation.** P2Rank predicts **zero pockets** in the REMOVED condition for 6 of 24 structures —
3N94, 3OAI, 5GPP, 5H7Q, 5JQE (MBP) and 5YQR (T4L) — all small target domains (86–182 modelled
residues). fpocket always returns at least one pocket.

**The trap.** Under the degenerate endpoint, "no pocket at all" satisfies "rank-1 is not
fusion-associated" vacuously and would be scored as a successful rescue. That would manufacture
support for the hypothesis out of a detector returning nothing.

**Resolution.** A structure with no REMOVED-condition prediction is assigned the explicit outcome
`NO_POCKET_PREDICTED`. It counts as **not** TARGET_DOMINATED, therefore **not** an improvement, and
it is never silently dropped from the denominator. This is the conservative direction: it works
**against** the study hypothesis. Counts are reported separately so the reader can see exactly how
many pairs are affected.

### C3 — Smoke-test output was seen during parser validation

During Phase 2C, detector output for **5IU7** (a structure that is *not* in the pilot cohort) was
inspected to validate the parsers. Protocol v1.1 was already frozen at that point, so no threshold,
class definition or gate could have been influenced. Recorded for completeness.

### C4 — fpocket pocket-file off-by-one, caught by tests before the pilot ran

The initial parser assumed fpocket wrote `pockets/pocket{N-1}_atm.pdb` for info-file block
`Pocket N`. fpocket 4.2.3 in fact writes **1-based** `pocket{N}_atm.pdb`. Left uncorrected, every
fpocket pocket's residue set would have been shifted by one rank and every fpocket classification
in the study would have been wrong. Caught by `tests/test_detector_parsers.py` check **D3**, which
verifies the file↔block correspondence independently using each pocket's alpha-sphere count.
**Fixed before any pilot run executed.**

---

## Amendment set D — protocol v1.1 → v1.2 (confirmatory freeze)

**Date:** 2026-09-15
**Outcome status:** **PILOT-INFORMED.** Every change below was made after seeing development-set
results and is labelled as such. The 24 pilot structures are permanently DEVELOPMENT data and may
never be reused as confirmatory observations. The confirmatory protocol (PROTOCOL.md Part II) is
frozen **before any detector has been run on any confirmatory structure** — no confirmatory
detector output exists at the time of writing.
**Version archive:** v1.0 in `PROTOCOL_v1.0_ARCHIVED.md`, v1.1 in `PROTOCOL_v1.1_ARCHIVED.md`.

### D1 — Target-pocket correspondence replaces the degenerate paired outcome (PILOT-INFORMED)

The pilot outcome "rank-1 becomes TARGET_DOMINATED" is retained **for audit only** and carries no
causal weight, for the reason given in C1 and restated by the reviewer: after FUSION/LINKER
deletion essentially every pocket is target-composed, so the transition cannot identify restoration
of a specific cavity. Protocol II.1 defines a deterministic correspondence rule using **only target
residues**, with rank and score excluded from the matching entirely. Endpoints B1 (displacement)
and B2 (paired rank change) are defined on matched cavities.

Threshold selection used only matching-quality metrics (match rate, mean Jaccard, mean centroid
distance, assignment ambiguity), computed and reported **before** any rescue quantity
(`results/correspondence_threshold_sweep.tsv`). The count of matched pairs is flat across
Jaccard 0.20-0.60, so the threshold has essentially no discretion.

### D2 — E9 REVISED (PILOT-INFORMED)

See PROTOCOL.md II.5. Three parts: fix the empty-segment fallback that made E9 appear to exclude
1,223-1,439 entities when its true marginal effect is 37; keep the 0.20 threshold on
measurement-validity grounds; add disorder as a pre-specified 0.35 sensitivity threshold and as a
model covariate.

**Correction to the Phase 2G report.** That report stated "the E9 disorder filter excluded 1,223 of
1,971 candidates — more than any other criterion" and treated it as the dominant population
determinant. That framing was **wrong**. Most of those flags were spurious (empty-segment fallback
on non-chimeric entities) or redundant with E7/E8. `results/PHASE2G_GO_NOGO_REPORT.md` carries a
correction note. No pilot result changes: eligibility was computed with all criteria applied, and
E9's marginal contribution to the eligible pool is +8.1%.

The decision was taken **before** any comparison of the fusion effect under alternative E9
settings. No such comparison was run.

### D3 — Isolated-partner probe states (PILOT-INFORMED)

The probe now emits `MATCHED_INTRINSIC_CAVITY`, `NO_MATCHING_CAVITY`, `DETECTOR_NO_PREDICTIONS`,
`REFERENCE_STRUCTURE_UNAVAILABLE`, `AMBIGUOUS`, using the same correspondence machinery applied to
FUSION residues. `DETECTOR_NO_PREDICTIONS` is never read as biological absence. Re-analysis changed
exactly one label: 6A73/P2Rank moves from `AMBIGUOUS` to `UNDETERMINED_DETECTOR_NO_PREDICTIONS`.
Mechanism counts are otherwise unchanged (21 intrinsic, 4 interface).

### D4 — Zero predictions are data (PILOT-INFORMED)

All six P2Rank zero-prediction cases audited and classified `GENUINE_DETECTOR_ABSTENTION`. No
structure may be removed from the study because a detector predicted nothing. See PROTOCOL.md II.6.

### D5 — Development / confirmatory separation (PILOT-INFORMED)

`data_manifest/development_set.json` (24 structures, permanent) and
`data_manifest/confirmatory_candidate_pool.json` (406 entities, 129 independent partner-target
units) created. The confirmatory pool excludes all development PDB entries and all entities sharing
a (partner, target accession) pair with the development set.

### D6 — Unequal strata are modelled, not equalised (PILOT-INFORMED)

Partner-stratified estimates are primary; any overall estimate is reported both raw-pooled and
partner-standardised. Balanced subsampling is rejected on design grounds — not by comparing effect
sizes under the two schemes. See PROTOCOL.md II.3.

### D7 — Confirmatory cohort uses the whole eligible fresh population (PILOT-INFORMED)

No sampling, no seed, no draw: all 129 independent units are used. This removes sampling discretion
entirely and maximises power. A consequence is that the confirmatory analysis must be run once and
reported in full; there is no second draw available if the result is unwelcome.

---

## Amendment set F — protocol v1.2 → v1.3 (FINAL CONFIRMATORY FREEZE)

**Date:** 2026-09-15
**Outcome status:** **PRE-CONFIRMATORY.** Recorded before any confirmatory detector was executed.
No confirmatory detector output existed at the time of writing; the only confirmatory data seen
were SIFTS annotations, entry metadata and cohort counts. Development results were of course
already known — every change below is therefore also DEVELOPMENT-INFORMED, which is the intended
design of a development/confirmation split.
**Archives:** v1.0, v1.1, v1.2 preserved verbatim as `PROTOCOL_v1.{0,1,2}_ARCHIVED.md`.
**Freeze record:** `environment/PROTOCOL_FREEZE.json` carries SHA-256 hashes of the protocol and
of every frozen manifest, script and test at freeze time.

### F1 — Target-disjoint confirmatory set (pre-confirmatory)

Any confirmatory candidate whose TARGET UniProt accession occurs anywhere in the development set is
removed, even under a different fusion partner. One entity was removed on this rule. Final
confirmatory **n = 128** (BRIL 29, T4L 38, MBP 61), target-accession overlap with development
**= 0**, verified by `tests/test_no_development_leakage.py` (9 checks, all passing).

### F2 — The structure is the inferential unit (pre-confirmatory)

**This withdraws an inferential claim made in the Phase 3A report.** That report presented a
Wilcoxon signed-rank test over 364 individual fpocket matched cavities (p = 4.0e-47) and 86 P2Rank
cavities (p = 5.3e-08). Those tests treated correlated cavities from the same structure as
independent observations — pseudoreplication — and the p-values are **not valid inferential
evidence**. They are retained in the audit trail as descriptive summaries only.

From the freeze onward: every confirmatory hypothesis test, CI and effect estimate uses
n = number of structures. Cavity-level statistics are descriptive only and always report the
number of structures alongside the number of cavities. Secondary cavity-level models must cluster
by target/structure.

### F3 — Non-degenerate causal rank-displacement endpoint (pre-confirmatory)

Primary Aim B is now `RANK_DISPLACEMENT = ORIGINAL_RANK - REMOVED_RANK` for the detector's
highest-ranked cavity in the fusion-removed structure, matched into the original by the frozen
correspondence algorithm. `DETECTOR_ABSTENTION` and `TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL` are
reported as explicit states and never assigned an invented numerical rank. Full definition in
PROTOCOL.md III.3.

The development-phase "rank-1 becomes TARGET_DOMINATED after fusion deletion" outcome is retained
in the audit trail only and is never cited as causal evidence.

### F4 — Within-confirmatory clustering, 128 structures / 123 clusters (pre-confirmatory)

Five targets appear twice, each as a cross-partner pair of the same receptor (beta-2 adrenergic,
M2 muscarinic, adenosine A2A, glucagon receptor, CRF receptor 1). No target is duplicated under the
same partner. They are retained because they are matched within-target partner contrasts and
because the statistical plan already clusters on target accession. Every reported quantity states
both n_structures = 128 and n_independent_clusters = 123. Detected by the leakage test during the
freeze, before any confirmatory run. See PROTOCOL.md III.8.

### F5 — Statistical additions (pre-confirmatory)

Separation check before any logistic regression, with pre-specified Firth penalised logistic
regression as the fallback, or an explicit statement that the conventional model is not estimable.
Only one interaction may be tested — detector x partner — and only if estimable. Modelling strategy
is never selected by which yields smaller p-values. See PROTOCOL.md III.5.

### F6 — Replication criteria, specified before any confirmatory result was inspected

Recorded while the confirmatory run was still executing and **before any confirmatory output file
had been opened or any confirmatory number viewed**. The confirmatory analysis scripts were written
but not yet run.

For each headline signal, comparing the confirmatory point estimate against the DEVELOPMENT 95%
Wilson interval for the same quantity:

- **REPLICATED** — the confirmatory point estimate lies inside the development 95% CI, and the
  direction of effect is the same.
- **PARTIALLY REPLICATED** — the two 95% CIs overlap but the confirmatory point estimate lies
  outside the development CI (same direction, materially different magnitude).
- **NOT REPLICATED** — the 95% CIs do not overlap, or the direction of effect differs, or the
  confirmatory interval is consistent with no effect where development indicated one.

For the displacement endpoint an additional requirement applies to REPLICATED: the confirmatory
cluster-bootstrap CI for "displacement > 0" must exclude zero.

Development and confirmatory cohorts are **never pooled** for confirmatory inference. A combined
descriptive estimate may be computed later for the manuscript, but only after the independent
replication has been reported.

---

## Amendment set G — technical deviations during confirmatory execution

**Outcome status:** these were made **AFTER primary confirmatory results were visible**. Each is an
implementation matter affecting only the sensitivity variants; none changes any scientific rule,
and none touches the frozen primary path.

### G1 — Chain-aware classifier for the assembly1 variant only (POST-RESULT-EXPOSURE, implementation)

**Problem.** The frozen classifier keys residues by `<resnum><icode>` with no chain identifier.
That is unambiguous for the single-chain primary input but collides in a multi-chain biological
assembly: residue 167 of chain B would inherit the class of residue 167 of chain A.

**Fix.** `scripts/21_classify_assembly1.py` parses the chain identifier from each detector's output
and looks up `<chain>:<resnum><icode>`. **The frozen path (scripts/05, scripts/15) is not
modified.** The variant classifier is used only for the assembly1 variant.

**Before/after test.** `tests/test_variant_classifier_equivalence.py` runs both parser pairs over
all 512 archived primary runs and asserts identical pocket counts, ranks, scores, centroids and
residue sets once the chain prefix is stripped: **46,448 assertions, all passing, over 7,656
pockets**. The primary result is therefore provably unaffected.

**Were results already exposed?** Yes. The primary confirmatory Aim A/B/C/E results were computed
and read before this fix was written. The fix cannot alter them — it is not in their code path, and
the equivalence test demonstrates the two parsers agree exactly on the primary inputs.

### G2 — gemmi assembly API correction (POST-RESULT-EXPOSURE, implementation)

The first assembly1 build called `gemmi.make_assembly(...)` with a wrong signature and silently
fell back to the asymmetric unit for all 128 structures; a second attempt used
`Structure.transform_to_assembly` correctly but reported a spurious `IndexError` note because
`st.assemblies` is cleared by the transform. Both were corrected before any assembly1 detector run.
The variant now builds **biological assembly 1** for 126/128 structures.

### G3 — Two structures not evaluable under the assembly1 variant (technical, recorded)

`5W0R` and `5EDU` (both MBP): the designated chain is absent from biological assembly 1, so the
variant cannot be constructed. They remain in the primary analysis and are reported as
not-evaluable for that sensitivity only. They are **not** removed from any other analysis.

### G4 — Chain-aware displacement for the assembly1 variant (POST-RESULT-EXPOSURE, implementation)

**Problem.** `scripts/16` computes correspondence with the frozen `corr.load_target_coords`, which
keys residues `<resnum><icode>`. For the assembly1 variant the pocket residue strings and the class
map are chain-qualified (`A:167`), so the two key spaces never met: every target centroid came back
`None`, every distance was infinite, and **all 126 assembly1 structures were misreported as
TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL with 0 evaluable**. That was an artifact of key formatting,
not a result, and it was obviously so — a 0/126 evaluable rate is not a plausible outcome.

**Fix.** `scripts/22_assembly1_displacement.py` supplies a chain-aware coordinate loader for this
variant only. The correspondence algorithm, its thresholds and the F3 endpoint definition are
unchanged; only the key format differs. Primary and ions variants are untouched and were not
recomputed.

**Were results already exposed?** Yes — the primary confirmatory results had been read. This fix
affects only the assembly1 sensitivity and cannot alter the primary numbers. The erroneous
0-evaluable output is superseded and is not reported anywhere as a finding.

---

## Amendment set H — Phase 5 dataset provenance audit

**Date:** 2026-09-15. This phase does not touch the frozen confirmatory study. No confirmatory
number was recomputed, and no confirmatory rule was consulted or changed.

### H1 — CORRECTION: "fpocket has no training set" was WRONG

The Phase 3A benchmark-feasibility note stated that "fpocket is not an ML method and has no
training set, so no training-contamination claim can attach to it." **That is incorrect and is
withdrawn.** Direct inspection of the built source at tag 4.2.3
(commit `4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066`) shows `src/pscoring.c` implements
`score_pocket()` as a **fitted linear model with hard-coded coefficients** (intercept
-0.03783394 plus five weighted descriptors), documented in-source as "determined using a logistic
regression based on an analysis of pocket descriptors", with adjacent comments recording PLS model
variants and training-set performance ("Train : 62/86 - 65/89"). fpocket's ranking therefore
**does** depend on data-derived parameters. `results/benchmark_audit_feasibility.md` carries a
correction note.

### H2 — Chain-identifier parsing gap in the JOINED membership (corrected)

The first membership parser did not recognise the DT198 naming form `1dug_A.pdb` and recorded
chain = None for 190 of 537 JOINED members. Corrected; JOINED now resolves 190 chain identifiers.
PDB entry identifiers were always parsed correctly, so the entry-level construct scan is
unaffected and was not re-run.

### H3 — sc-PDB membership not obtained: objective access limitation

sc-PDB v.2017 (16,034 entries / 4,782 proteins, frozen PDB data 2016-11) is distributed only as a
single 4,079,179,798-byte archive (`ressources/2016/scPDB.tar.gz`). A measured download-rate test
returned **~11 kB/s**, i.e. roughly four days for the full archive, and no lighter membership index
was located. sc-PDB is therefore recorded as **NOT AUDITED**, with the reason measured rather than
asserted. No prevalence figure is reported for it and none is inferred.

### H4 — LIGYSIS bulk membership not available by design

The LIGYSIS-web repository README states verbatim: **"We do not offer full LIGYSIS dataset
download."** Canonical bulk membership is therefore not obtainable, and LIGYSIS is recorded as
**NOT AUDITED** for prevalence. Its scale is reported from the authors' own documentation
(~65,000 binding sites across ~25,000 proteins, >100,000 PDBe structures).

### H5 — PLINDER audited by reverse intersection, not exhaustive scan

PLINDER 2024-06/v2 has 107,963 PDB entries; an exhaustive SIFTS scan was not attempted. Instead
the recognised-partner chimeric entry set was enumerated from RCSB/SIFTS and intersected with
PLINDER's frozen split membership. This is exact for the twelve partners enumerated and makes no
claim about partners outside that list — stated as a limitation rather than glossed.

---

## Amendment set I — Phase 6 publication package

**Date:** 2026-09-15. No scientific result was recomputed or changed. Phase 6 built tables,
figures and ledgers from already-frozen outputs and audited the numbers quoted in prose.

### I1 — Two denominator errors in the Phase 5 prose, CORRECTED

The canonical denominator audit (`scripts/28_publication_denominators.py`) found two wrong numbers
in the first version of `results/PHASE5_DATASET_AUDIT.md`:

1. **CHEN11 prevalence denominator.** The report used 251, which is the MEMBER count; CHEN11 has
   **241 unique PDB entries** (10 entries contribute two chains) and 242 unique chains. Corrected
   to 0/241 entries = 0.00% [0.00-1.57], with the member-based figure given alongside.
2. **Tier-1 total.** The report said "5,207 resolved entries". That figure omitted FPTRAIN and
   used CHEN11's member count. The correct **non-deduplicated sum is 5,419**; the **deduplicated
   union is 5,304**, since 115 entries appear in more than one Tier-1 dataset.

The PLINDER partner breakdown was also restated in both units (systems and unique PDB entries).
**No conclusion changes: the Tier-1 recognised-partner construct count remains 1 (1DUG).** The
corrected report carries a correction banner, and 11 automated assertions now check every
denominator quoted in that prose against `results/PUBLICATION_DENOMINATORS.tsv`.

### I2 — DECLARED GAP: the pre-specified E9 0.35 sensitivity was never executed

Protocol III.6 pre-specified a confirmatory sensitivity analysis at an E9 segment-disorder
threshold of 0.35, alongside the 0.20 primary. **That analysis was not run.** It would require
rebuilding the cohort (adding roughly 16 independent targets) and executing new detector runs,
which is a new experiment and was not authorised. It is therefore reported as an outstanding
pre-specified analysis, not silently dropped: supplementary figure S6 shows the E9 audit evidence
only and is labelled accordingly, and the claims ledger and readiness assessment both name the gap.
Segment disorder *was* included as a covariate in the secondary model, where it was not
significant.

### I3 — Ordering of the publication freeze

The freeze tag was created at the END of Phase 6 rather than the beginning, so that it captures
the complete package (results plus derived tables, figures and ledgers) in one state. All
scientific results predate Phase 6 and are unchanged; the diff between the confirmatory freeze
(`environment/PROTOCOL_FREEZE.json`) and the publication freeze
(`environment/PUBLICATION_FREEZE.json`) is derived artifacts plus the I1 prose corrections.

### I4 — Figure 5 renderings are C-alpha traces, not molecular-viewer cartoons

No molecular viewer is installed in this environment. Figure 5 panels are drawn from real
deposited coordinates as C-alpha traces with the rank-1 pocket centroid and deposited groups
marked; nothing is reconstructed by hand. `figures/FIG5_pymol.pml` reproduces each panel in PyMOL
for publication-quality rendering, and `figures/FIG5_case_validation.tsv` carries the per-case
validation facts (PDB, chain, boundaries, ligand identity and ownership, pocket membership,
detector and rank) required by 6H.

---

## Amendment set J — E9 = 0.35 disorder sensitivity (the pre-specified analysis declared missing in I2)

**Date:** 2026-09-16
**Outcome status:** the interpretation rule in J2 below was **recorded BEFORE the sensitivity cohort
was built and before any new detector was executed**. No E9 = 0.35 result existed at the time of
writing.

### J1 — Scope

Executes the single pre-specified sensitivity that amendment I2 declared unexecuted: protocol
III.6, segment disorder **0.35** instead of the **0.20** primary. Every other frozen rule is
unchanged — fusion-partner definitions, inclusion/exclusion E1-E8, zero target-UniProt overlap with
the development cohort, one-per-(partner, target) representative by highest resolution with
lowest-PDB-ID tie-break, target clustering, structure preparation, detector versions and
parameters, the fusion-associated definition (`FUSION_DOMINATED | INTERFACE | LINKER`), dominance
0.70, and the target-cavity correspondence rule (>= 3 target residues, Jaccard >= 0.40, centroid
<= 8.0 A, one-to-one, rank and score excluded from matching). The inferential unit remains one
structure / one independent target; no pocket-level pseudoreplication.

**All primary manuscript numbers remain the E9 = 0.20 confirmatory results. The 0.35 analysis is a
sensitivity analysis only and never replaces a primary value.**

### J2 — Interpretation rule, frozen before results

- **ROBUST** — the direction of every major conclusion is unchanged; MBP remains the strongest
  rank-1 partner for both detectors; construct-associated top-k exposure persists; displacement
  remains strongly directional and non-negative; no central claim must change.
- **QUANTITATIVELY SENSITIVE** — magnitudes change noticeably, but interpretation and the claim
  hierarchy are unchanged.
- **MATERIALLY CHANGED** — partner ordering changes materially, a central exposure conclusion
  reverses, displacement loses its directional pattern, or a locked paper claim becomes
  unsupported.

No numerical cutoff may be invented after seeing the result.

### J3 — Reuse of frozen detector output

Structures already present in the E9 = 0.20 primary cohort are **not re-run**; their frozen
detector outputs are reused unchanged. Only newly admitted structures are executed, at
2 conditions x 2 detectors.

### J4 — Watch condition

Relaxing a threshold should only add observations. If any primary structure is absent from the
0.35 cohort, that is investigated and reported before the analysis proceeds — it would indicate the
representative-selection rule swapping a representative rather than a true loss.

### J5 — Provenance preserved

`v1.0-paper-analysis-freeze` and `environment/PUBLICATION_FREEZE.json` are left unchanged as
historical provenance, including the known DEVIATIONS.md ordering mismatch documented in I3. The
final state is captured by a NEW tag and a NEW manifest.

### J6 — Watch condition J4 fired, investigated, benign: three representative UPGRADES

Building the 0.35 cohort triggered J4: three PDB entries present in the E9 = 0.20 primary cohort
(`4ZUD`, `5N2S`, `6C1Q`) are absent from the 0.35 cohort. Investigated before any analysis was run.

**Cause — not a loss of observations.** All three primary representatives remain **eligible** at
0.35. Each was displaced, under the unchanged frozen representative rule ("one structure per
(partner, target): highest resolution, tie-break lowest PDB ID"), by a higher-resolution entry for
the *same* target that had been excluded at 0.20 solely by `E9_fusion_disorder`:

| Unit (partner, target) | Primary rep | Challenger admitted at 0.35 | Challenger's 0.20 exclusion |
|---|---|---|---|
| BRIL, P21730 | 6C1Q, 2.90 A | **6C1R, 2.20 A** | `E9_fusion_disorder = 0.24` |
| BRIL, P30542 | 5N2S, 3.303 A | **5UEN, 3.20 A** | `E9_fusion_disorder = 0.26` |
| BRIL, P30556 | 4ZUD, 2.80 A | **6OS2, 2.70 A** | `E9_fusion_disorder = 0.24` |

Every swap moves to strictly better resolution. **No target cluster is lost**: the 0.35 cohort's
136 independent targets are a strict superset of the primary's 123. Zero entities eligible at 0.20
become ineligible at 0.35.

**Decision.** The frozen rule chain is applied exactly as written, including the representative
rule — that is what "the analysis as pre-specified, with disorder at 0.35" means, and freezing the
primary representatives instead would produce a hybrid cohort that no pre-specified rule generates.
Consequently the sensitivity cohort is not a strict superset of primary *structures* (it is of
*targets*), and 16 structures require new detector runs: 13 newly admitted units plus the 3
upgraded representatives. The other 125 primary structures reuse their frozen detector output
unchanged and are not re-run.

This is a consequence of applying the frozen rules, not a change to them.

### J7 — Result and classification: ROBUST

Judged against the rule frozen in J2 **before** the cohort was built.

| | primary 0.20 | sensitivity 0.35 | delta |
|---|---|---|---|
| cohort | 128 structures / 123 targets | 141 / 136 | +13 units |
| P2Rank rank-1 | 66/127 = 52.0% | 73/140 = 52.1% | +0.1 |
| fpocket rank-1 | 60/128 = 46.9% | 69/141 = 48.9% | +2.0 |
| P2Rank standardised | 39.7% | 40.3% | +0.6 |
| fpocket standardised | 36.2% | 38.8% | +2.6 |
| P2Rank MBP rank-1 | 93.4% | 94.0% | +0.6 |
| fpocket MBP rank-1 | 85.2% | 86.6% | +1.4 |
| P2Rank displacement > 0 | 42/103 = 40.8% | 47/114 = 41.2% | +0.5 |
| fpocket displacement > 0 | 65/128 = 50.8% | 76/141 = 53.9% | +3.1 |
| negative displacement | 0 and 0 | **0 and 0** | none |

**Verdict: ROBUST.** Every J2 criterion is met — direction unchanged; MBP remains the strongest
rank-1 partner in both detectors (94.0% and 86.6%, against 8.8-17.9% for BRIL and T4L);
construct-associated top-k exposure persists (top-5 80.0% and 83.0%); displacement remains strictly
directional with **zero** negative values across all 255 evaluable structure-detector pairs
(sign tests 47/0, p = 1.4e-14 and 76/0, p = 2.7e-23); and no central claim requires change.

**One observation reported rather than buried.** For fpocket rank-1 the BRIL and T4L strata swap
places: primary T4L 13.2% > BRIL 10.3%; sensitivity BRIL 17.1% > T4L 12.8%. This is **not** a
material ordering change under J2. The two intervals overlap almost completely in both analyses
(primary BRIL [3.6-26.4] vs T4L [5.8-27.3]; sensitivity BRIL [8.1-32.7] vs T4L [5.6-26.7]); the
material ordering — MBP far above both — is unchanged in both detectors; and no claim in the ledger
ever ranked BRIL against T4L. The secondary model had already reported T4L as indistinguishable
from BRIL (OR 2.40, p = 0.30). The largest single movement in the whole sensitivity is fpocket's
BRIL stratum (+6.8 points rank-1, +9.1 points displacement), driven by 6 newly admitted BRIL
structures, 4 of them cryo-EM.

**CENTRAL CLAIM CHANGED: NO.** Primary manuscript values remain the E9 = 0.20 confirmatory results.

### J8 — Artifacts updated for the sensitivity, and what was deliberately NOT touched

Updated: supplementary figure **S6** (now the primary-vs-sensitivity comparison, replacing the
audit-only version), **Table 6** and **Table 6b** (new), the sensitivity fields of
`CLAIMS_LEDGER.md` including new standing constraint 7, `results/READINESS.md` (new), and the two
manuscript-outline lines that referred to an unexecuted sensitivity.

Deliberately unchanged: every primary confirmatory value; Tables 1-5 and Figures 1-6; the Phase 5
dataset audit; the development and confirmatory analyses; `v1.0-paper-analysis-freeze` and
`environment/PUBLICATION_FREEZE.json`, preserved as historical provenance including the I3
ordering mismatch. Mechanism (Aim D) was not recomputed for the 0.35 cohort — outside the scope of
this sensitivity, and stated as such in the ledger.

### J9 — Final publication freeze `v1.1-paper-final-freeze`, and the ordering rule it obeys

**Type: PROCEDURAL. No scientific content.** This is the last entry written before the final
manifest is generated, and the last edit made to this file.

`v1.0-paper-analysis-freeze` recorded a manifest whose `DEVIATIONS.md` hash does not verify,
because amendment set I was appended *after* `environment/PUBLICATION_FREEZE.json` was written
(deviation **I3**). That manifest is **not** regenerated to repair the mismatch. Rewriting a freeze
to make it agree with itself would destroy the only evidence that the mismatch happened, and the
mismatch is itself part of the record. v1.0 is preserved exactly as it was published: same tag,
same commit `e0220f9276f69c0defc6a686e3696a5f2830205b`, same manifest bytes.

`environment/FINAL_FREEZE_v1.1.json` is therefore a **new** file, generated only after every
artifact it covers already existed — the E9 cohort and run files, J1-J9 of this log, the
regenerated S6, Tables 6 and 6b, the regenerated claims ledger, the corrected manuscript outline,
and `results/READINESS.md`. Its coverage is enumerated from `git ls-files` rather than from a
curated list, so no tracked file can be silently omitted; it is a strict superset of the 100 files
v1.0 hashed. Nothing it covers is edited afterwards.

A manifest cannot contain the hash of the commit that introduces it. `v1.1-paper-final-freeze`
points at the commit that adds the manifest; the manifest records the commit it was generated
from, and the only difference between the two is the manifest file itself. Every hashed file is
byte-identical across both.

**Scientific state at v1.1: identical to v1.0 for every primary value.** v1.1 adds one executed
sensitivity analysis and the artifacts describing it. It changes no primary result, no claim
direction, and no conclusion.

### J10 — Line-ending normalisation makes single-hash manifests platform-dependent

**Type: PROCEDURAL, found while generating the v1.1 manifest. No scientific content.**
**Supersedes one sentence of J9** — J9's claim that it was the last edit to this file. The
substantive ordering rule in J9 is unaffected: no manifest had been generated when this entry was
written, so the manifest still comes last and nothing it covers is edited after it.

**Finding.** This repository is configured with `core.autocrlf = true` and has no `.gitattributes`.
Of 224 tracked files, **128** are stored in git with LF and checked out on Windows with CRLF; 59
are LF in both; 33 are binary; 4 are empty. A byte-for-byte comparison confirms the divergence is
**entirely** line endings — zero files differ in content.

**Consequence.** A manifest that records one SHA-256 per file, computed on working-tree bytes,
verifies on Windows and fails on Linux or macOS for those 128 text files, for a reason that has
nothing to do with the data. `environment/PUBLICATION_FREEZE.json` (v1.0) has exactly this
property. It is not regenerated — see J9.

**Action.** No file is altered, no `.gitattributes` is added, and no working tree is renormalised;
doing any of those at freeze time would rewrite frozen artifacts to fix a cosmetic problem.
Instead `environment/FINAL_FREEZE_v1.1.json` records **two** hashes for every text file:

- `sha256` — the bytes on disk as frozen on this machine (Windows, CRLF where applicable);
- `sha256_lf` — the same content with CRLF normalised to LF, which is the platform-independent
  identity and is what a Linux or macOS checkout will reproduce.

Binary files carry a single hash. `sha256_lf` is the hash a third party should verify against.
This makes the freeze checkable on any platform without touching a single frozen byte, and it
declares a defect that v1.0 carried silently.
