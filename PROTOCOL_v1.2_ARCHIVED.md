# Frozen Study Protocol (v1.2)

**Project:** Fusion-Tag Hazard
**Frozen:** 2026-09-15, before any pocket detector has been run on any structure.
**Amended:** 2026-09-15 (v1.0 → v1.1), amendment set A in `DEVIATIONS.md`, recorded
**PRE-OUTCOME** — before any detector was installed or executed. Protocol v1.0 preserved verbatim
in `PROTOCOL_v1.0_ARCHIVED.md`. Changes in v1.1 are marked **[A1]**–**[A5]** inline.
**Amendment rule:** any change after this date is logged in `DEVIATIONS.md` with date, reason, and
whether the change was made before or after the analyst saw outcome data. Changes made after seeing
full-cohort results are permitted only as clearly labelled post-hoc analyses and may not be
reported as confirmatory.

---

## 1. Refined primary hypothesis

> In deposited structures of chimeric crystallization constructs, the fusion partner and the
> engineered target-fusion interface compete with, and frequently outrank, the intended target
> site in automated pocket prediction. Deleting the fusion residues from the same coordinate file
> causally restores the rank of the target site.

Two deliberate refinements to the hypothesis as originally posed:

**(a) The causal, within-structure comparison is primary; the between-structure comparison is
secondary.** A power calculation (Section 8.4) shows that the unpaired fusion-vs-native contrast
needs ~100-320 structures per arm for realistic effect sizes, and it confounds construct
composition with target identity, resolution, and experimental method. The paired
original-vs-fusion-removed design reaches adequate power at n = 60-150 and controls all of those by
construction. The paired design is the scientific core of the paper.

**(b) The claim is about pipelines that consume PDB coordinates without construct-aware
preprocessing** — large-scale pocket databases, ML training/benchmark sets, proteome-wide
druggability scans, and naive user workflows — not about careful practitioners who already delete
the fusion by hand. This scoping is what keeps the claim non-trivial and honest.

### Secondary questions (pre-specified)

- S1. Which fusion partner is most problematic (BRIL vs T4L vs MBP)?
- S2. Does the effect replicate across independent detectors?
- S3. Are artifact pockets the partner's own native ligand site, or newly created by the chimera?
- S4. Do high-ranking artifact pockets concentrate in the fusion body, the interface, or the linker?
- S5. Does the effect depend on fusion topology (internal insertion vs terminal fusion), fusion
  size, buried surface area at the target-fusion interface, resolution, or target class?
- S6. Does in-silico fusion removal restore target-pocket rank? (This is also the primary endpoint
  of the full study.)
- S7. **[A5] ADOPTED as a pre-specified second major study aim** (not a sub-question): do existing
  pocket-prediction datasets contain chimeric crystallization constructs? **Deferred to the phase
  after the pilot passes its gate; not performed during the pilot.** Binding language constraints:
  (i) no dataset is called "defective," "broken," or "invalid"; (ii) TRAINING, VALIDATION,
  TEST/BENCHMARK and REFERENCE datasets are distinguished throughout and never conflated; (iii) no
  claim that a predictor was *trained* on fusion structures unless that model's actual documented
  training set establishes it; (iv) mere **presence** of chimeric structures is reported as an
  observation and is never equated with **demonstrated bias** in reported performance.

---

## 2. Unit of analysis

One **(PDB entry, chimeric polymer entity, representative chain)** triple. One representative chain
per entity, chosen deterministically as the first `auth_asym_id` in ASCII sort order.

**Clustering unit for all inference:** the target UniProt accession. Structures sharing a target
accession are never treated as independent.

---

## 3. Cohort definition

### 3.1 Fusion partners (frozen list)

| Partner | UniProt | Rationale | Chimeric entities in PDB (measured 2026-09-15) |
|---|---|---|---|
| BRIL / apocytochrome b562 | P0ABE7 | Apo heme site is a vacated cavity by construction | **522** in 522 entries |
| MBP | P0AEX9 | Native maltose cleft; dominant soluble-protein fusion | **375** in 371 entries |
| T4 lysozyme | P00720 | Active-site cleft plus L99A-type hydrophobic cavity | **154** in 153 entries |

Union across the three: **1,041 unique PDB entries** (pairwise overlap is negligible: T4L/BRIL = 2,
BRIL/MBP = 3, T4L/MBP = 0).

**Deferred candidates**, measured but excluded from v1.0 to keep the partner set interpretable:
GST-Sj26 (P08515, 46 chimeric), rubredoxin (P00268, 21), thioredoxin (P0AA25, 10), GFP (P42212,
217 — but overwhelmingly a localisation tag rather than a crystallization chaperone, so a different
scientific question). Flavodoxin (P0ABY4) has 0 chimeric entities and is dropped. Any addition
after freeze goes in `DEVIATIONS.md`.

### 3.2 Inclusion criteria

1. Deposited X-ray or cryo-EM structure; not NMR, not a theoretical model; not obsolete/superseded.
2. Resolution <= 3.5 A (primary); sensitivity analysis at <= 3.0 A.
3. Exactly one chimeric polymer entity that maps to **exactly two** UniProt accessions via SIFTS,
   one of which is a frozen-list fusion partner.
4. SIFTS aligned regions cover >= 80% of the entity's modelled residues.
5. Fusion segment: >= 50 modelled residues (a domain small enough to be unmodelled cannot host a
   pocket; this also excludes vestigial tag remnants).
6. Target segment: >= 80 modelled residues.
7. <= 20% of either segment unmodelled/disordered within its SIFTS range.

### 3.3 Exclusion criteria

1. Entity maps to >= 3 UniProt accessions (fusion boundaries ambiguous).
2. Standalone structures of the partner itself (single-accession entities) — these form the
   partner-reference set for mechanistic analysis (S3), not the cohort.
3. Overlapping or inconsistent SIFTS segments for the two accessions.
4. Fusion partner is the *biological* subject of the study rather than a crystallization aid
   (flagged by manual review of the entry title and primary citation; logged with reason).

### 3.4 Fusion-boundary annotation — the evidence trail

Boundaries are taken from **SIFTS residue-level alignments** served by the RCSB Data API
(`rcsb_polymer_entity_align.aligned_regions`, `provenance_source = SIFTS`). Verified working on
representative entries:

```
5IU7_1  P29274 (A2A receptor)  entity 11-217, 324-420   |  P0ABE7 (BRIL) entity 218-322   -> internal insertion
3ODU_1  P61073 (CXCR4)         entity 12-239, 404-493   |  P00720 (T4L)  entity 242-401   -> internal insertion
3L2J_1  P0AEX9 (MBP)           entity 1-364             |  Q03431 (PTH1R) entity 371-529  -> N-terminal fusion, linker 365-370
```

Derived deterministically from these segments:

- **TARGET residues** = entity positions inside a target-accession aligned region.
- **FUSION residues** = entity positions inside a fusion-partner aligned region.
- **LINKER residues** = modelled entity positions in no aligned region and *flanked on both sides*
  by different accessions.
- **TAG residues** = modelled entity positions in no aligned region at either terminus
  (expression-tag remnants). Recorded separately; not counted as linker.
- **Fusion topology** = INTERNAL_INSERTION if the target has >= 2 aligned regions bracketing the
  fusion; TERMINAL_N or TERMINAL_C otherwise.

No fusion boundary is ever hand-drawn. Every annotation carries its SIFTS provenance and the
API-response hash. If SIFTS is missing or self-inconsistent for an entry, the entry is **excluded
and logged**, never manually patched.

### 3.5 Redundancy control

Measured redundancy in the candidate pool is severe and must be controlled:

- BRIL: 140 distinct target accessions; adenosine A2A receptor (P29274) alone accounts for
  128/522 = **25%** of the BRIL cohort.
- T4L: 62 distinct targets; beta2-adrenergic receptor (P07550) = 26/154 = 17%.
- MBP: 195 distinct targets; most even, top target 23/375 = 6%.

**Rule:** the primary cohort takes **one structure per (fusion partner, target accession)** pair —
the highest-resolution entry, ties broken by lowest PDB ID. All structures are still processed; the
multi-structure set is used only for a pre-specified sensitivity analysis with clustered inference.

---

## 4. Structure preparation (fixed; no per-structure tuning)

1. Fetch mmCIF from RCSB. **Never redistribute coordinates** — the repository stores IDs, SIFTS
   annotations, and derived results only.
2. Primary input = **the single designated chimeric chain**, first model, altloc A only.
   Rationale: isolates the fusion-vs-target question and makes the paired removal exact.
   Sensitivity analysis: biological assembly 1.
3. **Strip all non-polymer entities** — ligands, waters, ions, lipids, detergents, cryoprotectants.
   Otherwise fpocket and P2Rank output depends on which HETATMs the depositor happened to model,
   which is itself a crystallization artifact and would confound the exposure of interest.
   Sensitivity analysis: retain ions only.
4. No protonation, no minimisation, no repacking. Detectors receive coordinates as deposited.
5. **[A1] TAG residues are deleted from BOTH conditions** before pocket detection. The paired
   contrast is therefore exactly:

   | Condition | Residues present |
   |---|---|
   | ORIGINAL (fusion-present) | target + fusion + linker |
   | FUSION-REMOVED | target only |

   Terminal expression-tag remnants are absent from both arms, so the fusion domain is the single
   manipulated variable. Deleted TAG residue IDs are recorded per structure in the manifest field
   `tag_residues_deleted`.
6. **Fusion-removed counterpart:** identical file with all FUSION and LINKER residues deleted (TAG
   already absent from both per A1). Nothing else changes — no renumbering, no chain-break repair,
   no minimisation, and **no TARGET atom coordinate is altered** (asserted by automated test
   `test_target_coords_unchanged`). Any chain break introduced is left as-is and recorded.
7. **[A3] Deletion-boundary recording.** For every fusion-removed structure, record the residues
   flanking each excised span (`deletion_boundaries`: last target residue before, first target
   residue after, with CA coordinates). For every pocket predicted in the REMOVED condition,
   compute `dist_to_deletion_boundary` = Euclidean distance from pocket centroid to the nearest
   boundary CA.

---

## 5. Detectors

**Frozen for the pilot:**

- **P2Rank** — default model, `prank predict`, stock parameters. Java 11+.
- **fpocket** — stock defaults (fpocket 4.x).

Both are run identically on every structure in both conditions. No per-structure parameter changes,
ever. Versions, exact command lines, and output hashes are recorded per run.

**Third detector:** deferred to the full study, and added only if it strengthens the science rather
than the scope. Preferred candidate is an **ML method trained on PDB-derived data** (GrASP or
VN-EGNN), because that directly probes the training-contamination hypothesis (S7). DeepPocket and
PRANK are explicitly *not* independent third detectors — both rescore fpocket output.

**Declared installation risk:** fpocket has no supported native Windows build (bioconda ships
linux/osx only). It will be run under WSL2 or a Linux container. If it cannot be installed
reproducibly, that fact is **reported**, not silently worked around by substituting another tool.

---

## 6. Pocket classification rules (pre-specified, before any cohort results)

For each predicted pocket, take the detector's own reported pocket-residue set, map each residue to
TARGET / FUSION / LINKER / TAG via Section 3.4, and compute fractions over the residues that carry
a class (`f_target`, `f_fusion`, `f_linker`).

| Class | Rule |
|---|---|
| `TARGET_DOMINATED` | `f_target >= 0.70` |
| `FUSION_DOMINATED` | `f_fusion >= 0.70` |
| `LINKER` | `f_linker >= 0.50` |
| `INTERFACE` | neither dominance rule met **and** `f_target >= 0.20` **and** `f_fusion >= 0.20` |
| `MIXED` | any other assigned pocket |
| `UNASSIGNED` | fewer than 5 classifiable residues |

**Why 0.70 and not 0.50.** A simple majority rule makes "fusion-dominated" trivially easy to reach:
a pocket that is 51% fusion residues is a genuinely mixed object, and calling it fusion-dominated
inflates the headline rate — exactly the criticism a reviewer would raise. A 0.70 threshold
requires a clear plurality, and it leaves a real INTERFACE band (both sides >= 0.20) rather than
forcing every interface pocket into one camp. The 0.20 floor for INTERFACE prevents a pocket with
two incidental residues from the far side being labelled an interface pocket.

**Pre-specified sensitivity analyses on this choice:** thresholds 0.50 and 0.90; and an alternative
residue definition using a 5 A shell around pocket alpha-sphere centres instead of the detector's
reported residue list.

**"Fusion-associated"** = `FUSION_DOMINATED` OR `INTERFACE` OR `LINKER`. This is the adverse class
throughout.

---

## 7. Endpoints

### 7.1 Pilot primary endpoint (descriptive)

Proportion of pilot structures whose **rank-1** predicted pocket is fusion-associated, reported
separately for P2Rank and fpocket, with Wilson 95% CIs.

### 7.2 Full-study primary endpoint (confirmatory, causal) — **[A2] coverage-conditional**

Paired, within-structure, tested with an exact McNemar test, per detector. **Which quantity is
primary is decided by reference-site coverage, by a rule frozen here before any detector result was
viewed:**

- **Coverage ≥ 60%** of the confirmatory cohort → the primary endpoint **may** be
  **P(biological target reference site is the rank-1 pocket)**, removed vs original.
- **Coverage < 60%** → the primary endpoint **remains** the cohort-level fusion-capture / ranking
  endpoint, i.e. **P(rank-1 pocket is fusion-associated)**, removed vs original; biological
  target-site recovery drops to a **restricted-subset secondary** endpoint.

This supersedes the v1.0 wording, which named reference-site restoration as primary while §7.4
simultaneously forbade promoting it below 60% coverage. For the pilot, coverage is **reported** and
pilot selection is **not** altered on the basis of ligand availability.

### 7.3 Secondary endpoints

- Any fusion-associated pocket in top-3 and in top-5.
- Median rank of the first fusion-associated pocket.
- Score margin: best fusion-associated pocket score / best target-dominated pocket score.
- Change in rank of the best target pocket, original vs fusion-removed (paired, signed).
- Per-partner stratified rates (S1).
- Between-detector agreement on rank-1 class (Cohen's kappa) (S2).
- Fraction of FUSION_DOMINATED pockets overlapping the partner's own native ligand site, using
  standalone partner structures as reference (S3).
- Distribution across fusion body / interface / linker (S4).
- Covariate model: topology, fusion length, interface BSA, resolution, method, target class (S5).
- **[A3] Deletion-boundary sensitivity.** Primary paired analysis uses **all** predictions. The
  pre-specified sensitivity analysis disregards *newly appearing* REMOVED-condition pockets whose
  centroid lies **<= 8 A** from a deletion boundary, when interpreting apparent target-site rescue.
  The 8 A threshold is fixed before results and is not re-tuned. Purpose: separate genuine
  restoration of ranking from cavities generated directly by the in-silico cut.
- Benchmark-set construct presence (S7) — **deferred to the post-pilot phase**, not computed here.

### 7.4 Target reference site

Where a biologically relevant ligand is bound to the target segment, the reference site is the set
of target residues within 4.0 A of that ligand, with relevance judged by a pre-declared
BioLiP-style artifact exclusion list (buffers, cryoprotectants, detergents, additives, ions).
Ligand-recovery endpoints are computed **only** on the subset with a defined reference site.

**[A2] Reference sites are defined without inspecting pocket-prediction outcomes.** The definition
uses deposited ligand coordinates and the exclusion list only; no detector output is consulted.
Promotion of reference-site recovery to primary is governed solely by the 60% coverage rule in
§7.2.

---

## 8. Statistics

1. **Proportions:** Wilson 95% CIs for description; **clustered bootstrap by target UniProt
   accession** (10,000 resamples) for all inferential CIs.
2. **Paired binary outcomes:** exact McNemar. Effect reported as a risk difference with CI, not
   only a p-value.
3. **Paired rank changes:** Wilcoxon signed-rank, with the median paired rank change and a
   bootstrap CI.
4. **Covariate models:** mixed-effects logistic regression, random intercept for target accession.
5. **Multiplicity:** Holm correction within the family of secondary tests (3 partners x 2 detectors).
   The two primary endpoints are not corrected against each other; they are reported per detector
   and both must be stated regardless of outcome.
6. **Never** treat near-identical structures as independent; n is always reported alongside the
   number of independent target clusters.

### 8.4 Pre-specified power (computed 2026-09-15, simulation and closed form)

Paired McNemar, alpha = 0.05:

| n pairs | power, strong effect (p01=.25, p10=.05) | modest (.15/.05) | weak (.10/.05) |
|---|---|---|---|
| 24 | 0.44 | 0.18 | 0.08 |
| 60 | 0.86 | 0.40 | 0.15 |
| 120 | 0.99 | 0.72 | 0.30 |
| 150 | 1.00 | 0.81 | 0.37 |

Unpaired two-proportion comparison, 80% power, per group, before clustering inflation: 30% vs 10%
needs 62; 25% vs 10% needs 100; 20% vs 10% needs 199. With a design effect of ~1.6 for repeated
targets these become 100, 160, and 319 per group. **This is the quantitative reason the paired
design is primary.**

Single-rate precision (Wilson, observed p = 0.35): n = 24 gives a CI width of 0.35; n = 120 gives
0.17; n = 150 gives 0.15.

**Minimum sample size targets:** pilot n = 24 (8 per partner); full study n = 120-150 paired
structures, stratified by partner, subject to the one-structure-per-target rule. The measured pool
(1,041 entries, 397 distinct target accessions across the three partners) comfortably supports this.

---

## 9. Pilot design

- **n = 24**: 8 BRIL, 8 T4L, 8 MBP.
- **Selection (reproducible, non-cherry-picked):** enumerate all chimeric entities per partner from
  SIFTS; apply Section 3.2/3.3 filters; collapse to one structure per target accession (highest
  resolution, ties by lowest PDB ID); sort the surviving target accessions lexicographically; draw
  8 per partner with a fixed seed (`20260915`). No structure is inspected before selection, and no
  structure is swapped out after selection except for a documented technical failure.
- Both conditions (original, fusion-removed) x both detectors = 96 runs.
- The manifest (`data_manifest/pilot_manifest.tsv`) records PDB ID, chain, entity, target accession
  and name, target residue ranges, fusion partner, fusion residue range, linker residues, topology,
  resolution, method, bound ligands, and the SIFTS provenance hash.

---

## 10. Pilot go/no-go gate — **[A4]**

The numerical gates G1/G2/G3 below are **retained unchanged** from v1.0.

**The 10% figure is not an empirically established statistical null.** No study has measured the
rank-1 fusion-capture rate — that is the point of this work. Accordingly:

- **>= 6 of 24 is a pre-specified minimum practically interesting pilot signal and decision
  threshold**, not a rejection region.
- The binomial arithmetic against 10% (>= 6/24 → p = 0.028; >= 5/24 → p = 0.085) is reported as
  **descriptive context only**. It must never be presented as confirmatory hypothesis evidence, as
  a p-value for the scientific hypothesis, or as evidence against a real null.
- The pilot is a feasibility and signal-detection exercise and generates **no confirmatory
  inference**.

**GO** if any of:

- G1. Rank-1 fusion-associated rate >= 6/24 in one detector **and** >= 4/24 in the other.
- G2. Any single partner shows a rank-1 fusion-associated rate >= 50% (>= 4/8) in both detectors.
- G3. Paired removal flips the rank-1 pocket to the target site in >= 5/24 structures, with <= 1
  flip in the opposite direction.

**NO-GO** otherwise. Code already written is not a reason to continue.

If NO-GO, the pre-committed fallback question is the benchmark-contamination audit (S7) as a
standalone short paper: *"Engineered crystallization constructs in ligand-binding-site benchmark
sets"* — which does not depend on the artifact being highly ranked, only on it being present.

---

## 11. Claim boundaries (binding on the manuscript)

Will **not** claim: that all fusion tags cause artifacts; that pocket predictors are generally
unreliable; that any published docking result is invalid; that any biological ligand site is wrong;
or any causal claim beyond what the paired removal experiment supports.

Will claim, at most: *certain crystallographic fusion partners can introduce or promote
high-ranking predicted pockets that compete with the intended target site, producing a systematic
structural-bioinformatics failure mode that should be screened for before pocket-based CADD
workflows and before assembling pocket-prediction datasets.*

---

## 12. Immutability

Original mmCIF downloads are read-only and checksummed. Derived files are written to separate
paths. No script mutates an input. Detector versions, command lines, seeds, and output hashes are
recorded for every run.

---
---

# PART II — CONFIRMATORY PROTOCOL (v1.2)

**Frozen 2026-09-15, before any detector has been run on any confirmatory structure.**

Everything above (v1.1) governed the **development/exploratory pilot**. Part II governs the
**confirmatory study** and is informed by that pilot. The 24 pilot structures are permanently
development data (`data_manifest/development_set.json`) and may never appear as confirmatory
observations. Amendment set **D** in `DEVIATIONS.md` records every change and marks each as
PILOT-INFORMED.

## II.1 Target-pocket correspondence (replaces the degenerate paired outcome)

**Problem.** The pilot's paired outcome — "rank-1 becomes TARGET_DOMINATED" — cannot demonstrate
restoration of a *specific* cavity: once FUSION/LINKER residues are deleted, essentially every
predicted pocket is composed of TARGET residues. It is retained for audit only and carries no
causal weight.

**Frozen correspondence rule.** For each (structure, detector), pockets are matched across the
ORIGINAL and FUSION-REMOVED conditions by:

1. Restrict every pocket to its **TARGET residues** only.
2. Discard pockets with fewer than **3** target residues (not matchable).
3. Similarity = **Jaccard index on target-residue sets**.
4. Geometry = Euclidean distance between the **target-residue heavy-atom centroids**. Target
   coordinates are byte-identical between conditions (test I5), so this is directly comparable.
5. **One-to-one assignment** by `scipy.optimize.linear_sum_assignment` maximising total Jaccard.
6. A pair is MATCHED if **Jaccard >= 0.40 AND centroid distance <= 8.0 A**.
7. States: `MATCHED`, `UNMATCHED_ORIGINAL`, `NEW_AFTER_REMOVAL`.

**Detector rank and score are never inputs to correspondence.**

**Sensitivity rules (pre-specified):** J >= 0.25 with 12 A; J >= 0.40 with 8 A (primary);
J >= 0.60 with 5 A; plus an overlap-coefficient variant at >= 0.50.

**Why these thresholds, and why they cannot be accused of tuning.** On the development set the
matched pair count is *flat* across Jaccard 0.20-0.60 (P2Rank 86 pairs at every threshold;
fpocket 367 falling only to 363), and matched pairs have a **median Jaccard of 1.000 and a median
centroid distance of 0.00 A** — the same cavity is re-detected with an identical residue set. The
threshold therefore has almost no discretion to exercise; 0.40 sits in the middle of a plateau.

## II.2 Endpoint hierarchy (frozen)

### A — EXPOSURE / HAZARD (co-primary)

In the **ORIGINAL deposited construct**, the proportion of structures with a fusion-associated
pocket (FUSION_DOMINATED, INTERFACE or LINKER) at **rank 1**, in the **top 3**, and in the
**top 5**. Reported per detector and per partner. Rank-1 is the headline; top-3 and top-5 are
pre-specified, not post-hoc.

### B — PAIRED RANK INTERFERENCE (co-primary, causal)

Using II.1 correspondence:

- **B1 — DISPLACEMENT (primary binary endpoint).** Per structure, take the matched target cavity
  with the best rank in the REMOVED condition. Was it outranked in the ORIGINAL condition by at
  least one fusion-associated pocket? Analysed by exact McNemar against the removed condition and
  reported as a risk difference with a clustered-bootstrap CI.
- **B2 — PAIRED RANK CHANGE (primary continuous endpoint).** orig_rank minus removed_rank over all
  matched cavities; Wilcoxon signed-rank; median change with clustered-bootstrap CI.
- Structures where the detector predicts **no** pockets after removal are reported in an explicit
  `NO_POCKET_PREDICTED` stratum and are **never** counted as improvements.

### C — BIOLOGICAL TARGET-SITE SUBSET (secondary, restricted)

Restricted to structures with an outcome-independent biological reference site: the same B1/B2
quantities computed on the matched cavity that overlaps that site. Promotion of C to primary
requires at least 60% coverage (rule A2); development coverage was 54%, so **C is secondary**
unless the confirmatory cohort's coverage reaches 60%, measured before any detector is run on it.

### D — MECHANISM (secondary)

Each rank-1 fusion-associated pocket is assigned a probe state — `MATCHED_INTRINSIC_CAVITY`,
`NO_MATCHING_CAVITY`, `DETECTOR_NO_PREDICTIONS`, `REFERENCE_STRUCTURE_UNAVAILABLE`, `AMBIGUOUS` —
using the same correspondence rule applied to FUSION residues against the isolated-partner
structure. `DETECTOR_NO_PREDICTIONS` is **uninformative** and must never be read as biological
absence. Mechanism classes: `FUSION_INTRINSIC_CAVITY`, `CHIMERA_SPECIFIC_FUSION_CAVITY`,
`TARGET_FUSION_INTERFACE`, `LINKER_RELATED`, `MAPPING_ARTIFACT`, `UNDETERMINED_<state>`.

### E — DETECTOR DEPENDENCE (secondary)

Rank-1 and top-k rates, score margin, and Cohen's kappa on the rank-1 call, per partner and per
rank depth, for P2Rank and fpocket. No third detector.

## II.3 Statistical model (frozen)

**Unequal strata are modelled, not equalised.** MBP contributes 61 of 129 independent confirmatory
targets; BRIL 30; T4L 38.

1. **Partner-stratified estimates are the primary presentation.** BRIL, T4L and MBP are reported
   separately with Wilson and clustered-bootstrap CIs. There is no single headline number that
   hides them.
2. **Any overall estimate is reported twice:** the raw pooled proportion, and a **pre-specified
   partner-standardised estimate** giving each partner equal weight (1/3 each). Both are always
   shown together. Applying development rates to confirmatory stratum sizes, these differ by
   **12.5 points for P2Rank** (58.3% pooled vs 45.8% standardised) and 4.5 points for fpocket —
   large enough that reporting only one would be misleading.
3. **Model.** Mixed-effects logistic regression, outcome = fusion-associated rank-1; fixed effects
   = fusion partner, topology, target size, resolution, method, detector; random intercept for
   target UniProt accession. Paired endpoints use exact McNemar (binary) and Wilcoxon signed-rank
   (continuous), with clustered bootstrap by target accession for all CIs.
4. **Multiplicity.** Holm correction within the secondary family. Co-primary endpoints A and B are
   both reported unconditionally regardless of outcome.
5. Effect sizes with confidence intervals are reported everywhere; p-values never appear alone.

**Balanced subsampling is rejected.** Discarding valid MBP structures to force equal strata would
(i) throw away real data, (ii) reduce power in the stratum with the strongest and most
mechanistically informative effect, and (iii) make the result depend on an arbitrary subsample
draw. Standardisation achieves the same protection against MBP dominating the headline number
while retaining every observation and keeping the estimator transparent. This choice was made on
design grounds and **not** by comparing effect sizes under the two schemes.

## II.4 Cohort (frozen)

- **Source:** `data_manifest/confirmatory_candidate_pool.json` — 406 eligible entities,
  **129 independent (partner, target) units**: BRIL 30, T4L 38, MBP 61.
- Excludes all 24 development PDB entries and all entities sharing a (partner, target accession)
  pair with the development set.
- **Use all 129 independent units**, one structure per (partner, target): highest resolution, ties
  by lowest PDB ID. No further sampling, no seed, no draw — using the entire eligible fresh
  population removes sampling discretion altogether and maximises power.
- Inclusion/exclusion criteria are unchanged from v1.1 except E9 (see II.5).

## II.5 E9 decision — REVISE (PILOT-INFORMED)

Audit findings (`results/e9_audit/`):

- **The v1.1 implementation has a reporting defect.** When a segment is empty — which happens for
  every non-chimeric, single-UniProt entity — the fallback assigns a disorder of **1.0**, so E9
  fires on 988 entities that were never candidates and had already failed E1. This inflated the
  apparent reach of E9.
- **True marginal impact of E9 is small.** Among the 983 genuinely chimeric entities, E9 flags 451,
  but only **37** are excluded by E9 *alone*. Removing E9 entirely would take the eligible pool
  from 455 to 492 entities (**+8.1%**), adding 16 new independent targets.
- **The implementation is correct where it applies.** Recomputing disorder directly from deposited
  coordinates for 30 outcome-blind randomly sampled chimeric E9-excluded entities agreed with the
  pipeline in **30/30** cases. Missing target residues are predominantly **internal** (71%), not
  terminal, so this is genuine disorder rather than unmodelled termini. Author-numbering gaps were
  present in 20/30 entities but cannot be mistaken for disorder because the pipeline works in
  label_seq space.
- Among genuinely chimeric entities the target-disorder distribution is well behaved (median 0.083,
  p75 0.195, p90 0.324); the 0.20 threshold sits at roughly the 75th percentile.

**Decision, made on measurement-validity grounds and explicitly not on effect size** (no comparison
of the fusion effect under alternative E9 settings was performed before deciding):

1. **Fix the empty-segment defect.** Disorder is not-evaluable when a segment is absent; E9 is
   evaluated only on entities that have both segments. This changes no eligibility decision — the
   affected entities already fail E1 — but it stops the audit trail misattributing exclusions.
2. **Keep the 0.20 primary threshold.** Unmodelled residues inside a segment create artificial
   surface of exactly the kind the deletion-boundary control (A3) exists to guard against, so
   completeness is a genuine measurement-validity requirement, not a convenience. The marginal cost
   is 37 entities.
3. **Add completeness as a pre-specified sensitivity factor and covariate.** The confirmatory
   analysis is repeated at a 0.35 threshold (approximately p90 of the chimeric distribution), and
   segment disorder enters the mixed-effects model as a covariate.

## II.6 Zero predictions are data, not failures (PILOT-INFORMED)

P2Rank returned zero pockets in the REMOVED condition for 6 development structures. All six were
verified as **genuine detector abstentions**: return code 0, valid input (no NaN coordinates, no
chain breaks over 4.5 A, no residues lacking CA), and fpocket found cavities in the same files.
They are small, predominantly helical domains (86-182 target residues, 72-96% helix). No structure
is ever removed from the study because a detector predicted nothing; such cases are reported in the
`NO_POCKET_PREDICTED` stratum and counted against the hypothesis wherever direction matters.
