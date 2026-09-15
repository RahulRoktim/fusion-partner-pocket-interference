# Phase 1 — Frozen Study Protocol (v1.0)

**Project:** Fusion-Tag Hazard
**Frozen:** 2026-09-15, before any pocket detector has been run on any structure.
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
- S7. **[Recommended addition]** Are existing pocket-prediction benchmark and training sets
  (HOLO4K, COACH420, sc-PDB, LIGYSIS, PLINDER) contaminated with chimeric constructs?

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
5. **Fusion-removed counterpart:** identical file with all FUSION and LINKER residues deleted; TAG
   residues also deleted. Nothing else changes — no renumbering, no chain break repair, no
   minimisation. Any terminal chain break introduced is left as-is and recorded.

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

### 7.2 Full-study primary endpoint (confirmatory, causal)

Paired, within-structure: **P(target reference site is the rank-1 pocket)** in the fusion-removed
structure vs the original structure. Tested with an exact McNemar test, per detector. Restricted to
structures with a defined target reference site (Section 7.4).

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
- Benchmark-set contamination rate (S7).

### 7.4 Target reference site

Where a biologically relevant ligand is bound to the target segment, the reference site is the set
of target residues within 4.0 A of that ligand, with relevance judged by BioLiP-style artifact
exclusion lists (buffers, cryoprotectants, additives). Ligand-recovery endpoints are computed
**only** on the subset with a defined reference site; coverage is reported, and ligand recovery is
not promoted to a primary endpoint unless coverage exceeds 60% of the cohort.

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

## 10. Pilot go/no-go gate

Computed threshold: under a null rank-1 failure rate of 10%, observing **>= 6 of 24** failures has
p = 0.028; >= 5 has p = 0.085. So 6/24 is the smallest count that is inconsistent with an
uninteresting background rate.

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
