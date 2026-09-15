# Manuscript Outline

**Outline only — no manuscript prose is written here.** Every result slot points at a canonical
table, figure or claim ID. Nothing may be asserted in the manuscript that is not in
`results/CLAIMS_LEDGER.md`.

---

## TITLE OPTIONS

Ranked within each category. All avoid universal claims, "contamination", and "AI failure" framing.

**Technical / precise**

1. Crystallographic fusion partners introduce ligandable cavities that compete with target sites in automated pocket prediction
2. Partner-dependent ranking interference by crystallization fusion constructs in binding-site prediction
3. Construct-derived cavities and pocket-ranking interference in deposited protein structures
4. Maltose-binding protein fusions dominate rank-1 construct interference in automated binding-site prediction

**Broad / mechanism-forward**

5. A real cavity on the wrong protein: attribution failure in automated binding-site prediction of fusion constructs
6. When the crystallization tag outranks the target: construct-aware caveats for pocket-based drug discovery
7. Engineered crystallization constructs compete with intended targets in structure-based pocket detection

**Journal-of-Cheminformatics style**

8. Construct-aware evaluation of ligand-binding-site prediction: fusion partners, pocket ranking and the attribution problem
9. Quantifying fusion-partner interference in ligand-binding-site prediction: a target-disjoint confirmatory study
10. Assessing crystallographic fusion constructs as a source of pocket-ranking error in P2Rank and fpocket

**Preprint / discoverable**

11. Crystallization tags can outrank their targets in automated binding-pocket prediction
12. How often does a crystallization fusion partner capture the top-ranked predicted pocket?

### Recommendations

| Category | Choice | Why |
|---|---|---|
| **BEST TECHNICAL TITLE** | #2 — *Partner-dependent ranking interference by crystallization fusion constructs in binding-site prediction* | Names the mechanism (ranking interference), the exposure (crystallization fusion constructs), and the key qualifier (partner-dependent) without overreach. |
| **BEST BROAD TITLE** | #5 — *A real cavity on the wrong protein: attribution failure in automated binding-site prediction of fusion constructs* | Carries the paper's central intellectual point — attribution, not hallucination — in the title. |
| **BEST JOURNAL-OF-CHEMINFORMATICS STYLE** | #9 — *Quantifying fusion-partner interference in ligand-binding-site prediction: a target-disjoint confirmatory study* | Signals the quantitative contribution and the development/confirmation design, which is the methodological selling point. |
| **BEST PREPRINT TITLE** | #11 — *Crystallization tags can outrank their targets in automated binding-pocket prediction* | Discoverable and accurate; "can" is doing deliberate work. |

**RECOMMENDED WORKING TITLE:**
> *Partner-dependent ranking interference by crystallization fusion constructs in automated binding-site prediction*

It is the only option that states the finding, its principal qualifier, and its scope in one line,
and it will not need weakening at review.

---

## ABSTRACT LOGIC (structure, not prose)

- **Background.** Deposited structures frequently contain engineered crystallization constructs; automated pocket prediction increasingly runs at scale over deposited coordinates.
- **Question.** Can construct-derived cavities compete with cavities on the intended target, and if so how often, for which partners, and with what consequence?
- **Methods.** SIFTS-annotated fusion boundaries; paired original vs fusion-removed structures with unchanged target coordinates; two detectors at stock settings; frozen protocol; 24-structure exploratory development cohort, then an independent target-disjoint confirmatory cohort of 128 structures; frozen target-cavity correspondence for a non-degenerate paired endpoint.
- **Principal results.** Exposure rates by partner and rank depth [C1, C2, C3]; one-directional rank displacement of the same target cavity [C4]; mechanism dominated by cavities intrinsic to the partner [C5]; detector-dependent resilience where a characterised biological site exists [C6]; replication with one partner-level non-replication [C7].
- **Interpretation.** Attribution failure, strongly partner-dependent, most consequential for small or uncharacterised targets consumed without construct-aware preprocessing.

---

## INTRODUCTION

1. Deposited structures frequently contain engineered constructs — fusion partners are standard crystallization tools, and the resulting coordinates are deposited and reused as if they were the target alone. *(BACKGROUND)*
2. Automated pocket detection increasingly operates at scale — proteome-wide scans, pocket databases, ML training corpora. Human inspection does not scale with it. *(BACKGROUND)*
3. Pocket predictors score geometry and physicochemistry, not biological ownership: nothing in the input distinguishes "target" from "tag". *(BACKGROUND, motivates the attribution framing)*
4. The unresolved question: can construct-derived cavities compete with target cavities in ranked output, and does the construct's presence actively demote the target's own cavity?
5. Study aims and design: prespecified protocol, exploratory development cohort, independent target-disjoint confirmatory cohort, paired within-structure removal, two detectors, and a dataset-provenance audit.

*Prior-art positioning from `PHASE0_LITERATURE_AUDIT.md`. Closest precedents: Utgés & Barton 2024 (ranking-aware benchmarking), Bradford et al. 2021 (a crystallographic artifact propagating into computational predictions). Neither treats construct composition as a variable.*

---

## RESULTS

| § | Heading | Evidence | Claims |
|---|---|---|---|
| R1 | Study cohort and construct annotation | Table 1; Figure 1 | design |
| R2 | Fusion-associated pockets frequently enter ranked candidate lists | Table 2; Figures 2, 3 | C1, C3 |
| R3 | Rank-1 interference is strongly fusion-partner dependent | Table 2; Figure 2; secondary model | C2 |
| R4 | Fusion removal restores the ranking of the same target cavity | Table 3; Figure 4 | C4 |
| R5 | Intrinsic fusion-partner cavities explain most interference | Table 3b; Figure 5 | C5 |
| R6 | Biological target sites reveal detector-dependent resilience | Table 4; Figures S10, S7–S9 | C6, C11 |
| R7 | Findings replicate and are robust across structural, classification and disorder-threshold sensitivities | Figures S1–S6; Table 6 | C7 |
| R8 | Historical benchmark datasets contain few constructs, while modern PLINDER training data include construct-associated ligand systems | Table 5; Figure 6 | C8, C9, C10 |

**Order rule:** R2 and R3 must appear in that order — exposure first, then the partner decomposition — so the pooled number is never left standing alone. R8 is explicitly secondary and must open by reporting the Tier-1 **negative**.

---

## DISCUSSION

1. **Attribution failure, not pocket hallucination.** The detectors are right about the geometry and wrong about the owner. This reframes the problem from detector accuracy to input provenance. *(C5)*
2. **Why MBP is uniquely hazardous.** An occupied, high-quality native cleft; a tag routinely larger than the target; N-terminal topology; and targets that often lack a characterised site. *(C2, C6)*
3. **Detector dependence.** Score-margin and abstention behaviour differ; the interaction was not significant, so this is described, not modelled. *(C11)*
4. **Practical implications for CADD pipelines.** Where the hazard bites: small, uncharacterised targets; automated top-k inspection; pocket databases built from raw deposited coordinates.
5. **Implications for dataset curation.** Construct provenance is a curation dimension that ligand-artifact filtering does not cover. Report the Tier-1 negative honestly before the PLINDER observation. *(C8, C9)*
6. **Limits of the study.** Three partners; two detectors; association not causation for partner identity and target size; reference-site coverage below the promotion threshold in both the primary and the relaxed-threshold cohort.
7. **Preprocessing recommendation.** A construct-aware pre-flight check: SIFTS-annotate the entity, report pocket composition alongside rank, and treat a construct-dominated top pocket as a flag rather than a result.

---

## METHODS

1. Cohort construction: RCSB/SIFTS enumeration; inclusion/exclusion E1–E9; redundancy control; development/confirmatory separation and the leakage test.
2. Structure preparation: single designated chain; HETATM stripped; TAG deleted from both conditions; fusion-removed counterpart; deletion-boundary recording.
3. Detectors: P2Rank 2.5.1 and fpocket 4.2.3 with exact build provenance, stock parameters, command templates, hashes.
4. Pocket classification: dominance thresholds and the fusion-associated definition.
5. Target-cavity correspondence: the frozen matching rule and its validation.
6. Endpoints: Aim A exposure; Aim B F3 displacement with explicit non-numeric states; Aims C–E.
7. Statistics: structure as inferential unit; partner stratification and standardisation; clustered bootstrap; sign tests; the secondary model with its separation check.
8. Sensitivity analyses: all six pre-specified axes, executed. The E9 disorder axis (0.20 primary vs 0.35) is reported as sensitivity only; primary values are the 0.20 results.
9. Dataset provenance audit: frozen membership, construct detection, ground-truth localisation, PLINDER reverse intersection.
10. Reproducibility: protocol versions and hashes, `DEVIATIONS.md`, environment lock, `REPRODUCE.md`, data availability.

---

## FIGURE AND TABLE MANIFEST

| Item | Path |
|---|---|
| Figure 1 study design | `figures/FIG1_study_design.png` / `.pdf` |
| Figure 2 confirmatory exposure | `figures/FIG2_confirmatory_exposure.*` |
| Figure 3 rank-depth | `figures/FIG3_rank_depth.*` |
| Figure 4 displacement | `figures/FIG4_displacement.*` |
| Figure 5 structural cases | `figures/FIG5_structural_cases.*` (+ `FIG5_case_validation.tsv`, `FIG5_pymol.pml`) |
| Figure 6 dataset provenance | `figures/FIG6_dataset_provenance.*` |
| S1–S10 | `figures/supplementary/` (S6 = E9 0.20 vs 0.35 sensitivity) |
| Table 6 / 6b E9 sensitivity | `results/tables/TABLE6_e9_sensitivity.tsv`, `TABLE6b_e9_new_structures.tsv` |
| Table 1 cohort | `results/tables/TABLE1_cohort_composition.tsv` |
| Table 2 exposure | `results/tables/TABLE2_confirmatory_exposure.tsv` |
| Table 3 displacement | `results/tables/TABLE3_displacement.tsv` |
| Table 3b mechanism | `results/tables/TABLE3b_mechanism.tsv` |
| Table 4 reference site | `results/tables/TABLE4_reference_site.tsv` |
| Table 5 dataset provenance | `results/tables/TABLE5_dataset_provenance.tsv` |
| Denominators | `results/PUBLICATION_DENOMINATORS.tsv`, `results/PLINDER_DENOMINATORS.tsv` |
