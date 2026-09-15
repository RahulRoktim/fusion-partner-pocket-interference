# Table legends — MANUSCRIPT_DRAFT_v3

Main tables keep their v1/v2 numbering. The three tables that are supplementary in this version are renumbered S1–S3 so that they are cited in ascending order, as BMC/Journal of Cheminformatics format requires. **No frozen table file was modified or renamed.**

## Table-number mapping

| v3 number | Placement | Source file (unchanged on disk) | v2 number |
|---|---|---|---|
| Table 1 | Main | `results/tables/TABLE1_cohort_composition.tsv` | Table 1 |
| Table 2 | Main | `results/tables/TABLE2_confirmatory_exposure.tsv` | Table 2 |
| Table 3 | Main | `results/tables/TABLE3_displacement.tsv` | Table 3 |
| Table 3b | Main | `results/tables/TABLE3b_mechanism.tsv` | Table 3b |
| Table 4 | Main | `results/tables/TABLE4_reference_site.tsv` | Table 4 |
| **Table S1** | Additional file 1 | `results/tables/TABLE6_e9_sensitivity.tsv` | Table 6 |
| **Table S2** | Additional file 1 | `results/tables/TABLE6b_e9_new_structures.tsv` | Table 6b |
| **Table S3** | Additional file 1 | `results/tables/TABLE5_dataset_provenance.tsv` | Table 5 |

---

**Table 1. Cohort composition.**

Composition of the confirmatory cohort by fusion partner and overall, with the exploratory development cohort shown separately and labelled as such. Columns: structures; independent target clusters; median resolution and range; counts by experimental method; counts by fusion topology (internal insertion, N-terminal, C-terminal, complex); median modelled target and fusion segment lengths; and structures with an outcome-independent biological reference site. Confirmatory totals are 128 structures across 123 independent target clusters (BRIL 29, T4L 38, MBP 61). The two counts differ because five target accessions appear twice, each as a cross-partner pair of the same receptor solved with two different partners; all inference clusters on target accession. The development row is exploratory and is never pooled with the confirmatory rows. Target-accession overlap between the cohorts is zero.

**Table 2. Confirmatory exposure by rank depth.**

Proportion of confirmatory structures with a fusion-associated pocket at rank 1, within the top 3 and within the top 5, per detector, for each partner stratum and pooled. Each cell gives numerator, denominator, percentage and 95% Wilson confidence interval. Partner strata are listed before the pooled row. The final row for each detector is the pre-specified partner-standardised estimate weighting each partner one third, with a clustered-bootstrap interval on target accession; raw pooled and standardised estimates are always reported together. P2Rank denominators are 127 rather than 128 because one structure yielded no P2Rank prediction in the original condition; it is excluded from that detector's denominator and not imputed. The MBP top-3 and top-5 cells are at or near ceiling (61/61 and 61/61 for P2Rank; 60/61 and 60/61 for fpocket) and are not independent evidence of effect magnitude.

**Table 3. Target-cavity rank displacement.**

Structure-level displacement of the same target cavity between the original construct and the fusion-removed structure, by partner and detector. Columns: structures; evaluable structures; detector abstentions; unmatched structures; median, mean and maximum displacement; count of negative displacements; and the proportion with displacement greater than 0, at least 2 and at least 5, each with a 95% Wilson confidence interval. One structure contributes one observation. Detector abstention and target-cavity-not-recovered are reported as explicit states, never assigned an imputed rank and never counted as improvements. N_EVALUABLE is the denominator of every proportion in the same row and differs from N_STRUCTURES wherever abstentions or unmatched structures occur, most consequentially in the P2Rank MBP row, where 38 of 61 structures are evaluable (20 abstentions, 3 unmatched). The count of negative displacements is 0 in every partner stratum for both detectors, across all 231 evaluable structure-by-detector pairs. Structure-level exact sign tests give 42 up versus 0 down (p = 4.6 × 10⁻¹³) for P2Rank and 65 up versus 0 down (p = 5.4 × 10⁻²⁰) for fpocket; the direction and the absence of counter-examples, not the p-values, are the substantive result.

**Table 3b. Mechanism of rank-1 interference.**

Classification of every rank-1 fusion-associated result in the confirmatory cohort (126 across both detectors) by mechanism, using the frozen isolated-partner correspondence probe. Columns: rank-1 fusion-associated results; fusion-intrinsic cavity; target–fusion interface; chimera-specific fusion cavity; linker-related; mapping artifact; and the probe states (matched, no matching cavity, detector returned no predictions, ambiguous). Mechanism classes were fixed before the confirmatory analysis and none was added afterwards. Totals: 113 fusion-intrinsic, 11 interface, 2 chimera-specific, 0 linker-related and 0 mapping artifacts; of 119 probe matches, 101 correspond to the isolated partner's own rank-1 pocket. The probe state in which the isolated-partner structure yields no predictions is uninformative and is never read as biological absence; it did not occur. These are the counts plotted in Fig. 4A.

**Table 4. Biological reference-site subset (secondary analysis).**

Per detector: reference-site coverage; structures in which a site-corresponding pocket was found in the original construct; the number and proportion outranked by a fusion-associated pocket with a 95% Wilson confidence interval; rank changes after fusion removal (improved, worsened, unchanged, median); and top-1, top-3 and top-5 recovery of the reference site in the original construct. Coverage is 60/128 = 46.9%, below the 60% threshold fixed in advance as the condition for promoting this analysis to primary, so it remains secondary. Coverage is partner-skewed (BRIL 23/29, T4L 27/38, MBP 10/61). Reference sites were defined from deposited ligand coordinates and a pre-declared artifact-exclusion list without consulting any detector output.

The rank-change columns describe a different quantity from Table 3: here the tracked cavity is the pocket corresponding to the characterised biological site, whereas in Table 3 it is the cavity the detector itself ranks first in the target-only structure. The single fpocket worsening in this table is therefore not a counter-example to the zero negative displacements in Table 3. Rank change is computable in 58 of 58 structures for P2Rank (9 improved, 0 worsened, 49 unchanged) and in 59 of the 60 for fpocket (28 improved, 1 worsened, 30 unchanged).

---

## Supplementary tables (Additional file 1)

**Table S1. Segment-disorder sensitivity: primary 0.20 versus 0.35.**

Side-by-side comparison of every Aim A and Aim B quantity between the primary cohort (0.20 threshold, n = 128, 123 targets) and the sensitivity cohort (0.35, n = 141, 136 targets), per detector and partner, with numerators, denominators, percentages, 95% confidence intervals and the difference in percentage points. This is a sensitivity analysis: every primary manuscript value is the 0.20 column. The 0.35 cohort is a strict superset of the primary targets, not an independent replication cohort. The relaxation adds 13 independent units and, because the frozen selection rules are applied unchanged, selects a different, strictly higher-resolution representative structure for three BRIL targets. No conclusion changes.

**Table S2. Characteristics of the newly admitted structures.**

Descriptive characteristics of the 16 structures that required new detector runs for the 0.35 sensitivity cohort — 13 newly admitted independent units and 3 representative upgrades: PDB identifier, admission reason, partner, topology, experimental method, resolution, modelled target and fusion lengths, and target and fusion segment disorder. Descriptive only; no model was fitted to these structures and no post-hoc analysis was derived from them.

**Table S3. Dataset provenance audit (secondary analysis).**

Recognised crystallization fusion constructs in each audited dataset, with role, counting unit, canonical version pinned by SHA-256, denominators, construct counts, and — for PLINDER — the location of the labelled ligand site. Canonical denominators are in `results/PUBLICATION_DENOMINATORS.tsv` and the PLINDER split-level breakdown in `results/PLINDER_DENOMINATORS.tsv`. Counting units are stated per dataset and are not mixed: the classic datasets are counted in unique PDB entries; PLINDER is counted in systems and unique PDB entries, with both quoted together everywhere. Prevalence uses per-dataset denominators; the deduplicated union of entries examined (5,304; 5,419 before deduplication) is reported only as the number of distinct entries examined and is never used as a prevalence denominator. A standalone structure of a partner protein is not a construct and is counted separately. sc-PDB and LIGYSIS are recorded as not measured, with the access limitation stated, rather than assigned a value.
