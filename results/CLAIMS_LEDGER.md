# Claims Ledger

Every number below is injected from a canonical frozen table by `scripts/33_claims_ledger.py`; none is hand-typed. Sources are named per claim.

**Claim classes:** PRIMARY CONFIRMATORY · SECONDARY CONFIRMATORY · MECHANISTIC · EXTERNAL CORROBORATION · EXPLORATORY/DEVELOPMENT ONLY · BACKGROUND

---

## CENTRAL CLAIM (locked)

> Certain crystallographic fusion partners introduce genuine ligandable cavities that can compete with pockets on the intended target protein in automated binding-site prediction. The magnitude is strongly partner-dependent. MBP produces a particularly strong rank-1 hazard because its native maltose-binding cleft is frequently highly ranked, while BRIL and T4L generally show weaker rank-1 interference but can still contribute top-ranked candidate pockets. Paired fusion removal demonstrates ranking interference by following the same target cavity across unchanged target coordinates. The principal failure mode is **attribution rather than cavity detection**: predictors frequently detect a physically real cavity belonging to the engineered crystallization construct. The effect is most consequential for workflows consuming deposited structures without construct-aware preprocessing, particularly where the target lacks a dominant characterised pocket.

Scope is fixed at that wording. It may be refined stylistically but not widened.

---

### C1 — PRIMARY CONFIRMATORY

**EXACT CLAIM.** In deposited chimeric constructs, a fusion-associated pocket occupies rank 1 in a substantial minority of structures.

- **EVIDENCE:** Aim A, confirmatory cohort, ORIGINAL condition.
- **COHORT / DATASET:** CONFIRMATORY n = 128 structures / 123 independent targets; target-disjoint from development.
- **NUMERATOR / DENOMINATOR:** P2Rank 66/127 = 52.0% [43.3-60.5]; fpocket 60/128 = 46.9% [38.4-55.5]. Partner-standardised (1/3 each): P2Rank 39.7% [34.4-45.6 (cluster bootstrap)], fpocket 36.2% [30.5-42.4 (cluster bootstrap)].
- **DETECTOR:** P2Rank 2.5.1 and fpocket 4.2.3, stock parameters.
- **STATISTICAL SUPPORT:** Wilson intervals; cluster-bootstrap intervals on target accession.
- **SENSITIVITY SUPPORT:** Dominance 0.50/0.70/0.90 (S2), HETATM stripped vs ions (S5), single chain vs biological assembly 1 (S4), and E9 segment disorder 0.20 vs 0.35 (S6, Table 6; n=141, +13 units). All six pre-specified axes executed; no conclusion changes.
- **ALLOWED WORDING:** "A fusion-associated pocket occupied rank 1 in about half of deposited constructs overall, but the pooled figure is dominated by one partner and must be read together with the partner-standardised estimate."
- **PROHIBITED OVERSTATEMENT:** "All fusion tags strongly bias pocket prediction." Quoting the pooled rate without the partner strata or the standardised estimate.
- **FIGURE/TABLE SUPPORT:** Figure 2; Table 2.

### C2 — PRIMARY CONFIRMATORY

**EXACT CLAIM.** Rank-1 interference is strongly heterogeneous by fusion partner, and is dominated by MBP.

- **EVIDENCE:** Aim A, partner strata.
- **COHORT / DATASET:** CONFIRMATORY; BRIL n = 28/29, T4L n = 38, MBP n = 61.
- **NUMERATOR / DENOMINATOR:** P2Rank: BRIL 2/28 = 7.1% [2.0-22.6], T4L 7/38 = 18.4% [9.2-33.4], MBP 57/61 = 93.4% [84.3-97.4]. fpocket: BRIL 3/29 = 10.3%, T4L 5/38 = 13.2%, MBP 52/61 = 85.2%.
- **DETECTOR:** both
- **STATISTICAL SUPPORT:** Secondary mixed model, no separation detected (none); partner[MBP] odds ratio 49.5 [10.8-226.9], p = 5.2e-07 (n = 255 rows, 123 clusters). ASSOCIATION.
- **SENSITIVITY SUPPORT:** Stable across dominance thresholds, both structural variants, and the E9 0.35 cohort (MBP remains the strongest rank-1 partner in both detectors: 94.0% / 86.6%).
- **ALLOWED WORDING:** "MBP fusion was strongly associated with rank-1 construct-pocket capture." "BRIL and T4L showed weaker rank-1 interference."
- **PROHIBITED OVERSTATEMENT:** "MBP and small target size were proven causal determinants." "Fusion partners cause rank-1 capture." The model is associational.
- **FIGURE/TABLE SUPPORT:** Figure 2, Figure 3; Table 2.

### C3 — PRIMARY CONFIRMATORY

**EXACT CLAIM.** Lower-ranked construct exposure is broad across all three partners, even where rank-1 capture is rare.

- **EVIDENCE:** Aim A rank-depth.
- **COHORT / DATASET:** CONFIRMATORY.
- **NUMERATOR / DENOMINATOR:** top-5: P2Rank BRIL 50.0%, T4L 76.3%, MBP 100.0%; fpocket BRIL 69.0%, T4L 65.8%, MBP 98.4%.
- **DETECTOR:** both
- **STATISTICAL SUPPORT:** Wilson intervals per stratum.
- **SENSITIVITY SUPPORT:** Dominance thresholds (S2); E9 0.35 (top-5 80.0% / 83.0%).
- **ALLOWED WORDING:** "Any workflow that inspects more than the single best pocket is exposed regardless of tag."
- **PROHIBITED OVERSTATEMENT:** Presenting top-5 exposure as equivalent in consequence to rank-1 capture.
- **FIGURE/TABLE SUPPORT:** Figure 3; Table 2.

### C4 — PRIMARY CONFIRMATORY

**EXACT CLAIM.** Removing the fusion partner from the same coordinate file improves the rank of the SAME target cavity, and never worsens it.

- **EVIDENCE:** Aim B, F3 structure-level displacement with frozen target-cavity correspondence.
- **COHORT / DATASET:** CONFIRMATORY; P2Rank evaluable 103, abstention 22, unmatched 3; fpocket evaluable 128.
- **NUMERATOR / DENOMINATOR:** displacement > 0: P2Rank 42/103 = 40.8% [31.8-50.4]; fpocket 65/128 = 50.8% [42.2-59.3]. Negative displacement: P2Rank 0, fpocket 0. Median P2Rank 0, fpocket 1.0.
- **DETECTOR:** both
- **STATISTICAL SUPPORT:** Structure-level exact sign test, P2Rank p = None, fpocket p = None; cluster-bootstrap intervals. ONE STRUCTURE = ONE OBSERVATION.
- **SENSITIVITY SUPPORT:** Correspondence thresholds 0.25/12 A, 0.40/8 A, 0.60/5 A (S3); 8 A deletion-boundary control; assembly 1 (S4); E9 0.35 (S6) — zero negative displacements in all 255 evaluable structure-detector pairs of the 0.35 cohort.
- **ALLOWED WORDING:** "The presence of the fusion construct displaced the rank of the same target cavity." "Removal never improved on the target cavity's own preferred rank in any evaluable structure."
- **PROHIBITED OVERSTATEMENT:** Pocket-level Wilcoxon over individual cavities as inferential evidence (withdrawn as pseudoreplication, DEVIATIONS F2). Any use of the abandoned "rank-1 becomes TARGET_DOMINATED" endpoint as causal evidence.
- **FIGURE/TABLE SUPPORT:** Figure 4; Table 3.

### C5 — MECHANISTIC

**EXACT CLAIM.** Most rank-1 interference is caused by cavities intrinsic to the fusion partner, not by cavities created by the chimera; the failure is attribution, not hallucination.

- **EVIDENCE:** Aim D, isolated-partner correspondence probe with explicit states.
- **COHORT / DATASET:** All 126 rank-1 fusion-associated results in the confirmatory cohort (both detectors).
- **NUMERATOR / DENOMINATOR:** FUSION_INTRINSIC_CAVITY 113/126; TARGET_FUSION_INTERFACE 11/126; MAPPING_ARTIFACT 0/126.
- **DETECTOR:** both
- **STATISTICAL SUPPORT:** Descriptive counts; no inference.
- **SENSITIVITY SUPPORT:** Same correspondence thresholds as C4. Mechanism was not recomputed for the E9 0.35 cohort (out of scope for that sensitivity).
- **ALLOWED WORDING:** "Predictors detect a physically real cavity that belongs to the engineered construct." "DETECTOR_NO_PREDICTIONS is uninformative, not biological absence."
- **PROHIBITED OVERSTATEMENT:** "Predictors hallucinate pockets." "The cavity is an artifact of the prediction method."
- **FIGURE/TABLE SUPPORT:** Figure 5; Table 3b.

### C6 — SECONDARY CONFIRMATORY

**EXACT CLAIM.** Where the target has a characterised biological ligand site, detectors show markedly different resilience, and P2Rank usually still ranks that site first.

- **EVIDENCE:** Aim C, restricted subset.
- **COHORT / DATASET:** Reference-site coverage 60/128 = 46.9%, BELOW the pre-frozen 60% promotion rule, so this remains SECONDARY.
- **NUMERATOR / DENOMINATOR:** outranked by a fusion pocket: P2Rank 9/58 = 15.5% [8.4-26.9]; fpocket 27/60 = 45.0%. top-1 recovery P2Rank 81.0%, fpocket 40.0%.
- **DETECTOR:** both
- **STATISTICAL SUPPORT:** Wilson intervals; restricted subset.
- **SENSITIVITY SUPPORT:** Not promoted to primary; coverage rule applied exactly as frozen. E9 0.35 coverage 62/141 = 44.0%, also below 60% — status unchanged.
- **ALLOWED WORDING:** "The hazard concentrates in targets without a dominant characterised pocket." "Detector architecture materially changes susceptibility."
- **PROHIBITED OVERSTATEMENT:** "Fusion partners routinely outrank real drug-target sites." Generalising this subset to the whole cohort.
- **FIGURE/TABLE SUPPORT:** Figure S10; Table 4.

### C7 — PRIMARY CONFIRMATORY

**EXACT CLAIM.** The headline signals replicated on an independent, target-disjoint cohort; one development-stage partner signal did not.

- **EVIDENCE:** Development vs confirmatory replication, cohorts never pooled.
- **COHORT / DATASET:** DEVELOPMENT n = 24 (exploratory) vs CONFIRMATORY n = 128; target overlap 0.
- **NUMERATOR / DENOMINATOR:** P2Rank rank-1 dev 11/24 -> conf 66/127 : REPLICATED. fpocket rank-1 dev 15/24 -> conf 60/128 : REPLICATED. fpocket BRIL dev 5/8 -> conf 3/29 : NOT REPLICATED.
- **DETECTOR:** both
- **STATISTICAL SUPPORT:** Replication criteria frozen before any confirmatory result was inspected (DEVIATIONS F6).
- **SENSITIVITY SUPPORT:** n/a
- **ALLOWED WORDING:** "Both headline exposure signals and both displacement signals replicated. One partner-specific development signal (fpocket on BRIL) did not."
- **PROHIBITED OVERSTATEMENT:** Pooling development and confirmatory cohorts for any significance test. Reporting only the signals that replicated.
- **FIGURE/TABLE SUPPORT:** Figure S1.

### C8 — EXTERNAL CORROBORATION

**EXACT CLAIM.** The classic datasets that trained, developed and tested P2Rank and fpocket contain very few recognised crystallization fusion constructs.

- **EVIDENCE:** Phase 5 dataset provenance audit, SIFTS construct detection on frozen canonical membership lists.
- **COHORT / DATASET:** CHEN11 (train), JOINED (development), COACH420 and HOLO4K (test), FPTRAIN.
- **NUMERATOR / DENOMINATOR:** CHEN11 0/241 entries · JOINED 1/534 entries · COACH420 0/418 entries · HOLO4K 0/4004 entries · FPTRAIN 0/222 entries
- **DETECTOR:** n/a (dataset audit)
- **STATISTICAL SUPPORT:** Wilson intervals per dataset; descriptive only.
- **SENSITIVITY SUPPORT:** n/a
- **ALLOWED WORDING:** "Classic P2Rank training/development/test sets contained very few recognised crystallization fusion constructs."
- **PROHIBITED OVERSTATEMENT:** "Pocket-prediction benchmarks generally contain no fusion constructs." "The benchmarks are contaminated." Counting standalone MBP/T4L/BRIL structures as constructs.
- **FIGURE/TABLE SUPPORT:** Figure 6A; Table 5.

### C9 — EXTERNAL CORROBORATION

**EXACT CLAIM.** A modern protein-ligand training corpus contains training systems whose ligand-binding site lies on a crystallization fusion partner.

- **EVIDENCE:** PLINDER 2024-06/v2, reverse intersection with the recognised-partner chimera index.
- **COHORT / DATASET:** PLINDER train 309140 systems / 76901 entries; test 1036 systems / 1020 entries.
- **NUMERATOR / DENOMINATOR:** train: construct receptor 1385 systems / 454 unique entries; site labelled ON the fusion 342 systems / 158 unique entries; site on the target 939 systems. test: 0 construct-receptor systems.
- **DETECTOR:** n/a (dataset audit)
- **STATISTICAL SUPPORT:** Descriptive counts, both units always quoted together.
- **SENSITIVITY SUPPORT:** Partner recognition limited to 12 UniProt accessions.
- **ALLOWED WORDING:** "PLINDER contains training systems whose protein-ligand site lies on a crystallization fusion partner." Always give SYSTEM and UNIQUE-ENTRY counts together.
- **PROHIBITED OVERSTATEMENT:** "PLINDER is contaminated." "PLINDER's labels are wrong" — maltose genuinely binds MBP. Citing the system count without the entry count. Claiming any model was degraded by these systems.
- **FIGURE/TABLE SUPPORT:** Figure 6B; Table 5.

### C10 — EXTERNAL CORROBORATION

**EXACT CLAIM.** One benchmark member demonstrates a crystallization-carrier cavity being scored as a ground-truth success.

- **EVIDENCE:** 1DUG in JOINED/DT198; MOAD-2013 ground-truth ligand GSH.
- **COHORT / DATASET:** n = 1 of the JOINED development/validation set.
- **NUMERATOR / DENOMINATOR:** GSH contacts: 72 to the GST carrier, 0 to the fibrinogen target. P2Rank recovers the labelled site at rank 1 with a pocket of 15 carrier and 0 target residues; fpocket does not recover it at >=25% overlap.
- **DETECTOR:** both
- **STATISTICAL SUPPORT:** None. n = 1.
- **SENSITIVITY SUPPORT:** n/a
- **ALLOWED WORDING:** "A single benchmark member illustrates the failure mode: the labelled site belongs to the crystallization carrier."
- **PROHIBITED OVERSTATEMENT:** Any claim of measurable benchmark-score impact. Extrapolating from n = 1.
- **FIGURE/TABLE SUPPORT:** Figure 6C.

### C11 — SECONDARY CONFIRMATORY

**EXACT CLAIM.** Detector architecture materially changes susceptibility.

- **EVIDENCE:** Aim E.
- **COHORT / DATASET:** CONFIRMATORY.
- **NUMERATOR / DENOMINATOR:** rank-1 call agreement 88.2% (n = 127), Cohen's kappa 0.764. P2Rank abstains in the removed condition for 22 structures; fpocket for 0.
- **DETECTOR:** both
- **STATISTICAL SUPPORT:** Cohen's kappa; descriptive. Detector x partner interaction was estimable and NOT significant.
- **SENSITIVITY SUPPORT:** n/a
- **ALLOWED WORDING:** "Detector architecture materially changes susceptibility, most visibly in abstention behaviour and score margins."
- **PROHIBITED OVERSTATEMENT:** "P2Rank is more reliable than fpocket" (or the reverse) as a general claim. Asserting a significant detector x partner interaction.
- **FIGURE/TABLE SUPPORT:** Figures S7, S8, S9.

### C12 — EXPLORATORY / DEVELOPMENT ONLY

**EXACT CLAIM.** Development-cohort observations, including the pilot's mechanistic case studies and the pilot's partner rates.

- **EVIDENCE:** 24-structure development cohort.
- **COHORT / DATASET:** DEVELOPMENT n = 24, permanently exploratory.
- **NUMERATOR / DENOMINATOR:** See the pilot report; not restated as evidence.
- **DETECTOR:** both
- **STATISTICAL SUPPORT:** None permitted. All development p-values are descriptive; the cavity-level Wilcoxon results were withdrawn as pseudoreplication.
- **SENSITIVITY SUPPORT:** n/a
- **ALLOWED WORDING:** "In an exploratory development cohort we observed ..." clearly labelled.
- **PROHIBITED OVERSTATEMENT:** Pooling with confirmatory data. Any confirmatory significance claim. Reusing development structures as confirmatory observations.
- **FIGURE/TABLE SUPPORT:** Figure S1 (side-by-side only).

### C13 — BACKGROUND

**EXACT CLAIM.** fpocket's pocket ranking depends on a fitted, data-derived scoring function.

- **EVIDENCE:** Direct inspection of fpocket 4.2.3 source; fpocket 2009 paper.
- **COHORT / DATASET:** n/a
- **NUMERATOR / DENOMINATOR:** src/pscoring.c score_pocket(): intercept -0.03783394 plus five weighted descriptors; fitted by PLS on 307 proteins per the paper; descriptor normalisation recalibrated against the full PDB in 2017.
- **DETECTOR:** fpocket 4.2.3
- **STATISTICAL SUPPORT:** n/a
- **SENSITIVITY SUPPORT:** n/a
- **ALLOWED WORDING:** "fpocket's ranking depends on data-derived parameters."
- **PROHIBITED OVERSTATEMENT:** "fpocket has no training set" (explicitly withdrawn, DEVIATIONS H1). Claiming FPTRAIN is identical to fpocket's own training set.
- **FIGURE/TABLE SUPPORT:** n/a

---

## Global wording rules

| ALLOWED | NOT ALLOWED |
|---|---|
| "The presence of the fusion construct displaced the rank of the same target cavity." | "Fusion tags break pocket prediction." |
| "MBP fusion was strongly associated with rank-1 construct-pocket capture." | "MBP and small target size were proven causal determinants." |
| "PLINDER contains training systems whose protein-ligand site lies on a crystallization fusion partner." | "PLINDER is contaminated." |
| "Classic P2Rank training/development/test sets contained very few recognised crystallization fusion constructs." | "Pocket-prediction benchmarks generally contain no fusion constructs." |
| "Predictors detect a physically real cavity belonging to the construct." | "Predictors hallucinate pockets." |
| "Detector architecture materially changes susceptibility." | "P2Rank is more reliable than fpocket." |

## Standing constraints

1. The original broad hypothesis is **not** fully confirmed and must never be described as such. The confirmed hypothesis is the refined, partner-dependent one.
2. Partner strata are reported before any pooled number is interpreted; every overall exposure estimate appears as both raw pooled and partner-standardised.
3. Development and confirmatory cohorts are never pooled for inference.
4. One structure = one inferential observation. Cavity-level statistics are descriptive only.
5. `DETECTOR_ABSTENTION` and `TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL` are reported as states and never assigned an imputed rank.
6. MBP identity and target size are **associated** predictors, not independently proven causal factors.
7. All six pre-specified sensitivity axes have now been executed. The E9 = 0.35 analysis is SENSITIVITY ONLY: every primary manuscript value remains the E9 = 0.20 confirmatory result (n = 128). Verdict ROBUST against the rule frozen in DEVIATIONS J2.
