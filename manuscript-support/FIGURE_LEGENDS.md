# Figure legends — MANUSCRIPT_DRAFT_v3

**Figures are renumbered for v3 so that main and supplementary figures are cited in ascending order, as BMC/Journal of Cheminformatics format requires. No frozen figure file was modified, renamed or regenerated.** The mapping below is authoritative.

## Figure-number mapping

| v3 number | Source file (unchanged on disk) | v2 number | v1 number |
|---|---|---|---|
| **Figure 1** | `figures/FIG1_study_design.*` | Figure 1 | Figure 1 |
| **Figure 2** | `figures/FIG3_rank_depth.*` | Figure 2 | Figure 3 |
| **Figure 3** | `figures/FIG4_displacement.*` | Figure 3 | Figure 4 |
| **Figure 4A** | `figures/v2/FIG4B_mechanism_distribution.*` | Figure 4B | — (new in v2) |
| **Figure 4B** | `figures/FIG5_structural_cases.*` | Figure 4A | Figure 5 |
| **Figure 5** | `figures/supplementary/S10_reference_site.*` | Figure 5 | Figure S10 |
| **Fig. S1** | `figures/FIG2_confirmatory_exposure.*` | Figure S10 | Figure 2 |
| **Fig. S2** | `figures/supplementary/S7_abstention_by_target_size.*` | S7 | S7 |
| **Fig. S3** | `figures/supplementary/S8_score_margins.*` | S8 | S8 |
| **Fig. S4** | `figures/supplementary/S9_detector_agreement.*` | S9 | S9 |
| **Fig. S5** | `figures/supplementary/S1_replication.*` | S1 | S1 |
| **Fig. S6** | `figures/supplementary/S2_dominance_sensitivity.*` | S2 | S2 |
| **Fig. S7** | `figures/supplementary/S3_correspondence_sensitivity.*` | S3 | S3 |
| **Fig. S8** | `figures/supplementary/S5_hetatm_ions_sensitivity.*` | S5 | S5 |
| **Fig. S9** | `figures/supplementary/S4_assembly1_sensitivity.*` | S4 | S4 |
| **Fig. S10** | `figures/supplementary/S6_e9_disorder.*` | S6 | S6 |
| **Fig. S11** | `figures/FIG6_dataset_provenance.*` | Figure S11 | Figure 6 |

*Note on panel letters in Figure 4:* v2 placed the case studies as panel A and the distribution as panel B. v3 swaps them so that the panels are cited in letter order, since the systematic distribution is discussed before the illustrative cases. The underlying image files are unchanged.

---

## Main figures

**Figure 1. Study design.**

Schematic of the analysis pipeline. A deposited PDB entry containing a chimeric polymer entity is annotated at residue level from SIFTS alignments into target, fusion-partner, linker and terminal expression-tag classes. Terminal tag residues are deleted from both arms, so the paired contrast is exactly ORIGINAL (target + fusion + linker) against FUSION-REMOVED (target only), with no target atom coordinate altered between the two. Both conditions are processed by P2Rank 2.5.1 and fpocket 4.2.3 at stock parameters, and every predicted pocket is assigned a class from the composition of its residues. The lower track shows the paired endpoint: the same target cavity is matched across conditions by Jaccard index on target residues (≥ 0.40) and target-residue centroid distance (≤ 8 Å) under a one-to-one assignment in which detector rank and score are never used for matching, and rank displacement is the original rank minus the removed rank. The lower-left box records the two-stage design: a 24-structure exploratory development cohort, a frozen protocol, then an independent target-disjoint confirmatory cohort of 128 structures

**Figure 2. Exposure by rank depth and fusion partner.**

Proportion of confirmatory structures (n = 128, 123 independent targets, target-disjoint from development) with a fusion-associated pocket — fusion-dominated, interface or linker — within rank 1, the top 3 and the top 5, by partner, for P2Rank (left) and fpocket (right). Vertical bars are 95% Wilson confidence intervals. Rank-1 interference is concentrated in MBP (93.4% and 85.2%) and is far lower for BRIL (7.1% and 10.3%) and T4L (18.4% and 13.2%). Exposure broadens with rank depth for every partner: within the top 5, 50.0% and 69.0% of BRIL constructs and 76.3% and 65.8% of T4L constructs contain a fusion-associated pocket. The MBP top-5 values are at or near ceiling (61/61 and 60/61) and carry no additional information about effect magnitude in that stratum. Raw pooled and partner-standardised estimates are given in Table 2 and shown graphically in Fig. S1. Source data: `results/tables/TABLE2_confirmatory_exposure.tsv`

**Figure 3. Target-cavity rank displacement.**

Upper panels: for each evaluable structure, the rank of the same target cavity in the original construct (left) and in the fusion-removed structure (right), joined by a line coloured by fusion partner; the y-axis is inverted so that rank 1 is at the top. One structure contributes one line: the structure, not the pocket, is the inferential unit. No line descends from left to right, that is, no evaluable structure showed negative displacement in either detector, across 231 evaluable structure-by-detector pairs. Lower panels: the three mutually exclusive structure states, shown explicitly and never imputed — evaluable; detector abstention, where the removed structure yielded no prediction (22 structures for P2Rank, 20 of them MBP, and 0 for fpocket); and target cavity not recovered in the original (3 structures for P2Rank). Structures in the latter two states are never counted as improvements. Source data: `results/confirmatory/primary/displacement.tsv` and `results/tables/TABLE3_displacement.tsv`

**Figure 4. Mechanism of rank-1 interference.**

**A** Distribution of mechanism classes over all 126 rank-1 fusion-associated results in the confirmatory cohort, across both detectors, classified by the frozen isolated-partner correspondence probe: 113 cavities intrinsic to the fusion partner, 11 at the target–fusion interface, 2 chimera-specific fusion cavities, 0 linker-related and 0 mapping or parsing artifacts. Bars are stacked by detector (P2Rank 63/3/0/0/0; fpocket 50/8/2/0/0). The panel plots the counts in Table 3b and computes nothing.

**B** Four structural case studies illustrating the mechanism classes: 4EXK (MBP, N-terminal, 1.28 Å), 6ZX9 (T4L, N-terminal, 2.52 Å) and 6M97 (BRIL, internal insertion, 3.03 Å) as fusion-intrinsic cavities, and 4EPI (T4L, 1.74 Å) as a target–fusion interface cavity. Blue is the target Cα trace, red the fusion partner, grey the linker; the star marks the rank-1 pocket centroid, and gold marks deposited groups within 5 Å of it. Each panel title gives the rank-1 pocket class, its fusion-residue fraction, and the rank of the matched target cavity before and after fusion removal.

*Case selection for panel B.* Cases were chosen by a rule fixed in the figure-generation code before the figure was produced: for each of four (partner, mechanism) cells — MBP/fusion-intrinsic, T4L/fusion-intrinsic, BRIL/fusion-intrinsic and any-partner/target–fusion-interface — the confirmatory structure of best resolution among those whose rank-1 pocket is fusion-associated in both detectors, with ties broken by PDB identifier. The rule is deterministic, but its purpose is illustrative: these four structures are not a statistically representative sample, and the population-level evidence for the mechanism is panel A and Table 3b.

*Representation.* Panel B shows coordinate-derived Cα traces rather than molecular-surface renderings, drawn directly from the deposited coordinates so that every element is traceable to the underlying data. A PyMOL script (`figures/FIG5_pymol.pml`) reproduces the same four cases, selections and viewpoints for publication-quality cartoon or surface rendering; substituting those renderings changes presentation only and no scientific content. Per-case validation data are in `figures/FIG5_case_validation.tsv`

**Figure 5. Biological reference-site subset (secondary analysis).**

Restricted-subset analysis on the 60 of 128 structures (46.9%) with an outcome-independent biological reference site on the target. For each detector: whether the site-corresponding pocket was found in the original construct, its rank, whether a fusion-associated pocket outranked it, and top-1, top-3 and top-5 recovery. P2Rank found the site in 58 of 60 and ranked it first in 47 of 58 (81.0%), with a fusion-associated pocket outranking it in 9 of 58 (15.5%); fpocket found it in 60 of 60 and ranked it first in 24 of 60 (40.0%), outranked in 27 of 60 (45.0%). Coverage is below the 60% threshold fixed in advance as the condition for promotion, so this analysis remains secondary as pre-specified. Coverage is also strongly partner-skewed — BRIL 23/29, T4L 27/38, MBP 10/61 — so the subset under-represents the partner in which rank-1 interference is greatest. Source data: `results/tables/TABLE4_reference_site.tsv`

---

## Supplementary figures (Additional file 1)

**Figure S1. Rank-1 exposure with partner-standardised estimates.**
Forest plot of rank-1 exposure by partner and detector with 95% Wilson confidence intervals, numerators and denominators, together with the pre-specified partner-standardised estimates (P2Rank 39.7%, fpocket 36.2%) with clustered-bootstrap intervals on target accession. The same values are given numerically in Table 2.

**Figure S2. Detector abstention against target size.**
P2Rank returns no prediction for 22 of the 128 fusion-removed, target-only structures, 20 of them MBP, with abstention concentrated on small targets; fpocket never abstains. Abstention is reported as an explicit state and never counted as an improvement. This figure supports the statement that the P2Rank MBP displacement denominator is not a random subset of that stratum.

**Figure S3. Score margin by partner.**
Distribution of the ratio of the best fusion-associated pocket score to the best target-dominated pocket score in the original condition, by partner and detector, on a logarithmic axis with parity marked. Boxes show the median and interquartile range. Medians are 0.12 (BRIL, n = 21), 0.14 (T4L, n = 35) and 5.04 (MBP, n = 43) for P2Rank, and 0.30 (n = 29), 0.31 (n = 38) and 1.56 (n = 61) for fpocket. The ratio is defined only for structures in which both a fusion-associated and a target-dominated pocket are present, which is why n is smaller than the stratum size.

**Figure S4. Rank-1 call agreement between detectors.**
Cross-tabulation over the 127 structures evaluable in both detectors: both 55, neither 57, P2Rank only 11, fpocket only 4. Agreement 88.2%, Cohen's kappa 0.764.

**Figure S5. Development versus confirmatory replication.**
Side-by-side comparison of the exploratory development cohort (n = 24) and the independent confirmatory cohort (n = 128) under identical frozen endpoint definitions. The cohorts are displayed side by side and are never pooled; no combined estimate appears in this figure or anywhere in this paper. Of 20 pre-specified signals, 17 replicated, 1 replicated partially and 2 did not; both failures are the fpocket BRIL signal, 5/8 in development against 3/29 (rank-1) and 4/29 (displacement) in confirmation.

**Figure S6. Pocket-dominance threshold sensitivity.**
Confirmatory rank-1 exposure at dominance thresholds 0.50, 0.70 (primary) and 0.90. The pooled rate moves by at most 1.6 percentage points.

**Figure S7. Target-cavity correspondence threshold sensitivity.**
Displacement under the three executed correspondence settings: Jaccard ≥ 0.25 with 12 Å, ≥ 0.40 with 8 Å (primary) and ≥ 0.60 with 5 Å. The displacement rate moves by at most 0.4 percentage points, consistent with matched pairs having a median Jaccard of 1.000 and a median centroid distance of 0.00 Å.

**Figure S8. Heteroatom-handling sensitivity.**
All non-polymer entities stripped (primary) against ions retained; rank-1 exposure changes by at most 0.8 percentage points.

**Figure S9. Biological-assembly sensitivity.**
Exposure and displacement computed on biological assembly 1 rather than the single designated chain — the largest single sensitivity effect (rank-1 43.7% and 38.9%; displacement 31.5% and 46.4%) — with direction, partner ordering and the absence of negative displacement preserved. Two MBP structures (5W0R, 5EDU) are not evaluable under assembly 1 and are retained in every other analysis.

**Figure S10. Segment-disorder eligibility sensitivity.**
Primary cohort at a 0.20 segment-disorder threshold (n = 128) against the sensitivity cohort at 0.35 (n = 141), per detector. Rank-1 exposure moves from 52.0% to 52.1% (P2Rank) and from 46.9% to 48.9% (fpocket); MBP remains the strongest partner at 94.0% and 86.6%; displacement remains one-directional with no negative value. This is a sensitivity analysis: every primary value in the paper is the 0.20 result on the 128-structure cohort, and the relaxed cohort is not a second confirmatory sample.

**Figure S11. Dataset provenance audit (secondary analysis).**
**A** Recognised crystallization fusion constructs in the datasets connected to P2Rank and fpocket, with per-dataset denominators in unique PDB entries: CHEN11 0/241, JOINED 1/534, COACH420 0/418, HOLO4K 0/4,004, FPTRAIN 0/222. Standalone structures of MBP, T4L or BRIL are not constructs and are excluded from these counts. **B** PLINDER 2024-06/v2 by split, with both counting units shown — systems and unique PDB entries — for chimeric-receptor systems and for systems whose labelled site lies on the fusion partner; the test split contains none. **C** 1DUG, the single construct in any of the datasets in panel A: the benchmarked chain carries 217 GST-carrier residues and 10 fibrinogen target residues, and the ground-truth ligand, glutathione, taken from the Binding MOAD annotation used by that benchmark [16], makes 72 contacts with the carrier and 0 with the target. Panel C describes one member of 534 and supports no claim about aggregate benchmark scores, which were not recomputed.
