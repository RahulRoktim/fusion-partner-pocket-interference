# Statistical and figure provenance — manuscript branch v2

Three items raised in adversarial review required tracing a manuscript statement back to frozen source code or frozen output before v2 could state it. This note records what was found. **No frozen file was read-modified, no analysis was re-run, and no frozen tag was touched.** Everything below is an inspection result.

---

## 1. Score-margin summary statistic — **RESOLVED: median**

**Question.** `results/CONFIRMATORY_REPORT.md` reports score margins as `P2Rank BRIL 0.12 T4L 0.14 MBP 5.04` and `fpocket BRIL 0.30 T4L 0.31 MBP 1.56` without naming the summary statistic. Draft v1 asserted in a figure legend that these were medians; that assertion was not supported by the frozen text at the time it was made.

**Resolution.** The statistic is the **median**, confirmed from two independent frozen sources.

| Evidence | Location | Content |
|---|---|---|
| Generating code | `scripts/16_confirmatory_analysis.py`, Aim E block | The emitted line is `f"    {det:8s} {pt:5s}: n={len(vals):3d} median {median(vals):.2f}"`, where `median()` is defined at the top of the same file as the standard order statistic (mean of the two central values for even n). |
| Frozen output | `results/confirmatory/primary/confirmatory_summary.txt`, lines 128–136 | The frozen text itself carries the word `median` and the per-stratum n, which the prose summary in `CONFIRMATORY_REPORT.md` omits. |

**Definition.** For each structure in the ORIGINAL condition, the margin is the score of the highest-scoring fusion-associated pocket divided by the score of the highest-scoring TARGET_DOMINATED pocket. **The ratio is defined only where both pocket classes are present in the same structure**, which is why n is smaller than the stratum size.

| Detector | BRIL | T4L | MBP | ALL |
|---|---|---|---|---|
| P2Rank | 0.12 (n = 21) | 0.14 (n = 35) | 5.04 (n = 43) | 0.96 (n = 99) |
| fpocket | 0.30 (n = 29) | 0.31 (n = 38) | 1.56 (n = 61) | 0.76 (n = 128) |

**Action taken in v2.** Manuscript text and the Figure S8 legend now state "median", give the per-stratum n, and state the condition under which the ratio is defined.

---

## 2. Multiplicity / Holm correction — **AMBIGUOUS: family not reconstructible; no correction applied**

**Question.** The protocol states a Holm correction that no results section ever applied. The revision brief authorised computing Holm-adjusted p-values from already-frozen raw p-values **if and only if** the pre-specified family is unambiguous, and required this item to stop and report otherwise.

**What the protocol specifies.**

| Source | Text |
|---|---|
| `PROTOCOL.md` §8.5 | "**Multiplicity:** Holm correction within the family of secondary tests (3 partners x 2 detectors). The two primary endpoints are not corrected against each other..." |
| `PROTOCOL.md` §III.5.4 | "**Multiplicity.** Holm correction within the secondary family. Co-primary endpoints A and B are both reported unconditionally regardless of outcome." |

The family is therefore defined by its **shape** — six cells, one per (partner, detector) — but the protocol never names the *test* to be performed in each cell.

**What exists in frozen output.** Every p-value produced by the confirmatory analysis was enumerated:

| p-value | Value | Frozen location | Role |
|---|---|---|---|
| Structure-level exact sign test, P2Rank displacement | 4.55e-13 | `results/confirmatory/primary/confirmatory_summary.txt` line 73 | **Co-primary Aim B** — protocol excludes it from correction |
| Structure-level exact sign test, fpocket displacement | 5.42e-20 | same file, line 91 | **Co-primary Aim B** — excluded |
| Secondary model, partner[MBP] | 5.2e-07 | `results/CONFIRMATORY_REPORT.md`, model table | Regression coefficient |
| Secondary model, partner[T4L] | 0.30 | same | Regression coefficient |
| Secondary model, target_size / 100 residues | 0.0016 | same | Regression coefficient |
| Secondary model, resolution | 0.12 | same | Regression coefficient |
| Secondary model, detector[fpocket] | 0.12 | same | Regression coefficient |
| Secondary model, target_disorder | 0.38 | same | Regression coefficient |
| Interaction, fpocket × T4L | 0.41 | same | Interaction term |
| Interaction, fpocket × MBP | 0.13 | same | Interaction term |

**Finding.** **No family of six partner-by-detector hypothesis tests was ever instantiated.** The partner strata (Aim A) are reported as proportions with Wilson and clustered-bootstrap intervals and are not tested against anything. The secondary model is a single regression; its six non-intercept coefficients happen to number six but are covariate terms, not a partner × detector grid, and mapping them onto the pre-specified family would be a post-hoc reinterpretation. The only two genuine hypothesis tests in the confirmatory analysis are the co-primary sign tests, which the protocol explicitly exempts from correction.

**Decision.** The family cannot be reconstructed unambiguously from the pre-specified protocol. Per the brief, this item **stops here**: no Holm-adjusted values were computed and **`manuscript/HOLM_SECONDARY_AUDIT.tsv` was deliberately not created**, because producing it would have required inventing a family the protocol does not define.

**Action taken in v2.** The multiplicity statement is **retained, not deleted**, and is stated accurately in Methods: the protocol pre-specified Holm correction within a three-partner-by-two-detector secondary family; no such family of tests was performed; no adjustment was therefore applicable or applied; model p-values are reported unadjusted alongside effect sizes and intervals; and no secondary conclusion rests on a significance threshold.

**Recommendation for any future analysis.** If a partner-stratified hypothesis test is ever wanted, it must be pre-specified as a test (for example, six per-cell exact tests of the rank-1 rate against a stated reference), and the Holm correction applied to that named family. Retrofitting one now would be post-hoc and is not done.

---

## 3. Figure 4A (formerly Figure 5) case selection — **RESOLVED: a genuine deterministic rule, illustrative in purpose**

**Question.** The v1 legend asserted a "stated deterministic rule" without giving it, which review flagged as unverifiable and potentially cherry-picking.

**Resolution.** A deterministic rule exists and is written into the frozen generating script's docstring and implemented in its `pick()` function, both fixed before the figure was produced.

**Source:** `scripts/31_figure5_structural.py`, module docstring and `pick(partner, mechanism)`.

**The rule, as implemented.** Four (partner, mechanism) cells are requested in fixed order:

| Case | Partner | Mechanism class |
|---|---|---|
| C1 | MBP | FUSION_INTRINSIC_CAVITY |
| C2 | T4L | FUSION_INTRINSIC_CAVITY |
| C3 | BRIL | FUSION_INTRINSIC_CAVITY |
| C4 | any | TARGET_FUSION_INTERFACE |

For each cell, candidate structures are those with that mechanism assignment in the confirmatory cohort. Candidates whose rank-1 pocket is fusion-associated in **both** detectors are preferred; if none exists, all candidates for that cell are used. Among the resulting pool the structure with the **best (lowest) resolution** is taken, with ties broken by the **lowest PDB identifier**. The selected structures are 4EXK, 6ZX9, 6M97 and 4EPI.

**Purpose.** The same frozen docstring states: *"These are illustrative, NOT statistically representative."* Both facts are true simultaneously and v2 states both: the selection is deterministic and pre-coded, and its purpose is illustration. The systematic evidence for the mechanism is the full 126-result distribution in Table 3b, now also plotted as Figure 4B.

**Action taken in v2.** The Figure 4A legend and the Results text give the rule in full and state plainly that the four cases are illustrative rather than a representative sample. **The selection was not upgraded into a claim of deterministic sampling from the cohort.**

---

## Files inspected

`PROTOCOL.md` · `scripts/16_confirmatory_analysis.py` · `scripts/31_figure5_structural.py` · `scripts/32_figures_supp.py` · `results/confirmatory/primary/confirmatory_summary.txt` · `results/confirmatory/primary/secondary_model.json` · `results/CONFIRMATORY_REPORT.md` · `results/tables/TABLE3b_mechanism.tsv`

All were opened read-only. `environment/FINAL_FREEZE_v1.1.json` re-verifies unchanged after this work.
