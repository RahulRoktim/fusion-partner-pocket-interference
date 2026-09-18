# Submission checklist — Journal of Cheminformatics

> ## STATUS: SUBMITTED — 2026-09-17
>
> The manuscript was submitted to *Journal of Cheminformatics* (Springer Nature /
> BMC) as a Research article on **17 September 2026** and was in technical check
> at the time of writing. Peer review is single anonymous. The submission
> identifier and portal tracking link are deliberately not published here.
>
> The archived artefact cited in the manuscript's data availability statement is
> the Zenodo **version** DOI [`10.5281/zenodo.22793220`](https://doi.org/10.5281/zenodo.22793220),
> published 2026-09-16 from tag `v1.0.0-submission`. The concept DOI
> `10.5281/zenodo.22793219` is deliberately not cited in the manuscript.
>
> **This checklist below is a pre-submission artefact and is retained as a
> record of preparation, not as a statement of current status.** Where it reads
> as though submission is still pending, that describes the state when it was
> written. It has not been rewritten, because editing a preparation record after
> the fact would destroy its value as provenance.

Prepared during Phase 8 (editorial and consistency QA). Requirements were read from the journal's
public submission guidance in September 2026 and are separated from recommendations, as instructed.

> **READINESS AT TIME OF WRITING: READY FOR AUTHOR METADATA + SUBMISSION PACKAGE.**
> **Declared protocol deviation:** one incompletely specified, unexecuted correspondence sensitivity (`POST_FREEZE_DEVIATIONS.md` PF1). Scientific analysis is CLOSED; the remaining items below are author metadata and packaging, not science.

Sources: [Submission guidelines](https://jcheminf.biomedcentral.com/submission-guidelines) · [Research article guidance](https://jcheminf.biomedcentral.com/submission-guidelines/preparing-your-manuscript/research) · [BMC submission guidelines PDF](https://cibb2023.dei.unipd.it/submission/BMC/authors/BMC_submission_guidelines.pdf)

> **Access note.** The journal's guidance pages on `link.springer.com` redirect to an authenticated
> endpoint and were not fetched. The requirements below come from the BioMed Central mirror of the
> same guidance and from the BMC submission-guidelines document. **Every item marked REQUIREMENT
> must be re-read on the journal site by a human before submission**, because journal instructions
> change and because some details (graphical abstract, article-processing charge, current editorial
> policies) could not be confirmed from the accessible pages.

---

## A. Article type

| | |
|---|---|
| **REQUIREMENT** | Research article |
| Status | **MET.** The manuscript reports original primary research with a pre-specified confirmatory design |
| RECOMMENDATION | Not a Methodology article: the contribution is a measured phenomenon and a workflow implication, not a new algorithm. Not a Data Note: the datasets are derived, not the primary deliverable |

## B. Section structure

| | |
|---|---|
| **REQUIREMENT** | Abstract · Keywords · Background · Methods · Results · Discussion · Conclusions · List of abbreviations · Declarations · References · Figure legends · Tables · Additional files |
| Status | **MET in v3.** v2 used Introduction/Results/Discussion/Methods order; v3 renames Introduction to Background, moves Methods before Results, splits the former Discussion subsection "Interpretation" into a standalone Conclusions section, and adds List of abbreviations and Declarations |
| RECOMMENDATION | Keep Results and Discussion separate rather than combining them; the Discussion carries substantial limitation material that would be diluted by merging |

## C. Abstract

| | |
|---|---|
| **REQUIREMENT** | Must not exceed **350 words**. Unstructured for Research articles |
| Status | **MET.** 238 words, unstructured |
| RECOMMENDATION | The 238-word length is a deliberate self-imposed limit, not a journal constraint. There is room to expand to ~300 words if a reviewer asks for the pooled exposure rates or the sensitivity summary to appear in the abstract |

## D. Graphical abstract

| | |
|---|---|
| **REQUIREMENT** | Optional. If supplied: **920 × 300 pixels, maximum 150 KB, JPEG/PNG/SVG** |
| Status | **NOT PREPARED** |
| RECOMMENDATION | A cropped panel of Figure 3 (the paired displacement lines, none descending) would make an effective graphical abstract and requires no new analysis. Author decision |

## E. Keywords

| | |
|---|---|
| **REQUIREMENT** | Keywords required |
| Status | **MET.** Nine keywords supplied on the title page |

## F. References

| | |
|---|---|
| **REQUIREMENT** | Numbered citation style, numbered in order of first appearance in the text |
| Status | **MET.** 19 references, all cited, all defined, numbered in order of first citation; verified in the Phase 8 cross-reference check |
| RECOMMENDATION | Convert to the journal's exact reference template at submission (author list truncation rules, journal-abbreviation style). Reference 19 (SciPy) has 30+ authors and will need the journal's "et al." rule applied |

## G. Figures

| | |
|---|---|
| **REQUIREMENT** | Figures numbered in order of first citation; legends supplied in the manuscript, not embedded in the image files |
| Status | **MET in v3.** Main figures 1–5 and supplementary figures S1–S11 are cited in ascending order; all legends are in the manuscript; a v2→v3 renumbering map is in `FIGURE_LEGENDS_v3.md` |
| **REQUIREMENT** | Figure files in an accepted format at publication resolution |
| Status | **PARTIAL.** PNG at 350–400 dpi and vector PDF exist for every figure |
| RECOMMENDATION | Submit the PDF versions. Before submission, regenerate **Figure 4B** (structural cases) from `figures/FIG5_pymol.pml` for publication-quality rendering; the current panel is a coordinate-derived Cα trace, which is scientifically accurate but visually plain. **This changes presentation only and no scientific content** |

## H. Tables

| | |
|---|---|
| **REQUIREMENT** | Tables cited in order, self-contained, with descriptive titles and footnotes defining abbreviations |
| Status | **MET in v3.** Main tables 1, 2, 3, 3b, 4; supplementary tables S1–S3, all cited in ascending order |
| RECOMMENDATION | "Table 3b" is non-standard for BMC. If the editor objects, renumber the main tables 1–5 and update the five main-text callouts; the underlying files need not be renamed |

## I. Additional files

| | |
|---|---|
| **REQUIREMENT** | Named "Additional file 1", "Additional file 2", … with descriptive titles |
| Status | **MET in v3.** Two additional files are declared and described |
| RECOMMENDATION | Assemble Additional file 1 as a single PDF (supplementary figures S1–S11, supplementary tables S1–S3, the artifact-exclusion list, the power table, the 1DUG worked example, the environment specification) before submission |

## J. Reproducibility policy

| | |
|---|---|
| **REQUIREMENT** | The journal publishes only work that is **entirely reproducible by third parties**. Datasets, software and algorithms needed to reach the stated conclusions must be provided as supplemental material or otherwise accessible **without registration, login, or licence terms other than Creative Commons (data/text) and OSI-approved open-source licences (software)**. For software, **source code must be provided** |
| Status | **SUBSTANTIALLY MET, ONE AUTHOR DECISION OUTSTANDING.** All analysis source code, protocol versions, deviation log, cohort manifests, detector run records, classified pocket data and canonical result tables are provided. Atomic coordinates are not redistributed but are retrievable from the PDB without registration. Raw detector output trees (~900 MB) are not deposited but are regenerable from the manifests, with the SHA-256 of every run and output tree recorded |
| **ACTION REQUIRED** | The repository must be public and under an OSI-approved licence for the code and a Creative Commons licence for the data **before submission**, and the URL must replace the placeholder. Confirm that not depositing the ~900 MB raw output trees satisfies the editor, given that they are byte-verifiably regenerable |

## K. Declarations

| Subsection | Status |
|---|---|
| Ethics approval and consent to participate | **MET** — "Not applicable" (no human or animal subjects) |
| Consent for publication | **MET** — "Not applicable" |
| Availability of data and materials | **PRESENT, contains one placeholder** (repository URL) |
| Competing interests | **PLACEHOLDER — author action required** |
| Funding | **PLACEHOLDER — author action required** |
| Authors' contributions | **PLACEHOLDER — author action required** |
| Acknowledgements | **PLACEHOLDER — author action required** |

## L. ORCID

| | |
|---|---|
| **REQUIREMENT** | Could not be confirmed from the accessible guidance pages. Springer Nature journals generally require or strongly encourage ORCID for the corresponding author |
| Status | **NOT SUPPLIED** — no author identities exist yet |
| **ACTION REQUIRED** | Collect ORCID identifiers for all authors; verify the journal's current requirement at submission |

## M. Preprint policy

| | |
|---|---|
| **REQUIREMENT** | Could not be confirmed from the accessible guidance pages. BMC journals generally permit preprint deposition and require it to be declared |
| Status | **NOT APPLICABLE YET** — no preprint has been posted |
| RECOMMENDATION | If a preprint is posted, declare it at submission and link it to the final article. Note separately that **reference 17 (PLINDER) is itself a preprint** and is labelled as such in the reference list |

## N. Word and page limits

| | |
|---|---|
| **REQUIREMENT** | No article word limit was found for Research articles; the 350-word abstract limit is the only length constraint identified |
| Status | **MET.** Body 10,649 words (Abstract 238, Background 937, Methods 3,420, Results 3,526, Discussion 2,311, Conclusions 201) |
| RECOMMENDATION | Length is at the upper end of typical for this journal but is justified by the Methods, which carry the reproducibility argument. If the editor requests compression, the first candidates are the sensitivity-analysis Methods paragraph and the dataset-provenance Methods paragraph, both of which are fully documented in the Supplement |

---

## Outstanding placeholders (must be resolved before submission)

| # | Placeholder | Location in `MANUSCRIPT_DRAFT_v3.md` |
|---|---|---|
| 1 | `[AUTHOR NAMES — PLACEHOLDER]` | Title page, "Authors" |
| 2 | `[AFFILIATIONS — PLACEHOLDER]` | Title page, "Affiliations" |
| 3 | `[CORRESPONDING AUTHOR — PLACEHOLDER]` | Title page, "Corresponding author" |
| 4 | `[REPOSITORY URL — PLACEHOLDER]` | Declarations → Availability of data and materials |
| 5 | `[COMPETING INTERESTS — PLACEHOLDER]` | Declarations → Competing interests |
| 6 | `[FUNDING — PLACEHOLDER]` | Declarations → Funding |
| 7 | `[AUTHORS' CONTRIBUTIONS — PLACEHOLDER]` | Declarations → Authors' contributions |
| 8 | `[ACKNOWLEDGEMENTS — PLACEHOLDER]` | Declarations → Acknowledgements |
| 9 | Corresponding-author e-mail | Not present; required by the submission system |
| 10 | ORCID identifiers | Not present |
| 11 | Repository DOI or archive DOI | Not present; required if the repository is archived |

**No placeholder value has been invented.**

---

## Outstanding scientific-documentation items for the authors

These are not manuscript defects and none blocks submission, but each should be resolved in the
manuscript branch. **No frozen file may be edited to fix them.**

1. **RESOLVED as a declared protocol deviation — no longer outstanding.** The protocol listed four
   correspondence settings; three were executed. The fourth, an overlap-coefficient threshold of
   0.50, was **incompletely specified** — the centroid-distance criterion, the assignment objective
   and the relationship to the Jaccard acceptance rule were never defined, and the frozen
   machine-readable correspondence specification contains no overlap rule. Because the omission was
   found after confirmatory results were known, the implementation choices were **not** defined
   retrospectively and the sensitivity **remains unexecuted, with an unknown result**. Two frozen
   documents state otherwise (`results/READINESS.md` line 24;
   `environment/FINAL_FREEZE_v1.1.json`, `sensitivity_axes_executed`); neither was edited, and both
   are superseded by `manuscript/POST_FREEZE_DEVIATIONS.md` entries PF1 and PF2. Robustness claims
   in the manuscript are restricted to the three executed Jaccard/distance settings.
2. **Claims-ledger p-value gap.** The frozen `results/CLAIMS_LEDGER.md` renders the structure-level
   sign-test p-values for claim C4 as `None`. The correct values are in
   `results/CONFIRMATORY_REPORT.md` and in the frozen confirmatory summary, and are used in the
   manuscript. `manuscript/CLAIMS_LEDGER_v2.md` carries the corrected field.
3. **Reference 17 status.** PLINDER remains a bioRxiv preprint. Re-check for a peer-reviewed
   successor immediately before submission.
4. **Phase 0 bibliography.** Two author attributions in `PHASE0_LITERATURE_AUDIT.md` are wrong and
   are corrected in the manuscript reference list only; the frozen audit is unmodified.
