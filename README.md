# Partner-dependent interference by crystallization fusion partners in ligand-binding-site prediction

Reproducibility package: analysis code, pre-specified protocol, cohort manifests and canonical
results for a target-disjoint confirmatory study of how crystallization fusion partners affect the
ranking of predicted ligand-binding sites.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22793220.svg)](https://doi.org/10.5281/zenodo.22793220)

| | |
|---|---|
| **Title** | Partner-dependent interference by crystallization fusion partners in ligand-binding-site prediction |
| **Author** | Md. Rahul Reza Roktim |
| **ORCID** | [0009-0003-6518-0495](https://orcid.org/0009-0003-6518-0495) |
| **Affiliation** | Department of Pharmacy, Daffodil International University, Dhaka, Bangladesh |
| **Contact** | roktim2311091058@diu.edu.bd |
| **Status** | Submission / preprint reproducibility repository. Not yet submitted; no preprint posted |
| **Scientific freeze** | `v1.1-paper-final-freeze` — all 224 tracked files hash-verified in `environment/FINAL_FREEZE_v1.1.json` |
| **Archived release** | [`10.5281/zenodo.22793220`](https://doi.org/10.5281/zenodo.22793220) — Zenodo snapshot of `v1.0.0-submission`. **Cite this version DOI.** |
| **All versions** | [`10.5281/zenodo.22793219`](https://doi.org/10.5281/zenodo.22793219) — concept DOI, always resolves to the newest version |
| **Licences** | Code MIT; data, tables, figures and documentation CC BY 4.0; third-party resources retain their own terms |

> **Current archive note (2026-09-16):** the version archive above is published. The
> pre-archive plan in `ARCHIVE_STRATEGY.md` is retained as historical provenance and is
> superseded by this table. Post-release maintenance may improve verification and portability,
> but it does not rewrite the immutable `v1.0.0-submission` tag or the archived scientific record.

---

## 1. Scientific question

Deposited protein structures are frequently engineered constructs. Crystallization fusion partners —
maltose-binding protein (MBP), T4 lysozyme (T4L) and thermostabilised apocytochrome b562RIL (BRIL) —
contribute large ordered domains that carry their own cavities, and the resulting coordinates are
reused as though they described the target alone. Automated binding-site predictors rank cavities by
geometry, not by which polypeptide a cavity belongs to.

This study asks how often a construct-associated cavity enters the highest ranks of an automated
prediction, how that varies by fusion partner and by detector, whether removing the fusion improves
the rank of the target's own preferred cavity, and whether the responsible cavities are intrinsic to
the partner or created by the chimeric junction.

## 2. Repository scope

This repository contains everything needed to verify or reproduce the analysis: the pre-specified
protocol and its amendment log, the cohort manifests with residue-level construct annotations, the
full analysis pipeline, the classified detector output, every canonical result table, every figure,
and four automated test suites.

It does **not** contain atomic coordinates, detector distributions, or third-party benchmark files.
Those are referenced by identifier, version and checksum and retrieved at run time — see
`THIRD_PARTY.md`.

## 3. Manuscript

> **Partner-dependent interference by crystallization fusion partners in ligand-binding-site
> prediction**
> Md. Rahul Reza Roktim, Department of Pharmacy, Daffodil International University, Dhaka,
> Bangladesh. ORCID [0009-0003-6518-0495](https://orcid.org/0009-0003-6518-0495)

Status: prepared for submission. Not yet submitted, and not yet published as a preprint.

## 4. Frozen scientific version

The scientific analysis is closed and frozen. Two tags mark the frozen states, and a third marks
this public release:

| Tag | Meaning |
|---|---|
| `v1.0-paper-analysis-freeze` | Publication analysis freeze |
| `v1.1-paper-final-freeze` | Final scientific state, after the last pre-specified sensitivity. `environment/FINAL_FREEZE_v1.1.json` records the SHA-256 of all 224 tracked files at that point |
| `v1.0.0-submission` | This public release, at commit `7ce70ea`. Adds packaging, licensing and manuscript files only; changes no scientific result. Archived at Zenodo as [`10.5281/zenodo.22793220`](https://doi.org/10.5281/zenodo.22793220) |

**No scientific result has changed since `v1.1-paper-final-freeze`.** Everything added after it is
manuscript text, packaging and documentation.

## 5. Directory structure

```
├── PROTOCOL.md                  pre-specified protocol v1.3 (+ three archived versions)
├── DEVIATIONS.md                append-only amendment log, pre- and post-outcome marked
├── POST_FREEZE_DEVIATIONS.md    deviations recorded after the final freeze
├── REPRODUCE.md                 clean third-party reproduction workflow
├── THIRD_PARTY.md               every upstream dependency, version, source and licence
├── LICENSE                      MIT, for code
├── LICENSE_DATA.md              CC BY 4.0, for data/figures/documentation
├── CITATION.cff
├── environment/                 environment lock and freeze manifests
├── scripts/                     numbered analysis pipeline (38 files)
├── scripts_release/             retrieval helpers for third-party objects
├── tests/                       four automated test suites
├── data_manifest/               cohort enumeration, eligibility, per-structure manifests
├── annotations/                 residue-level chimera indexes for the three partners
├── results/                     classified pockets, endpoints, sensitivities, reports, tables
├── figures/                     Figures 1–5 and S1–S11 (PNG + PDF)
├── dataset_audit/               dataset provenance audit outputs
├── logs/                        execution logs from the scientific runs
└── manuscript-support/          claims ledger, number audit, reference verification, readiness
```

## 6. Quick reproduction path — verify the reported numbers (minutes, no detectors)

Every headline value in the manuscript is derived from files included here. You do **not** need to
run P2Rank or fpocket to check them.

```bash
git clone https://github.com/RahulRoktim/fusion-partner-pocket-interference && cd fusion-tag-hazard
python -m pip install numpy scipy statsmodels matplotlib gemmi
python scripts/16_confirmatory_analysis.py      # rebuilds the confirmatory analysis from
                                                # results/confirmatory/primary/pockets_classified.json
python scripts/29_canonical_tables.py           # rebuilds Tables 1-6b
python scripts/30_figures_main.py               # rebuilds the main figures
python tests/test_no_development_leakage.py     # cohort disjointness
```

Storage: ~120 MB. Runtime: a few minutes.

## 7. Full reproduction path — from PDB identifiers (hours + ~2.7 GB)

See `REPRODUCE.md` for the complete step-by-step workflow, including detector installation. In
outline: retrieve third-party objects → enumerate and freeze the cohort → prepare structures →
run both detectors in both conditions → classify pockets → compute correspondence and endpoints →
run the sensitivity analyses → run the dataset audit → regenerate tables and figures.

## 8. Software requirements

| | |
|---|---|
| Python | 3.12.10 (3.11+ expected to work) |
| Python packages | gemmi 0.7.5, numpy 2.5.2, scipy 1.18.1, statsmodels 0.15.0, matplotlib 3.11.1 |
| Java | Temurin OpenJDK 17.0.19+10 (required by P2Rank) |
| P2Rank | 2.5.1, stock default model |
| fpocket | source tag 4.2.3, commit `4bb0d844`, built with a **serial** `make` |
| OS | Developed on Windows 11 with WSL2 Ubuntu 24.04 for fpocket. The analysis scripts are OS-independent; fpocket requires a Linux or macOS build |

Exact versions, hashes and two environment traps are in `environment/ENVIRONMENT.md`.

## 9. Storage requirements

| Item | Size |
|---|---|
| This repository | ≈ 112 MB |
| Downloaded mmCIF structure cache | ≈ 1.5 GB |
| Prepared detector inputs | ≈ 120 MB |
| Raw detector output, all four variants | ≈ 900 MB |
| P2Rank + fpocket distributions | ≈ 200 MB |
| **Total for a full reproduction** | **≈ 2.8 GB** |

## 10. Runtime

Measured on the original hardware, single-threaded, from the frozen run records:

| Stage | Runs | Measured total |
|---|---|---|
| Confirmatory primary, P2Rank | 256 | 21.0 min (mean 4.9 s/run) |
| Confirmatory primary, fpocket | 256 | 4.0 min (mean 0.9 s/run) |
| All variants combined (primary, ions, assembly 1, disorder sensitivity) | 2,092 | **1.75 h** |

Structure download and preparation dominate wall-clock time on a first run and depend on network
speed. The analysis and figure steps take minutes.

## 11. External datasets required

| For | Resource | Obtained from |
|---|---|---|
| The main study | PDB coordinates (by identifier), SIFTS alignments | RCSB (`files.rcsb.org`, `data.rcsb.org`) |
| The dataset provenance audit | CHEN11, JOINED, COACH420, HOLO4K, FPTRAIN membership lists | `github.com/rdk/p2rank-datasets` |
| The dataset provenance audit | PLINDER 2024-06/v2 split file | `storage.googleapis.com/plinder` |

`scripts_release/fetch_third_party.sh` retrieves all of these and verifies each SHA-256 against
`EXCLUDED_THIRD_PARTY.tsv`.

## 12. Licences

- **Code** (`scripts/`, `scripts_release/`, `tests/`): **MIT** — see `LICENSE`
- **Data, derived tables, figures, documentation**: **CC BY 4.0** — see `LICENSE_DATA.md`
- **Third-party software and data**: retain their own licences and are **not** relicensed here. See
  `THIRD_PARTY.md`

## 13. Citation

Cite the **version DOI**, which resolves to the exact archived state analysed in the manuscript:

> Roktim, Md. Rahul Reza (2026). *Partner-dependent interference by crystallization fusion partners
> in ligand-binding-site prediction* (v1.0.0-submission) \[Software\]. Zenodo. https://doi.org/10.5281/zenodo.22793220

```bibtex
@software{roktim_2026_fusion_partner_interference,
  author    = {Roktim, Md. Rahul Reza},
  title     = {Partner-dependent interference by crystallization fusion
               partners in ligand-binding-site prediction},
  version   = {v1.0.0-submission},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22793220},
  url       = {https://doi.org/10.5281/zenodo.22793220}
}
```

To cite the project as a whole rather than this snapshot, use the concept DOI [`10.5281/zenodo.22793219`](https://doi.org/10.5281/zenodo.22793219),
which always resolves to the newest version. Machine-readable metadata is in `CITATION.cff`.

Please also cite the originating depositors' publications for any PDB entry you use, and the
P2Rank and fpocket method papers if you run those tools.

## 14. Known deviations

Two are declared in `POST_FREEZE_DEVIATIONS.md`, both with **no effect on any result**:

1. **One pre-specified sensitivity was never executed.** A correspondence sensitivity using an
   overlap coefficient of 0.50 was listed in the protocol but was incompletely specified — the
   centroid-distance criterion, the assignment objective and its relationship to the Jaccard
   acceptance rule were never defined. The omission was found during submission quality control,
   after the confirmatory results were known, so those choices were **not** defined retrospectively
   and the sensitivity **remains unexecuted with an unknown result**. Three fully specified
   correspondence settings were executed and are reported. Robustness claims are restricted to
   those three.
2. **Two frozen documents describe that sensitivity as executed.** `results/READINESS.md` and
   `environment/FINAL_FREEZE_v1.1.json` list it among executed alternatives. Neither has been
   edited: the frozen record is preserved as historical provenance and superseded by
   `POST_FREEZE_DEVIATIONS.md`, not rewritten. No result value is affected; the manifest's file
   hashes remain correct.

The full amendment history, including every pre-outcome and post-outcome change during the study,
is in `DEVIATIONS.md`.

## 15. Provenance note: local paths and execution identifiers in frozen records

The frozen execution records in `results/**/detector_runs.json` and four scripts under `scripts/`
contain **local scratch paths** of the form `C:\Users\user\AppData\Local\Temp\claude\...` and **opaque execution/session identifiers**
such as `5c4a9bb4-b9dd-40f9-b5cc-b8b58139430b`.

These are retained verbatim, deliberately:

- They are part of the run provenance — the exact command and input path for each of the 512
  confirmatory detector runs and every sensitivity run.
- The Windows account name is the generic `user`. No personal name, home directory, machine
  identifier or institutional path is exposed.
- **They contain no password, token, key or credential of any kind.** The execution identifiers are
  opaque local directory names: they authenticate nothing, grant no access to any external service,
  and are useful only for linking records within this repository.
- Every one of these files is hashed in `environment/FINAL_FREEZE_v1.1.json`. Rewriting them would
  break the reproducibility hashes and destroy the exact provenance, for no privacy benefit.

A full secret scan of the package found zero API keys, tokens, passwords, private keys or cloud
credentials; see `RELEASE_QA.md`. If any identifier in these records is ever found to function as an
authentication credential, it should be reported rather than published.

## 16. Contact

Md. Rahul Reza Roktim — <roktim2311091058@diu.edu.bd>
Department of Pharmacy, Daffodil International University, Dhaka, Bangladesh
ORCID [0009-0003-6518-0495](https://orcid.org/0009-0003-6518-0495)

Issues and questions about reproduction are best raised on the repository's issue tracker.
