# Release, archive and preprint strategy

**Recommendation only. Nothing has been published, released, archived or submitted.** Each step
below is a separate decision for the author, and the order matters because the identifiers created
at each step are referenced by the next.

---

## Recommended order

### Step 1 — Public GitHub repository

| | |
|---|---|
| **What becomes public** | The assembled package: 234 files, ≈ 112 MB. Code (MIT), derived data, tables, figures and documentation (CC BY 4.0). No atomic coordinates, no detector distributions, no third-party benchmark files |
| **Prerequisite** | The QA in `RELEASE_QA.md` passes, and the author has reviewed the two items it flags |
| **Tag to create** | `v1.0.0-submission`, created **only after** the repository content is final |
| **What it does not have yet** | A DOI. GitHub URLs are mutable and are not citable identifiers |

The repository must be public before the journal is approached, because *Journal of Cheminformatics*
publishes only work that third parties can reproduce, with code and necessary data accessible
without registration or restrictive licence terms.

### Step 2 — Zenodo archive

| | |
|---|---|
| **What becomes public** | An immutable snapshot of the tagged GitHub release |
| **Why** | GitHub is not an archive. A Zenodo deposit gives a permanent DOI that survives repository renaming, transfer or deletion, which is what a reference list needs |
| **How** | Enable the Zenodo–GitHub integration, then publish the `v1.0.0-submission` tag as a release; Zenodo captures it automatically and mints the DOI |
| **DOI relationship** | Zenodo issues a **concept DOI** that always resolves to the newest version, and a **version DOI** for this specific snapshot. **Cite the version DOI in the manuscript** so a reader reaches exactly the analysed state; put the concept DOI in the README |
| **Metadata** | Populate from `CITATION.cff`: title, sole author, ORCID, affiliation, MIT + CC BY 4.0, keywords. Set the publication date to the archive date, not the analysis date |

### Step 3 — Preprint (optional; author's decision)

| | |
|---|---|
| **Venue** | ChemRxiv fits the subject matter. bioRxiv is a reasonable alternative given the structural-biology content |
| **What becomes public** | The manuscript and Additional file 1 |
| **Arguments for** | Establishes priority and a citable date; makes the work readable while under review; *Journal of Cheminformatics* is an open-access BMC title and BMC journals generally permit preprint deposition |
| **Arguments against** | Once posted it cannot be withdrawn, only marked; and if the manuscript changes substantially at review, the preprint and the final article will differ visibly |
| **DOI relationship** | The preprint gets its own DOI, distinct from the Zenodo software/data DOI and from the eventual article DOI. The preprint should cite the Zenodo version DOI for the code and data |
| **If posted** | Declare it at submission and link it to the final article. Note that reference 17 in this manuscript is itself a preprint and is labelled as such |

### Step 4 — Journal submission

| | |
|---|---|
| **Prerequisite** | Steps 1 and 2 complete, so that `[REPOSITORY URL — PLACEHOLDER]` and `[ARCHIVE DOI — PLACEHOLDER]` can be replaced with real values in the manuscript **before** the files are uploaded |
| **Files** | `manuscript/submission/` — manuscript DOCX and PDF, Additional file 1, six figure files in PNG and vector PDF, cover letter, metadata sheet |
| **Not yet done** | Author instructions must be re-read on the journal's own site; the accessible mirror was used for this package because the Springer-hosted pages redirect to an authenticated endpoint |

---

## DOI relationships at a glance

```
GitHub repo  ──tag v1.0.0-submission──►  Zenodo concept DOI   (always latest)
                                          └─ Zenodo version DOI   ◄── cited in the manuscript
                                                    ▲
                     ChemRxiv preprint DOI ─────────┘ (cites the version DOI)
                                 │
                                 └──────────►  Journal article DOI (links the preprint)
```

---

## Tag naming

| Tag | Meaning | Create when |
|---|---|---|
| `v1.0-paper-analysis-freeze` | **Exists.** Publication analysis freeze | — |
| `v1.1-paper-final-freeze` | **Exists.** Final scientific state | — |
| `v1.0.0-submission` | Public release: packaging, licensing and manuscript only | After `RELEASE_QA.md` passes and the author has resolved its two flagged items |

The scientific tags use the project's internal `v1.x-<phase>` scheme; the public release uses
semantic versioning because that is what Zenodo and package tooling expect. Keeping both is
deliberate: the scientific freezes are historical and must not move, and the release tag is what a
citation points at.

---

## What must be true before any of this happens

1. `RELEASE_QA.md` passes, and its two flagged items have been read and accepted by the author.
2. The generative-AI disclosure decision has been taken (`AI_DISCLOSURE_NOTE.md`).
3. Reference 17 (PLINDER) has been re-checked for a peer-reviewed successor.
4. The author is satisfied that publishing the repository under their own name, ORCID and
   institutional e-mail is what they intend.
