# Third-party software and data

Every upstream resource this study depends on, with the information needed to obtain the exact
version used. **None of these is redistributed by this repository**, and none is relicensed by it.
Licences are as stated by their owners; where a licence could not be verified from an authoritative
source it is marked *unverified* rather than guessed.

---

## Software

| Name | Version used | Source | Licence | Retrieval |
|---|---|---|---|---|
| **P2Rank** | 2.5.1 | `https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz` | MIT | Download the pinned release tarball. Verify `SHA-256(p2rank_2.5.1.tar.gz) = d243f2d9036ac053fefb9407b5fe1c85f4fe077c519fd975ac585e995feab274` and `SHA-256(bin/p2rank.jar) = 4d73a85b796bd5ec5563d840abb5b1005f37b4651fd7f8aaad4b01a936ea1ece`. Stock default model; no retraining |
| **fpocket** | source tag `4.2.3`, commit `4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066` | `https://github.com/Discngine/fpocket.git` | GPL-3.0 (bundles qhull under its own licence) | `git clone`, `git checkout 4.2.3`, build with a **serial** `make` (a parallel build races on the bundled qhull subbuild). Verify `SHA-256(fpocket) = 90e2a3fc83ede6b06a06f779eb556596e6627cc2669fb0e2b42d9d46ac4ff271`. Note: the built binary's banner reports "fpocket 4.0"; the authoritative version is the git tag and commit |
| **Java (Temurin OpenJDK)** | 17.0.19+10 | Eclipse Adoptium | GPL-2.0 with Classpath Exception | Required by P2Rank |
| **Python** | 3.12.10 | python.org | PSF | |
| **gemmi** | 0.7.5 | PyPI | MPL-2.0 | mmCIF parsing, numbering correspondence, coordinate writing |
| **NumPy** | 2.5.2 | PyPI | BSD-3-Clause | |
| **SciPy** | 1.18.1 | PyPI | BSD-3-Clause | Wilson intervals, exact tests, linear-sum assignment |
| **statsmodels** | 0.15.0 | PyPI | BSD-3-Clause | |
| **Matplotlib** | 3.11.1 | PyPI | Matplotlib licence (BSD-style) | Figures |

### Citations for the software

- P2Rank: Krivák R, Hoksza D. *J Cheminform* 2018;10(1):39. [10.1186/s13321-018-0285-8](https://doi.org/10.1186/s13321-018-0285-8)
- fpocket: Le Guilloux V, Schmidtke P, Tuffery P. *BMC Bioinformatics* 2009;10:168. [10.1186/1471-2105-10-168](https://doi.org/10.1186/1471-2105-10-168)
- SciPy: Virtanen P, et al. *Nat Methods* 2020;17(3):261–272. [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2)

---

## Structural data

| Resource | Version / date | Source | Licence | Redistributed? |
|---|---|---|---|---|
| **Protein Data Bank** coordinates (mmCIF) | as deposited, retrieved 2026-09-15 | `https://files.rcsb.org/download/<PDB_ID>.cif` | wwPDB, CC0 1.0 | **No.** Every structure is referenced by PDB identifier in `data_manifest/`; coordinates are fetched at run time into a git-ignored cache |
| **SIFTS** residue-level UniProt↔PDB alignments | as served 2026-09-15 | RCSB Data API, `rcsb_polymer_entity_align` with `provenance_source = SIFTS` | PDBe / EMBL-EBI, CC BY 4.0 | **No.** Derived residue ranges and the SHA-256 of each API response are stored in `data_manifest/` and `annotations/` |
| **RCSB Search and Data APIs** | queried 2026-09-15 | `https://search.rcsb.org`, `https://data.rcsb.org` | RCSB terms | **No** |

### Citations

- Protein Data Bank: Berman HM, et al. *Nucleic Acids Res* 2000;28(1):235–242. [10.1093/nar/28.1.235](https://doi.org/10.1093/nar/28.1.235)
- SIFTS: Velankar S, et al. *Nucleic Acids Res* 2013;41(D1):D483–D489. [10.1093/nar/gks1258](https://doi.org/10.1093/nar/gks1258)

---

## Benchmark and corpus data (dataset provenance audit)

These were analysed but are **not** redistributed. `EXCLUDED_THIRD_PARTY.tsv` gives the exact byte
size and SHA-256 of each object as used, so a reproducer can confirm they obtained the identical
file. `scripts/23_dataset_membership.py` consumes them from `dataset_audit/membership/`.

| Resource | Version | Source | Licence | SHA-256 (first 16 hex) |
|---|---|---|---|---|
| CHEN11 membership (`chen11.ds`) | p2rank-datasets, master, accessed 2026-09-15 | `https://github.com/rdk/p2rank-datasets` | **No licence declared upstream** — see note | `91caf227eb213f5a` |
| JOINED (`joined.ds`, `joined(mlig).ds`) | same | same | **No licence declared upstream** | `749f23f5da3bac82`, `237286874af9dcff` |
| COACH420 (`coach420.ds`, `coach420(mlig).ds`) | same | same | **No licence declared upstream** | `4c33519448ae2df7`, `87d6e1724637adfa` |
| HOLO4K (`holo4k.ds`, `holo4k(mlig).ds`) | same | same | **No licence declared upstream** | `442a1abb07ee5ff9`, `28730da34f04848b` |
| FPTRAIN (`fptrain.ds`) | same | same | **No licence declared upstream** | `af16cd4819b5947c` |
| **PLINDER** split file | release 2024-06, schema v2 | `https://storage.googleapis.com/plinder/2024-06/v2/splits/split.parquet` | Apache-2.0 for PLINDER-curated data, per the upstream README — see note | `2959fb4b32f8c5cc` |
| **1DUG** coordinates (worked example) | as deposited | `https://files.rcsb.org/download/1DUG.cif` | wwPDB, CC0 1.0 | `c75dbb952653fa91` |

> **Note on licensing, re-checked 2026-09-16.** Nothing below is redistributed by this repository,
> and no licence is inferred for material this project does not own.
>
> **`rdk/p2rank-datasets` declares no licence.** The repository contains no `LICENSE`, `LICENSE.txt`
> or `COPYING` file, and the GitHub API reports no detected licence. That absence is a finding, not
> an omission by this project: it means no redistribution right is granted, so the membership files
> are referenced by source URL and SHA-256 and are **not** copied here. They are publicly readable
> and clonable without registration or login, which is what the reproducibility requirement needs.
>
> **PLINDER.** The upstream README states that data curated by PLINDER are made available under the
> Apache License 2.0; that is recorded here as the authoritative statement for the split file. An
> inconsistency exists upstream which this project neither resolves nor has authority to resolve:
> the README and badge say Apache-2.0, while `LICENSE.txt` on the default branch is GPL-2.0 and
> GitHub detects GPL-2.0 for the code. The README statement is the one addressing the *data*; the
> code licence is a separate question for the upstream project. The split file is served anonymously
> from a public bucket and needs no account.
>
> A reproducer should consult the upstream projects for their current terms before redistributing
> anything derived from them.

### Citations

- CHEN11 dataset: Chen K, Mizianty MJ, Gao J, Kurgan L. *Structure* 2011;19(5):613–621. [10.1016/j.str.2011.02.015](https://doi.org/10.1016/j.str.2011.02.015)
- COACH420, HOLO4K, FPTRAIN as distributed with P2Rank: Krivák & Hoksza 2018 (above)
- Binding MOAD (ground-truth annotation referenced for 1DUG): Hu L, Benson ML, Smith RD, Lerner MG, Carlson HA. *Proteins* 2005;60(3):333–340. [10.1002/prot.20512](https://doi.org/10.1002/prot.20512)
- PLINDER: Durairaj J, et al. *bioRxiv* 2024. [10.1101/2024.07.17.603955](https://doi.org/10.1101/2024.07.17.603955) — **preprint, not peer-reviewed**
- LIGYSIS and the benchmarking framework: Utgés JS, Barton GJ. *J Cheminform* 2024;16(1):126. [10.1186/s13321-024-00923-z](https://doi.org/10.1186/s13321-024-00923-z)
- sc-PDB (referenced in Limitations; **not measured** by this study): Desaphy J, Bret G, Rognan D, Kellenberger E. *Nucleic Acids Res* 2015;43(D1):D399–D404. [10.1093/nar/gku928](https://doi.org/10.1093/nar/gku928)

---

## Retrieval helper

`scripts_release/fetch_third_party.sh` in this package downloads the P2Rank release, clones and
builds fpocket at the pinned tag, and fetches the `.ds` membership files and the PLINDER split
file, verifying every SHA-256 against `EXCLUDED_THIRD_PARTY.tsv` and the values above. It does not
bypass any access control and does not redistribute anything; it simply automates lawful retrieval
from the sources named here.
