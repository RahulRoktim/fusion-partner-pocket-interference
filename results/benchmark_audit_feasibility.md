# Benchmark / Training-Set Construct Audit — FEASIBILITY ONLY

> **CORRECTION (2026-09-15, Phase 5).** The statement below that "fpocket is not an ML method and
> has no training set, so no training-contamination claim can attach to it" is **WRONG and
> withdrawn**. Inspection of fpocket 4.2.3 source (`src/pscoring.c`, commit `4bb0d84`) shows
> `score_pocket()` is a fitted linear model with hard-coded coefficients, documented in-source as
> derived by logistic regression / PLS from a training set. The fpocket paper reports the scoring
> function was fitted by PLS on **307 proteins**. See `results/PHASE5_DATASET_AUDIT.md` and
> `DEVIATIONS.md` H1.


**No contamination claim is made here, and none is permitted yet.** This document answers one
question: *could* the audit be done, and with what evidence. Nothing below reports how many
chimeric constructs any dataset contains — that measurement has not been run.

Three distinct statements must never be collapsed into one another:

| Statement | What it requires |
|---|---|
| **PRESENT IN DATASET** | the dataset's membership list contains a PDB entry that our SIFTS method annotates as chimeric |
| **USED TO TRAIN MODEL** | that model's *own documented training set* contains it — dataset name similarity is not evidence |
| **BIAS DEMONSTRATED** | a measured change in a reported performance figure attributable to those entries |

Only the first is testable with the method we already have. The second requires per-model
documentation. The third requires a separate experiment that has not been designed.

---

## Evidence table

| Dataset | Exact membership obtainable? | PDB IDs / chains available? | Role | Chimeric membership testable by our SIFTS method? | Documented use by pocket predictors |
|---|---|---|---|---|---|
| **HOLO4K** | **Yes.** `holo4k.ds` from `github.com/rdk/p2rank-datasets`, fetched 2026-09-15, HTTP 200, **4,009 lines** | PDB **ID only** (`holo4k/121p.pdb`); whole entries, multi-chain, no chain selector | **TEST / BENCHMARK** | **Yes**, at entry level. Chain-level attribution needs a rule because the list is entry-level and HOLO4K is explicitly multi-chain | Test set for **P2Rank**; widely reused as a benchmark by later methods |
| **COACH420** | **Yes.** `coach420.ds`, HTTP 200, **420 lines**; `coach420(mlig).ds` variant also present | PDB ID **+ chain** (`coach420/148lE.pdb` → 148L chain E) | **TEST / BENCHMARK** | **Yes**, at the exact chain level our pipeline already works in — the best-matched dataset of the five | Test set for **P2Rank**; widely reused |
| **CHEN11** | **Yes.** `chen11.ds`, HTTP 200, 258 lines incl. header/params | PDB ID + chain | **TRAINING** (P2Rank's own) | **Yes** | **Documented training set of the distributed default P2Rank model.** The only dataset here for which a "used to train" claim could be supported for P2Rank |
| **LIGYSIS** | **Partly.** Pipeline and web code public (`bartongroup/LIGYSIS*`); LIGYSIS-web hosts ~65,000 sites over ~25,000 proteins. A single flat membership file was not located in this feasibility pass | Per-protein/site records via the web resource; bulk export needs checking | **REFERENCE** dataset, used as an evaluation substrate in Utgés & Barton 2024 | **Likely yes**, conditional on obtaining a bulk membership export | Evaluation substrate for the 13-method benchmark; **not** a training set for any method here |
| **sc-PDB** | **Partly.** Distributed from `bioinfo-pharma.u-strasbg.fr/scPDB`; 2015 release registers 9,283 sites / 3,678 proteins. Access terms and current availability need checking | Entry-level with site definitions | **REFERENCE** resource, in practice the de-facto **TRAINING** source for several CNN predictors | **Yes** if the membership list can be obtained | Documented training data for **Kalasanty**, **PUResNet** (5,020 structures), **DeepSite**, **DeepSurf**, **EquiPocket**, **DeepPocket** |
| **PLINDER** | **Yes, and best-structured of the five.** `plinder-org/plinder`, versioned releases (2024-04/v0), explicit **train / validation / test** splits | System-level identifiers derived from PDB entries and chains | **TRAINING + VALIDATION + TEST**, explicitly separated by the authors | **Yes**, and the split labels allow each role to be reported separately — which is exactly the distinction the language constraints demand | Primarily co-folding / docking models rather than the two detectors in this study |

---

## Feasibility verdict

**Technically feasible for COACH420, HOLO4K, CHEN11 and PLINDER**; feasible for sc-PDB and LIGYSIS
**conditional on obtaining a bulk membership list**, which was not resolved in this pass.

Membership retrieval was verified live: the three P2Rank `.ds` files fetched successfully and
parse into PDB identifiers with no further tooling. Our existing SIFTS annotation path consumes
exactly these identifiers, so no new method is required — only the membership lists.

## Three methodological cautions to build into the audit design

1. **Standalone partner entries are not chimeric constructs.** These datasets contain structures of
   T4 lysozyme, MBP and cytochrome b562 *in their own right*. Our pipeline already separates these
   (single-UniProt entities fail criterion E1) from genuine chimeras (exactly two UniProt
   accessions). Conflating them would manufacture a contamination figure out of ordinary entries.
   Any audit must report the two categories separately.
2. **Entry-level vs chain-level lists differ.** COACH420 names a chain; HOLO4K names an entry. A
   chimeric entity in a HOLO4K entry does not imply the benchmarked chain was chimeric. The audit
   needs a pre-specified attribution rule, and entry-level and chain-level results must be reported
   apart.
3. **"P2Rank was trained on CHEN11" is a supportable training claim** for the detectors in this
   study. ~~fpocket is not a machine-learning method at all — it has no training set, so no
   training-contamination claim can apply to it.~~ **[CORRECTED — see the note at the top of this
   file. fpocket 4.2.3 carries a fitted, data-derived scoring function; its scoring function was
   fitted by PLS on 307 proteins per the fpocket paper, and `FPTRAIN` (222 entries) is a
   third-party characterisation of that training set whose exact identity is not established.]**
   Any statement about sc-PDB-trained CNNs concerns *other* predictors that this study does not run.

## What is still required before any claim

- A pre-specified attribution rule for entry-level lists.
- Confirmed bulk membership access for sc-PDB and LIGYSIS, or explicit exclusion of those two.
- A separate, designed experiment before the word "bias" is used at all.
