# Fusion-Tag Hazard

Do crystallographic fusion partners systematically bias automated binding-pocket prediction by
introducing high-ranking pockets on the fusion partner or at the fusion interface, rather than on
the intended target protein?

**Status: Confirmatory analysis complete (n = 128, target-disjoint). PARTIAL CONFIRMATION.**
Awaiting external scientific review. No manuscript written; no benchmark audit started; no third
detector added.

The 24 pilot structures are permanently DEVELOPMENT/EXPLORATORY data and were never pooled into
confirmatory inference. Confirmatory cohort: 128 structures / 123 independent target clusters,
zero target-accession overlap with development.
Not a software product — a computational research project intended for a preprint.

## Where things are

| Path | Contents |
|---|---|
| `PHASE0_LITERATURE_AUDIT.md` | Prior-art audit and novelty verdict (MODERATE) |
| `PROTOCOL.md` | Part I = development protocol v1.1; **Part II = confirmatory protocol v1.2** |
| `PROTOCOL_v1.1_ARCHIVED.md` | v1.1, preserved verbatim |
| `PROTOCOL_v1.0_ARCHIVED.md` | Original v1.0, preserved verbatim |
| `DEVIATIONS.md` | Append-only amendment log; each entry marked PRE- or POST-OUTCOME |
| `data_manifest/` | Candidate pool, eligibility verdicts with per-entity exclusion reasons, seeded selection, pilot manifest |
| `scripts/` | Numbered pipeline, 01 → 07 |
| `tests/` | Mapping and parser invariants; must pass before any pilot run |
| `results/` | Detector run records, classified pockets, summary statistics |
| `results/qc/` | Phase 2F mechanistic QC + revised isolated-partner probe |
| `results/CONFIRMATORY_REPORT.md` | **Confirmatory results and final decision** |
| `results/PHASE5_DATASET_AUDIT.md` | Dataset provenance / training / benchmark construct audit |
| `dataset_audit/` | Frozen dataset membership, construct scan, PLINDER audit |
| `results/confirmatory/` | Confirmatory runs and analyses (primary, ions, assembly1) |
| `results/PHASE3A_REPORT.md` | Phase 3A methods freeze and confirmatory cohort |
| `environment/PROTOCOL_FREEZE.json` | Hash-tagged frozen protocol state |
| `results/e9_audit/` | E9 filter audit with independent coordinate validation |
| `results/benchmark_audit_feasibility.md` | Dataset audit feasibility (no contamination claim) |
| `environment/` | Exact versions, hashes, command templates, known traps |
| `REPRODUCE.md` | Step-by-step rerun instructions |
| `DATA_LICENSES.md` | Sources and licensing; coordinates are never redistributed |

## Pilot headline

24 structures (8 BRIL, 8 T4L, 8 MBP), drawn by seeded outcome-blind sampling from 455 eligible
chimeric entities; 96 detector runs, all successful.

| | P2Rank 2.5.1 | fpocket 4.2.3 |
|---|---|---|
| rank-1 pocket is fusion-associated | 11/24 = 46% [28–65] | 15/24 = 63% [43–79] |
| any fusion-associated pocket in top-5 | 20/24 = 83% | 22/24 = 92% |
| rank-1 becomes target-dominated after in-silico fusion removal | 5/24, 0 reverse | 15/24, 0 reverse |

Mechanism: 21 of 26 rank-1 failures are cavities **intrinsic to the fusion partner** — in six MBP
structures the top-ranked pocket contains deposited maltose, MBP's own natural ligand. Four are
target–fusion interface cavities. None was a mapping or parser artifact.

## What this does not claim

Not that all fusion tags cause artifacts; not that pocket predictors are generally unreliable; not
that any published docking result is invalid; not that any biological ligand site is wrong. Claim
boundaries are binding and stated in `PROTOCOL.md` §11.
