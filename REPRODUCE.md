# Reproduction guide

Written for a third party starting from a fresh clone, with no access to the original machine.
Every step names the script that performs it. Where the original run recorded a measured value —
a version, a hash, a runtime — that value is given so you can confirm you obtained the same thing.

> **Two levels.** §A verifies every number in the manuscript in minutes without installing a
> detector. §B–§L reproduce the study from PDB identifiers upward. Most readers want §A.

---

## 0. Assumptions

| | |
|---|---|
| OS | The analysis scripts are OS-independent (Python 3.11+). **fpocket requires Linux or macOS**; the original study built it under WSL2 Ubuntu 24.04 on Windows 11 |
| Shell | POSIX shell for the commands below |
| Network | Required for §B (RCSB, GitHub, Google Storage) |
| Disk | ≈ 2.8 GB for a full reproduction |
| Privileges | None. fpocket is built into a user prefix; no `sudo` is needed |

---

## A. Verify the reported numbers without running detectors

```bash
git clone https://github.com/RahulRoktim/fusion-partner-pocket-interference fusion-tag-hazard
cd fusion-tag-hazard
python -m pip install "numpy==2.5.2" "scipy==1.18.1" "statsmodels==0.15.0" \
                      "matplotlib==3.11.1" "gemmi==0.7.5"

python scripts/16_confirmatory_analysis.py     # Aims A–E from the frozen classified pockets
python scripts/17_confirmatory_mechanism.py    # mechanism classification
python scripts/18_replication.py               # development vs confirmatory, never pooled
python scripts/19_secondary_model.py           # secondary multivariable model
python scripts/29_canonical_tables.py          # Tables 1, 2, 3, 3b, 4, S1, S2, S3
python scripts/30_figures_main.py              # main figures
python scripts/32_figures_supp.py              # supplementary figures
```

Then run the test suites:

```bash
python tests/test_mapping_invariants.py            # residue-mapping invariants
python tests/test_detector_parsers.py              # detector output parsing
python tests/test_no_development_leakage.py        # cohort disjointness, target overlap must be 0
python tests/test_variant_classifier_equivalence.py
```

All four must pass. Runtime for §A: a few minutes. Storage: ~120 MB.

---

## B. Retrieve third-party objects

Nothing in §B is redistributed by this repository. The helper verifies every checksum.

```bash
bash scripts_release/fetch_third_party.sh
```

It performs, and you can equally do by hand:

**P2Rank 2.5.1**

```bash
mkdir -p tools && cd tools
curl -sSL -O https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz
# expect SHA-256 d243f2d9036ac053fefb9407b5fe1c85f4fe077c519fd975ac585e995feab274
tar xzf p2rank_2.5.1.tar.gz
# expect SHA-256(p2rank_2.5.1/bin/p2rank.jar)
#   = 4d73a85b796bd5ec5563d840abb5b1005f37b4651fd7f8aaad4b01a936ea1ece
cd ..
```

**fpocket 4.2.3**

```bash
git clone https://github.com/Discngine/fpocket.git
cd fpocket && git checkout 4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066
make            # SERIAL. `make -j` races on the bundled qhull subbuild and fails with
                # "libqhull/libqhull.h: No such file or directory"
make install PREFIX=$HOME/opt/fpocket
# expect SHA-256($HOME/opt/fpocket/bin/fpocket)
#   = 90e2a3fc83ede6b06a06f779eb556596e6627cc2669fb0e2b42d9d46ac4ff271
cd ..
```

> The built binary prints the banner "fpocket 4.0". That string is stale upstream; the authoritative
> version is the git tag and commit above.

**Benchmark membership files and the PLINDER split file**

```bash
# .ds membership lists
git clone --depth 1 https://github.com/rdk/p2rank-datasets
cp p2rank-datasets/{chen11,joined,coach420,holo4k,fptrain}*.ds dataset_audit/membership/

# PLINDER 2024-06/v2 split
curl -sSL -o dataset_audit/membership/plinder_2024-06_v2_split.parquet \
  https://storage.googleapis.com/plinder/2024-06/v2/splits/split.parquet
# expect SHA-256 prefix 2959fb4b32f8c5cc...

# 1DUG coordinates for the worked example
curl -sSL -o /tmp/1DUG.cif https://files.rcsb.org/download/1DUG.cif
```

Verify every hash against `EXCLUDED_THIRD_PARTY.tsv`. After this step,
`environment/FINAL_FREEZE_v1.1.json` verifies completely.

## C. Java and environment

```bash
export JAVA_HOME=/path/to/temurin-jdk-17     # P2Rank's launcher honours JAVA_HOME and will fail
java -version                                # with "No such file or directory" if it points at an
                                             # absent JDK. Confirm 17.x before continuing.
export FUSIONTAG_WORK=/abs/path/to/scratch   # working directory for detector runs
```

> **Known portability defect.** `scripts/04_run_detectors.py` reads `FUSIONTAG_WORK` from the
> environment, but `scripts/15_run_confirmatory.py`, `scripts/17_confirmatory_mechanism.py` and
> `scripts/35_e9_run_new.py` contain a hard-coded Windows scratch path in their `WORK =` assignment
> (lines 23, 18 and 22 respectively). These scripts are part of the frozen scientific record and
> have deliberately **not** been edited. To run them elsewhere, change that one assignment in each
> to a writable local path. Doing so changes no scientific behaviour: `WORK` is a scratch directory
> for detector output only, and every output is hash-recorded.

## D. Cohort construction

```bash
python scripts/01_build_pool.py          # RCSB/SIFTS enumeration + eligibility (E1-E9)
python scripts/02_select_pilot.py        # 24-structure development cohort, fixed seed 20260915
python scripts/13_freeze_sets.py
python scripts/14_final_confirmatory_set.py   # target-disjoint confirmatory set, n = 128
```

The distributed `data_manifest/` already contains the frozen outputs of these steps. Re-running
queries RCSB live; because SIFTS is re-released weekly and PDB entries can be superseded, **a later
run is not guaranteed to reproduce the archived manifests byte for byte**. The archived manifests
are the authoritative cohort record; treat any divergence as an upstream change, not an error.

## E. Structure preparation

```bash
python scripts/03_prepare_structures.py       # development cohort
python scripts/20_prepare_variants.py         # ions and biological-assembly-1 variants
```

Single designated chimeric chain, first model, altloc A, all non-polymer entities stripped,
terminal tag remnants deleted from both arms. The fusion-removed counterpart deletes fusion and
linker residues and nothing else; no target atom coordinate changes.

## F. Detector execution

```bash
python scripts/04_run_detectors.py                       # development cohort
FUSIONTAG_VARIANT=primary   python scripts/15_run_confirmatory.py
FUSIONTAG_VARIANT=ions      python scripts/15_run_confirmatory.py
FUSIONTAG_VARIANT=assembly1 python scripts/15_run_confirmatory.py
```

512 runs for the primary variant (128 structures × 2 conditions × 2 detectors). Measured runtime on
the original hardware: 21.0 min for P2Rank and 4.0 min for fpocket, single-threaded. Each run
records its command, versions, input and output SHA-256, stdout/stderr, return code and runtime to
`results/confirmatory/<variant>/detector_runs.json`. Compare against `results/RAW_OUTPUT_HASHES.tsv`.

## G. Pocket classification

```bash
python scripts/05_classify_pockets.py
python scripts/21_classify_assembly1.py       # chain-aware classifier for the assembly variant
```

Dominance 0.70 primary (0.50 and 0.90 as sensitivities); interface floor 0.20;
*fusion-associated* = fusion-dominated ∨ interface ∨ linker.

## H. Correspondence and endpoints

```bash
python scripts/08_pocket_correspondence.py
python scripts/22_assembly1_displacement.py
python scripts/16_confirmatory_analysis.py
```

Correspondence: target residues only; minimum 3 target residues; Jaccard on target-residue sets;
target-residue heavy-atom centroid distance; one-to-one assignment by `scipy.optimize.
linear_sum_assignment` maximising total Jaccard; matched at Jaccard ≥ 0.40 **and** centroid
distance ≤ 8.0 Å. Detector rank and score are never inputs to the matching.

## I. Sensitivity analyses

Executed by `scripts/16_confirmatory_analysis.py` and, for the disorder axis,
`scripts/34_e9_sensitivity_cohort.py` → `35_e9_run_new.py` → `36_e9_compare.py`.

| Axis | Settings | Status |
|---|---|---|
| Pocket dominance | 0.50 / **0.70 primary** / 0.90 | executed |
| Correspondence | J ≥ 0.25 & 12 Å / **J ≥ 0.40 & 8 Å primary** / J ≥ 0.60 & 5 Å | executed |
| Correspondence | **overlap coefficient ≥ 0.50** | **PRE-SPECIFIED BUT NOT EXECUTED DUE TO INCOMPLETE SPECIFICATION** |
| Deletion boundary | exclude removed-condition cavities within 8 Å of a cut | executed |
| Structure representation | **single designated chain primary** / biological assembly 1 | executed |
| HETATM handling | **stripped primary** / ions retained | executed |
| Segment disorder | **0.20 primary** / 0.35 | executed |

> **The overlap-coefficient variant is deliberately not implemented in this repository.** The
> protocol names the metric and the 0.50 threshold but never defines the accompanying
> centroid-distance criterion, the assignment objective, or whether the overlap test replaces or
> supplements the Jaccard acceptance rule. The omission was discovered after the confirmatory
> results were known, so supplying those choices retrospectively would introduce post-outcome
> analytical discretion. Its result is unknown, and no implementation is provided here. See
> `POST_FREEZE_DEVIATIONS.md` PF1. If you implement it yourself, please report it as your own
> post-hoc analysis, not as this study's pre-specified sensitivity.

## J. Dataset provenance audit

```bash
python scripts/23_dataset_membership.py
python scripts/24_dataset_construct_scan.py
python scripts/25_dataset_groundtruth.py
python scripts/26_plinder_audit.py
python scripts/27_benchmark_impact.py
python scripts/28_publication_denominators.py     # asserts 11 prose denominators; must print 11/11
```

Requires the third-party membership files from §B.

## K. Figures and tables

```bash
python scripts/29_canonical_tables.py
python scripts/30_figures_main.py
python scripts/31_figure5_structural.py     # needs the structure cache from §D-E
python scripts/32_figures_supp.py
python scripts/33_claims_ledger.py
```

`figures/FIG5_pymol.pml` reproduces the four structural cases in PyMOL at publication quality;
substituting those renderings changes presentation only.

## L. Verify the freeze

```bash
python - <<'EOF'
import json, hashlib
m = json.load(open("environment/FINAL_FREEZE_v1.1.json"))
bad = [p for g in m["sha256"].values() for p, h in g.items()
       if hashlib.sha256(open(p, "rb").read()).hexdigest() != h["sha256"]]
print("files:", sum(len(g) for g in m["sha256"].values()), "mismatches:", bad or "none")
EOF
```

Expect 224 files and no mismatches once §B has restored the third-party objects. Each text file also
carries `sha256_lf`, a line-ending-normalised hash that verifies identically on Linux, macOS and
Windows; prefer it for cross-platform checks.

---

## Figure and table numbering

Manuscript figure numbers differ from the historical file names, because figures were renumbered for
submission so that they are cited in ascending order. The authoritative mapping is in
`manuscript-support/FIGURE_LEGENDS.md`. No image file was renamed.
