# Reproduction guide

Written for a third party starting from a fresh clone, with no access to the original machine.
Every step names the script that performs it. Where the original run recorded a measured value —
a version, a hash, a runtime — that value is given so you can confirm you obtained the same thing.

> **Three routes.** §A verifies the archived/current package without detectors. §B recomputes
> analyses from committed classified/intermediate results without detector binaries. §C performs
> full raw reproduction and requires P2Rank, fpocket, prepared structures and third-party inputs.
> The commands are separated so “detector-free” never silently launches a detector.

---

## 0. Assumptions

| | |
|---|---|
| OS | The analysis scripts are OS-independent (Python 3.11+). **fpocket requires Linux or macOS**; the original study built it under WSL2 Ubuntu 24.04 on Windows 11 |
| Shell | POSIX shell for the commands below |
| Network | Required for full reproduction §C (RCSB, GitHub, Google Storage) |
| Disk | ≈ 2.8 GB for a full reproduction |
| Privileges | None. fpocket is built into a user prefix; no `sudo` is needed |

---

## A. Archive and package integrity — no detectors

```bash
git clone https://github.com/RahulRoktim/fusion-partner-pocket-interference fusion-tag-hazard
cd fusion-tag-hazard

python scripts_release/verify_freeze.py --package
python scripts_release/verify_freeze.py
python tests/test_no_development_leakage.py
```

These commands use the Python standard library, read existing files, and launch no detector. The
first verifies the current checkout against `PACKAGE_MANIFEST.tsv`. The second verifies the frozen
scientific state; post-release maintenance code is checked against the immutable
`v1.0.0-submission` tag, while public-facing replacements and non-redistributed third-party inputs
are reported separately. The leakage check verifies the exploratory/confirmatory cohort boundary.

`tests/test_variant_classifier_equivalence.py` is intentionally **not** in this route. It requires
the archived raw detector trees. Without them it exits 2 with `UNAVAILABLE`; it can never report
PASS from empty parser outputs.

---

## B. Recompute from committed classified/intermediate results — no detectors

Use a disposable clone or expect regenerated tracked artifacts to appear in `git diff`. Install the
analysis-only dependencies, then run only the scripts below:

```bash
python -m pip install "numpy==2.5.2" "scipy==1.18.1" "statsmodels==0.15.0" \
                      "matplotlib==3.11.1" "gemmi==0.7.5"

python scripts/16_confirmatory_analysis.py     # Aims A–E from frozen classified pockets
python scripts/18_replication.py               # development vs confirmatory, never pooled
python scripts/19_secondary_model.py           # secondary multivariable model
python scripts/29_canonical_tables.py          # canonical tables
python scripts/30_figures_main.py              # main figures
python scripts/32_figures_supp.py              # supplementary figures
```

None of those six scripts invokes P2Rank or fpocket. This route does **not** include
`scripts/17_confirmatory_mechanism.py`: that script creates isolated-partner structures and launches
both detectors when their outputs are absent. It also excludes `tests/test_mapping_invariants.py`,
`tests/test_detector_parsers.py` and `tests/test_variant_classifier_equivalence.py`, which need
structures, detector installations or archived raw detector outputs.

---

## C. Full raw detector reproduction

The remainder of this guide rebuilds from PDB identifiers upward. It requires network retrieval,
prepared structures, P2Rank 2.5.1, fpocket 4.2.3, Java 17 and approximately 2.8 GB of disk. After
the detector stages have restored the hash-pinned raw output trees, run:

```bash
python scripts/17_confirmatory_mechanism.py
python tests/test_mapping_invariants.py            # residue-mapping invariants
python tests/test_detector_parsers.py              # detector output parsing
python tests/test_variant_classifier_equivalence.py
```

An equivalence PASS requires verified raw trees and a nonzero pocket comparison count.

---

### C1. Retrieve third-party objects

Nothing in §C1 is redistributed by this repository. The helper verifies every checksum.

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

### C2. Java and environment

```bash
export JAVA_HOME=/path/to/temurin-jdk-17     # P2Rank's launcher honours JAVA_HOME and will fail
java -version                                # with "No such file or directory" if it points at an
                                             # absent JDK. Confirm 17.x before continuing.
export FUSIONTAG_WORK=/abs/path/to/scratch   # working directory for detector runs
```

`FUSIONTAG_WORK` is optional. If unset, detector scripts use the platform temporary directory under
`fusiontag_hazard/`; if set, it is expanded and resolved to an absolute path. The confirmatory,
mechanism and E9 scripts create separate subdirectories beneath it. This changes only scratch-file
placement, never scientific parameters, prepared inputs or hash-recorded outputs.

### C3. Cohort construction

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

### C4. Structure preparation

```bash
python scripts/03_prepare_structures.py       # development cohort
python scripts/20_prepare_variants.py         # ions and biological-assembly-1 variants
```

Single designated chimeric chain, first model, altloc A, all non-polymer entities stripped,
terminal tag remnants deleted from both arms. The fusion-removed counterpart deletes fusion and
linker residues and nothing else; no target atom coordinate changes.

### C5. Detector execution

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

### C6. Pocket classification

```bash
python scripts/05_classify_pockets.py
python scripts/21_classify_assembly1.py       # chain-aware classifier for the assembly variant
```

Dominance 0.70 primary (0.50 and 0.90 as sensitivities); interface floor 0.20;
*fusion-associated* = fusion-dominated ∨ interface ∨ linker.

### C7. Correspondence and endpoints

```bash
python scripts/08_pocket_correspondence.py
python scripts/22_assembly1_displacement.py
python scripts/16_confirmatory_analysis.py
```

Correspondence: target residues only; minimum 3 target residues; Jaccard on target-residue sets;
target-residue heavy-atom centroid distance; one-to-one assignment by `scipy.optimize.
linear_sum_assignment` maximising total Jaccard; matched at Jaccard ≥ 0.40 **and** centroid
distance ≤ 8.0 Å. Detector rank and score are never inputs to the matching.

### C8. Sensitivity analyses

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

### C9. Dataset provenance audit

```bash
python scripts/23_dataset_membership.py
python scripts/24_dataset_construct_scan.py
python scripts/25_dataset_groundtruth.py
python scripts/26_plinder_audit.py
python scripts/27_benchmark_impact.py
python scripts/28_publication_denominators.py     # asserts 11 prose denominators; must print 11/11
```

Requires the third-party membership files from §C1.

### C10. Figures and tables

```bash
python scripts/29_canonical_tables.py
python scripts/30_figures_main.py
python scripts/31_figure5_structural.py     # needs the structure cache from §C3-C4
python scripts/32_figures_supp.py
python scripts/33_claims_ledger.py
```

`figures/FIG5_pymol.pml` reproduces the four structural cases in PyMOL at publication quality;
substituting those renderings changes presentation only.

### C11. Verify the freeze

```bash
python - <<'EOF'
import json, hashlib
m = json.load(open("environment/FINAL_FREEZE_v1.1.json"))
bad = [p for g in m["sha256"].values() for p, h in g.items()
       if hashlib.sha256(open(p, "rb").read()).hexdigest() != h["sha256"]]
print("files:", sum(len(g) for g in m["sha256"].values()), "mismatches:", bad or "none")
EOF
```

Expect 224 files and no mismatches once §C1 has restored the third-party objects. Each text file also
carries `sha256_lf`, a line-ending-normalised hash that verifies identically on Linux, macOS and
Windows; prefer it for cross-platform checks.

---

## Figure and table numbering

Manuscript figure numbers differ from the historical file names, because figures were renumbered for
submission so that they are cited in ascending order. The authoritative mapping is in
`manuscript-support/FIGURE_LEGENDS.md`. No image file was renamed.
