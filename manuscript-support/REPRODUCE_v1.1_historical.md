# Reproducing the Pilot

Total runtime on the reference machine: roughly 2 hours, of which ~50 minutes is the one-off RCSB
pool enumeration (cached afterwards in `data_manifest/pool_raw.json`).

## 0. Prerequisites

- Windows with **WSL2** (Ubuntu 24.04 used here) — required for fpocket, which has no supported
  native Windows build. A plain Linux host works too; drop the `wsl.exe` prefix in
  `scripts/common_runner.py`.
- **JDK 17+** for P2Rank.
- **Python 3.12** with `gemmi==0.7.5 numpy scipy`.

Full versions, hashes and caveats: `environment/ENVIRONMENT.md`.

> **JAVA_HOME trap.** P2Rank's launcher honours `JAVA_HOME`. On the reference machine the inherited
> value pointed at a JDK that is not installed. `scripts/common_runner.py` invokes the JVM directly
> and sets `JAVA_HOME` explicitly; override with `FUSIONTAG_JAVA_HOME` if your JDK lives elsewhere.

## 1. Install the detectors

```bash
mkdir -p tools && cd tools
curl -sSL -O https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz
# expect sha256 d243f2d9036ac053fefb9407b5fe1c85f4fe077c519fd975ac585e995feab274
tar xzf p2rank_2.5.1.tar.gz
```

fpocket, inside WSL — **serial `make` only**; `make -j` hits a race in fpocket's makefile between
the bundled-qhull subbuild and the top-level compile:

```bash
git clone https://github.com/Discngine/fpocket.git ~/opt/fpocket
cd ~/opt/fpocket && git checkout tags/4.2.3 && make
# expect bin/fpocket sha256 90e2a3fc83ede6b06a06f779eb556596e6627cc2669fb0e2b42d9d46ac4ff271
```

## 2. Validate before running anything

Both suites must pass. They are not optional: check **D3** catches an fpocket pocket-file
off-by-one that would silently shift every fpocket classification by one rank.

```bash
python tests/test_mapping_invariants.py   # 51 checks
python tests/test_detector_parsers.py     # 21 checks
```

## 3. Run the pilot

```bash
python scripts/01_build_pool.py        # ~50 min first run; cached thereafter
python scripts/02_select_pilot.py      # seeded draw, deterministic
python scripts/03_prepare_structures.py
python scripts/04_run_detectors.py     # 96 runs
python scripts/05_classify_pockets.py
python scripts/06_analyze_pilot.py
python scripts/07_qc_mechanism.py      # Phase 2F, QC only
```

## 4. Determinism

The cohort draw is `random.Random(f"20260915-{PARTNER}").sample(sorted(targets), 8)` — stable
across machines and Python builds. Neither detector randomises at prediction time. Every detector
run records input hash, output hash, command, return code and runtime in
`results/detector_runs.json`.

**One caveat on exact reproducibility:** `scripts/01_build_pool.py` queries the *live* PDB. As new
chimeric structures are deposited the eligible pool grows, which can change the draw. To reproduce
this exact cohort, use the committed `data_manifest/pool_raw.json`
(sha256 `0f6989e3c0ddd3c53501833318f4cca04ee8e60952c8e622babc7c3fcdf32e1c`, snapshot 2026-09-15),
which the script prefers over refetching when present.
