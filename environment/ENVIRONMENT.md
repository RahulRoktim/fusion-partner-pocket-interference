# Environment Specification (pilot, protocol v1.1)

Captured 2026-09-15. Machine-readable copy: `environment/env_lock.json`.

## Host

| Component | Value |
|---|---|
| OS | Windows 11 Pro 10.0.26200 |
| Shell used for orchestration | Git Bash (MINGW64_NT-10.0-26200, msys 3.6.4) |
| Linux subsystem | WSL2, Ubuntu 24.04.2 LTS, kernel 6.18.33.2-microsoft-standard-WSL2 |
| Container runtime | Docker 29.7.2 present but **daemon not running**; not used |

## Language runtimes

| Component | Version | Path |
|---|---|---|
| Python | 3.12.10 | `C:\Users\user\AppData\Local\Programs\Python\Python312\python.exe` |
| Java (JDK) | Temurin OpenJDK 17.0.19+10 | `C:\Program Files\Eclipse Adoptium\jdk-17.0.19.10-hotspot` |
| gcc (WSL) | 13.3.0 (Ubuntu 13.3.0-6ubuntu2~24.04.1) | `/usr/bin/gcc` |

**JAVA_HOME caveat.** The inherited `JAVA_HOME` points at
`C:\Users\user\AppData\Local\Programs\Eclipse Adoptium\jdk-25.0.2.10-hotspot\`, which **does not
exist** on this machine. P2Rank's `prank` launcher honours `JAVA_HOME` and fails with
"No such file or directory" unless it is overridden. All run scripts explicitly export
`JAVA_HOME="/c/Program Files/Eclipse Adoptium/jdk-17.0.19.10-hotspot"`. Anyone reproducing this
must do the same or fix their `JAVA_HOME`.

## Python libraries

| Package | Version | Role |
|---|---|---|
| gemmi | 0.7.5 | mmCIF parsing, entity/auth numbering correspondence, coordinate writing |
| numpy | 2.5.2 | numerics |
| scipy | 1.18.1 | statistics (Wilson CIs, exact tests) |

## Detectors

### P2Rank — **2.5.1** (pinned)

| Item | Value |
|---|---|
| Source | `https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz` |
| Tarball SHA-256 | `d243f2d9036ac053fefb9407b5fe1c85f4fe077c519fd975ac585e995feab274` |
| `bin/p2rank.jar` SHA-256 | `4d73a85b796bd5ec5563d840abb5b1005f37b4651fd7f8aaad4b01a936ea1ece` |
| Model | stock/default (`default.model`); no alternative model, no retraining |
| Self-reported version | `P2Rank 2.5.1` |

2.6-alpha exists (2026-05-20) and was **not** used: pre-release software is inappropriate for a
reproducibility-critical study.

### fpocket — source tag **4.2.3** (pinned), built from source

| Item | Value |
|---|---|
| Source | `https://github.com/Discngine/fpocket.git` |
| Tag | `4.2.3` |
| Commit | `4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066` |
| Build | WSL2 Ubuntu 24.04, gcc 13.3.0, **serial `make`** |
| Install | user prefix `$HOME/opt/fpocket/bin/fpocket`; **no root, no sudo** |
| Binary SHA-256 | `90e2a3fc83ede6b06a06f779eb556596e6627cc2669fb0e2b42d9d46ac4ff271` |
| Parameters | stock defaults; no flags beyond `-f <input.pdb>` |

**Two build facts worth recording.**

1. **Parallel make fails.** `make -j4` dies with
   `src/qhull/src/qvoronoi/qvoronoi.c:20:10: fatal error: libqhull/libqhull.h: No such file or
   directory` — a race in fpocket's makefile between the bundled-qhull subbuild and the top-level
   compile. **Serial `make` succeeds.** Reproducers must not use `-j`.
2. **The binary under-reports its own version.** `fpocket -h` prints the banner `fpocket 4.0`
   although it is built from tag 4.2.3. The banner string is stale upstream. The authoritative
   version for this study is the **git tag and commit above**, not the banner.

### Why source-build rather than package or container

`sudo` in WSL requires an interactive password (unavailable to an automated run), so `apt install`
was not possible; and the Docker daemon was not running, so a container build was not possible
without user action. Building from a pinned tag into a user prefix is fully reproducible and
required neither. Recorded as amendment **B1** in `DEVIATIONS.md`. **fpocket was not substituted,
downgraded, or dropped.**

## Command templates (exact, as executed)

P2Rank:

```bash
export JAVA_HOME="/c/Program Files/Eclipse Adoptium/jdk-17.0.19.10-hotspot"
"<tools>/p2rank_2.5.1/prank" predict -f <input.pdb> -o <outdir> -threads 1
```

fpocket (invoked from Windows through WSL; input path translated to `/mnt/c/...`):

```bash
wsl.exe -d Ubuntu -- bash -lc '$HOME/opt/fpocket/bin/fpocket -f <input_wsl_path.pdb>'
```

`-threads 1` is set for P2Rank so that results are deterministic and runtimes comparable; it is not
a tuning parameter and is identical for every structure and both conditions.

## Determinism

Neither detector uses a random seed at prediction time (P2Rank applies a pre-trained model;
fpocket is deterministic geometry). Cohort **sampling** uses Python `random.Random(20260915)`.
Every detector run records input SHA-256, output SHA-256, full command, return code and runtime.
