#!/usr/bin/env bash
# Retrieve the third-party objects this repository does not redistribute, and verify each against
# the checksum recorded when the study was run. Nothing here bypasses access control; every source
# is a public URL named in THIRD_PARTY.md. No file fetched here is relicensed by this repository.
#
# Usage:  bash scripts_release/fetch_third_party.sh [--skip-fpocket]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
MEM="dataset_audit/membership"
mkdir -p tools "$MEM"

sha() { python -c "import hashlib,sys;print(hashlib.sha256(open(sys.argv[1],'rb').read()).hexdigest())" "$1"; }
check() {  # check <file> <expected-sha256> <label>
  local got; got="$(sha "$1")"
  if [ "$got" = "$2" ]; then
    echo "  OK    $3"
  else
    echo "  FAIL  $3"; echo "        expected $2"; echo "        got      $got"; exit 1
  fi
}

# ----------------------------------------------------------------- detectors
echo "== P2Rank 2.5.1 =="
[ -f tools/p2rank_2.5.1.tar.gz ] || curl -fsSL -o tools/p2rank_2.5.1.tar.gz \
  https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz
check tools/p2rank_2.5.1.tar.gz \
  d243f2d9036ac053fefb9407b5fe1c85f4fe077c519fd975ac585e995feab274 "p2rank_2.5.1.tar.gz"
tar xzf tools/p2rank_2.5.1.tar.gz -C tools
check tools/p2rank_2.5.1/bin/p2rank.jar \
  4d73a85b796bd5ec5563d840abb5b1005f37b4651fd7f8aaad4b01a936ea1ece "p2rank.jar"

echo "== fpocket 4.2.3 =="
if [ "${1:-}" = "--skip-fpocket" ]; then
  echo "  skipped on request"
else
  [ -d tools/fpocket ] || git clone -q https://github.com/Discngine/fpocket.git tools/fpocket
  (
    cd tools/fpocket
    git checkout -q 4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066
    make                        # SERIAL. `make -j` races on the bundled qhull subbuild.
    make install PREFIX="$HOME/opt/fpocket"
  )
  check "$HOME/opt/fpocket/bin/fpocket" \
    90e2a3fc83ede6b06a06f779eb556596e6627cc2669fb0e2b42d9d46ac4ff271 "fpocket binary"
fi

# ----------------------------------------------------------------- datasets
echo "== benchmark membership lists =="
# -c core.autocrlf=false is REQUIRED, not cosmetic. The recorded SHA-256 values are of the
# upstream bytes, which are LF. With core.autocrlf=true (the Windows default) git checks the
# .ds files out as CRLF and every checksum below fails, wrongly reporting upstream change.
[ -d tools/p2rank-datasets ] || \
  git -c core.autocrlf=false clone -q --depth 1 \
    https://github.com/rdk/p2rank-datasets tools/p2rank-datasets
for f in chen11.ds joined.ds "joined(mlig).ds" coach420.ds "coach420(mlig).ds" \
         holo4k.ds "holo4k(mlig).ds" fptrain.ds; do
  cp "tools/p2rank-datasets/$f" "$MEM/$f"
done

echo "== PLINDER 2024-06/v2 split =="
[ -f "$MEM/plinder_2024-06_v2_split.parquet" ] || curl -fsSL \
  -o "$MEM/plinder_2024-06_v2_split.parquet" \
  https://storage.googleapis.com/plinder/2024-06/v2/splits/split.parquet

echo "== 1DUG coordinates (worked example) =="
# The study used the benchmarked chain extracted from PDB entry 1DUG. Fetch the deposited entry;
# scripts/25_dataset_groundtruth.py derives the benchmarked chain from it.
[ -f "$MEM/1DUG.cif" ] || curl -fsSL -o "$MEM/1DUG.cif" https://files.rcsb.org/download/1DUG.cif

# ----------------------------------------------------------------- verify
echo
echo "== verifying every retrieved object against EXCLUDED_THIRD_PARTY.tsv =="
python - <<'PY'
import csv, hashlib, os, sys
missing, bad = [], []
for r in csv.DictReader(open("EXCLUDED_THIRD_PARTY.tsv", encoding="utf-8"), delimiter="\t"):
    p = r["REPO_PATH"]
    if not os.path.exists(p):
        print(f"  MISSING  {p}")
        missing.append(p)
        continue
    got = hashlib.sha256(open(p, "rb").read()).hexdigest()
    ok = got == r["SHA256"]
    print(f"  {'OK   ' if ok else 'FAIL '}   {p}")
    if not ok:
        bad.append(p)
if bad:
    print("\nCHECKSUM MISMATCH.")
    print("If a .ds file is listed above, check line endings first: the recorded hashes")
    print("are of the upstream LF bytes; a clone made with core.autocrlf=true yields CRLF.")
    print("This script clones with core.autocrlf=false to avoid exactly that. Otherwise")
    print("the upstream object changed and the analysis used a different file.")
if missing:
    print("\nNot retrieved by this script:")
    for p in missing:
        print(f"  {p}")
    print("1dug_benchmark.pdb is derived, not downloaded; see scripts/25_dataset_groundtruth.py.")
sys.exit(1 if bad else 0)
PY

echo
echo "Done. Now run the freeze check in REPRODUCE.md section L."
