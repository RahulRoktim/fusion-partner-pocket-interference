#!/usr/bin/env python
"""Verify this package against the frozen manifest environment/FINAL_FREEZE_v1.1.json.

The manifest hashes the full working repository. The public package differs from it in two
documented, deliberate ways, and this script reports each separately instead of failing opaquely:

  * 10 third-party files are referenced rather than redistributed. Run
    `bash scripts_release/fetch_third_party.sh` to retrieve them, after which they verify.
  * 2 root documents (README.md, REPRODUCE.md) were replaced with public-facing versions.
    The originals are preserved verbatim at manuscript-support/README_v1.1_historical.md and
    manuscript-support/REPRODUCE_v1.1_historical.md, and are verified here in their stead.

Exit status is non-zero only for a genuine mismatch.
"""
import csv
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
SUBSTITUTED = {"README.md": "manuscript-support/README_v1.1_historical.md",
               "REPRODUCE.md": "manuscript-support/REPRODUCE_v1.1_historical.md"}


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    man = json.load(open("environment/FINAL_FREEZE_v1.1.json", encoding="utf-8"))
    excluded = {r["REPO_PATH"] for r in csv.DictReader(
        open("EXCLUDED_THIRD_PARTY.tsv", encoding="utf-8"), delimiter="\t")}

    ok = mismatch = []
    ok, mismatch, pending, substituted = 0, [], [], []
    for group in man["sha256"].values():
        for path, h in group.items():
            target = SUBSTITUTED.get(path, path)
            if not os.path.exists(target):
                (pending if path in excluded else mismatch).append(path)
                continue
            if sha(target) == h["sha256"]:
                ok += 1
                if path in SUBSTITUTED:
                    substituted.append(path)
            else:
                mismatch.append(path)

    total = sum(len(g) for g in man["sha256"].values())
    print(f"frozen manifest: {total} files")
    print(f"  verified                     {ok}")
    print(f"  substituted but verified     {len(substituted)}  {substituted}")
    print(f"  pending third-party retrieval {len(pending)}")
    for p in sorted(pending):
        print(f"      {p}")
    if pending:
        print("      -> run: bash scripts_release/fetch_third_party.sh")
    print(f"  GENUINE MISMATCHES           {len(mismatch)}")
    for p in sorted(mismatch):
        print(f"      {p}")
    print("\nEach text file also carries sha256_lf, a line-ending-normalised hash that verifies")
    print("identically on Linux, macOS and Windows; prefer it for cross-platform checks.")
    return 1 if mismatch else 0


if __name__ == "__main__":
    sys.exit(main())
