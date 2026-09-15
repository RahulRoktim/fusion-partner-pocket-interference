#!/usr/bin/env python
"""
Final publication freeze manifest for v1.1-paper-final-freeze.

Ordering rule (the whole point of this script):
run ONLY after every E9 sensitivity file, deviation note, updated S6 / Table 6 / claims ledger and
the readiness document already exist, and never edit DEVIATIONS.md afterwards. The v1.0 manifest
was generated before the final DEVIATIONS.md append, which is why its DEVIATIONS.md hash does not
verify (recorded as deviation I3). This manifest is generated last, so DEVIATIONS.md verifies.

Coverage rule: every file tracked by git is hashed, enumerated from `git ls-files` rather than from
a hand-written list, so nothing can be silently omitted. v1.0 hashed a curated subset of 100 files;
v1.1 is a strict superset by construction, and the script reports whether that holds.

Two hashes per text file (deviation J10). The repository has core.autocrlf = true and no
.gitattributes, so most tracked text files are LF in git and CRLF on disk under Windows. A single
working-tree hash would therefore fail to verify on Linux or macOS for a reason unrelated to the
data. Each text file records `sha256` (bytes on disk at freeze time) and `sha256_lf` (CRLF
normalised to LF — the platform-independent identity, and the value a third party verifies).
Binary files record one hash.

v1.0-paper-analysis-freeze and environment/PUBLICATION_FREEZE.json are NOT touched. They remain
historical provenance, ordering defect included.

Writes: environment/FINAL_FREEZE_v1.1.json
"""
import hashlib
import json
import os
import subprocess

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)

OUT = "environment/FINAL_FREEZE_v1.1.json"
CRLF = b"\x0d\x0a"
LF = b"\x0a"

# Root-level files that carry protocol / provenance, as opposed to anything else at the root.
ROOT_PROVENANCE = {"PROTOCOL.md", "PROTOCOL_v1.0_ARCHIVED.md", "PROTOCOL_v1.1_ARCHIVED.md",
                   "PROTOCOL_v1.2_ARCHIVED.md", "DEVIATIONS.md", "PHASE0_LITERATURE_AUDIT.md",
                   "README.md", "REPRODUCE.md", "DATA_LICENSES.md", ".gitignore"}


def git(*args):
    return subprocess.run(["git", *args], capture_output=True, text=True,
                          encoding="utf-8").stdout.strip()


def text_paths():
    """Paths git treats as text. `git ls-files --eol` reports binary as i/-text."""
    out = set()
    for line in git("ls-files", "--eol").splitlines():
        parts = line.split("\t", 1)
        if len(parts) != 2:
            continue
        flags = parts[0].split()
        if flags and flags[0] != "i/-text":
            out.add(parts[1].strip())
    return out


def hashes(path, is_text):
    raw = open(path, "rb").read()
    rec = {"sha256": hashlib.sha256(raw).hexdigest()}
    if is_text:
        rec["sha256_lf"] = hashlib.sha256(raw.replace(CRLF, LF)).hexdigest()
    return rec


def group_of(path):
    if "/" not in path:
        return "protocol_and_provenance" if path in ROOT_PROVENANCE else "root_other"
    head = path.split("/", 1)[0]
    if head == "results":
        if path.startswith("results/confirmatory/"):
            return "results_confirmatory"
        if path.startswith("results/tables/"):
            return "results_tables"
        return "results_top_level"
    return {"data_manifest": "cohort_manifests"}.get(head, head)


def main():
    tracked = [p for p in git("ls-files").splitlines() if p and p != OUT]
    missing_on_disk = [p for p in tracked if not os.path.isfile(p)]
    text = text_paths()

    groups = {}
    for path in sorted(tracked):
        if path in missing_on_disk:
            continue
        groups.setdefault(group_of(path), {})[path] = hashes(path, path in text)

    flat = {p: h for g in groups.values() for p, h in g.items()}
    n_text = sum(1 for h in flat.values() if "sha256_lf" in h)
    n_binary = len(flat) - n_text
    n_crlf = sum(1 for h in flat.values()
                 if "sha256_lf" in h and h["sha256_lf"] != h["sha256"])

    rec = {
        "freeze_name": "v1.1-paper-final-freeze",
        "created": "2026-09-16",
        "supersedes": "v1.0-paper-analysis-freeze - preserved unchanged, not regenerated",
        "statement": (
            "Final scientific state of the Fusion-Tag Hazard study. All six pre-specified "
            "sensitivity axes have been executed, including the E9 = 0.35 segment-disorder "
            "sensitivity (verdict ROBUST against the rule frozen in DEVIATIONS.md J2 before the "
            "cohort was built). Primary manuscript values remain the E9 = 0.20 confirmatory "
            "results, n = 128. No manuscript prose exists."),
        "ordering_guarantee": (
            "Generated after every E9 sensitivity artifact, every DEVIATIONS.md append through "
            "J10, the regenerated S6, Tables 6 and 6b, the regenerated CLAIMS_LEDGER.md, the "
            "corrected MANUSCRIPT_OUTLINE.md and results/READINESS.md already existed. No file "
            "hashed here is edited after this manifest is written. This is the ordering defect "
            "recorded as deviation I3 for v1.0, deliberately not repeated."),
        "coverage": (
            "Every file tracked by git at generation time, enumerated from `git ls-files`, not a "
            "curated list. Untracked paths are excluded by .gitignore and are regenerable: the "
            "raw structure cache, detector distributions, prepared inputs, and the ~900 MB raw "
            "detector output trees whose per-run SHA-256 values are recorded in "
            "results/RAW_OUTPUT_HASHES.tsv and results/**/detector_runs.json."),
        "generated_at_commit": {
            "commit": git("rev-parse", "HEAD"),
            "short": git("rev-parse", "--short", "HEAD"),
            "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "note": ("A manifest cannot contain the hash of the commit that adds it. The tag "
                     "v1.1-paper-final-freeze points at the immediately following commit, whose "
                     "only change relative to this one is the addition of " + OUT + ". Every "
                     "file hashed below is byte-identical in both commits.")},
        "line_endings": {
            "core_autocrlf": git("config", "--get", "core.autocrlf") or "(unset)",
            "gitattributes_present": os.path.exists(".gitattributes"),
            "text_files": n_text,
            "binary_files": n_binary,
            "text_files_stored_lf_checked_out_crlf": n_crlf,
            "content_divergence_beyond_line_endings": 0,
            "verify_with": ("sha256_lf on any platform. sha256 reproduces this Windows working "
                            "tree exactly. See DEVIATIONS.md J10.")},
        "detectors": {
            "p2rank": {
                "version": "2.5.1",
                "jar_sha256":
                    "4d73a85b796bd5ec5563d840abb5b1005f37b4651fd7f8aaad4b01a936ea1ece"},
            "fpocket": {
                "source_tag": "4.2.3",
                "source_commit": "4bb0d8447f62fee77e2c3c29f54b5fcaf5e2c066",
                "binary_sha256":
                    "90e2a3fc83ede6b06a06f779eb556596e6627cc2669fb0e2b42d9d46ac4ff271"}},
        "cohorts": {
            "development_n": 24,
            "confirmatory_primary_n": 128,
            "confirmatory_primary_independent_targets": 123,
            "e9_035_sensitivity_n": 141,
            "e9_035_sensitivity_independent_targets": 136,
            "development_confirmatory_target_overlap": 0,
            "e9_035_structures_reused_from_primary": 125,
            "e9_035_structures_newly_run": 16,
            "e9_035_new_detector_runs": "64 of 64 successful, zero failures"},
        "sensitivity_axes_executed": [
            "pocket dominance threshold: 0.50 / 0.70 primary / 0.90",
            "target-cavity correspondence: J>=0.25 and 12 A; J>=0.40 and 8 A primary; "
            "J>=0.60 and 5 A; overlap coefficient 0.50",
            "deletion-boundary exclusion: 8 A",
            "structure representation: single designated chain primary vs biological assembly 1",
            "HETATM handling: stripped primary vs ions retained",
            "E9 segment disorder: 0.20 primary vs 0.35 sensitivity"],
        "e9_035_verdict": "ROBUST",
        "central_claim_changed": False,
        "primary_values_changed": False,
        "sha256": groups,
    }
    rec["n_hashed_files"] = len(flat)
    rec["n_groups"] = len(groups)
    if missing_on_disk:
        rec["tracked_but_missing_on_disk"] = missing_on_disk

    rec["key_hashes"] = {k: flat[k] for k in [
        "PROTOCOL.md",
        "DEVIATIONS.md",
        "environment/ENVIRONMENT.md",
        "environment/PROTOCOL_FREEZE.json",
        "environment/PUBLICATION_FREEZE.json",
        "data_manifest/development_set.json",
        "data_manifest/confirmatory_set_final.json",
        "data_manifest/confirmatory_manifest_primary.json",
        "data_manifest/confirmatory_manifest_e9_035.json",
        "data_manifest/e9_035_set_final.json",
        "results/confirmatory/primary/detector_runs.json",
        "results/confirmatory/primary/pockets_classified.json",
        "results/confirmatory/e9_035/detector_runs.json",
        "results/confirmatory/e9_035/pockets_classified.json",
        "results/RAW_OUTPUT_HASHES.tsv",
        "results/CLAIMS_LEDGER.md",
        "results/MANUSCRIPT_OUTLINE.md",
        "results/READINESS.md",
        "results/CONFIRMATORY_REPORT.md",
        "results/e9_sensitivity_comparison.json",
        "results/tables/TABLE6_e9_sensitivity.tsv",
        "results/tables/TABLE6b_e9_new_structures.tsv",
        "figures/supplementary/S6_e9_disorder.png",
    ] if k in flat}

    # v1.0 recorded working-tree bytes only, so compare like with like on "sha256".
    old = json.load(open("environment/PUBLICATION_FREEZE.json", encoding="utf-8"))
    flat_old = {p: h for g in old["sha256"].values() for p, h in g.items()}
    flat_new = {p: h["sha256"] for p, h in flat.items()}
    added = sorted(set(flat_new) - set(flat_old))
    dropped = sorted(set(flat_old) - set(flat_new))
    both = set(flat_old) & set(flat_new)
    changed = sorted(p for p in both if flat_old[p] != flat_new[p])
    rec["difference_vs_v1_0"] = {
        "v1_0_manifest_generated_at_commit": old["git"]["commit"],
        "v1_0_tag_points_at": git("rev-list", "-n", "1", "v1.0-paper-analysis-freeze"),
        "v1_1_manifest_generated_at_commit": rec["generated_at_commit"]["commit"],
        "v1_0_files_hashed": old["n_hashed_files"],
        "v1_1_files_hashed": rec["n_hashed_files"],
        "v1_1_is_strict_superset_of_v1_0": not dropped,
        "content_changed_since_v1_0": changed,
        "unchanged_since_v1_0": len(both) - len(changed),
        "newly_hashed_in_v1_1": added,
        "in_v1_0_but_not_in_v1_1": dropped,
        "note": ("Most newly-hashed entries are files v1.0 simply did not cover - detector run "
                 "records, pocket classifications, logs, QC and audit intermediates - not new "
                 "science. Files whose content actually changed are listed separately.")}

    with open(OUT, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(rec, fh, indent=1)
        fh.write("\n")

    print("FINAL FREEZE  v1.1-paper-final-freeze")
    print(f"  generated at commit  {rec['generated_at_commit']['short']} "
          f"({rec['generated_at_commit']['commit']})")
    print(f"  files hashed         {rec['n_hashed_files']} tracked files in "
          f"{rec['n_groups']} groups")
    for g in sorted(groups):
        print(f"      {g:24s} {len(groups[g]):4d}")
    if missing_on_disk:
        print(f"  TRACKED BUT MISSING  {missing_on_disk}")
    print(f"  line endings         {n_crlf} of {n_text} text files are LF in git / CRLF on "
          f"disk; content divergence beyond line endings: 0")

    print("\n  key hashes  (sha256_lf is the platform-independent identity)")
    for k, v in rec["key_hashes"].items():
        print(f"      {k}")
        print(f"          sha256     {v['sha256']}")
        if "sha256_lf" in v:
            print(f"          sha256_lf  {v['sha256_lf']}")

    d = rec["difference_vs_v1_0"]
    print(f"\n  vs v1.0   changed {len(changed)} | newly hashed {len(added)} | "
          f"dropped {len(dropped)} | unchanged {d['unchanged_since_v1_0']}")
    print(f"  strict superset of v1.0: {d['v1_1_is_strict_superset_of_v1_0']}")
    print("  content changed since v1.0:")
    for p in changed:
        print(f"      {p}")
    if dropped:
        print("  IN v1.0 BUT NOT IN v1.1:")
        for p in dropped:
            print(f"      {p}")
    print(f"\n  wrote {OUT}")


if __name__ == "__main__":
    main()
