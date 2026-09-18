#!/usr/bin/env python
"""
Equivalence test for the chain-aware variant classifier (scripts/21).

Asserts that on SINGLE-CHAIN input the chain-aware parsers recover exactly the same pockets,
ranks, scores, centroids and residue sets as the frozen single-chain parsers, once the chain
prefix is stripped. This is what licenses using the chain-aware path for the multi-chain
assembly1 variant while leaving the frozen path untouched for everything else.

Runs against the already-archived PRIMARY confirmatory detector output — no detector is re-run,
and no primary result is modified.

Run:  python tests/test_variant_classifier_equivalence.py
"""
import csv, hashlib, importlib.util, json, math, os, sys, tempfile
from pathlib import Path

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRIM = os.path.join(ROOT, "results", "confirmatory", "primary")
RAW_HASHES = os.path.join(ROOT, "results", "RAW_OUTPUT_HASHES.tsv")

EXIT_FAILED = 1
EXIT_UNAVAILABLE = 2

FAILURES, CHECKS = [], [0]


def check(cond, label, detail=""):
    CHECKS[0] += 1
    if not cond:
        FAILURES.append(f"{label} {detail}")
        print(f"  FAIL  {label}  {detail}")


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(ROOT, "scripts", f))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


frozen = _load("frozen", "05_classify_pockets.py")
variant = _load("variant", "21_classify_assembly1.py")


def strip_chain(keys):
    return [k.split(":", 1)[1] if ":" in k else k for k in keys]


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_tree(path):
    digest = hashlib.sha256()
    for current, dirs, files in os.walk(path):
        dirs.sort()
        for filename in sorted(files):
            file_path = os.path.join(current, filename)
            digest.update(os.path.relpath(file_path, path).replace("\\", "/").encode())
            digest.update(sha256_file(file_path).encode())
    return digest.hexdigest()


def load_declared_hashes(path=RAW_HASHES):
    with open(path, newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {(r["TREE"], r["DETECTOR"], r["RUN_KEY"]): r["SHA256"] for r in rows}


def preflight_archived_outputs(runs, root=ROOT, declared_hashes=None):
    """Require every successful archived detector tree and verify its frozen digest."""
    declared_hashes = load_declared_hashes() if declared_hashes is None else declared_hashes
    missing, mismatched, undeclared = [], [], []
    successful = [r for r in runs if r.get("success")]
    if not successful:
        return "FAILED", ["detector_runs.json contains no successful runs"]

    for run in successful:
        relative = run.get("archived_output", "").replace("\\", "/").strip("/")
        parts = relative.split("/") if relative else []
        if len(parts) < 3:
            undeclared.append(f"{run.get('detector')}:{run.get('key')} has no valid archived_output")
            continue
        tree = "/".join(parts[:-2])
        key = (tree, run["detector"], run["key"])
        recorded = run.get("output_sha256")
        declared = declared_hashes.get(key)
        if not recorded or declared != recorded:
            undeclared.append(
                f"{run['detector']}:{run['key']} run={recorded or '<missing>'} "
                f"manifest={declared or '<missing>'}"
            )
            continue
        output_dir = Path(root, *parts)
        if not output_dir.is_dir():
            missing.append(relative)
            continue
        actual = sha256_tree(output_dir)
        if actual != recorded:
            mismatched.append(f"{relative}: expected {recorded}, got {actual}")

    if undeclared or mismatched:
        return "FAILED", undeclared + mismatched
    if missing:
        return "UNAVAILABLE", missing
    return "READY", []


def equivalence_exit_code(n_runs, n_pockets, failures):
    if failures or n_runs <= 0:
        return EXIT_FAILED
    if n_pockets <= 0:
        return EXIT_UNAVAILABLE
    return 0


def self_test():
    """Regression guards for the former empty-vs-empty false PASS."""
    empty_hash = hashlib.sha256(b"").hexdigest()
    run = {
        "success": True,
        "detector": "fpocket",
        "key": "EMPTY_CASE",
        "archived_output": "raw/fpocket/EMPTY_CASE",
        "output_sha256": empty_hash,
    }
    declared = {("raw", "fpocket", "EMPTY_CASE"): empty_hash}
    with tempfile.TemporaryDirectory(prefix="fusion-equivalence-guard-") as temp_root:
        status, _ = preflight_archived_outputs([run], temp_root, declared)
        assert status == "UNAVAILABLE", status
        Path(temp_root, "raw", "fpocket", "EMPTY_CASE").mkdir(parents=True)
        status, detail = preflight_archived_outputs([run], temp_root, declared)
        assert status == "READY", detail
        assert equivalence_exit_code(1, 0, []) == EXIT_UNAVAILABLE
        assert equivalence_exit_code(1, 1, []) == 0
        Path(temp_root, "raw", "fpocket", "EMPTY_CASE", "unexpected.txt").write_text(
            "changed", encoding="utf-8"
        )
        status, _ = preflight_archived_outputs([run], temp_root, declared)
        assert status == "FAILED", status
    print("SELF-TEST PASS: missing, empty and hash-mismatched inputs cannot produce PASS")
    return 0


def main():
    print("=" * 78)
    print("VARIANT CLASSIFIER EQUIVALENCE (single-chain inputs)")
    print("=" * 78)
    run_manifest = json.load(open(os.path.join(PRIM, "detector_runs.json")))
    runs = run_manifest["runs"]
    expected_runs = run_manifest.get("expected_runs")
    actual_runs = run_manifest.get("actual_runs")
    if expected_runs != len(runs) or actual_runs != len(runs):
        print(
            "FAILED: detector run manifest count mismatch: "
            f"expected_runs={expected_runs}, actual_runs={actual_runs}, rows={len(runs)}"
        )
        return EXIT_FAILED
    ok = [r for r in runs if r["success"]]
    print(f"  comparing parsers on {len(ok)} archived single-chain runs\n")

    availability, detail = preflight_archived_outputs(ok)
    if availability == "UNAVAILABLE":
        print(
            "UNAVAILABLE: required archived raw detector outputs are absent "
            f"({len(detail)} of {len(ok)} successful runs)."
        )
        print("Restore or regenerate the hash-pinned raw_detector trees before running equivalence.")
        return EXIT_UNAVAILABLE
    if availability == "FAILED":
        print(f"FAILED: archived detector output integrity failed ({len(detail)} problem(s)).")
        for problem in detail[:10]:
            print("  -", problem)
        return EXIT_FAILED

    n_pockets = 0
    for r in ok:
        od = os.path.join(ROOT, *r["archived_output"].replace("\\", "/").split("/"))
        base = r["key"]
        if r["detector"] == "p2rank":
            a = frozen.parse_p2rank(od, base + ".pdb")
            b = variant.parse_p2rank_chainaware(od, base + ".pdb")
        else:
            a = frozen.parse_fpocket(od, base)
            b = variant.parse_fpocket_chainaware(od, base)
        check(len(a) == len(b), f"E1 pocket count {r['detector']} {base}",
              f"{len(a)} vs {len(b)}")
        for pa, pb in zip(a, b):
            n_pockets += 1
            check(pa["rank"] == pb["rank"], f"E2 rank {base}")
            check((math.isnan(pa["score"]) and math.isnan(pb["score"]))
                  or abs(pa["score"] - pb["score"]) < 1e-9, f"E3 score {base} rank{pa['rank']}")
            for i in range(3):
                ca, cb = pa["centroid"][i], pb["centroid"][i]
                check((math.isnan(ca) and math.isnan(cb)) or abs(ca - cb) < 1e-6,
                      f"E4 centroid {base} rank{pa['rank']}")
            check(list(pa["residues"]) == strip_chain(pb["residues"]),
                  f"E5 residue set {r['detector']} {base} rank{pa['rank']}",
                  f"{pa['residues'][:4]} vs {strip_chain(pb['residues'])[:4]}")

    print(f"  {CHECKS[0]} assertions over {n_pockets} pockets")
    print("\n" + "=" * 78)
    result = equivalence_exit_code(len(ok), n_pockets, FAILURES)
    if result == EXIT_FAILED:
        print(f"RESULT: {len(FAILURES)} FAILED of {CHECKS[0]}")
        for f in FAILURES[:10]:
            print("  -", f)
        return EXIT_FAILED
    if result == EXIT_UNAVAILABLE:
        print("RESULT: UNAVAILABLE — zero pockets were compared; zero observations cannot PASS")
        return EXIT_UNAVAILABLE
    print(f"RESULT: all {CHECKS[0]} assertions PASSED — chain-aware parsers are equivalent "
          f"to the frozen parsers on single-chain input")
    return 0


if __name__ == "__main__":
    raise SystemExit(self_test() if "--self-test" in sys.argv[1:] else main())
