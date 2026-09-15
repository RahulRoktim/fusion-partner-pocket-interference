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
import importlib.util, json, math, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRIM = os.path.join(ROOT, "results", "confirmatory", "primary")

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


def main():
    print("=" * 78)
    print("VARIANT CLASSIFIER EQUIVALENCE (single-chain inputs)")
    print("=" * 78)
    runs = json.load(open(os.path.join(PRIM, "detector_runs.json")))["runs"]
    ok = [r for r in runs if r["success"]]
    print(f"  comparing parsers on {len(ok)} archived single-chain runs\n")

    n_pockets = 0
    for r in ok:
        od = os.path.join(ROOT, r["archived_output"])
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
    if FAILURES:
        print(f"RESULT: {len(FAILURES)} FAILED of {CHECKS[0]}")
        for f in FAILURES[:10]:
            print("  -", f)
        sys.exit(1)
    print(f"RESULT: all {CHECKS[0]} assertions PASSED — chain-aware parsers are equivalent "
          f"to the frozen parsers on single-chain input")


if __name__ == "__main__":
    main()
