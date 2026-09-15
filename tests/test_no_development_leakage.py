#!/usr/bin/env python
"""
Amendment F1 leakage test. FAILS if any target UniProt accession, PDB entry, or polymer entity
occurs in both the development set and the final confirmatory set.

This test must pass before any confirmatory detector is run, and again before the confirmatory
results are reported.

Run:  python tests/test_no_development_leakage.py
"""
import json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")

FAILURES, CHECKS = [], [0]


def check(cond, label, detail=""):
    CHECKS[0] += 1
    print(("  PASS  " if cond else "  FAIL  ") + label + (f"  {detail}" if detail and not cond else ""))
    if not cond:
        FAILURES.append(f"{label} {detail}")


def main():
    print("=" * 78)
    print("F1 DEVELOPMENT / CONFIRMATORY LEAKAGE TEST")
    print("=" * 78)

    dev = json.load(open(os.path.join(MANI, "development_set.json")))["structures"]
    conf = json.load(open(os.path.join(MANI, "confirmatory_set_final.json")))["structures"]

    dev_t = {d["target_accession"] for d in dev}
    dev_p = {d["pdb_id"] for d in dev}
    dev_e = {d["entity_id"] for d in dev}
    con_t = {c["target_accession"] for c in conf}
    con_p = {c["pdb_id"] for c in conf}
    con_e = {c["entity_id"] for c in conf}

    print(f"\n  development: {len(dev)} structures, {len(dev_t)} targets")
    print(f"  confirmatory: {len(conf)} structures, {len(con_t)} targets\n")

    check(not (dev_t & con_t), "L1 zero TARGET ACCESSION overlap",
          f"shared: {sorted(dev_t & con_t)}")
    check(not (dev_p & con_p), "L2 zero PDB ENTRY overlap", f"shared: {sorted(dev_p & con_p)}")
    check(not (dev_e & con_e), "L3 zero POLYMER ENTITY overlap",
          f"shared: {sorted(dev_e & con_e)}")

    # a (partner, target) pair must not recur either
    dev_pair = {(d["fusion_partner"], d["target_accession"]) for d in dev}
    con_pair = {(c["fusion_partner"], c["target_accession"]) for c in conf}
    check(not (dev_pair & con_pair), "L4 zero (partner, target) pair overlap",
          f"shared: {sorted(dev_pair & con_pair)}")

    # the confirmatory set must itself hold one structure per (partner, target)
    check(len(con_pair) == len(conf),
          "L5 confirmatory set holds exactly one structure per (partner, target)",
          f"{len(conf)} structures vs {len(con_pair)} pairs")

    # Targets are the CLUSTERING unit. A target may legitimately appear twice under two
    # different fusion partners (a within-target partner contrast). What must never happen is an
    # accidental duplicate under the SAME partner. [F4]
    import collections
    dupes = {t: [c for c in conf if c["target_accession"] == t]
             for t, n in collections.Counter(c["target_accession"] for c in conf).items() if n > 1}
    bad = {t: [c["fusion_partner"] for c in rows] for t, rows in dupes.items()
           if len({c["fusion_partner"] for c in rows}) != len(rows)}
    check(not bad, "L6 every duplicated target is a legitimate CROSS-PARTNER pair", f"{bad}")
    print(f"        {len(conf)} structures over {len(con_t)} independent target clusters; "
          f"{len(dupes)} targets appear under two partners: "
          f"{ {t: [c['fusion_partner'] for c in r] for t, r in dupes.items()} }")
    check(len(con_t) >= 0.9 * len(conf),
          "L6b independent clusters are at least 90% of structures",
          f"{len(con_t)}/{len(conf)}")

    check(len(dev) == 24, "L7 development set is still exactly 24 structures", f"{len(dev)}")
    check(all(d.get("status") == "DEVELOPMENT_EXPLORATORY_PERMANENT" for d in dev),
          "L8 every development record is flagged permanently exploratory")

    print("\n" + "=" * 78)
    if FAILURES:
        print(f"RESULT: {len(FAILURES)} FAILED of {CHECKS[0]} checks — LEAKAGE PRESENT")
        for f in FAILURES:
            print("  -", f)
        sys.exit(1)
    print(f"RESULT: all {CHECKS[0]} checks PASSED — sets are disjoint")


if __name__ == "__main__":
    main()
