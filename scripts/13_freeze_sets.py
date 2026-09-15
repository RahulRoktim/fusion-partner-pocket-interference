#!/usr/bin/env python
"""
Phase 3A.5 — permanently separate DEVELOPMENT from CONFIRMATORY data.

The 24 pilot structures are permanently DEVELOPMENT/EXPLORATORY observations and may never be
reused as confirmatory observations. The confirmatory candidate pool additionally excludes every
entity sharing a (fusion partner, target accession) pair with the development set, so that a
closely related alternative entry of a development target cannot leak in.

No detector is run here.

Outputs:
  data_manifest/development_set.json / .tsv
  data_manifest/confirmatory_candidate_pool.json / .tsv
  data_manifest/confirmatory_pool_summary.json
"""
import collections, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")


def write_tsv(path, rows, cols):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(
                "" if r.get(c) is None else str(r.get(c)).replace("\t", " ") for c in cols) + "\n")


def main():
    pool = json.load(open(os.path.join(MANI, "pool_eligibility.json")))["records"]
    sel = json.load(open(os.path.join(MANI, "pilot_selection.json")))
    prep = {p["pdb_id"]: p for p in
            json.load(open(os.path.join(MANI, "pilot_manifest.json")))["prepared"]
            if "PREPARATION_FAILED" not in p}

    # ---------------- development set
    dev = []
    for s in sel["selection"]:
        p = prep[s["pdb_id"]]
        dev.append({
            "pdb_id": s["pdb_id"], "entity_id": s["entity_id"],
            "chain": p["auth_chain"], "fusion_partner": s["fusion_partner"],
            "fusion_accession": s["fusion_accession"],
            "target_accession": s["target_accession"], "topology": s["topology"],
            "resolution": s["resolution"], "method": s["method"],
            "description": s["description"],
            "target_residues_modelled": p["residue_class_counts"].get("TARGET"),
            "fusion_residues_modelled": p["residue_class_counts"].get("FUSION"),
            "reference_site_available": p["reference_site_available"],
            "status": "DEVELOPMENT_EXPLORATORY_PERMANENT",
        })
    dev_pdb = {d["pdb_id"] for d in dev}
    dev_pairs = {(d["fusion_partner"], d["target_accession"]) for d in dev}
    dev_targets = {d["target_accession"] for d in dev}

    json.dump({"frozen": "2026-09-15", "n": len(dev),
               "policy": "PERMANENTLY DEVELOPMENT. These observations may never be used as "
                         "confirmatory observations, and no confirmatory threshold, rule or "
                         "endpoint may be re-tuned on them after the confirmatory freeze.",
               "seed_base": sel["seed_base"], "structures": dev},
              open(os.path.join(MANI, "development_set.json"), "w"), indent=1)
    write_tsv(os.path.join(MANI, "development_set.tsv"), dev, list(dev[0].keys()))

    # ---------------- confirmatory candidate pool
    eligible = [r for r in pool if r["eligible"]]
    conf, excl = [], collections.Counter()
    for r in eligible:
        if r["pdb_id"] in dev_pdb:
            excl["same_pdb_entry_as_development"] += 1
            continue
        if (r["fusion_partner"], r["target_accession"]) in dev_pairs:
            excl["same_partner_target_pair_as_development"] += 1
            continue
        conf.append(r)

    # cross-partner leakage check: target used by a DIFFERENT partner in development
    cross = [r for r in conf if r["target_accession"] in dev_targets]
    for r in conf:
        r["target_also_in_development_under_other_partner"] = \
            r["target_accession"] in dev_targets

    print("=" * 96)
    print("DEVELOPMENT / CONFIRMATORY SPLIT")
    print("=" * 96)
    print(f"\nDEVELOPMENT SET: {len(dev)} structures, frozen, permanently exploratory")
    for pt in ("BRIL", "T4L", "MBP"):
        d = [x for x in dev if x["fusion_partner"] == pt]
        print(f"  {pt:5s}: {len(d)} structures / {len({x['target_accession'] for x in d})} targets")

    print(f"\neligible entities in the full pool: {len(eligible)}")
    for k, v in excl.most_common():
        print(f"  removed, {k}: {v}")
    print(f"CONFIRMATORY CANDIDATE POOL: {len(conf)} entities")

    print("\nindependent targets remaining per partner (after development exclusion):")
    summary = {}
    total_targets = set()
    for pt in ("BRIL", "T4L", "MBP"):
        sub = [r for r in conf if r["fusion_partner"] == pt]
        tg = sorted({r["target_accession"] for r in sub})
        total_targets |= {(pt, t) for t in tg}
        # collapsed: one structure per (partner,target), highest resolution then lowest PDB ID
        collapsed = {}
        for r in sorted(sub, key=lambda r: (r["resolution"], r["pdb_id"])):
            collapsed.setdefault(r["target_accession"], r)
        topo = collections.Counter(r["topology"] for r in collapsed.values())
        res = sorted(r["resolution"] for r in collapsed.values())
        summary[pt] = {"n_entities": len(sub), "n_independent_targets": len(tg),
                       "collapsed_units": len(collapsed),
                       "topology": dict(topo),
                       "median_resolution": res[len(res)//2] if res else None,
                       "targets": tg}
        print(f"  {pt:5s}: {len(sub):4d} entities -> {len(tg):3d} independent targets "
              f"(median res {res[len(res)//2] if res else float('nan'):.2f} A)  "
              f"topology {dict(topo)}")
    print(f"  TOTAL: {sum(v['n_independent_targets'] for v in summary.values())} "
          f"independent (partner, target) units")
    print(f"\n  targets appearing in development under a DIFFERENT partner "
          f"(retained, flagged): {len(cross)} entities")

    cols = ["entity_id", "pdb_id", "fusion_partner", "fusion_accession", "target_accession",
            "topology", "method", "resolution", "representative_chain", "entity_length",
            "target_modelled", "fusion_modelled", "target_disorder", "fusion_disorder",
            "align_coverage_of_modelled",
            "target_also_in_development_under_other_partner", "description"]
    write_tsv(os.path.join(MANI, "confirmatory_candidate_pool.tsv"), conf, cols)
    json.dump({"frozen": "2026-09-15",
               "policy": "Confirmatory candidates. Detectors have NOT been run on these. "
                         "Excludes all development PDB entries and all entities sharing a "
                         "(partner, target accession) pair with the development set.",
               "n_entities": len(conf),
               "excluded_counts": dict(excl),
               "per_partner": summary,
               "entities": conf},
              open(os.path.join(MANI, "confirmatory_candidate_pool.json"), "w"), indent=1)
    json.dump(summary, open(os.path.join(MANI, "confirmatory_pool_summary.json"), "w"), indent=1)

    # ---------------- what a relaxed-E9 pool would add (reported, not adopted)
    sole_e9 = [r for r in pool
               if r["exclusions"] and all(e.startswith("E9") for e in r["exclusions"])
               and r["pdb_id"] not in dev_pdb
               and (r["fusion_partner"], r["target_accession"]) not in dev_pairs]
    print(f"\n  for reference: removing E9 entirely would add {len(sole_e9)} entities to the "
          f"confirmatory pool")
    add = collections.Counter()
    for r in sole_e9:
        add[r["fusion_partner"]] += 1
    newt = {pt: len({r["target_accession"] for r in sole_e9 if r["fusion_partner"] == pt}
                    - set(summary[pt]["targets"])) for pt in ("BRIL", "T4L", "MBP")}
    print(f"    by partner: {dict(add)};  NEW independent targets added: {newt}")


if __name__ == "__main__":
    main()
