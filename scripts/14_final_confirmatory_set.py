#!/usr/bin/env python
"""
Amendment F1 — build the TARGET-DISJOINT final confirmatory set.

Removes any confirmatory candidate whose TARGET UniProt accession occurs anywhere in the
24-structure development set, even under a different fusion partner. Requirement: ZERO
target-accession overlap between development and confirmation.

One structure per (partner, target): highest resolution, ties by lowest PDB ID. No sampling.

Outputs:
  data_manifest/confirmatory_set_final.json / .tsv
  data_manifest/confirmatory_set_final_excluded.tsv
"""
import collections, json, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
PARTNERS = ["BRIL", "T4L", "MBP"]


def write_tsv(path, rows, cols):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else
                               str(r.get(c)).replace("\t", " ") for c in cols) + "\n")


def main():
    dev = json.load(open(os.path.join(MANI, "development_set.json")))["structures"]
    pool = json.load(open(os.path.join(MANI, "confirmatory_candidate_pool.json")))["entities"]

    dev_targets = {d["target_accession"] for d in dev}
    dev_pdb = {d["pdb_id"] for d in dev}
    print(f"development targets: {len(dev_targets)}  development entries: {len(dev_pdb)}")

    kept, dropped = [], []
    for r in pool:
        if r["target_accession"] in dev_targets:
            r = dict(r)
            r["drop_reason"] = "F1_target_accession_in_development_set"
            dropped.append(r)
        elif r["pdb_id"] in dev_pdb:
            r = dict(r)
            r["drop_reason"] = "F1_pdb_entry_in_development_set"
            dropped.append(r)
        else:
            kept.append(r)
    print(f"candidate entities {len(pool)} -> {len(kept)} after F1 "
          f"({len(dropped)} removed)")
    for k, v in collections.Counter(d["drop_reason"] for d in dropped).items():
        print(f"  {k}: {v}")

    # one structure per (partner, target): highest resolution, ties by lowest PDB ID
    collapsed = {}
    for r in sorted(kept, key=lambda r: (r["resolution"], r["pdb_id"])):
        key = (r["fusion_partner"], r["target_accession"])
        if key not in collapsed:
            r = dict(r)
            r["n_alternative_entries"] = 0
            collapsed[key] = r
        else:
            collapsed[key]["n_alternative_entries"] += 1

    final = [collapsed[k] for k in sorted(collapsed)]
    for f in final:
        f["status"] = "CONFIRMATORY_FINAL"

    # ---- hard leakage assertion
    conf_targets = {f["target_accession"] for f in final}
    overlap = conf_targets & dev_targets
    assert not overlap, f"LEAKAGE: {overlap}"
    assert not ({f['pdb_id'] for f in final} & dev_pdb), "LEAKAGE: shared PDB entry"

    print(f"\nFINAL CONFIRMATORY SET: {len(final)} structures")
    for pt in PARTNERS:
        sub = [f for f in final if f["fusion_partner"] == pt]
        res = sorted(f["resolution"] for f in sub)
        topo = collections.Counter(f["topology"] for f in sub)
        meth = collections.Counter(f["method"] for f in sub)
        ntg = len({f["target_accession"] for f in sub})
        print(f"  {pt:5s}: {len(sub):3d} structures / {ntg} targets"
              f"  median res {res[len(res)//2]:.2f} A  topology {dict(topo)}")
        print(f"         method {dict(meth)}")
    print(f"  target-accession overlap with development: {len(overlap)}  (required 0)")

    cols = ["entity_id", "pdb_id", "fusion_partner", "fusion_accession", "target_accession",
            "topology", "method", "resolution", "representative_chain", "entity_length",
            "target_modelled", "fusion_modelled", "target_disorder", "fusion_disorder",
            "align_coverage_of_modelled", "n_alternative_entries", "status", "description"]
    write_tsv(os.path.join(MANI, "confirmatory_set_final.tsv"), final, cols)
    write_tsv(os.path.join(MANI, "confirmatory_set_final_excluded.tsv"), dropped,
              cols[:-2] + ["drop_reason"])
    json.dump({"frozen": "2026-09-15", "amendment": "F1",
               "policy": "Target-disjoint from the development set. One structure per "
                         "(partner, target): highest resolution, ties by lowest PDB ID. "
                         "No sampling, no seed.",
               "n": len(final),
               "per_partner": {pt: sum(1 for f in final if f["fusion_partner"] == pt)
                               for pt in PARTNERS},
               "development_target_overlap": 0,
               "structures": final},
              open(os.path.join(MANI, "confirmatory_set_final.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
