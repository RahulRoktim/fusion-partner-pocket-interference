#!/usr/bin/env python
"""
Amendment J — build the E9 = 0.35 sensitivity cohort.

Re-runs the FROZEN eligibility logic from scripts/01_build_pool.py against the FROZEN cached pool
(data_manifest/pool_raw.json), changing exactly one constant: MAX_SEGMENT_DISORDER 0.20 -> 0.35.
Every other rule is imported, not reimplemented, so the two cohorts cannot drift by accident.

Then applies the same frozen chain as the primary:
  eligible -> drop development PDB entries -> drop (partner, target) pairs seen in development
           -> [F1] drop any target UniProt accession seen in development
           -> one structure per (partner, target): highest resolution, tie-break lowest PDB ID

Outputs:
  data_manifest/e9_035_pool_eligibility.json
  data_manifest/e9_035_set_final.json / .tsv
  data_manifest/e9_035_cohort_diff.json
"""
import collections, importlib.util, json, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
PARTNERS = ["BRIL", "T4L", "MBP"]
SENS_DISORDER = 0.35


def _load(name, fn):
    s = importlib.util.spec_from_file_location(name, os.path.join(ROOT, "scripts", fn))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


build = _load("build_pool", "01_build_pool.py")


def main():
    print(f"frozen primary threshold : {build.MAX_SEGMENT_DISORDER}")
    build.MAX_SEGMENT_DISORDER = SENS_DISORDER          # the ONLY change
    print(f"sensitivity threshold    : {build.MAX_SEGMENT_DISORDER}")
    for k in ("MAX_RESOLUTION", "MIN_FUSION_MODELLED", "MIN_TARGET_MODELLED",
              "MIN_ALIGN_COVERAGE"):
        print(f"  unchanged {k:22s} = {getattr(build, k)}")

    raw = json.load(open(os.path.join(MANI, "pool_raw.json")))
    print(f"\nfrozen cached pool: {len(raw)} entity records")
    evaluated = [build.evaluate(r) for r in raw]
    eligible = [e for e in evaluated if e["eligible"]]
    json.dump({"threshold": SENS_DISORDER, "n_records": len(evaluated),
               "records": evaluated},
              open(os.path.join(MANI, "e9_035_pool_eligibility.json"), "w"))
    print(f"eligible at 0.35: {len(eligible)}")

    # sanity: primary eligibility recomputed identically at 0.20?
    prim = json.load(open(os.path.join(MANI, "pool_eligibility.json")))["records"]
    prim_elig = {r["entity_id"] for r in prim if r["eligible"]}
    sens_elig = {e["entity_id"] for e in eligible}
    print(f"eligible at 0.20 (frozen file): {len(prim_elig)}")
    lost = prim_elig - sens_elig
    print(f"entities eligible at 0.20 but NOT at 0.35: {len(lost)}"
          + (f"  {sorted(lost)[:10]}" if lost else "  (as expected: relaxing adds only)"))

    # ---- same frozen chain as the primary
    dev = json.load(open(os.path.join(MANI, "development_set.json")))["structures"]
    dev_pdb = {d["pdb_id"] for d in dev}
    dev_pairs = {(d["fusion_partner"], d["target_accession"]) for d in dev}
    dev_targets = {d["target_accession"] for d in dev}

    kept, dropped = [], collections.Counter()
    for r in eligible:
        if r["pdb_id"] in dev_pdb:
            dropped["pdb_entry_in_development"] += 1
        elif (r["fusion_partner"], r["target_accession"]) in dev_pairs:
            dropped["partner_target_pair_in_development"] += 1
        elif r["target_accession"] in dev_targets:
            dropped["F1_target_accession_in_development"] += 1
        else:
            kept.append(r)
    print(f"\nafter development exclusions: {len(kept)} entities  (dropped {dict(dropped)})")

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
        f["status"] = "E9_035_SENSITIVITY"

    conf_targets = {f["target_accession"] for f in final}
    assert not (conf_targets & dev_targets), "LEAKAGE: development target in sensitivity cohort"
    assert not ({f["pdb_id"] for f in final} & dev_pdb), "LEAKAGE: development entry"

    # ---- diff against the frozen primary confirmatory set
    prim_set = json.load(open(os.path.join(MANI,
                                           "confirmatory_set_final.json")))["structures"]
    prim_by_pair = {(p["fusion_partner"], p["target_accession"]): p for p in prim_set}
    sens_by_pair = {(f["fusion_partner"], f["target_accession"]): f for f in final}

    added_pairs = sorted(set(sens_by_pair) - set(prim_by_pair))
    lost_pairs = sorted(set(prim_by_pair) - set(sens_by_pair))
    shared = sorted(set(prim_by_pair) & set(sens_by_pair))
    swapped = [(k, prim_by_pair[k]["pdb_id"], sens_by_pair[k]["pdb_id"])
               for k in shared if prim_by_pair[k]["pdb_id"] != sens_by_pair[k]["pdb_id"]]
    prim_pdbs = {p["pdb_id"] for p in prim_set}
    sens_pdbs = {f["pdb_id"] for f in final}
    missing_prim_structures = sorted(prim_pdbs - sens_pdbs)

    print("\n" + "=" * 92)
    print("COHORT COMPARISON")
    print("=" * 92)
    def summarise(name, st):
        print(f"\n{name}: n = {len(st)} structures, "
              f"{len({s['target_accession'] for s in st})} independent target clusters")
        for pt in PARTNERS:
            sub = [s for s in st if s["fusion_partner"] == pt]
            print(f"  {pt:5s}: {len(sub):3d} structures / "
                  f"{len({s['target_accession'] for s in sub})} targets")
    summarise("E9 = 0.20 PRIMARY (frozen)", prim_set)
    summarise("E9 = 0.35 SENSITIVITY", final)

    print(f"\nNEW (partner, target) units added : {len(added_pairs)}")
    print(f"units present in primary but absent in sensitivity : {len(lost_pairs)}")
    print(f"units whose REPRESENTATIVE STRUCTURE changed       : {len(swapped)}")
    print(f"primary PDB structures absent from the 0.35 cohort : "
          f"{len(missing_prim_structures)}"
          + (f"  {missing_prim_structures}" if missing_prim_structures else ""))
    if swapped:
        print("  representative swaps (unit, primary -> sensitivity):")
        for k, a, b in swapped:
            pa, pb = prim_by_pair[k], sens_by_pair[k]
            print(f"    {k[0]:5s} {k[1]:10s}  {a} ({pa['resolution']} A) -> "
                  f"{b} ({pb['resolution']} A)")

    if added_pairs:
        print(f"\nNEW STRUCTURES ({len(added_pairs)}):")
        print(f"  {'pdb':6s} {'partner':7s} {'target':10s} {'res':>6s} {'topology':18s} "
              f"{'tgt_dis':>8s} {'fus_dis':>8s}  description")
        for k in added_pairs:
            f = sens_by_pair[k]
            print(f"  {f['pdb_id']:6s} {f['fusion_partner']:7s} {f['target_accession']:10s} "
                  f"{f['resolution']:6.2f} {str(f['topology']):18s} "
                  f"{f['target_disorder']:8.3f} {f['fusion_disorder']:8.3f}  "
                  f"{(f['description'] or '')[:46]}")
    by_p = collections.Counter(sens_by_pair[k]["fusion_partner"] for k in added_pairs)
    print(f"\n  added by partner: {dict(by_p)}")

    json.dump({"threshold": SENS_DISORDER, "amendment": "J",
               "policy": "Sensitivity cohort only. Primary manuscript numbers remain E9 = 0.20.",
               "n": len(final),
               "per_partner": {pt: sum(1 for f in final if f["fusion_partner"] == pt)
                               for pt in PARTNERS},
               "independent_targets": len(conf_targets),
               "development_target_overlap": 0,
               "structures": final},
              open(os.path.join(MANI, "e9_035_set_final.json"), "w"), indent=1)

    cols = ["entity_id", "pdb_id", "fusion_partner", "target_accession", "topology", "method",
            "resolution", "representative_chain", "target_modelled", "fusion_modelled",
            "target_disorder", "fusion_disorder", "align_coverage_of_modelled", "description"]
    with open(os.path.join(MANI, "e9_035_set_final.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for f in final:
            fh.write("\t".join("" if f.get(c) is None else
                               str(f.get(c)).replace("\t", " ") for c in cols) + "\n")

    json.dump({"added_units": [list(k) for k in added_pairs],
               "added_pdb_ids": sorted(sens_by_pair[k]["pdb_id"] for k in added_pairs),
               "lost_units": [list(k) for k in lost_pairs],
               "representative_swaps": [[list(k), a, b] for k, a, b in swapped],
               "primary_structures_absent_from_sensitivity": missing_prim_structures,
               "entities_eligible_at_020_not_at_035": sorted(lost)},
              open(os.path.join(MANI, "e9_035_cohort_diff.json"), "w"), indent=1)

    if missing_prim_structures or lost_pairs:
        print("\n*** WATCH CONDITION J4 TRIGGERED — investigate before proceeding ***")
    else:
        print("\nWatch condition J4: clear. Relaxing the threshold only added observations.")


if __name__ == "__main__":
    main()
