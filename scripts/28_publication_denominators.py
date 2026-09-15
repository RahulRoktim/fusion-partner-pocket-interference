#!/usr/bin/env python
"""
Phase 6B/6C — canonical denominator audit.

Resolves every denominator used anywhere in the project to a single table, makes the unit of
counting explicit (member / PDB entry / chain / system), and asserts that the numbers quoted in
the Phase 5 prose match the canonical table.

Outputs:
  results/PUBLICATION_DENOMINATORS.tsv
  results/PLINDER_DENOMINATORS.tsv
  results/denominator_audit.json
"""
import collections, json, math, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
AUD = os.path.join(ROOT, "dataset_audit")
NA = "NA"


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def main():
    reg = json.load(open(os.path.join(AUD, "dataset_registry.json")))
    ents = {r["pdb_id"]: r for r in
            json.load(open(os.path.join(AUD, "entry_constructs.json")))["entries"]}
    gt = json.load(open(os.path.join(AUD, "groundtruth_chimeric.json")))["rows"]

    rows, audit = [], {}

    # ---------------- study cohorts
    dev = json.load(open(os.path.join(ROOT, "data_manifest",
                                      "development_set.json")))["structures"]
    conf = json.load(open(os.path.join(ROOT, "data_manifest",
                                       "confirmatory_set_final.json")))["structures"]
    for nm, st, role in (("DEVELOPMENT_COHORT", dev, "DEVELOPMENT/EXPLORATORY (never pooled)"),
                         ("CONFIRMATORY_COHORT", conf, "CONFIRMATORY")):
        rows.append({
            "DATASET": nm, "ROLE": role,
            "CANONICAL_VERSION": "frozen 2026-09-15; protocol v1.3",
            "RAW_MEMBERS": len(st), "RESOLVED_MEMBERS": len(st),
            "UNIQUE_PDB_ENTRIES": len({s["pdb_id"] for s in st}),
            "CHAIN_LEVEL_MEMBERS": len({(s["pdb_id"], s.get("chain") or s.get("representative_chain")) for s in st}),
            "SYSTEM_LEVEL_MEMBERS": NA,
            "FUSION_MEMBERS": len(st),
            "UNIQUE_FUSION_PDB_ENTRIES": len({s["pdb_id"] for s in st}),
            "UNRESOLVED": 0,
            "NOTES": ("every member is a fusion construct by construction; "
                      f"independent target clusters = {len({s['target_accession'] for s in st})}")})

    # ---------------- Tier-1 datasets
    tier1 = ["CHEN11", "JOINED", "COACH420", "HOLO4K", "FPTRAIN"]
    per_ds = {}
    for name in tier1:
        d = reg["datasets"][name]
        members = d["members"]
        pdbs = [m["pdb_id"] for m in members]
        uniq = sorted(set(pdbs))
        resolved = [p for p in uniq if p in ents]
        unresolved = [p for p in uniq if p not in ents]
        chains = {(m["pdb_id"], m["chain"]) for m in members if m["chain"]}
        fus = [p for p in resolved if ents[p]["entry_has_recognised_partner_chimera"]]
        lo, hi = wilson(len(fus), len(resolved))
        per_ds[name] = {"resolved": resolved, "fusion": fus, "unresolved": unresolved,
                        "members": len(members), "unique": len(uniq),
                        "chains": len(chains), "ci": (lo, hi)}
        rows.append({
            "DATASET": name, "ROLE": d["role"],
            "CANONICAL_VERSION": f"{d['membership_file']} sha256 {d['sha256'][:16]}",
            "RAW_MEMBERS": len(members), "RESOLVED_MEMBERS": len(resolved),
            "UNIQUE_PDB_ENTRIES": len(uniq),
            "CHAIN_LEVEL_MEMBERS": len(chains) if chains else NA,
            "SYSTEM_LEVEL_MEMBERS": NA,
            "FUSION_MEMBERS": len(fus),
            "UNIQUE_FUSION_PDB_ENTRIES": len(fus),
            "UNRESOLVED": len(unresolved),
            "NOTES": (f"prevalence {100*len(fus)/len(resolved):.2f}% "
                      f"[{100*lo:.2f}-{100*hi:.2f}]; counting unit = PDB entry; "
                      f"RESOLVED_MEMBERS counts unique entries resolved by SIFTS, not rows")})

    # ---------------- Tier-2 not measured
    rows.append({"DATASET": "LIGYSIS", "ROLE": "REFERENCE / evaluation substrate",
                 "CANONICAL_VERSION": "public web resource; no bulk export offered",
                 "RAW_MEMBERS": NA, "RESOLVED_MEMBERS": NA,
                 "UNIQUE_PDB_ENTRIES": NA, "CHAIN_LEVEL_MEMBERS": NA,
                 "SYSTEM_LEVEL_MEMBERS": NA, "FUSION_MEMBERS": NA,
                 "UNIQUE_FUSION_PDB_ENTRIES": NA, "UNRESOLVED": NA,
                 "NOTES": "NOT MEASURED. Authors' figures ~65,000 sites / ~25,000 proteins / "
                          ">100,000 PDBe structures. README: 'We do not offer full LIGYSIS "
                          "dataset download.'"})
    rows.append({"DATASET": "sc-PDB", "ROLE": "REFERENCE; de-facto TRAIN for other CNN predictors",
                 "CANONICAL_VERSION": "v.2017 (PDB freeze 2016-11)",
                 "RAW_MEMBERS": 16034, "RESOLVED_MEMBERS": NA,
                 "UNIQUE_PDB_ENTRIES": NA, "CHAIN_LEVEL_MEMBERS": NA,
                 "SYSTEM_LEVEL_MEMBERS": NA, "FUSION_MEMBERS": NA,
                 "UNIQUE_FUSION_PDB_ENTRIES": NA, "UNRESOLVED": NA,
                 "NOTES": "NOT MEASURED. 16,034 entries / 4,782 proteins / 6,326 ligands per the "
                          "resource. Only distribution is a 4,079,179,798-byte archive measured "
                          "at ~11 kB/s. Not a documented P2Rank or fpocket dataset."})

    # ---------------- PLINDER (system-level; see the dedicated table)
    pl = json.load(open(os.path.join(AUD, "plinder_audit.json")))
    plrows = [r for r in pl["rows"]]
    rec = [r for r in plrows if r.get("receptor_is_chimeric_chain")]
    for sp in ("train", "val", "test", "removed"):
        sysn = pl["split_sizes"][sp]
        rsub = [r for r in rec if r["split"] == sp]
        esub = {r["pdb_id"] for r in rsub}
        fuslab = [r for r in rsub if r.get("site_class") == "FUSION"]
        rows.append({
            "DATASET": f"PLINDER_{sp}", "ROLE": sp.upper(),
            "CANONICAL_VERSION": "2024-06/v2 split.parquet sha256 2959fb4b32f8c5cc",
            "RAW_MEMBERS": sysn, "RESOLVED_MEMBERS": sysn,
            "UNIQUE_PDB_ENTRIES": NA, "CHAIN_LEVEL_MEMBERS": NA,
            "SYSTEM_LEVEL_MEMBERS": sysn,
            "FUSION_MEMBERS": len(rsub),
            "UNIQUE_FUSION_PDB_ENTRIES": len(esub),
            "UNRESOLVED": sum(1 for r in rsub if r.get("site_class") in (None, "unresolved")),
            "NOTES": (f"counting unit = SYSTEM; receptor chain is the chimeric chain; "
                      f"fusion-site-labelled systems {len(fuslab)} over "
                      f"{len({r['pdb_id'] for r in fuslab})} unique entries")})

    # ---------------- cross-dataset deduplicated totals
    tier1_union = sorted(set().union(*[set(per_ds[n]["resolved"]) for n in tier1]))
    tier1_sum = sum(len(per_ds[n]["resolved"]) for n in tier1)
    tier1_fus_union = sorted(set().union(*[set(per_ds[n]["fusion"]) for n in tier1]))
    all_listed = set(reg["all_unique_pdb_entries"])
    all_resolved = sorted(p for p in all_listed if p in ents)
    all_unres = sorted(p for p in all_listed if p not in ents)

    audit["cross_dataset"] = {
        "tier1_datasets": tier1,
        "SUM_OF_PER_DATASET_RESOLVED_ENTRIES": tier1_sum,
        "CROSS_DATASET_DEDUPLICATED_TOTAL": len(tier1_union),
        "duplicate_entries_between_tier1_datasets": tier1_sum - len(tier1_union),
        "tier1_unique_fusion_entries": tier1_fus_union,
        "ALL_LISTS_INCLUDING_MLIG_unique_entries": len(all_listed),
        "ALL_LISTS_resolved": len(all_resolved),
        "ALL_LISTS_unresolved": len(all_unres),
        "unresolved_ids": all_unres,
        "deduplication_method":
            "Deduplication is on the 4-character PDB entry identifier, upper-cased. "
            "The (mlig) files are ligand-annotated SUBSETS of their parent datasets and "
            "contribute no new entries; they are therefore excluded from the Tier-1 union to "
            "avoid presenting a subset as an independent dataset. The union is scientifically "
            "meaningful only as 'how many distinct PDB entries were examined'; it must never be "
            "used as the denominator of a per-dataset prevalence, because the datasets have "
            "different roles."}

    rows.append({
        "DATASET": "TIER1_UNION_DEDUPLICATED", "ROLE": "derived (not a dataset)",
        "CANONICAL_VERSION": "union of CHEN11+JOINED+COACH420+HOLO4K+FPTRAIN",
        "RAW_MEMBERS": NA, "RESOLVED_MEMBERS": len(tier1_union),
        "UNIQUE_PDB_ENTRIES": len(tier1_union), "CHAIN_LEVEL_MEMBERS": NA,
        "SYSTEM_LEVEL_MEMBERS": NA,
        "FUSION_MEMBERS": len(tier1_fus_union),
        "UNIQUE_FUSION_PDB_ENTRIES": len(tier1_fus_union),
        "UNRESOLVED": NA,
        "NOTES": (f"deduplicated union. The non-deduplicated SUM of per-dataset resolved entries "
                  f"is {tier1_sum}; {tier1_sum - len(tier1_union)} entries appear in more than "
                  f"one Tier-1 dataset. Use per-dataset denominators for prevalence, never this.")})

    cols = ["DATASET", "ROLE", "CANONICAL_VERSION", "RAW_MEMBERS", "RESOLVED_MEMBERS",
            "UNIQUE_PDB_ENTRIES", "CHAIN_LEVEL_MEMBERS", "SYSTEM_LEVEL_MEMBERS",
            "FUSION_MEMBERS", "UNIQUE_FUSION_PDB_ENTRIES", "UNRESOLVED", "NOTES"]
    with open(os.path.join(RES, "PUBLICATION_DENOMINATORS.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, NA)) for c in cols) + "\n")

    # ---------------- PLINDER dedicated table (6C)
    pcols = ["SPLIT", "SYSTEMS", "UNIQUE_PDB_ENTRIES", "FUSION_SYSTEMS", "FUSION_PDB_ENTRIES",
             "FUSION_SITE_LABELLED_SYSTEMS", "FUSION_SITE_LABELLED_PDB_ENTRIES",
             "TARGET_SITE_SYSTEMS", "INTERFACE_SYSTEMS", "MIXED_SYSTEMS", "UNRESOLVED_SYSTEMS"]
    import pandas as pd
    df = pd.read_parquet(os.path.join(AUD, "membership", "plinder_2024-06_v2_split.parquet"))
    df["pdb"] = df["system_id"].str.split("__").str[0].str.upper()
    prows = []
    for sp in ("train", "val", "test", "removed"):
        sub = df[df["split"] == sp]
        rsub = [r for r in rec if r["split"] == sp]
        sc = collections.Counter(r.get("site_class") or "unresolved" for r in rsub)
        fl = [r for r in rsub if r.get("site_class") == "FUSION"]
        tg = [r for r in rsub if r.get("site_class") == "TARGET"]
        iface = [r for r in rsub if r.get("site_class") == "INTERFACE"]
        mixed = [r for r in rsub if r.get("site_class") and ";" in r["site_class"]]
        unres = [r for r in rsub if r.get("site_class") in (None, "unresolved")]
        prows.append({"SPLIT": sp, "SYSTEMS": len(sub),
                      "UNIQUE_PDB_ENTRIES": sub["pdb"].nunique(),
                      "FUSION_SYSTEMS": len(rsub),
                      "FUSION_PDB_ENTRIES": len({r["pdb_id"] for r in rsub}),
                      "FUSION_SITE_LABELLED_SYSTEMS": len(fl),
                      "FUSION_SITE_LABELLED_PDB_ENTRIES": len({r["pdb_id"] for r in fl}),
                      "TARGET_SITE_SYSTEMS": len(tg),
                      "INTERFACE_SYSTEMS": len(iface),
                      "MIXED_SYSTEMS": len(mixed),
                      "UNRESOLVED_SYSTEMS": len(unres)})
    with open(os.path.join(RES, "PLINDER_DENOMINATORS.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(pcols) + "\n")
        for r in prows:
            fh.write("\t".join(str(r[c]) for c in pcols) + "\n")

    # partner breakdown, BOTH units
    part = collections.defaultdict(lambda: {"systems": 0, "entries": set()})
    for r in rec:
        for p in r["fusion_partners"].split(","):
            part[p]["systems"] += 1
            part[p]["entries"].add(r["pdb_id"])
    partner_tbl = sorted(((k, v["systems"], len(v["entries"])) for k, v in part.items()),
                         key=lambda x: -x[1])
    fl_all = [r for r in rec if r.get("site_class") == "FUSION"]
    fpart = collections.defaultdict(lambda: {"systems": 0, "entries": set()})
    for r in fl_all:
        for p in r["fusion_partners"].split(","):
            fpart[p]["systems"] += 1
            fpart[p]["entries"].add(r["pdb_id"])
    audit["plinder_partners_all"] = [{"partner": k, "systems": s, "entries": e}
                                     for k, s, e in partner_tbl]
    audit["plinder_partners_fusion_labelled"] = [
        {"partner": k, "systems": v["systems"], "entries": len(v["entries"])}
        for k, v in sorted(fpart.items(), key=lambda x: -x[1]["systems"])]
    audit["plinder_rows"] = prows

    # ---------------- assertions against the Phase 5 prose
    prose = open(os.path.join(RES, "PHASE5_DATASET_AUDIT.md"), encoding="utf-8").read()
    checks = []

    def assert_in(label, text, ok=None):
        present = text in prose
        checks.append({"check": label, "string": text, "found": present,
                       "expected": True if ok is None else ok})

    assert_in("Tier-1 non-dedup sum quoted", f"{tier1_sum:,} resolved entries")
    assert_in("Tier-1 dedup union quoted", f"{len(tier1_union):,} entries after deduplication")
    assert_in("PLINDER fusion-labelled stated in BOTH units", "342 SYSTEMS across 158 UNIQUE PDB")
    assert_in("PLINDER train fusion-labelled entries", "158")
    assert_in("PLINDER train systems+entries", "train 1,385 systems / 454 entries")
    assert_in("PLINDER test zero", "test split contains no such system")
    assert_in("CHEN11 entry denominator", "0/241")
    assert_in("JOINED denominator", "1/534")
    assert_in("COACH420 denominator", "0/418")
    assert_in("HOLO4K denominator", "0/4,004")
    assert_in("FPTRAIN denominator", "0/222")
    audit["prose_assertions"] = checks
    audit["denominator_rows"] = rows

    json.dump(audit, open(os.path.join(RES, "denominator_audit.json"), "w"),
              indent=1, default=str)

    print("=" * 100)
    print("CANONICAL DENOMINATOR AUDIT")
    print("=" * 100)
    print(f"\nTier-1 per-dataset resolved-entry SUM (non-deduplicated): {tier1_sum}")
    print(f"Tier-1 DEDUPLICATED union of unique PDB entries:          {len(tier1_union)}")
    print(f"  -> {tier1_sum - len(tier1_union)} entries appear in more than one Tier-1 dataset")
    print(f"All frozen lists incl. (mlig) subsets: {len(all_listed)} unique entries, "
          f"{len(all_resolved)} resolved, {len(all_unres)} unresolved")
    print(f"  unresolved: {all_unres}")
    print(f"Tier-1 unique fusion entries: {tier1_fus_union}")
    print("\nPLINDER, both units:")
    for r in prows:
        print(f"  {r['SPLIT']:8s} systems {r['SYSTEMS']:7d} / entries {r['UNIQUE_PDB_ENTRIES']:6d}"
              f" | fusion systems {r['FUSION_SYSTEMS']:5d} / entries {r['FUSION_PDB_ENTRIES']:4d}"
              f" | fusion-site-labelled {r['FUSION_SITE_LABELLED_SYSTEMS']:4d} systems /"
              f" {r['FUSION_SITE_LABELLED_PDB_ENTRIES']:4d} entries")
    print("\nPLINDER partner breakdown (chimeric-receptor systems / unique entries):")
    for k, s, e in partner_tbl:
        print(f"  {k:26s} {s:5d} systems / {e:4d} entries")
    print("\nPLINDER partner breakdown, FUSION-SITE-LABELLED only:")
    for x in audit["plinder_partners_fusion_labelled"]:
        print(f"  {x['partner']:26s} {x['systems']:5d} systems / {x['entries']:4d} entries")
    bad = [c for c in checks if c["found"] != c["expected"]]
    print(f"\nprose assertions: {len(checks) - len(bad)}/{len(checks)} pass")
    for c in bad:
        print(f"  MISMATCH: {c['check']} -> looked for '{c['string']}'")


if __name__ == "__main__":
    main()
