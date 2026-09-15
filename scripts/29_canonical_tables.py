#!/usr/bin/env python
"""
Phase 6G — manuscript-ready canonical tables, generated from frozen result files only.
No headline number is typed by hand anywhere in this script.

Outputs (results/tables/):
  TABLE1_cohort_composition.tsv
  TABLE2_confirmatory_exposure.tsv
  TABLE3_displacement.tsv
  TABLE4_reference_site.tsv
  TABLE5_dataset_provenance.tsv
"""
import collections, json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
CONF = os.path.join(RES, "confirmatory", "primary")
MANI = os.path.join(ROOT, "data_manifest")
TAB = os.path.join(RES, "tables")
os.makedirs(TAB, exist_ok=True)
PARTNERS = ["BRIL", "T4L", "MBP"]
DETECTORS = ["p2rank", "fpocket"]


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def write(name, cols, rows):
    with open(os.path.join(TAB, name), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("NA" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    print(f"  wrote {name}  ({len(rows)} rows)")


def med(v):
    s = sorted(v)
    n = len(s)
    return float("nan") if not n else (s[n // 2] if n % 2 else (s[n//2-1]+s[n//2])/2)


def main():
    conf_set = json.load(open(os.path.join(MANI, "confirmatory_set_final.json")))
    dev_set = json.load(open(os.path.join(MANI, "development_set.json")))
    cman = {p["pdb_id"]: p for p in
            json.load(open(os.path.join(MANI,
                                        "confirmatory_manifest_primary.json")))["prepared"]
            if "PREPARATION_FAILED" not in p}
    summ = json.load(open(os.path.join(CONF, "confirmatory_summary.json")))
    disp = [r for r in open(os.path.join(CONF, "displacement.tsv"), encoding="utf-8")]
    hdr = disp[0].rstrip("\n").split("\t")
    drows = [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in disp[1:]]
    mech = json.load(open(os.path.join(CONF, "mechanism.json")))["rows"]
    print("building canonical tables ...")

    # ---------------- TABLE 1 : cohort composition
    t1 = []
    for label, st, man in (("CONFIRMATORY", conf_set["structures"], cman),
                           ("DEVELOPMENT (exploratory, never pooled)", dev_set["structures"], None)):
        for pt in PARTNERS + ["ALL"]:
            sub = [s for s in st if pt == "ALL" or s["fusion_partner"] == pt]
            if not sub:
                continue
            res = sorted(s["resolution"] for s in sub if s.get("resolution"))
            topo = collections.Counter(s["topology"] for s in sub)
            meth = collections.Counter(s["method"] for s in sub)
            if man:
                tsz = [man[s["pdb_id"]]["residue_class_counts"].get("TARGET", 0)
                       for s in sub if s["pdb_id"] in man]
                fsz = [man[s["pdb_id"]]["residue_class_counts"].get("FUSION", 0)
                       for s in sub if s["pdb_id"] in man]
                refs = sum(1 for s in sub if s["pdb_id"] in man
                           and man[s["pdb_id"]]["reference_site_available"])
            else:
                tsz = [s.get("target_residues_modelled") or 0 for s in sub]
                fsz = [s.get("fusion_residues_modelled") or 0 for s in sub]
                refs = sum(1 for s in sub if s.get("reference_site_available"))
            t1.append({
                "COHORT": label, "PARTNER": pt, "N_STRUCTURES": len(sub),
                "N_INDEPENDENT_TARGETS": len({s["target_accession"] for s in sub}),
                "MEDIAN_RESOLUTION_A": round(med(res), 2) if res else None,
                "RESOLUTION_RANGE_A": f"{min(res):.2f}-{max(res):.2f}" if res else None,
                "N_XRAY": meth.get("X-RAY DIFFRACTION", 0),
                "N_CRYOEM": meth.get("ELECTRON MICROSCOPY", 0),
                "TOPOLOGY_INTERNAL": topo.get("INTERNAL_INSERTION", 0),
                "TOPOLOGY_TERMINAL_N": topo.get("TERMINAL_N", 0),
                "TOPOLOGY_TERMINAL_C": topo.get("TERMINAL_C", 0),
                "TOPOLOGY_COMPLEX": topo.get("COMPLEX", 0),
                "MEDIAN_TARGET_RESIDUES": int(med(tsz)) if tsz else None,
                "MEDIAN_FUSION_RESIDUES": int(med(fsz)) if fsz else None,
                "N_WITH_BIOLOGICAL_REFERENCE_SITE": refs})
    write("TABLE1_cohort_composition.tsv", list(t1[0].keys()), t1)

    # ---------------- TABLE 2 : confirmatory exposure
    t2 = []
    for det in DETECTORS:
        D = summ["detectors"][det]["A"]
        for pt in PARTNERS:
            bp = D["by_partner"][pt]
            row = {"DETECTOR": det, "STRATUM": pt, "N_STRUCTURES": bp["n"]}
            for depth, key in (("RANK1", None), ("TOP3", "top-3"), ("TOP5", "top-5")):
                if key is None:
                    k, n = bp["k"], bp["n"]
                else:
                    k, n = bp[key]["k"], bp[key]["n"]
                lo, hi = wilson(k, n)
                row[f"{depth}_K"] = k
                row[f"{depth}_N"] = n
                row[f"{depth}_PCT"] = round(100 * k / n, 1) if n else None
                row[f"{depth}_CI95"] = f"{100*lo:.1f}-{100*hi:.1f}" if n else None
            t2.append(row)
        # pooled + standardised
        a = D["depths"]
        row = {"DETECTOR": det, "STRATUM": "ALL (raw pooled)", "N_STRUCTURES": a["rank-1"]["n"]}
        for depth, key in (("RANK1", "rank-1"), ("TOP3", "top-3"), ("TOP5", "top-5")):
            k, n = a[key]["k"], a[key]["n"]
            lo, hi = wilson(k, n)
            row[f"{depth}_K"], row[f"{depth}_N"] = k, n
            row[f"{depth}_PCT"] = round(100 * k / n, 1)
            row[f"{depth}_CI95"] = f"{100*lo:.1f}-{100*hi:.1f}"
        t2.append(row)
        sl, sh = D["standardised_ci"]
        srow = {"DETECTOR": det, "STRATUM": "ALL (partner-standardised 1/3 each)",
                "N_STRUCTURES": a["rank-1"]["n"],
                "RANK1_K": "NA", "RANK1_N": "NA",
                "RANK1_PCT": round(100 * D["standardised"], 1),
                "RANK1_CI95": f"{100*float(sl):.1f}-{100*float(sh):.1f} (cluster bootstrap)"}
        for depth in ("TOP3", "TOP5"):
            rates = [D["by_partner"][p][depth.replace("TOP", "top-")]["k"]
                     / D["by_partner"][p][depth.replace("TOP", "top-")]["n"] for p in PARTNERS]
            srow[f"{depth}_K"] = "NA"
            srow[f"{depth}_N"] = "NA"
            srow[f"{depth}_PCT"] = round(100 * sum(rates) / 3, 1)
            srow[f"{depth}_CI95"] = "NA"
        t2.append(srow)
    write("TABLE2_confirmatory_exposure.tsv", list(t2[0].keys()), t2)

    # ---------------- TABLE 3 : displacement
    t3 = []
    for det in DETECTORS:
        B = summ["detectors"][det]["B"]
        sub = [r for r in drows if r["detector"] == det]
        ev = [r for r in sub if r["state"] == "EVALUABLE"]
        for pt in PARTNERS + ["ALL"]:
            e = [r for r in ev if pt == "ALL" or r["partner"] == pt]
            s = [r for r in sub if pt == "ALL" or r["partner"] == pt]
            d = [int(r["displacement"]) for r in e]
            row = {"DETECTOR": det, "STRATUM": pt,
                   "N_STRUCTURES": len(s), "N_EVALUABLE": len(e),
                   "N_DETECTOR_ABSTENTION": sum(1 for r in s
                                                if r["state"] == "DETECTOR_ABSTENTION"),
                   "N_UNMATCHED": sum(1 for r in s
                                      if r["state"] == "TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL"),
                   "MEDIAN_DISPLACEMENT": med(d) if d else None,
                   "MEAN_DISPLACEMENT": round(sum(d) / len(d), 2) if d else None,
                   "MAX_DISPLACEMENT": max(d) if d else None,
                   "N_NEGATIVE_DISPLACEMENT": sum(1 for x in d if x < 0)}
            for thr, lab in ((1, "GT0"), (2, "GE2"), (5, "GE5")):
                k = sum(1 for x in d if x >= thr)
                lo, hi = wilson(k, len(d))
                row[f"{lab}_K"] = k
                row[f"{lab}_N"] = len(d)
                row[f"{lab}_PCT"] = round(100 * k / len(d), 1) if d else None
                row[f"{lab}_CI95"] = f"{100*lo:.1f}-{100*hi:.1f}" if d else None
            if pt == "ALL":
                row["SIGN_TEST_P"] = f"{float(B['sign_test_p']):.3g}"
                row["MEAN_CI95_CLUSTER_BOOT"] = (f"{float(B['mean_ci'][0]):+.2f} to "
                                                 f"{float(B['mean_ci'][1]):+.2f}")
            t3.append(row)
    write("TABLE3_displacement.tsv", list(t3[0].keys()), t3)

    # ---------------- TABLE 4 : biological reference-site subset
    t4 = []
    cov = summ["reference_site_coverage"]
    for det in DETECTORS:
        C = summ["detectors"][det].get("C")
        if not C:
            continue
        found = C["found"]
        lo, hi = wilson(C["outranked"], found)
        t4.append({
            "DETECTOR": det,
            "COVERAGE_K": cov["k"], "COVERAGE_N": cov["n"],
            "COVERAGE_PCT": round(100 * cov["fraction"], 1),
            "SITE_POCKET_FOUND_K": found, "SITE_POCKET_FOUND_N": C["coverage"],
            "OUTRANKED_BY_FUSION_K": C["outranked"], "OUTRANKED_BY_FUSION_N": found,
            "OUTRANKED_PCT": round(100 * C["outranked"] / found, 1),
            "OUTRANKED_CI95": f"{100*lo:.1f}-{100*hi:.1f}",
            "RANK_IMPROVED_AFTER_REMOVAL": C["improved"],
            "RANK_WORSENED": C["worsened"], "RANK_UNCHANGED": C["unchanged"],
            "MEDIAN_RANK_CHANGE": C["median_change"],
            "RECOVERY_TOP1_K": C["top1"], "RECOVERY_TOP1_PCT": round(100 * C["top1"] / found, 1),
            "RECOVERY_TOP3_K": C["top3"], "RECOVERY_TOP3_PCT": round(100 * C["top3"] / found, 1),
            "RECOVERY_TOP5_K": C["top5"], "RECOVERY_TOP5_PCT": round(100 * C["top5"] / found, 1)})
    write("TABLE4_reference_site.tsv", list(t4[0].keys()), t4)

    # ---------------- TABLE 5 : dataset provenance
    den = [l.rstrip("\n").split("\t") for l in
           open(os.path.join(RES, "PUBLICATION_DENOMINATORS.tsv"), encoding="utf-8")]
    dh, drow = den[0], [dict(zip(den[0], r)) for r in den[1:]]
    keep = ["CHEN11", "JOINED", "COACH420", "HOLO4K", "FPTRAIN", "LIGYSIS", "sc-PDB",
            "PLINDER_train", "PLINDER_val", "PLINDER_test", "PLINDER_removed"]
    pl = [l.rstrip("\n").split("\t") for l in
          open(os.path.join(RES, "PLINDER_DENOMINATORS.tsv"), encoding="utf-8")]
    plr = {r[0]: dict(zip(pl[0], r)) for r in pl[1:]}
    t5 = []
    for r in drow:
        if r["DATASET"] not in keep:
            continue
        sp = r["DATASET"].replace("PLINDER_", "") if r["DATASET"].startswith("PLINDER") else None
        p = plr.get(sp, {})
        t5.append({
            "DATASET": r["DATASET"], "ROLE": r["ROLE"],
            "COUNTING_UNIT": ("SYSTEM" if sp else
                              ("NOT MEASURED" if r["RESOLVED_MEMBERS"] == "NA" else "PDB ENTRY")),
            "CANONICAL_VERSION": r["CANONICAL_VERSION"],
            "N_MEMBERS": r["RAW_MEMBERS"], "N_RESOLVED": r["RESOLVED_MEMBERS"],
            "N_UNIQUE_PDB_ENTRIES": (p.get("UNIQUE_PDB_ENTRIES") if sp
                                     else r["UNIQUE_PDB_ENTRIES"]),
            "N_CHAINS": r["CHAIN_LEVEL_MEMBERS"],
            "FUSION_CONSTRUCTS": r["FUSION_MEMBERS"],
            "FUSION_UNIQUE_PDB_ENTRIES": r["UNIQUE_FUSION_PDB_ENTRIES"],
            "FUSION_SITE_LABELLED_SYSTEMS": p.get("FUSION_SITE_LABELLED_SYSTEMS", "NA"),
            "FUSION_SITE_LABELLED_ENTRIES": p.get("FUSION_SITE_LABELLED_PDB_ENTRIES", "NA"),
            "TARGET_SITE_SYSTEMS": p.get("TARGET_SITE_SYSTEMS", "NA"),
            "UNRESOLVED": r["UNRESOLVED"], "NOTES": r["NOTES"]})
    write("TABLE5_dataset_provenance.tsv", list(t5[0].keys()), t5)

    # mechanism summary appended as Table 3b
    t3b = []
    for det in DETECTORS:
        for pt in PARTNERS + ["ALL"]:
            sub = [m for m in mech if m["detector"] == det and (pt == "ALL" or m["partner"] == pt)]
            if not sub:
                continue
            c = collections.Counter(m["mechanism"] for m in sub)
            s = collections.Counter(m["probe_state"] for m in sub)
            t3b.append({"DETECTOR": det, "STRATUM": pt, "N_RANK1_FAILURES": len(sub),
                        "FUSION_INTRINSIC_CAVITY": c.get("FUSION_INTRINSIC_CAVITY", 0),
                        "TARGET_FUSION_INTERFACE": c.get("TARGET_FUSION_INTERFACE", 0),
                        "CHIMERA_SPECIFIC_FUSION_CAVITY":
                            c.get("CHIMERA_SPECIFIC_FUSION_CAVITY", 0),
                        "LINKER_RELATED": c.get("LINKER_RELATED", 0),
                        "MAPPING_ARTIFACT": c.get("MAPPING_ARTIFACT", 0),
                        "PROBE_MATCHED": s.get("MATCHED_INTRINSIC_CAVITY", 0),
                        "PROBE_NO_MATCHING_CAVITY": s.get("NO_MATCHING_CAVITY", 0),
                        "PROBE_DETECTOR_NO_PREDICTIONS": s.get("DETECTOR_NO_PREDICTIONS", 0),
                        "PROBE_AMBIGUOUS": s.get("AMBIGUOUS", 0)})
    write("TABLE3b_mechanism.tsv", list(t3b[0].keys()), t3b)


if __name__ == "__main__":
    main()
