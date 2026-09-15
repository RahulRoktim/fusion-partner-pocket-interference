#!/usr/bin/env python
"""
Phase 3A.4 — audit of the six structures for which P2Rank returned zero pockets in the
FUSION-REMOVED condition: 3N94, 3OAI, 5GPP, 5H7Q, 5JQE, 5YQR.

These are NOT assumed to be technical failures. The question is whether each is a genuine
detector abstention (P2Rank ran correctly on valid input and predicted nothing) or an
implementation/input problem. No structure is removed from the study on account of a zero
prediction.
"""
import collections, json, math, os
import gemmi

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
RESULTS = os.path.join(ROOT, "results")
CACHE = os.path.join(ROOT, ".structure_cache")
CASES = ["3N94", "3OAI", "5GPP", "5H7Q", "5JQE", "5YQR"]


def secondary_structure(pdb_id, chain, target_auth_keys):
    """Depositor-assigned SS from the mmCIF, restricted to TARGET residues."""
    try:
        blk = gemmi.cif.read(os.path.join(CACHE, f"{pdb_id.lower()}.cif")).sole_block()
    except Exception:
        return None
    want = {int(''.join(c for c in k if c.isdigit() or c == '-')) for k in target_auth_keys}
    helix = strand = 0
    tab = blk.find("_struct_conf.", ["conf_type_id", "beg_auth_asym_id", "beg_auth_seq_id",
                                     "end_auth_seq_id"])
    for row in tab:
        try:
            if row[1].strip('"\'') != chain:
                continue
            b, e = int(row[2]), int(row[3])
            n = len([p for p in want if b <= p <= e])
            if row[0].upper().startswith("HELX"):
                helix += n
        except Exception:
            continue
    tab = blk.find("_struct_sheet_range.", ["beg_auth_asym_id", "beg_auth_seq_id",
                                            "end_auth_seq_id"])
    for row in tab:
        try:
            if row[0].strip('"\'') != chain:
                continue
            b, e = int(row[1]), int(row[2])
            strand += len([p for p in want if b <= p <= e])
        except Exception:
            continue
    n = len(want)
    return {"n_target": n, "helix": helix, "strand": strand,
            "coil": n - helix - strand,
            "pct_helix": round(100 * helix / n, 1) if n else None,
            "pct_strand": round(100 * strand / n, 1) if n else None}


def validate_input(path):
    st = gemmi.read_structure(path)
    res = [r for ch in st[0] for r in ch]
    atoms = [a for r in res for a in r]
    no_ca = [f"{r.seqid.num}" for r in res if r.find_atom("CA", "*") is None]
    xs = [a.pos.x for a in atoms]; ys = [a.pos.y for a in atoms]; zs = [a.pos.z for a in atoms]
    cx, cy, cz = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
    rg = math.sqrt(sum((a.pos.x-cx)**2 + (a.pos.y-cy)**2 + (a.pos.z-cz)**2
                       for a in atoms) / len(atoms))
    # chain breaks: consecutive CA-CA distance > 4.5 A
    cas = [r.find_atom("CA", "*") for r in res]
    cas = [c for c in cas if c is not None]
    breaks = sum(1 for a, b in zip(cas, cas[1:]) if a.pos.dist(b.pos) > 4.5)
    nan = sum(1 for a in atoms if any(math.isnan(v) for v in (a.pos.x, a.pos.y, a.pos.z)))
    return {"n_residues": len(res), "n_atoms": len(atoms), "residues_without_CA": len(no_ca),
            "chain_breaks_gt_4.5A": breaks, "nan_coords": nan,
            "radius_of_gyration": round(rg, 2),
            "max_dim": round(max(max(xs)-min(xs), max(ys)-min(ys), max(zs)-min(zs)), 1)}


def main():
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI, "pilot_manifest.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    P = json.load(open(os.path.join(RESULTS, "pockets_classified.json")))["pockets"]
    runs = json.load(open(os.path.join(RESULTS, "detector_runs.json")))["runs"]
    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    runmap = {(r["detector"], r["pdb_id"], r["condition"]): r for r in runs}

    rows = []
    print("=" * 100)
    print("ZERO-PREDICTION AUDIT — P2Rank, FUSION-REMOVED condition")
    print("=" * 100)
    for s in CASES:
        rec = man[s]
        removed = os.path.join(ROOT, rec["removed_pdb"])
        cmap = rec["auth_class_map"]
        tkeys = [k for k, v in cmap.items() if v == "TARGET"]
        v = validate_input(removed)
        ss = secondary_structure(s, rec["auth_chain"], tkeys)
        run = runmap[("p2rank", s, "REMOVED")]
        fp_rem = len(idx.get(("fpocket", s, "REMOVED"), []))
        p2_orig = len(idx.get(("p2rank", s, "ORIGINAL"), []))
        fp_orig = len(idx.get(("fpocket", s, "ORIGINAL"), []))
        row = {
            "pdb_id": s, "partner": rec["fusion_partner"], "topology": rec["topology"],
            "resolution": rec["resolution"],
            "removed_target_residues": rec["residue_class_counts"].get("TARGET"),
            "p2rank_returncode": run["returncode"], "p2rank_runtime_s": run["runtime_s"],
            "p2rank_completed_normally": run["success"],
            "p2rank_pockets_ORIGINAL": p2_orig, "p2rank_pockets_REMOVED": 0,
            "fpocket_pockets_ORIGINAL": fp_orig, "fpocket_pockets_REMOVED": fp_rem,
            **{f"input_{k}": val for k, val in v.items()},
            **({f"ss_{k}": val for k, val in ss.items()} if ss else {}),
        }
        # fpocket's best REMOVED score, as a cross-check on "is there anything there?"
        fpr = sorted(idx.get(("fpocket", s, "REMOVED"), []), key=lambda r: r["rank"])
        row["fpocket_best_REMOVED_score"] = fpr[0]["score"] if fpr else None
        row["fpocket_REMOVED_scores_top3"] = [round(p["score"], 3) for p in fpr[:3]]
        # classification
        problems = []
        if not run["success"] or run["returncode"] != 0:
            problems.append("p2rank_did_not_complete")
        if v["nan_coords"]:
            problems.append("nan_coordinates")
        if v["n_residues"] < 1 or v["n_atoms"] < 10:
            problems.append("empty_or_tiny_input")
        if v["residues_without_CA"] > 0.2 * v["n_residues"]:
            problems.append("many_residues_without_CA")
        row["implementation_problems"] = ";".join(problems) if problems else "none"
        row["classification"] = ("IMPLEMENTATION_OR_INPUT_PROBLEM" if problems
                                 else "GENUINE_DETECTOR_ABSTENTION")
        rows.append(row)

        print(f"\n{s}  {rec['fusion_partner']}  {rec['topology']}  {rec['resolution']} A")
        print(f"  removed target residues        {row['removed_target_residues']}")
        print(f"  input: {v['n_residues']} residues, {v['n_atoms']} atoms, "
              f"Rg {v['radius_of_gyration']} A, max dim {v['max_dim']} A, "
              f"chain breaks {v['chain_breaks_gt_4.5A']}, no-CA {v['residues_without_CA']}, "
              f"NaN {v['nan_coords']}")
        if ss:
            print(f"  target secondary structure     helix {ss['pct_helix']}%, "
                  f"strand {ss['pct_strand']}%, coil "
                  f"{round(100*ss['coil']/ss['n_target'],1) if ss['n_target'] else 0}%")
        print(f"  P2Rank  rc={run['returncode']} runtime={run['runtime_s']}s "
              f"completed={run['success']}   ORIGINAL {p2_orig} pockets -> REMOVED 0")
        print(f"  fpocket ORIGINAL {fp_orig} pockets -> REMOVED {fp_rem} pockets "
              f"(best score {row['fpocket_best_REMOVED_score']}, "
              f"top3 {row['fpocket_REMOVED_scores_top3']})")
        print(f"  -> {row['classification']}  ({row['implementation_problems']})")

    cols = sorted({k for r in rows for k in r})
    with open(os.path.join(RESULTS, "zero_prediction_audit.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    json.dump(rows, open(os.path.join(RESULTS, "zero_prediction_audit.json"), "w"), indent=1)

    print("\n" + "=" * 100)
    print("SUMMARY:", dict(collections.Counter(r["classification"] for r in rows)))
    # context: what does P2Rank do on small targets generally?
    small = []
    for s in sorted(man):
        n = man[s]["residue_class_counts"].get("TARGET", 0)
        small.append((n, s, len(idx.get(("p2rank", s, "REMOVED"), []))))
    small.sort()
    print("\nP2Rank REMOVED-condition pocket count vs target size (all 24, ascending size):")
    for n, s, k in small:
        print(f"  {s:6s} {man[s]['fusion_partner']:5s} target={n:4d} residues -> {k} pockets")


if __name__ == "__main__":
    main()
