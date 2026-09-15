#!/usr/bin/env python
"""
Phase 3A.1 — deterministic TARGET-POCKET CORRESPONDENCE across the ORIGINAL and
FUSION-REMOVED conditions, developed on the 24 DEVELOPMENT structures.

Rationale. The pilot's paired outcome ("rank-1 becomes TARGET_DOMINATED") is structurally
degenerate: after FUSION/LINKER residues are deleted, essentially every pocket is composed of
TARGET residues, so the transition cannot demonstrate restoration of any *specific* cavity.
Correspondence fixes this by identifying the SAME cavity in both conditions and then asking what
happened to its rank.

Method constraints (all satisfied below):
  * only TARGET residues are used;
  * matching is deterministic geometry / residue overlap;
  * one-to-one assignment;
  * detector score and rank are NEVER inputs to correspondence;
  * thresholds are chosen on matching-quality metrics that are independent of the rescue outcome.

TARGET coordinates are byte-identical between conditions (asserted by test I5), so a target-residue
centroid is directly comparable across conditions.

Outputs:
  results/correspondence_threshold_sweep.tsv   method-selection evidence (no endpoint numbers)
  results/matched_target_pockets.tsv           per matched pair
  results/correspondence_summary.json
"""
import collections, json, math, os, sys
import gemmi
from scipy.optimize import linear_sum_assignment

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
RESULTS = os.path.join(ROOT, "results")

MIN_TARGET_RESIDUES = 3          # a pocket needs this many TARGET residues to be matchable

# ---- FROZEN PRIMARY RULE (selected below on matching-quality evidence, see report) ----
PRIMARY_JACCARD = 0.40
PRIMARY_MAX_CENTROID_DIST = 8.0
# ---- pre-specified sensitivity rules ----
SENSITIVITY = [("J0.25/12A", 0.25, 12.0), ("J0.40/8A", 0.40, 8.0), ("J0.60/5A", 0.60, 5.0)]
SWEEP_J = [0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70]
SWEEP_D = [5.0, 8.0, 12.0, 999.0]


def load_target_coords(manrec):
    """auth residue key -> list of heavy-atom coords, from the prepared ORIGINAL file.
    Target coordinates are identical in the REMOVED file (test I5)."""
    st = gemmi.read_structure(os.path.join(ROOT, manrec["original_pdb"]))
    cmap = manrec["auth_class_map"]
    out = {}
    for ch in st[0]:
        for res in ch:
            key = f"{res.seqid.num}{res.seqid.icode.strip()}"
            if cmap.get(key) == "TARGET":
                out[key] = [(a.pos.x, a.pos.y, a.pos.z) for a in res]
    return out


def target_set(pocket, cmap):
    res = pocket["residues"].split(";") if pocket["residues"] else []
    return {r for r in res if cmap.get(r) == "TARGET"}


def centroid(keys, coords):
    pts = [p for k in keys for p in coords.get(k, [])]
    if not pts:
        return None
    n = len(pts)
    return (sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n, sum(p[2] for p in pts) / n)


def jaccard(a, b):
    u = len(a | b)
    return len(a & b) / u if u else 0.0


def overlap_coef(a, b):
    m = min(len(a), len(b))
    return len(a & b) / m if m else 0.0


def build_pairs(orig, rem, coords, cmap):
    """Returns (candidates, O, R) where candidates[i][j] holds the metrics for orig i vs rem j."""
    O = [p for p in orig if len(target_set(p, cmap)) >= MIN_TARGET_RESIDUES]
    R = [p for p in rem if len(target_set(p, cmap)) >= MIN_TARGET_RESIDUES]
    cand = []
    for o in O:
        to = target_set(o, cmap)
        co = centroid(to, coords)
        row = []
        for r in R:
            tr = target_set(r, cmap)
            cr = centroid(tr, coords)
            d = math.dist(co, cr) if (co and cr) else float("inf")
            row.append({"jaccard": jaccard(to, tr), "overlap": overlap_coef(to, tr),
                        "dist": d, "n_shared": len(to & tr),
                        "n_orig": len(to), "n_rem": len(tr)})
        cand.append(row)
    return cand, O, R


def assign(cand, O, R, jmin, dmax):
    """Optimal one-to-one assignment maximising total Jaccard, then threshold filtering.
    Rank and score play no part."""
    if not O or not R:
        return []
    cost = [[-(c["jaccard"] if c["dist"] <= dmax else 0.0) for c in row] for row in cand]
    ri, ci = linear_sum_assignment(cost)
    out = []
    for i, j in zip(ri, ci):
        m = cand[i][j]
        if m["jaccard"] >= jmin and m["dist"] <= dmax:
            out.append((i, j, m))
    return out


def main():
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI, "pilot_manifest.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    P = json.load(open(os.path.join(RESULTS, "pockets_classified.json")))["pockets"]
    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    for v in idx.values():
        v.sort(key=lambda r: r["rank"])

    coords = {s: load_target_coords(man[s]) for s in man}
    structures = sorted(man)
    detectors = ["p2rank", "fpocket"]

    # cache candidate matrices once
    cache = {}
    for det in detectors:
        for s in structures:
            o = idx.get((det, s, "ORIGINAL"), [])
            r = idx.get((det, s, "REMOVED"), [])
            cache[(det, s)] = build_pairs(o, r, coords[s], man[s]["auth_class_map"])

    # ---------------------------------------------------------------- method selection
    # Metrics here are all matching-quality metrics. None of them involves rank change,
    # rescue, or any endpoint quantity.
    print("=" * 104)
    print("CORRESPONDENCE METHOD SELECTION — matching-quality metrics only")
    print("(no endpoint / rescue quantity is computed in this section)")
    print("=" * 104)
    sweep_rows = []
    print(f"{'Jmin':>5s} {'dmax':>6s} {'det':8s} {'matched':>8s} {'%remMatched':>12s} "
          f"{'meanJ':>7s} {'meanD':>7s} {'ambig':>7s} {'medShared':>9s}")
    for jmin in SWEEP_J:
        for dmax in SWEEP_D:
            for det in detectors:
                nm = nrem = 0
                js, ds, amb, shared = [], [], [], []
                for s in structures:
                    cand, O, R = cache[(det, s)]
                    nrem += len(R)
                    pairs = assign(cand, O, R, jmin, dmax)
                    nm += len(pairs)
                    for i, j, m in pairs:
                        js.append(m["jaccard"]); ds.append(m["dist"])
                        shared.append(m["n_shared"])
                        # ambiguity: 2nd-best Jaccard for this removed pocket / best
                        col = sorted((cand[k][j]["jaccard"] for k in range(len(O))),
                                     reverse=True)
                        amb.append(col[1] / col[0] if len(col) > 1 and col[0] > 0 else 0.0)
                row = {"jmin": jmin, "dmax": dmax, "detector": det, "n_matched": nm,
                       "n_removed_matchable": nrem,
                       "pct_removed_matched": round(100 * nm / nrem, 1) if nrem else None,
                       "mean_jaccard": round(sum(js) / len(js), 3) if js else None,
                       "mean_dist": round(sum(ds) / len(ds), 2) if ds else None,
                       "mean_ambiguity": round(sum(amb) / len(amb), 3) if amb else None,
                       "median_shared_residues": sorted(shared)[len(shared)//2] if shared else None}
                sweep_rows.append(row)
                if dmax in (8.0, 999.0):
                    print(f"{jmin:5.2f} {dmax:6.0f} {det:8s} {nm:8d} "
                          f"{row['pct_removed_matched'] or 0:12.1f} "
                          f"{row['mean_jaccard'] or 0:7.3f} {row['mean_dist'] or 0:7.2f} "
                          f"{row['mean_ambiguity'] or 0:7.3f} "
                          f"{str(row['median_shared_residues']):>9s}")
    with open(os.path.join(RESULTS, "correspondence_threshold_sweep.tsv"), "w",
              encoding="utf-8") as fh:
        cols = list(sweep_rows[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in sweep_rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---------------------------------------------------------------- apply frozen rule
    print("\n" + "=" * 104)
    print(f"FROZEN PRIMARY RULE: Jaccard(target residues) >= {PRIMARY_JACCARD} "
          f"AND target-centroid distance <= {PRIMARY_MAX_CENTROID_DIST} A, "
          f"one-to-one optimal assignment")
    print("=" * 104)

    rows, per_struct = [], []
    for det in detectors:
        for s in structures:
            cand, O, R = cache[(det, s)]
            pairs = assign(cand, O, R, PRIMARY_JACCARD, PRIMARY_MAX_CENTROID_DIST)
            matched_o = {i for i, _, _ in pairs}
            matched_r = {j for _, j, _ in pairs}
            cmap = man[s]["auth_class_map"]
            orig_all = idx.get((det, s, "ORIGINAL"), [])
            for i, j, m in pairs:
                o, r = O[i], R[j]
                above = [p for p in orig_all
                         if p["rank"] < o["rank"] and p["fusion_associated"]]
                rows.append({
                    "detector": det, "pdb_id": s, "partner": man[s]["fusion_partner"],
                    "topology": man[s]["topology"],
                    "orig_rank": o["rank"], "removed_rank": r["rank"],
                    "rank_change": o["rank"] - r["rank"],
                    "orig_score": o["score"], "removed_score": r["score"],
                    "score_change": r["score"] - o["score"],
                    "jaccard": round(m["jaccard"], 3), "overlap_coef": round(m["overlap"], 3),
                    "centroid_dist": round(m["dist"], 2), "n_shared_target_res": m["n_shared"],
                    "n_target_res_orig": m["n_orig"], "n_target_res_removed": m["n_rem"],
                    "orig_class": o["pocket_class"],
                    "fusion_assoc_above_in_original": len(above),
                    "orig_recovers_ref": o["recovers_reference_site"],
                    "removed_recovers_ref": r["recovers_reference_site"],
                })
            # per-structure state accounting
            per_struct.append({
                "detector": det, "pdb_id": s, "partner": man[s]["fusion_partner"],
                "n_orig_matchable": len(O), "n_removed_matchable": len(R),
                "n_matched": len(pairs),
                "n_unmatched_original": len(O) - len(matched_o),
                "n_new_after_removal": len(R) - len(matched_r),
                "n_removed_pockets_total": len(idx.get((det, s, "REMOVED"), [])),
            })

    cols = list(rows[0].keys())
    with open(os.path.join(RESULTS, "matched_target_pockets.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    for det in detectors:
        sub = [r for r in rows if r["detector"] == det]
        ps = [p for p in per_struct if p["detector"] == det]
        print(f"\n  {det}: {len(sub)} MATCHED pairs over "
              f"{len({r['pdb_id'] for r in sub})} structures")
        print(f"       UNMATCHED_ORIGINAL   {sum(p['n_unmatched_original'] for p in ps)}")
        print(f"       NEW_AFTER_REMOVAL    {sum(p['n_new_after_removal'] for p in ps)}")
        if sub:
            js = sorted(r["jaccard"] for r in sub)
            ds = sorted(r["centroid_dist"] for r in sub)
            print(f"       median Jaccard {js[len(js)//2]:.3f}, "
                  f"median centroid dist {ds[len(ds)//2]:.2f} A, "
                  f"median shared target residues "
                  f"{sorted(r['n_shared_target_res'] for r in sub)[len(sub)//2]}")

    # ---------------------------------------------------------------- validation examples
    print("\n" + "=" * 104)
    print("VALIDATION — highest-confidence matches, showing the shared target residues")
    print("=" * 104)
    for det in detectors:
        sub = sorted([r for r in rows if r["detector"] == det],
                     key=lambda r: (-r["jaccard"], r["centroid_dist"]))[:4]
        for r in sub:
            cmap = man[r["pdb_id"]]["auth_class_map"]
            o = next(p for p in idx[(det, r["pdb_id"], "ORIGINAL")] if p["rank"] == r["orig_rank"])
            rm = next(p for p in idx[(det, r["pdb_id"], "REMOVED")]
                      if p["rank"] == r["removed_rank"])
            to = sorted(target_set(o, cmap), key=lambda x: int(''.join(c for c in x if c.isdigit() or c=='-')))
            tr = sorted(target_set(rm, cmap), key=lambda x: int(''.join(c for c in x if c.isdigit() or c=='-')))
            print(f"\n  {det} {r['pdb_id']} rank {r['orig_rank']} -> {r['removed_rank']}  "
                  f"J={r['jaccard']:.3f} d={r['centroid_dist']:.2f}A")
            print(f"    ORIGINAL target residues ({len(to)}): {','.join(to[:18])}"
                  f"{' ...' if len(to) > 18 else ''}")
            print(f"    REMOVED  target residues ({len(tr)}): {','.join(tr[:18])}"
                  f"{' ...' if len(tr) > 18 else ''}")
            print(f"    shared {r['n_shared_target_res']}; "
                  f"only-in-original {sorted(set(to)-set(tr))[:6]}; "
                  f"only-in-removed {sorted(set(tr)-set(to))[:6]}")

    json.dump({"primary_rule": {"jaccard_min": PRIMARY_JACCARD,
                                "max_centroid_dist": PRIMARY_MAX_CENTROID_DIST,
                                "min_target_residues": MIN_TARGET_RESIDUES,
                                "assignment": "scipy linear_sum_assignment maximising total "
                                              "Jaccard on TARGET residues; rank/score unused"},
               "sensitivity_rules": SENSITIVITY,
               "matched": rows, "per_structure": per_struct},
              open(os.path.join(RESULTS, "correspondence_summary.json"), "w"), indent=1)
    print(f"\nwrote {len(rows)} matched pairs")


if __name__ == "__main__":
    main()
