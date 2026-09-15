#!/usr/bin/env python
"""
Phase 2E — pilot analysis. Results are produced independently for P2Rank and fpocket.

Nothing is filtered, hidden, or excluded here: every structure that produced a successful
detector run is reported, including those that contradict the hypothesis.

The pilot is a feasibility and signal-detection exercise (protocol v1.1 A4). No confirmatory
inference is drawn and no p-value is presented as evidence for the scientific hypothesis.

Outputs:
  results/pilot_summary.json
  results/pilot_summary.txt
  results/per_structure_rank1.tsv
"""
import json, math, os, collections

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
RESULTS = os.path.join(ROOT, "results")

DETECTORS = ["p2rank", "fpocket"]
PARTNERS = ["BRIL", "T4L", "MBP"]
BOUNDARY_CUTOFF = 8.0          # protocol v1.1 A3, fixed before results
FUSION_CLASSES = {"FUSION_DOMINATED", "INTERFACE", "LINKER"}


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def fmt_rate(k, n):
    if n == 0:
        return "n/a"
    lo, hi = wilson(k, n)
    return f"{k}/{n} = {100*k/n:4.1f}%  [95% CI {100*lo:.1f}-{100*hi:.1f}]"


def exact_mcnemar(b, c):
    """Two-sided exact binomial test on discordant pairs. Descriptive at pilot scale."""
    n = b + c
    if n == 0:
        return 1.0
    from math import comb
    obs = min(b, c)
    tail = sum(comb(n, i) for i in range(0, obs + 1)) / (2 ** n)
    return min(1.0, 2 * tail)


def main():
    P = json.load(open(os.path.join(RESULTS, "pockets_classified.json")))["pockets"]
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI, "pilot_manifest.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    runs = json.load(open(os.path.join(RESULTS, "detector_runs.json")))["runs"]

    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    for v in idx.values():
        v.sort(key=lambda r: r["rank"])

    structures = sorted(man)
    out = {"protocol_version": "1.1", "n_structures": len(structures),
           "boundary_cutoff_A": BOUNDARY_CUTOFF, "detectors": {}}
    lines = []

    def emit(s=""):
        lines.append(s)
        print(s)

    emit("=" * 96)
    emit("PILOT ANALYSIS — protocol v1.1 — signal detection only, no confirmatory inference")
    emit("=" * 96)

    # ---------- run accounting
    emit("\nDETECTOR RUN ACCOUNTING")
    for det in DETECTORS:
        sub = [r for r in runs if r["detector"] == det]
        ok = [r for r in sub if r["success"]]
        emit(f"  {det:8s}: {len(ok)}/{len(sub)} successful runs")
        for f in [r for r in sub if not r["success"]]:
            emit(f"      FAILURE {f['key']} rc={f['returncode']}")

    # ---------- reference-site coverage (A2)
    cov = [s for s in structures if man[s]["reference_site_available"]]
    emit(f"\nTARGET REFERENCE-SITE COVERAGE: {len(cov)}/{len(structures)} = "
         f"{100*len(cov)/max(len(structures),1):.0f}%  "
         f"(A2 threshold 60% -> "
         f"{'reference-site endpoint ELIGIBLE for primary' if len(cov)/max(len(structures),1) >= 0.6 else 'fusion-capture endpoint REMAINS primary'})")
    out["reference_site_coverage"] = {"n_with": len(cov), "n_total": len(structures),
                                      "fraction": len(cov) / max(len(structures), 1),
                                      "pdb_ids": cov}

    per_struct_rows = []

    for det in DETECTORS:
        emit("\n" + "=" * 96)
        emit(f"DETECTOR: {det.upper()}")
        emit("=" * 96)
        D = {}

        evaluable = [s for s in structures if idx.get((det, s, "ORIGINAL"))]
        n = len(evaluable)

        # ---- 1-3: rank-1 / top-3 / top-5 fusion-associated (ORIGINAL condition)
        r1 = [s for s in evaluable if idx[(det, s, "ORIGINAL")][0]["fusion_associated"]]
        t3 = [s for s in evaluable
              if any(p["fusion_associated"] for p in idx[(det, s, "ORIGINAL")][:3])]
        t5 = [s for s in evaluable
              if any(p["fusion_associated"] for p in idx[(det, s, "ORIGINAL")][:5])]
        emit(f"\n  rank-1 fusion-associated : {fmt_rate(len(r1), n)}")
        emit(f"  top-3  fusion-associated : {fmt_rate(len(t3), n)}")
        emit(f"  top-5  fusion-associated : {fmt_rate(len(t5), n)}")
        D["rank1"] = {"k": len(r1), "n": n, "ci": wilson(len(r1), n), "pdb_ids": r1}
        D["top3"] = {"k": len(t3), "n": n, "ci": wilson(len(t3), n), "pdb_ids": t3}
        D["top5"] = {"k": len(t5), "n": n, "ci": wilson(len(t5), n), "pdb_ids": t5}

        # rank-1 class breakdown
        cls_counts = collections.Counter(
            idx[(det, s, "ORIGINAL")][0]["pocket_class"] for s in evaluable)
        emit(f"  rank-1 class breakdown   : {dict(cls_counts)}")
        D["rank1_class_counts"] = dict(cls_counts)

        # ---- 4: partner-stratified rank-1
        emit("\n  partner-stratified rank-1 fusion-associated:")
        D["by_partner"] = {}
        for pt in PARTNERS:
            sub = [s for s in evaluable if man[s]["fusion_partner"] == pt]
            k = sum(1 for s in sub if idx[(det, s, "ORIGINAL")][0]["fusion_associated"])
            emit(f"    {pt:5s}: {fmt_rate(k, len(sub))}")
            D["by_partner"][pt] = {"k": k, "n": len(sub), "ci": wilson(k, len(sub)),
                                   "pdb_ids": [s for s in sub
                                               if idx[(det, s, 'ORIGINAL')][0]['fusion_associated']]}

        # ---- 5: median rank of first fusion-associated pocket
        first_ranks = []
        for s in evaluable:
            hit = next((p["rank"] for p in idx[(det, s, "ORIGINAL")] if p["fusion_associated"]),
                       None)
            first_ranks.append(hit)
        present = sorted(r for r in first_ranks if r is not None)
        med = (present[len(present)//2] if len(present) % 2 else
               (present[len(present)//2 - 1] + present[len(present)//2]) / 2) if present else None
        emit(f"\n  first fusion-associated pocket present in {len(present)}/{n} structures; "
             f"median rank among those = {med}")
        D["median_rank_first_fusion"] = med
        D["n_with_any_fusion_pocket"] = len(present)

        # ---- 6: score margin
        margins = []
        for s in evaluable:
            pk = idx[(det, s, "ORIGINAL")]
            bf = max((p["score"] for p in pk if p["fusion_associated"]), default=None)
            bt = max((p["score"] for p in pk if p["pocket_class"] == "TARGET_DOMINATED"),
                     default=None)
            if bf is not None and bt not in (None, 0):
                margins.append(bf / bt)
        margins.sort()
        mm = (margins[len(margins)//2] if margins else None)
        emit(f"  score margin (best fusion-assoc / best target-dominated): "
             f"n={len(margins)} median={mm:.2f}" if margins else
             "  score margin: not computable")
        D["score_margin_median"] = mm
        D["score_margin_n"] = len(margins)
        D["score_margin_values"] = [round(x, 3) for x in margins]

        # ---- 9: paired original vs removed
        # [C1] The fusion-association version of this contrast is STRUCTURALLY DEGENERATE:
        # the REMOVED condition contains only TARGET residues, so no REMOVED pocket can ever be
        # fusion-associated and the "not-fusion -> fusion" cell is fixed at 0 by construction.
        # The meaningful outcome variable is "rank-1 pocket is TARGET_DOMINATED" (the literal
        # reading of "flips the rank-1 pocket to the target site"). Both are reported.
        # [C2] A structure with NO predicted pocket in REMOVED is NOT a rescue; it is recorded
        # as NO_POCKET_PREDICTED and counted as not-TARGET_DOMINATED (conservative).
        emit("\n  PAIRED REMOVAL (original -> fusion-removed):")
        paired = list(evaluable)                       # no structure is dropped
        b = c = 0                       # on TARGET_DOMINATED (non-degenerate, [C1])
        b_degen = c_degen = 0           # on fusion-association (degenerate, reported for audit)
        n_no_pocket = 0
        rank_changes, rescue_all, rescue_filtered = [], [], []
        ref_b = ref_c = 0
        for s in paired:
            o1 = idx[(det, s, "ORIGINAL")][0]
            rem = idx.get((det, s, "REMOVED"))
            r1p = rem[0] if rem else None
            if r1p is None:
                n_no_pocket += 1
            o_td = o1["pocket_class"] == "TARGET_DOMINATED"
            r_td = (r1p is not None and r1p["pocket_class"] == "TARGET_DOMINATED")
            if not o_td and r_td:
                b += 1
            elif o_td and not r_td:
                c += 1
            if o1["fusion_associated"] and (r1p is None or not r1p["fusion_associated"]):
                b_degen += 1
            elif not o1["fusion_associated"] and r1p is not None and r1p["fusion_associated"]:
                c_degen += 1
            if r1p is None:
                continue
            # best target-dominated rank in each condition
            def best_target_rank(cond):
                return next((p["rank"] for p in idx[(det, s, cond)]
                             if p["pocket_class"] == "TARGET_DOMINATED"), None)
            ro, rr = best_target_rank("ORIGINAL"), best_target_rank("REMOVED")
            if ro is not None and rr is not None:
                rank_changes.append((s, ro, rr, ro - rr))
            # reference-site recovery at rank 1
            if man[s]["reference_site_available"]:
                oref = o1["recovers_reference_site"]
                rref = r1p["recovers_reference_site"]
                if oref and not rref:
                    ref_b += 1
                elif rref and not oref:
                    ref_c += 1
            # [A3] rescue with and without the 8 A boundary filter
            if not o_td and r_td:
                rescue_all.append(s)
                d = r1p["dist_to_deletion_boundary"]
                if d is None or d > BOUNDARY_CUTOFF:
                    rescue_filtered.append(s)

        emit(f"    pairs evaluable: {len(paired)} "
             f"(of which {n_no_pocket} produced NO pocket at all after removal)")
        emit(f"    [C1] PRIMARY (non-degenerate), outcome = 'rank-1 is TARGET_DOMINATED':")
        emit(f"      not-target -> target  (improved): {b}")
        emit(f"      target -> not-target  (worsened): {c}")
        emit(f"      exact McNemar on discordant pairs (DESCRIPTIVE ONLY): "
             f"p = {exact_mcnemar(b, c):.4f}")
        emit(f"    [C1] degenerate version, outcome = 'rank-1 is fusion-associated':")
        emit(f"      fusion-assoc -> not: {b_degen};  not -> fusion-assoc: {c_degen} "
             f"(structurally fixed at 0 — no fusion residues exist in REMOVED)")
        improved = [x for x in rank_changes if x[3] > 0]
        worse = [x for x in rank_changes if x[3] < 0]
        emit(f"    best target-dominated pocket rank: improved in {len(improved)}, "
             f"worsened in {len(worse)}, unchanged in "
             f"{len(rank_changes)-len(improved)-len(worse)} of {len(rank_changes)}")
        emit(f"    reference-site recovery at rank1: lost {ref_b}, gained {ref_c} "
             f"(of {sum(1 for s in paired if man[s]['reference_site_available'])} with a site)")
        D["paired"] = {"n_pairs": len(paired), "b_improved": b, "c_worsened": c,
                       "mcnemar_p_descriptive": exact_mcnemar(b, c),
                       "n_no_pocket_after_removal": n_no_pocket,
                       "degenerate_b": b_degen, "degenerate_c": c_degen,
                       "rank_changes": rank_changes,
                       "ref_site_lost": ref_b, "ref_site_gained": ref_c}

        # ---- 10: deletion-boundary sensitivity
        emit(f"\n  [A3] DELETION-BOUNDARY SENSITIVITY (cutoff {BOUNDARY_CUTOFF} A):")
        emit(f"    apparent rank-1 rescues, all predictions            : {len(rescue_all)}")
        emit(f"    rank-1 rescues surviving the {BOUNDARY_CUTOFF} A filter        : "
             f"{len(rescue_filtered)}")
        near = [(s, idx[(det, s, 'REMOVED')][0]['dist_to_deletion_boundary'])
                for s in rescue_all if s not in rescue_filtered]
        if near:
            emit(f"    rescues discounted as cut-adjacent: {near}")
        newpk = [p for s in paired for p in idx[(det, s, "REMOVED")]
                 if p["dist_to_deletion_boundary"] is not None
                 and p["dist_to_deletion_boundary"] <= BOUNDARY_CUTOFF]
        emit(f"    REMOVED-condition pockets within {BOUNDARY_CUTOFF} A of a cut: "
             f"{len(newpk)} of {sum(len(idx[(det,s,'REMOVED')]) for s in paired)}")
        D["boundary_sensitivity"] = {"rescues_all": rescue_all,
                                     "rescues_after_filter": rescue_filtered,
                                     "discounted": near,
                                     "n_pockets_near_cut": len(newpk)}

        # ---- 11: classification threshold sensitivity
        emit("\n  CLASSIFICATION SENSITIVITY (rank-1 fusion-associated rate):")
        D["threshold_sensitivity"] = {}
        for t in ("0.50", "0.70", "0.90"):
            k = sum(1 for s in evaluable
                    if idx[(det, s, "ORIGINAL")][0][f"class_at_{t}"] in FUSION_CLASSES)
            emit(f"    dominance {t}: {fmt_rate(k, n)}")
            D["threshold_sensitivity"][t] = {"k": k, "n": n}

        # per-structure rows
        for s in evaluable:
            o1 = idx[(det, s, "ORIGINAL")][0]
            r1p = idx[(det, s, "REMOVED")][0] if idx.get((det, s, "REMOVED")) else None
            per_struct_rows.append({
                "detector": det, "pdb_id": s, "partner": man[s]["fusion_partner"],
                "target": man[s]["target_accession"], "topology": man[s]["topology"],
                "resolution": man[s]["resolution"],
                "orig_rank1_class": o1["pocket_class"],
                "orig_rank1_score": o1["score"],
                "orig_rank1_f_target": o1["f_target"], "orig_rank1_f_fusion": o1["f_fusion"],
                "orig_rank1_f_linker": o1["f_linker"],
                "orig_rank1_fusion_assoc": o1["fusion_associated"],
                "orig_rank1_recovers_ref": o1["recovers_reference_site"],
                "removed_rank1_class": r1p["pocket_class"] if r1p else None,
                "removed_rank1_fusion_assoc": r1p["fusion_associated"] if r1p else None,
                "removed_rank1_recovers_ref": r1p["recovers_reference_site"] if r1p else None,
                "removed_rank1_dist_to_cut": r1p["dist_to_deletion_boundary"] if r1p else None,
                "ref_site_available": man[s]["reference_site_available"],
                "n_pockets_original": len(idx[(det, s, "ORIGINAL")]),
            })

        out["detectors"][det] = D

    # ---- 7: between-detector agreement on rank-1 class
    emit("\n" + "=" * 96)
    emit("BETWEEN-DETECTOR AGREEMENT (rank-1, ORIGINAL condition)")
    both = [s for s in structures
            if idx.get(("p2rank", s, "ORIGINAL")) and idx.get(("fpocket", s, "ORIGINAL"))]
    a = [idx[("p2rank", s, "ORIGINAL")][0]["fusion_associated"] for s in both]
    bl = [idx[("fpocket", s, "ORIGINAL")][0]["fusion_associated"] for s in both]
    agree = sum(1 for x, y in zip(a, bl) if x == y)
    po = agree / len(both) if both else float("nan")
    pa = ((sum(a) / len(both)) * (sum(bl) / len(both)) +
          (1 - sum(a) / len(both)) * (1 - sum(bl) / len(both))) if both else float("nan")
    kappa = (po - pa) / (1 - pa) if both and pa < 1 else float("nan")
    emit(f"  structures scored by both: {len(both)}")
    emit(f"  agreement on 'rank-1 is fusion-associated': {agree}/{len(both)} = {100*po:.0f}%")
    emit(f"  Cohen's kappa: {kappa:.3f}")
    emit(f"  both call fusion-associated: {sum(1 for x,y in zip(a,bl) if x and y)}; "
         f"P2Rank only: {sum(1 for x,y in zip(a,bl) if x and not y)}; "
         f"fpocket only: {sum(1 for x,y in zip(a,bl) if y and not x)}; "
         f"neither: {sum(1 for x,y in zip(a,bl) if not x and not y)}")
    out["detector_agreement"] = {"n": len(both), "agreement": po, "kappa": kappa,
                                 "both": sum(1 for x, y in zip(a, bl) if x and y),
                                 "p2rank_only": sum(1 for x, y in zip(a, bl) if x and not y),
                                 "fpocket_only": sum(1 for x, y in zip(a, bl) if y and not x),
                                 "neither": sum(1 for x, y in zip(a, bl) if not x and not y)}

    # ---- gates
    emit("\n" + "=" * 96)
    emit("GO / NO-GO GATES (protocol v1.1 section 10; thresholds are decision thresholds,")
    emit("not statistical rejection regions — see amendment A4)")
    emit("=" * 96)
    k_p2 = out["detectors"]["p2rank"]["rank1"]["k"]
    k_fp = out["detectors"]["fpocket"]["rank1"]["k"]
    g1 = (k_p2 >= 6 and k_fp >= 4) or (k_fp >= 6 and k_p2 >= 4)
    emit(f"  G1  >=6/24 in one detector AND >=4/24 in the other: "
         f"P2Rank={k_p2}, fpocket={k_fp}  ->  {'PASS' if g1 else 'FAIL'}")
    g2_hits = []
    for pt in PARTNERS:
        a2 = out["detectors"]["p2rank"]["by_partner"][pt]
        b2 = out["detectors"]["fpocket"]["by_partner"][pt]
        if a2["n"] and b2["n"] and a2["k"] / a2["n"] >= 0.5 and b2["k"] / b2["n"] >= 0.5:
            g2_hits.append(pt)
    emit(f"  G2  any partner >=50% in BOTH detectors: {g2_hits or 'none'}  ->  "
         f"{'PASS' if g2_hits else 'FAIL'}")
    g3_det = {}
    for det in DETECTORS:
        p = out["detectors"][det]["paired"]
        g3_det[det] = (p["b_improved"] >= 5 and p["c_worsened"] <= 1)
    g3 = any(g3_det.values())
    emit(f"  G3  paired removal makes rank-1 TARGET_DOMINATED in >=5/24 with <=1 reverse "
         f"[C1 non-degenerate reading]: "
         f"{ {d: (out['detectors'][d]['paired']['b_improved'], out['detectors'][d]['paired']['c_worsened']) for d in DETECTORS} }"
         f"  ->  {'PASS' if g3 else 'FAIL'}")
    overall = g1 or bool(g2_hits) or g3
    emit(f"\n  OVERALL: {'GO' if overall else 'NO-GO'}")
    out["gates"] = {"G1": g1, "G2": bool(g2_hits), "G2_partners": g2_hits,
                    "G3": g3, "G3_by_detector": g3_det, "overall_GO": overall}

    json.dump(out, open(os.path.join(RESULTS, "pilot_summary.json"), "w"), indent=1,
              default=str)
    open(os.path.join(RESULTS, "pilot_summary.txt"), "w", encoding="utf-8").write(
        "\n".join(lines) + "\n")

    cols = list(per_struct_rows[0].keys()) if per_struct_rows else []
    with open(os.path.join(RESULTS, "per_structure_rank1.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in per_struct_rows:
            fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
