#!/usr/bin/env python
"""
Phase 3A.1b — the non-degenerate paired endpoint, computed on DEVELOPMENT data.

Uses the frozen correspondence rule from 08 to follow the SAME target cavity across conditions
and asks what fusion removal did to its rank.

Endpoint definitions evaluated here (candidates for the confirmatory freeze):

  E-B1  DISPLACEMENT (per structure, binary).
        Take the matched target cavity with the best rank in the REMOVED condition ("the target's
        best cavity"). Was it outranked by >=1 fusion-associated pocket in the ORIGINAL condition?
        Deterministic; needs no score; non-degenerate.

  E-B2  PAIRED RANK CHANGE (per matched cavity, continuous).
        orig_rank - removed_rank over all matched pairs; Wilcoxon signed-rank.

  E-C   BIOLOGICAL REFERENCE SITE (restricted subset).
        Same two quantities restricted to the matched cavity that overlaps the outcome-independent
        biological reference site.
"""
import collections, json, math, os
from scipy.stats import wilcoxon

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
RESULTS = os.path.join(ROOT, "results")
DETECTORS = ["p2rank", "fpocket"]
PARTNERS = ["BRIL", "T4L", "MBP"]


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"),) * 2
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def rate(k, n):
    if n == 0:
        return "n/a"
    lo, hi = wilson(k, n)
    return f"{k}/{n} = {100*k/n:5.1f}% [{100*lo:.1f}-{100*hi:.1f}]"


def main():
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI, "pilot_manifest.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    corr = json.load(open(os.path.join(RESULTS, "correspondence_summary.json")))
    matched = corr["matched"]
    P = json.load(open(os.path.join(RESULTS, "pockets_classified.json")))["pockets"]
    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    for v in idx.values():
        v.sort(key=lambda r: r["rank"])

    by = collections.defaultdict(list)
    for m in matched:
        by[(m["detector"], m["pdb_id"])].append(m)

    out, lines = {}, []

    def emit(s=""):
        lines.append(s); print(s)

    emit("=" * 100)
    emit("PAIRED TARGET-CAVITY ENDPOINT — DEVELOPMENT DATA (n=24), exploratory")
    emit("=" * 100)

    for det in DETECTORS:
        emit(f"\n{'='*100}\nDETECTOR: {det.upper()}\n{'='*100}")
        D = {}

        evaluable = [s for s in sorted(man) if by.get((det, s))]
        emit(f"\n  structures with >=1 matched target cavity: {len(evaluable)}/24")
        no_match = [s for s in sorted(man) if not by.get((det, s))]
        if no_match:
            emit(f"  structures with NO matched target cavity: {no_match}")

        # ---------------- E-B1 displacement
        disp, rows = [], []
        for s in evaluable:
            best = min(by[(det, s)], key=lambda m: m["removed_rank"])
            displaced = best["fusion_assoc_above_in_original"] > 0
            disp.append((s, displaced, best))
            rows.append({"detector": det, "pdb_id": s, "partner": man[s]["fusion_partner"],
                         "orig_rank": best["orig_rank"], "removed_rank": best["removed_rank"],
                         "rank_change": best["rank_change"],
                         "n_fusion_above": best["fusion_assoc_above_in_original"],
                         "displaced": displaced, "jaccard": best["jaccard"],
                         "orig_score": best["orig_score"],
                         "removed_score": best["removed_score"]})
        k = sum(1 for _, d, _ in disp if d)
        emit(f"\n  [E-B1] DISPLACEMENT — the target's best cavity was outranked by a")
        emit(f"         fusion-associated pocket in the deposited construct:")
        emit(f"         {rate(k, len(disp))}")
        D["displacement"] = {"k": k, "n": len(disp), "ci": wilson(k, len(disp)),
                             "pdb_ids": [s for s, d, _ in disp if d]}
        emit("         by partner:")
        D["displacement_by_partner"] = {}
        for pt in PARTNERS:
            sub = [(s, d) for s, d, _ in disp if man[s]["fusion_partner"] == pt]
            kk = sum(1 for _, d in sub if d)
            emit(f"           {pt:5s}: {rate(kk, len(sub))}")
            D["displacement_by_partner"][pt] = {"k": kk, "n": len(sub)}

        # ---------------- E-B2 paired rank change
        allm = [m for m in matched if m["detector"] == det]
        deltas = [m["rank_change"] for m in allm]
        imp = sum(1 for d in deltas if d > 0)
        wor = sum(1 for d in deltas if d < 0)
        same = sum(1 for d in deltas if d == 0)
        emit(f"\n  [E-B2] PAIRED RANK CHANGE over all {len(allm)} matched cavities:")
        emit(f"         improved {imp}, worsened {wor}, unchanged {same}")
        sd = sorted(deltas)
        emit(f"         median change {sd[len(sd)//2]}, mean {sum(deltas)/len(deltas):+.2f}, "
             f"range [{min(deltas)}, {max(deltas)}]")
        try:
            w = wilcoxon([d for d in deltas if d != 0], alternative="two-sided")
            emit(f"         Wilcoxon signed-rank on non-zero changes (DESCRIPTIVE): "
                 f"p = {w.pvalue:.3g}")
            wp = float(w.pvalue)
        except Exception as exc:
            wp = None
            emit(f"         Wilcoxon not computable: {exc}")
        D["rank_change"] = {"n": len(allm), "improved": imp, "worsened": wor, "unchanged": same,
                            "median": sd[len(sd)//2], "mean": sum(deltas)/len(deltas),
                            "wilcoxon_p_descriptive": wp}

        # restricted to the displaced cavities
        dcav = [m for m in allm if m["fusion_assoc_above_in_original"] > 0]
        if dcav:
            dd = [m["rank_change"] for m in dcav]
            emit(f"         restricted to cavities with a fusion pocket above them "
                 f"(n={len(dcav)}): median change {sorted(dd)[len(dd)//2]}, "
                 f"improved {sum(1 for x in dd if x>0)}, worsened {sum(1 for x in dd if x<0)}")
            D["rank_change_displaced"] = {"n": len(dcav),
                                          "median": sorted(dd)[len(dd)//2],
                                          "improved": sum(1 for x in dd if x > 0),
                                          "worsened": sum(1 for x in dd if x < 0)}

        # ---------------- E-C biological reference site
        emit(f"\n  [E-C] BIOLOGICAL REFERENCE-SITE SUBSET")
        refs = [s for s in sorted(man) if man[s]["reference_site_available"]]
        emit(f"        coverage: {len(refs)}/24 structures")
        found = []
        for s in refs:
            ms = by.get((det, s), [])
            cand = [m for m in ms if m["orig_recovers_ref"] or m["removed_recovers_ref"]]
            if not cand:
                found.append((s, None))
                continue
            best = min(cand, key=lambda m: m["removed_rank"])
            found.append((s, best))
        have = [(s, m) for s, m in found if m]
        emit(f"        matched cavity overlapping the reference site found in "
             f"{len(have)}/{len(refs)}")
        if have:
            kk = sum(1 for _, m in have if m["fusion_assoc_above_in_original"] > 0)
            emit(f"        biological site outranked by a fusion pocket in the construct: "
                 f"{rate(kk, len(have))}")
            dl = [m["rank_change"] for _, m in have]
            emit(f"        paired rank change of the biological site: "
                 f"median {sorted(dl)[len(dl)//2]}, improved {sum(1 for x in dl if x>0)}, "
                 f"worsened {sum(1 for x in dl if x<0)}, unchanged {sum(1 for x in dl if x==0)}")
            for s, m in have:
                emit(f"          {s:6s} {man[s]['fusion_partner']:5s} "
                     f"rank {m['orig_rank']} -> {m['removed_rank']} "
                     f"(fusion pockets above in original: {m['fusion_assoc_above_in_original']}) "
                     f"J={m['jaccard']}")
            D["reference_site"] = {"coverage": len(refs), "matched": len(have),
                                   "displaced": kk,
                                   "rank_changes": [(s, m["orig_rank"], m["removed_rank"])
                                                    for s, m in have]}
        missing = [s for s, m in found if not m]
        if missing:
            emit(f"        no matched cavity overlapping the reference site: {missing}")

        out[det] = D
        if det == DETECTORS[0]:
            allrows = rows
        else:
            allrows += rows

    # ---------------- comparison with the degenerate pilot endpoint
    emit("\n" + "=" * 100)
    emit("COMPARISON WITH THE DEGENERATE PILOT ENDPOINT (audit only)")
    emit("=" * 100)
    pilot = json.load(open(os.path.join(RESULTS, "pilot_summary.json")))
    for det in DETECTORS:
        p = pilot["detectors"][det]["paired"]
        emit(f"  {det:8s} degenerate 'rank-1 becomes TARGET_DOMINATED': "
             f"{p['b_improved']} improved / {p['c_worsened']} worsened")
        emit(f"  {det:8s} non-degenerate displacement of the target's own best cavity: "
             f"{out[det]['displacement']['k']}/{out[det]['displacement']['n']}")

    cols = list(allrows[0].keys())
    with open(os.path.join(RESULTS, "paired_endpoint_per_structure.tsv"), "w",
              encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in allrows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    json.dump(out, open(os.path.join(RESULTS, "paired_endpoint_summary.json"), "w"),
              indent=1, default=str)
    open(os.path.join(RESULTS, "paired_endpoint_summary.txt"), "w",
         encoding="utf-8").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
