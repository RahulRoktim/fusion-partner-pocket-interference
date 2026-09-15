#!/usr/bin/env python
"""
Confirmatory analysis — protocol v1.3. Aims A, B, C, E and all pre-specified sensitivities.

[F2] THE STRUCTURE IS THE INFERENTIAL UNIT. Every CI and test below uses n = structures, with
clustering on target accession. Cavity-level numbers appear only as labelled descriptives.

[F3] Aim B is RANK_DISPLACEMENT = ORIGINAL_RANK - REMOVED_RANK for the detector's highest-ranked
cavity in the FUSION-REMOVED structure, matched into ORIGINAL by the frozen correspondence rule.
DETECTOR_ABSTENTION and TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL carry no numerical rank.
"""
import collections, importlib.util, json, math, os, random, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
VARIANT = os.environ.get("FUSIONTAG_VARIANT", "primary")
RES = os.path.join(ROOT, "results", "confirmatory", VARIANT)

DETECTORS = ["p2rank", "fpocket"]
PARTNERS = ["BRIL", "T4L", "MBP"]
BOOT = 10000
BOOT_SEED = 20260915


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(ROOT, "scripts", f))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


corr = _load("corr", "08_pocket_correspondence.py")   # frozen correspondence algorithm


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def cluster_bootstrap(units, stat, reps=BOOT, seed=BOOT_SEED):
    """units: list of (cluster_key, value). Resamples CLUSTERS with replacement."""
    by = collections.defaultdict(list)
    for k, v in units:
        by[k].append(v)
    keys = sorted(by)
    if not keys:
        return (float("nan"), float("nan"))
    rng = random.Random(seed)
    out = []
    for _ in range(reps):
        vals = []
        for _ in range(len(keys)):
            vals += by[keys[rng.randrange(len(keys))]]
        try:
            out.append(stat(vals))
        except Exception:
            pass
    out.sort()
    if not out:
        return (float("nan"), float("nan"))
    return out[int(0.025 * len(out))], out[min(len(out) - 1, int(0.975 * len(out)))]


def mean(v):
    return sum(v) / len(v) if v else float("nan")


def median(v):
    s = sorted(v)
    n = len(s)
    if not n:
        return float("nan")
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def fmt(k, n, lo=None, hi=None):
    if n == 0:
        return "n/a"
    s = f"{k}/{n} = {100*k/n:5.1f}%"
    a, b = wilson(k, n)
    s += f"  Wilson[{100*a:.1f}-{100*b:.1f}]"
    if lo is not None and not math.isnan(lo):
        s += f"  boot[{100*lo:.1f}-{100*hi:.1f}]"
    return s


def main():
    frz = json.load(open(os.path.join(ROOT, "environment", "PROTOCOL_FREEZE.json")))
    man = {p["pdb_id"]: p for p in
           json.load(open(os.path.join(MANI,
                                       f"confirmatory_manifest_{VARIANT}.json")))["prepared"]
           if "PREPARATION_FAILED" not in p}
    P = json.load(open(os.path.join(RES, "pockets_classified.json")))["pockets"]
    runs = json.load(open(os.path.join(RES, "detector_runs.json")))["runs"]

    idx = collections.defaultdict(list)
    for r in P:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    for v in idx.values():
        v.sort(key=lambda r: r["rank"])

    structures = sorted(man)
    clus = {s: man[s]["target_accession"] for s in structures}
    n_clusters = len({clus[s] for s in structures})

    lines = []

    def emit(s=""):
        lines.append(s)
        print(s)

    emit("=" * 104)
    emit(f"CONFIRMATORY ANALYSIS — protocol v1.3, variant '{VARIANT}'")
    emit(f"freeze sha256(PROTOCOL.md) = {frz['sha256']['PROTOCOL.md'][:32]}...")
    emit(f"n_structures = {len(structures)}   n_independent_target_clusters = {n_clusters}")
    emit("[F2] every CI and test below uses n = STRUCTURES, clustered on target accession")
    emit("=" * 104)

    out = {"variant": VARIANT, "n_structures": len(structures), "n_clusters": n_clusters,
           "detectors": {}}

    # ---------------- technical accounting
    emit("\nTECHNICAL SUCCESS")
    for det in DETECTORS:
        for cond in ("ORIGINAL", "REMOVED"):
            sub = [r for r in runs if r["detector"] == det and r["condition"] == cond]
            ok = sum(1 for r in sub if r["success"])
            emit(f"  {det:8s} {cond:9s}: {ok}/{len(sub)}")
    fails = [r for r in runs if not r["success"]]
    emit(f"  technical failures: {len(fails)}" +
         (f" -> {[f'{f['detector']}:{f['key']}' for f in fails][:6]}" if fails else ""))

    # ---------------- AIM A : exposure
    emit("\n" + "=" * 104)
    emit("PRIMARY AIM A — FUSION-POCKET EXPOSURE (ORIGINAL deposited construct)")
    emit("=" * 104)
    for det in DETECTORS:
        D = out["detectors"].setdefault(det, {})
        ev = [s for s in structures if idx.get((det, s, "ORIGINAL"))]
        emit(f"\n{det.upper()}  (evaluable structures {len(ev)})")
        A = {}
        for depth, label in ((1, "rank-1"), (3, "top-3"), (5, "top-5")):
            hit = {s: any(p["fusion_associated"] for p in idx[(det, s, "ORIGINAL")][:depth])
                   for s in ev}
            k = sum(hit.values())
            lo, hi = cluster_bootstrap([(clus[s], 1 if hit[s] else 0) for s in ev], mean)
            emit(f"  {label:7s} : {fmt(k, len(ev), lo, hi)}")
            A[label] = {"k": k, "n": len(ev), "wilson": wilson(k, len(ev)),
                        "cluster_boot": (lo, hi),
                        "pdb_ids": [s for s in ev if hit[s]]}
        emit("  partner strata (reported before any pooled interpretation):")
        by_p = {}
        for pt in PARTNERS:
            sub = [s for s in ev if man[s]["fusion_partner"] == pt]
            k = sum(1 for s in sub if idx[(det, s, "ORIGINAL")][0]["fusion_associated"])
            lo, hi = cluster_bootstrap(
                [(clus[s], 1 if idx[(det, s, 'ORIGINAL')][0]["fusion_associated"] else 0)
                 for s in sub], mean)
            emit(f"    {pt:5s} rank-1 : {fmt(k, len(sub), lo, hi)}")
            by_p[pt] = {"k": k, "n": len(sub), "rate": k / len(sub) if sub else float("nan"),
                        "wilson": wilson(k, len(sub)), "cluster_boot": (lo, hi)}
            for depth, label in ((3, "top-3"), (5, "top-5")):
                kk = sum(1 for s in sub
                         if any(p["fusion_associated"] for p in idx[(det, s, "ORIGINAL")][:depth]))
                by_p[pt][label] = {"k": kk, "n": len(sub)}
                emit(f"          {label} : {fmt(kk, len(sub))}")
        rates = [by_p[pt]["rate"] for pt in PARTNERS if by_p[pt]["n"]]
        std = sum(rates) / len(rates) if rates else float("nan")
        pooled = A["rank-1"]["k"] / A["rank-1"]["n"]
        emit(f"  RAW POOLED rank-1          : {100*pooled:.1f}%")
        emit(f"  PARTNER-STANDARDISED (1/3) : {100*std:.1f}%   "
             f"(difference {100*(pooled-std):+.1f} points)")
        # standardised CI by clustered bootstrap of the stratified estimator
        def std_stat(vals):
            acc = collections.defaultdict(list)
            for pt, v in vals:
                acc[pt].append(v)
            rr = [mean(acc[pt]) for pt in PARTNERS if acc[pt]]
            return sum(rr) / len(rr) if rr else float("nan")
        slo, shi = cluster_bootstrap(
            [(clus[s], (man[s]["fusion_partner"],
                        1 if idx[(det, s, 'ORIGINAL')][0]["fusion_associated"] else 0))
             for s in ev], std_stat)
        emit(f"  standardised cluster-bootstrap 95% CI: [{100*slo:.1f}-{100*shi:.1f}]")
        cc = collections.Counter(idx[(det, s, "ORIGINAL")][0]["pocket_class"] for s in ev)
        emit(f"  rank-1 class breakdown: {dict(cc)}")
        D["A"] = {"depths": A, "by_partner": by_p, "pooled": pooled, "standardised": std,
                  "standardised_ci": (slo, shi), "rank1_classes": dict(cc)}

    # ---------------- AIM B : F3 displacement
    emit("\n" + "=" * 104)
    emit("PRIMARY AIM B — TARGET-CAVITY RANK DISPLACEMENT  [F3, structure-level]")
    emit("=" * 104)
    coords = {s: corr.load_target_coords(man[s]) for s in structures}
    disp_rows = []
    for det in DETECTORS:
        D = out["detectors"][det]
        recs = []
        for s in structures:
            cmap = man[s]["auth_class_map"]
            orig = idx.get((det, s, "ORIGINAL"), [])
            rem = idx.get((det, s, "REMOVED"), [])
            if not rem:
                recs.append({"pdb_id": s, "partner": man[s]["fusion_partner"],
                             "target": clus[s], "state": "DETECTOR_ABSTENTION",
                             "removed_rank": None, "original_rank": None,
                             "displacement": None})
                continue
            top = rem[0]
            if len(corr.target_set(top, cmap)) < corr.MIN_TARGET_RESIDUES:
                recs.append({"pdb_id": s, "partner": man[s]["fusion_partner"],
                             "target": clus[s],
                             "state": "TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL",
                             "removed_rank": top["rank"], "original_rank": None,
                             "displacement": None,
                             "note": "removed rank-1 cavity has <3 target residues"})
                continue
            cand, O, R = corr.build_pairs(orig, rem, coords[s], cmap)
            pairs = corr.assign(cand, O, R, corr.PRIMARY_JACCARD, corr.PRIMARY_MAX_CENTROID_DIST)
            jrem = next((j for j, r in enumerate(R) if r["rank"] == top["rank"]), None)
            hit = next(((i, j, m) for i, j, m in pairs if j == jrem), None)
            if hit is None:
                recs.append({"pdb_id": s, "partner": man[s]["fusion_partner"],
                             "target": clus[s],
                             "state": "TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL",
                             "removed_rank": top["rank"], "original_rank": None,
                             "displacement": None})
                continue
            i, j, m = hit
            orank = O[i]["rank"]
            above = sum(1 for p in orig if p["rank"] < orank and p["fusion_associated"])
            recs.append({"pdb_id": s, "partner": man[s]["fusion_partner"], "target": clus[s],
                         "state": "EVALUABLE", "removed_rank": top["rank"],
                         "original_rank": orank, "displacement": orank - top["rank"],
                         "jaccard": round(m["jaccard"], 3),
                         "centroid_dist": round(m["dist"], 2),
                         "fusion_pockets_above_in_original": above,
                         "orig_score": O[i]["score"], "removed_score": top["score"]})
        for r in recs:
            r["detector"] = det
        disp_rows += recs

        ev = [r for r in recs if r["state"] == "EVALUABLE"]
        ab = [r for r in recs if r["state"] == "DETECTOR_ABSTENTION"]
        um = [r for r in recs if r["state"] == "TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL"]
        d = [r["displacement"] for r in ev]
        emit(f"\n{det.upper()}")
        emit(f"  evaluable                            {len(ev)}")
        emit(f"  DETECTOR_ABSTENTION                  {len(ab)}"
             + (f"  {[r['pdb_id'] for r in ab][:8]}" if ab else ""))
        emit(f"  TARGET_CAVITY_NOT_RECOVERED_IN_ORIG  {len(um)}"
             + (f"  {[r['pdb_id'] for r in um][:8]}" if um else ""))
        if ev:
            lo_m, hi_m = cluster_bootstrap([(r["target"], r["displacement"]) for r in ev], median)
            lo_a, hi_a = cluster_bootstrap([(r["target"], r["displacement"]) for r in ev], mean)
            emit(f"  median displacement                  {median(d):g}   "
                 f"cluster-boot 95% CI [{lo_m:g}, {hi_m:g}]")
            emit(f"  mean displacement                    {mean(d):+.2f}  "
                 f"cluster-boot 95% CI [{lo_a:+.2f}, {hi_a:+.2f}]")
            emit(f"  range                                [{min(d)}, {max(d)}]")
            props = {}
            for thr, lab in ((1, ">0"), (2, ">=2"), (5, ">=5")):
                k = sum(1 for x in d if x >= thr)
                lo, hi = cluster_bootstrap(
                    [(r["target"], 1 if r["displacement"] >= thr else 0) for r in ev], mean)
                emit(f"  displacement {lab:4s}                    {fmt(k, len(ev), lo, hi)}")
                props[lab] = {"k": k, "n": len(ev), "wilson": wilson(k, len(ev)),
                              "cluster_boot": (lo, hi)}
            neg = sum(1 for x in d if x < 0)
            emit(f"  displacement < 0 (fusion HELPED)     {neg}/{len(ev)}")
            emit(f"  displacement == 0 (no displacement)  {sum(1 for x in d if x == 0)}/{len(ev)}")
            # structure-level sign test (exact binomial on non-zero structures)
            pos = sum(1 for x in d if x > 0)
            nz = pos + neg
            from math import comb
            pv = (min(1.0, 2 * sum(comb(nz, i) for i in range(0, min(pos, neg) + 1)) / 2 ** nz)
                  if nz else float("nan"))
            emit(f"  structure-level sign test on {nz} non-zero structures: "
                 f"{pos} up / {neg} down, exact p = {pv:.3g}")
            emit("  partner-stratified displacement:")
            bp = {}
            for pt in PARTNERS:
                sub = [r for r in ev if r["partner"] == pt]
                if not sub:
                    continue
                dd = [r["displacement"] for r in sub]
                kk = sum(1 for x in dd if x > 0)
                emit(f"    {pt:5s} n={len(sub):3d}  median {median(dd):g}  mean {mean(dd):+.2f}  "
                     f">0 in {fmt(kk, len(sub))}")
                bp[pt] = {"n": len(sub), "median": median(dd), "mean": mean(dd),
                          "k_gt0": kk, "wilson": wilson(kk, len(sub))}
            D["B"] = {"evaluable": len(ev), "abstention": len(ab), "unmatched": len(um),
                      "median": median(d), "mean": mean(d),
                      "median_ci": (lo_m, hi_m), "mean_ci": (lo_a, hi_a),
                      "props": props, "sign_test_p": pv, "by_partner": bp,
                      "n_negative": neg}

    cols = sorted({k for r in disp_rows for k in r})
    with open(os.path.join(RES, "displacement.tsv"), "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in disp_rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")

    # ---------------- AIM C : biological reference site
    emit("\n" + "=" * 104)
    emit("SECONDARY AIM C — BIOLOGICAL REFERENCE SITE")
    emit("=" * 104)
    refs = [s for s in structures if man[s]["reference_site_available"]]
    cov = len(refs) / len(structures)
    emit(f"  coverage overall: {len(refs)}/{len(structures)} = {100*cov:.1f}%  -> "
         f"{'>=60%: PROMOTED per the pre-frozen rule' if cov >= 0.60 else '<60%: remains a restricted-subset SECONDARY analysis'}")
    for pt in PARTNERS:
        sub = [s for s in structures if man[s]["fusion_partner"] == pt]
        k = sum(1 for s in sub if man[s]["reference_site_available"])
        emit(f"    {pt:5s}: {k}/{len(sub)}")
    out["reference_site_coverage"] = {"k": len(refs), "n": len(structures), "fraction": cov}
    for det in DETECTORS:
        emit(f"\n  {det.upper()}")
        rows = []
        for s in refs:
            cmap = man[s]["auth_class_map"]
            orig = idx.get((det, s, "ORIGINAL"), [])
            rem = idx.get((det, s, "REMOVED"), [])
            o = next((p for p in orig if p["recovers_reference_site"]), None)
            r = next((p for p in rem if p["recovers_reference_site"]), None)
            above = (sum(1 for p in orig if p["rank"] < o["rank"] and p["fusion_associated"])
                     if o else None)
            rows.append({"pdb_id": s, "partner": man[s]["fusion_partner"], "target": clus[s],
                         "orig_rank": o["rank"] if o else None,
                         "removed_rank": r["rank"] if r else None,
                         "rank_change": (o["rank"] - r["rank"]) if (o and r) else None,
                         "fusion_above": above,
                         "top1": bool(o and o["rank"] == 1),
                         "top3": bool(o and o["rank"] <= 3),
                         "top5": bool(o and o["rank"] <= 5)})
        found = [r for r in rows if r["orig_rank"] is not None]
        emit(f"    site-corresponding pocket found in ORIGINAL: {len(found)}/{len(refs)}")
        if found:
            ch = [r["rank_change"] for r in found if r["rank_change"] is not None]
            k = sum(1 for r in found if (r["fusion_above"] or 0) > 0)
            lo, hi = cluster_bootstrap(
                [(r["target"], 1 if (r["fusion_above"] or 0) > 0 else 0) for r in found], mean)
            emit(f"    outranked by a fusion-associated pocket: {fmt(k, len(found), lo, hi)}")
            emit(f"    rank change (orig->removed): improved {sum(1 for x in ch if x>0)}, "
                 f"worsened {sum(1 for x in ch if x<0)}, unchanged {sum(1 for x in ch if x==0)} "
                 f"(median {median(ch):g})")
            for lab in ("top1", "top3", "top5"):
                kk = sum(1 for r in found if r[lab])
                emit(f"    recovery in ORIGINAL {lab}: {fmt(kk, len(found))}")
            out["detectors"][det]["C"] = {
                "coverage": len(refs), "found": len(found), "outranked": k,
                "improved": sum(1 for x in ch if x > 0),
                "worsened": sum(1 for x in ch if x < 0),
                "unchanged": sum(1 for x in ch if x == 0),
                "median_change": median(ch),
                "top1": sum(1 for r in found if r["top1"]),
                "top3": sum(1 for r in found if r["top3"]),
                "top5": sum(1 for r in found if r["top5"])}
        with open(os.path.join(RES, f"reference_site_{det}.tsv"), "w", encoding="utf-8") as fh:
            cc = list(rows[0].keys()) if rows else []
            fh.write("\t".join(cc) + "\n")
            for r in rows:
                fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cc) + "\n")

    # ---------------- AIM E : detector dependence
    emit("\n" + "=" * 104)
    emit("SECONDARY AIM E — DETECTOR DEPENDENCE")
    emit("=" * 104)
    both = [s for s in structures
            if idx.get(("p2rank", s, "ORIGINAL")) and idx.get(("fpocket", s, "ORIGINAL"))]
    a = [idx[("p2rank", s, "ORIGINAL")][0]["fusion_associated"] for s in both]
    b = [idx[("fpocket", s, "ORIGINAL")][0]["fusion_associated"] for s in both]
    agree = sum(1 for x, y in zip(a, b) if x == y)
    po = agree / len(both)
    pa = (sum(a)/len(both))*(sum(b)/len(both)) + (1-sum(a)/len(both))*(1-sum(b)/len(both))
    kappa = (po - pa) / (1 - pa) if pa < 1 else float("nan")
    emit(f"  structures scored by both: {len(both)}")
    emit(f"  agreement on rank-1 fusion-associated: {agree}/{len(both)} = {100*po:.1f}%")
    emit(f"  Cohen's kappa: {kappa:.3f}")
    emit(f"  both {sum(1 for x,y in zip(a,b) if x and y)} | "
         f"P2Rank only {sum(1 for x,y in zip(a,b) if x and not y)} | "
         f"fpocket only {sum(1 for x,y in zip(a,b) if y and not x)} | "
         f"neither {sum(1 for x,y in zip(a,b) if not x and not y)}")
    emit("  score margin (best fusion-associated / best target-dominated), ORIGINAL:")
    for det in DETECTORS:
        for pt in PARTNERS + ["ALL"]:
            vals = []
            for s in structures:
                if pt != "ALL" and man[s]["fusion_partner"] != pt:
                    continue
                pk = idx.get((det, s, "ORIGINAL"), [])
                bf = max((p["score"] for p in pk if p["fusion_associated"]), default=None)
                bt = max((p["score"] for p in pk if p["pocket_class"] == "TARGET_DOMINATED"),
                         default=None)
                if bf is not None and bt:
                    vals.append(bf / bt)
            emit(f"    {det:8s} {pt:5s}: n={len(vals):3d} median "
                 f"{median(vals):.2f}" if vals else f"    {det:8s} {pt:5s}: n=0")
    emit("  abstention pattern (REMOVED condition, zero predictions):")
    for det in DETECTORS:
        ab = [s for s in structures if not idx.get((det, s, "REMOVED"))]
        bypt = collections.Counter(man[s]["fusion_partner"] for s in ab)
        emit(f"    {det:8s}: {len(ab)}/{len(structures)}  {dict(bypt)}")
    out["E"] = {"n_both": len(both), "agreement": po, "kappa": kappa}

    # ---------------- SENSITIVITIES (recomputable without re-running detectors)
    emit("\n" + "=" * 104)
    emit("PRE-SPECIFIED SENSITIVITY ANALYSES (this variant)")
    emit("=" * 104)
    emit("\n  Pocket dominance threshold — Aim A rank-1 exposure:")
    sens = {}
    for det in DETECTORS:
        ev = [s for s in structures if idx.get((det, s, "ORIGINAL"))]
        row = {}
        for t in ("0.50", "0.70", "0.90"):
            k = sum(1 for s in ev if idx[(det, s, "ORIGINAL")][0][f"class_at_{t}"]
                    in ("FUSION_DOMINATED", "INTERFACE", "LINKER"))
            row[t] = {"k": k, "n": len(ev)}
            emit(f"    {det:8s} dominance {t}: {fmt(k, len(ev))}")
        sens.setdefault("dominance", {})[det] = row

    emit("\n  Correspondence rule — Aim B displacement:")
    for det in DETECTORS:
        for lab, jm, dm in corr.SENSITIVITY:
            ev, dsp = 0, []
            for s in structures:
                cmap = man[s]["auth_class_map"]
                orig = idx.get((det, s, "ORIGINAL"), [])
                rem = idx.get((det, s, "REMOVED"), [])
                if not rem:
                    continue
                top = rem[0]
                if len(corr.target_set(top, cmap)) < corr.MIN_TARGET_RESIDUES:
                    continue
                cand, O, R = corr.build_pairs(orig, rem, coords[s], cmap)
                pairs = corr.assign(cand, O, R, jm, dm)
                jrem = next((j for j, r in enumerate(R) if r["rank"] == top["rank"]), None)
                hit = next(((i, j, m) for i, j, m in pairs if j == jrem), None)
                if hit:
                    ev += 1
                    dsp.append(O[hit[0]]["rank"] - top["rank"])
            k = sum(1 for x in dsp if x > 0)
            emit(f"    {det:8s} {lab:12s}: evaluable {ev:3d}  median {median(dsp):g}  "
                 f">0 {fmt(k, len(dsp))}")
            sens.setdefault("correspondence", {}).setdefault(det, {})[lab] = {
                "evaluable": ev, "median": median(dsp), "k_gt0": k, "n": len(dsp)}

    emit("\n  Deletion-boundary 8 A control (REMOVED-condition rank-1 cavities near a cut):")
    for det in DETECTORS:
        near = [s for s in structures
                if idx.get((det, s, "REMOVED"))
                and idx[(det, s, "REMOVED")][0]["dist_to_deletion_boundary"] is not None
                and idx[(det, s, "REMOVED")][0]["dist_to_deletion_boundary"] <= 8.0]
        ev = [r for r in disp_rows if r["detector"] == det and r["state"] == "EVALUABLE"]
        filt = [r for r in ev if r["pdb_id"] not in set(near)]
        k = sum(1 for r in filt if r["displacement"] > 0)
        emit(f"    {det:8s}: {len(near)} structures have a rank-1 REMOVED cavity within 8 A of a "
             f"cut; excluding them, displacement >0 in {fmt(k, len(filt))} "
             f"(median {median([r['displacement'] for r in filt]):g})")
        sens.setdefault("boundary8A", {})[det] = {"n_near": len(near), "k_gt0": k,
                                                  "n": len(filt)}
    out["sensitivity"] = sens

    json.dump(out, open(os.path.join(RES, "confirmatory_summary.json"), "w"),
              indent=1, default=str)
    open(os.path.join(RES, "confirmatory_summary.txt"), "w",
         encoding="utf-8").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
