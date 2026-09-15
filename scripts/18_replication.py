#!/usr/bin/env python
"""
DEVELOPMENT vs CONFIRMATORY replication, protocol v1.3 / amendment F6.

The development cohort is re-described with the FROZEN CONFIRMATORY endpoints (in particular the
F3 displacement endpoint) so that the two cohorts are compared like for like. Re-describing
development data with the confirmatory endpoint changes no confirmatory rule.

Cohorts are NEVER pooled for confirmatory inference.

Criteria (frozen in DEVIATIONS F6, before any confirmatory result was inspected):
  REPLICATED           confirmatory point estimate inside the development 95% Wilson CI,
                       same direction  (+ for displacement: confirmatory boot CI excludes 0)
  PARTIALLY REPLICATED CIs overlap but the confirmatory estimate lies outside the development CI
  NOT REPLICATED       CIs disjoint, or direction differs, or confirmatory CI consistent with no
                       effect where development indicated one
"""
import collections, importlib.util, json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANI = os.path.join(ROOT, "data_manifest")
RES = os.path.join(ROOT, "results")
CONF = os.path.join(RES, "confirmatory", "primary")
DETECTORS = ["p2rank", "fpocket"]
PARTNERS = ["BRIL", "T4L", "MBP"]


def _load(n, f):
    s = importlib.util.spec_from_file_location(n, os.path.join(ROOT, "scripts", f))
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


corr = _load("corr", "08_pocket_correspondence.py")
ana = _load("ana", "16_confirmatory_analysis.py")


def wilson(k, n, z=1.96):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def verdict(dev_k, dev_n, con_k, con_n, extra_ok=True):
    if dev_n == 0 or con_n == 0:
        return "NOT EVALUABLE"
    dlo, dhi = wilson(dev_k, dev_n)
    clo, chi = wilson(con_k, con_n)
    p = con_k / con_n
    if dlo <= p <= dhi and extra_ok:
        return "REPLICATED"
    if not (chi < dlo or clo > dhi):
        return "PARTIALLY REPLICATED"
    return "NOT REPLICATED"


def f3_displacement(man, idx, det, structures):
    """Frozen F3 endpoint, applied to any cohort."""
    coords = {s: corr.load_target_coords(man[s]) for s in structures}
    out = []
    for s in structures:
        cmap = man[s]["auth_class_map"]
        orig = idx.get((det, s, "ORIGINAL"), [])
        rem = idx.get((det, s, "REMOVED"), [])
        if not rem:
            out.append({"pdb_id": s, "partner": man[s]["fusion_partner"],
                        "state": "DETECTOR_ABSTENTION", "displacement": None})
            continue
        top = rem[0]
        if len(corr.target_set(top, cmap)) < corr.MIN_TARGET_RESIDUES:
            out.append({"pdb_id": s, "partner": man[s]["fusion_partner"],
                        "state": "TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL",
                        "displacement": None})
            continue
        cand, O, R = corr.build_pairs(orig, rem, coords[s], cmap)
        pairs = corr.assign(cand, O, R, corr.PRIMARY_JACCARD, corr.PRIMARY_MAX_CENTROID_DIST)
        jrem = next((j for j, r in enumerate(R) if r["rank"] == top["rank"]), None)
        hit = next(((i, j, m) for i, j, m in pairs if j == jrem), None)
        if hit is None:
            out.append({"pdb_id": s, "partner": man[s]["fusion_partner"],
                        "state": "TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL",
                        "displacement": None})
            continue
        i, j, m = hit
        out.append({"pdb_id": s, "partner": man[s]["fusion_partner"], "state": "EVALUABLE",
                    "removed_rank": top["rank"], "original_rank": O[i]["rank"],
                    "displacement": O[i]["rank"] - top["rank"]})
    return out


def load(manifest, pockets):
    man = {p["pdb_id"]: p for p in json.load(open(manifest))["prepared"]
           if "PREPARATION_FAILED" not in p}
    idx = collections.defaultdict(list)
    for r in json.load(open(pockets))["pockets"]:
        idx[(r["detector"], r["pdb_id"], r["condition"])].append(r)
    for v in idx.values():
        v.sort(key=lambda r: r["rank"])
    return man, idx


def main():
    dman, didx = load(os.path.join(MANI, "pilot_manifest.json"),
                      os.path.join(RES, "pockets_classified.json"))
    cman, cidx = load(os.path.join(MANI, "confirmatory_manifest_primary.json"),
                      os.path.join(CONF, "pockets_classified.json"))
    dst, cst = sorted(dman), sorted(cman)

    lines = []

    def emit(s=""):
        lines.append(s)
        print(s)

    emit("=" * 112)
    emit("DEVELOPMENT vs CONFIRMATORY REPLICATION")
    emit(f"development n = {len(dst)} (exploratory, never pooled)   "
         f"confirmatory n = {len(cst)}   target overlap = 0")
    emit("=" * 112)

    results = {}

    # ---- Aim A rank-1 / top-3 / top-5
    emit("\nPRIMARY AIM A — RANK-1 / TOP-3 / TOP-5 EXPOSURE")
    emit(f"{'signal':34s} {'development':>22s} {'confirmatory':>22s}  verdict")
    for det in DETECTORS:
        for depth, lab in ((1, "rank-1"), (3, "top-3"), (5, "top-5")):
            dev = [s for s in dst if didx.get((det, s, "ORIGINAL"))]
            con = [s for s in cst if cidx.get((det, s, "ORIGINAL"))]
            dk = sum(1 for s in dev
                     if any(p["fusion_associated"] for p in didx[(det, s, "ORIGINAL")][:depth]))
            ck = sum(1 for s in con
                     if any(p["fusion_associated"] for p in cidx[(det, s, "ORIGINAL")][:depth]))
            v = verdict(dk, len(dev), ck, len(con))
            dl, dh = wilson(dk, len(dev)); cl, ch = wilson(ck, len(con))
            emit(f"{det+' '+lab:34s} {dk:3d}/{len(dev):3d}={100*dk/len(dev):5.1f}% "
                 f"{'':1s}{ck:3d}/{len(con):3d}={100*ck/len(con):5.1f}%   {v}")
            emit(f"{'':34s} [{100*dl:4.1f}-{100*dh:4.1f}]        [{100*cl:4.1f}-{100*ch:4.1f}]")
            results[f"A_{det}_{lab}"] = {"dev": (dk, len(dev)), "con": (ck, len(con)),
                                         "verdict": v}

    emit("\nAIM A BY PARTNER (rank-1)")
    emit(f"{'signal':34s} {'development':>22s} {'confirmatory':>22s}  verdict")
    for det in DETECTORS:
        for pt in PARTNERS:
            dev = [s for s in dst if dman[s]["fusion_partner"] == pt
                   and didx.get((det, s, "ORIGINAL"))]
            con = [s for s in cst if cman[s]["fusion_partner"] == pt
                   and cidx.get((det, s, "ORIGINAL"))]
            dk = sum(1 for s in dev if didx[(det, s, "ORIGINAL")][0]["fusion_associated"])
            ck = sum(1 for s in con if cidx[(det, s, "ORIGINAL")][0]["fusion_associated"])
            v = verdict(dk, len(dev), ck, len(con))
            emit(f"{det+' '+pt:34s} {dk:3d}/{len(dev):3d}={100*dk/max(len(dev),1):5.1f}% "
                 f"{'':1s}{ck:3d}/{len(con):3d}={100*ck/max(len(con),1):5.1f}%   {v}")
            results[f"A_{det}_{pt}"] = {"dev": (dk, len(dev)), "con": (ck, len(con)),
                                         "verdict": v}

    # ---- Aim B displacement, both cohorts under the frozen F3 endpoint
    emit("\nPRIMARY AIM B — F3 RANK DISPLACEMENT (development re-described with the frozen endpoint)")
    for det in DETECTORS:
        dd = f3_displacement(dman, didx, det, dst)
        cc = f3_displacement(cman, cidx, det, cst)
        de = [r for r in dd if r["state"] == "EVALUABLE"]
        ce = [r for r in cc if r["state"] == "EVALUABLE"]
        dv = [r["displacement"] for r in de]
        cv = [r["displacement"] for r in ce]
        dk = sum(1 for x in dv if x > 0)
        ck = sum(1 for x in cv if x > 0)
        clo, chi = ana.cluster_bootstrap(
            [(cman[r["pdb_id"]]["target_accession"], 1 if r["displacement"] > 0 else 0)
             for r in ce], ana.mean)
        excl0 = clo > 0
        v = verdict(dk, len(de), ck, len(ce), extra_ok=excl0)
        emit(f"\n  {det.upper()}")
        emit(f"    development : evaluable {len(de):3d}/{len(dd)}  "
             f"abstention {sum(1 for r in dd if r['state']=='DETECTOR_ABSTENTION'):2d}  "
             f"unmatched {sum(1 for r in dd if r['state']=='TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL'):2d}  "
             f"median {ana.median(dv):g}  >0 in {dk}/{len(de)} = {100*dk/max(len(de),1):.1f}%")
        emit(f"    confirmatory: evaluable {len(ce):3d}/{len(cc)}  "
             f"abstention {sum(1 for r in cc if r['state']=='DETECTOR_ABSTENTION'):2d}  "
             f"unmatched {sum(1 for r in cc if r['state']=='TARGET_CAVITY_NOT_RECOVERED_IN_ORIGINAL'):2d}  "
             f"median {ana.median(cv):g}  >0 in {ck}/{len(ce)} = {100*ck/max(len(ce),1):.1f}%")
        emit(f"    confirmatory cluster-bootstrap CI for '>0': "
             f"[{100*clo:.1f}-{100*chi:.1f}]  excludes zero: {excl0}")
        emit(f"    VERDICT: {v}")
        results[f"B_{det}"] = {"dev": (dk, len(de)), "con": (ck, len(ce)),
                               "dev_median": ana.median(dv), "con_median": ana.median(cv),
                               "con_ci": (clo, chi), "verdict": v}
        emit("    by partner (displacement > 0):")
        for pt in PARTNERS:
            ds = [r for r in de if r["partner"] == pt]
            cs = [r for r in ce if r["partner"] == pt]
            dkk = sum(1 for r in ds if r["displacement"] > 0)
            ckk = sum(1 for r in cs if r["displacement"] > 0)
            vv = verdict(dkk, len(ds), ckk, len(cs))
            emit(f"      {pt:5s} dev {dkk}/{len(ds)}  conf {ckk}/{len(cs)}   {vv}")
            results[f"B_{det}_{pt}"] = {"dev": (dkk, len(ds)), "con": (ckk, len(cs)),
                                        "verdict": vv}

    # ---- summary
    emit("\n" + "=" * 112)
    emit("REPLICATION SUMMARY")
    emit("=" * 112)
    tally = collections.Counter(v["verdict"] for v in results.values())
    for k, n in tally.most_common():
        emit(f"  {k:22s} {n}")
    emit("\n  headline signals:")
    for key in ("A_p2rank_rank-1", "A_fpocket_rank-1", "B_p2rank", "B_fpocket"):
        if key in results:
            emit(f"    {key:22s} {results[key]['verdict']}")

    json.dump(results, open(os.path.join(RES, "replication.json"), "w"), indent=1, default=str)
    open(os.path.join(RES, "replication.txt"), "w", encoding="utf-8").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
